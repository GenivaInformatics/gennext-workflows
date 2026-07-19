#!/usr/bin/env python3
"""cnv-cohort-builder: lazy per-kit CNVkit pooled reference builder.

Adapts the offline build_kit_reference.py from the cnv-ref-dataset project to
run inside an Argo pod. Differences from the original:

  - cnvkit.py is invoked directly (the bioconda binary is on PATH in this image),
    not via nested Docker.
  - Multi-track D4 reads use the `d4tools` CLI (also on PATH) instead of pyd4,
    because pyd4 only ships sdists and needs Rust to build.
  - All filesystem paths are supplied via flags. No hardcoded cohort layout.
  - Builds are flock-serialised on <kits-dir>/<kit-name>/.build.lock so that
    concurrent first-time-on-this-kit invocations don't double-build.
  - Idempotent: if pooled_reference.cnn already exists, exit 0 immediately.

Usage:
  cnv-cohort-builder ensure-kit \\
      --kit-bed   /ref/regions/agilent_cgp.bed \\
      --cohort-d4 /cohort/cnv_cohort.d4 \\
      --ref-fa    /ref/genome/GRCh38.fa \\
      --ref-flat  /ref/refflat/refFlat.txt \\
      --kits-dir  /kits \\
      [--kit-name <name>]              # default: <kit-bed basename minus .bed>
      [--mode hybrid|amplicon|auto]    # default: auto
      [--workers 12]
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

# ----------------------------------------------------------------------------
# cnvkit / d4tools shell-outs
# ----------------------------------------------------------------------------


def cnvkit(args: list[str], log: Path | None = None, timeout: int = 14400) -> None:
    """Run a cnvkit.py subcommand. Captures and logs output, raises on non-zero."""
    cmd = ["cnvkit.py", *args]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            f"# cmd: {' '.join(cmd)}\n# rc: {res.returncode}\n"
            f"# stdout:\n{res.stdout}\n# stderr:\n{res.stderr}\n"
        )
    if res.returncode != 0:
        last = res.stderr.strip().splitlines()[-1] if res.stderr.strip() else "rc!=0"
        raise RuntimeError(f"cnvkit {args[0]} failed: {last}")


def d4_list_tracks(d4: Path) -> list[str]:
    """Return the bare list of track names in a multi-track D4 file.

    `d4tools ls-track <file>` emits lines of the form `<file-path>:<track-name>`
    (one per line). Strip the leading `<file-path>:` prefix so callers get just
    the track names — required because `d4tools stat -F '^<track>$' <file>`
    matches against the bare track-label, not the file:track string.
    """
    res = subprocess.run(
        ["d4tools", "ls-track", str(d4)],
        capture_output=True, text=True, check=True,
    )
    file_prefix = f"{d4}:"
    tracks: list[str] = []
    for line in res.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        token = line.split()[0]
        if token.startswith(file_prefix):
            token = token[len(file_prefix):]
        # Older d4tools versions just print the bare track name with no prefix;
        # also tolerate a leading `:` or square-bracket framing.
        tracks.append(token.lstrip(":").strip("[]"))
    if not tracks:
        raise RuntimeError(f"d4tools ls-track produced no tracks for {d4}")
    return tracks


def d4_region_means(d4: Path, track: str, bed: Path) -> dict[tuple[str, int, int], float]:
    """Per-BED-region mean depth for one track of a multi-track D4 file.

    Uses `d4tools stat -F '^<track>$' -r <bed> -s mean` so only the requested
    track is processed. Output is parsed as a tab-separated stream with
    chrom/start/end and a mean column; the position of the mean column
    varies between d4tools versions when a track-name column is included,
    so we sniff it from the row width.
    """
    cmd = [
        "d4tools", "stat",
        "-F", f"^{re.escape(track)}$",
        "-r", str(bed),
        "-s", "mean",
        str(d4),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    if res.returncode != 0:
        last = res.stderr.strip().splitlines()[-1] if res.stderr.strip() else "rc!=0"
        raise RuntimeError(f"d4tools stat (track={track}) failed: {last}")

    out: dict[tuple[str, int, int], float] = {}
    for line in res.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 4:
            continue
        try:
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
        except ValueError:
            # Header row from -H, or unexpected line — skip.
            continue
        # Mean is the last numeric column; len-1 covers both 4-col
        # (chrom/start/end/mean) and 5-col (chrom/start/end/track/mean).
        try:
            mean = float(cols[-1])
        except ValueError:
            continue
        out[(chrom, start, end)] = mean
    return out


# ----------------------------------------------------------------------------
# BED + CNN helpers (lifted near-verbatim from build_kit_reference.py)
# ----------------------------------------------------------------------------


def detect_panel_mode(bed: Path, mb_threshold: float = 1.0) -> tuple[str, dict]:
    """<1 Mb merged-BED footprint ⇒ amplicon, else hybrid. See cnv-ref-dataset
    docs/gotchas for why this single-signal split is reliable in practice."""
    rows: list[tuple[str, int, int]] = []
    with bed.open() as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                rows.append((parts[0], int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    rows.sort()
    if not rows:
        raise RuntimeError(f"empty BED: {bed}")
    merged_bp = 0
    cur_chrom, cur_s, cur_e = rows[0]
    for chrom, s, e in rows[1:]:
        if chrom == cur_chrom and s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged_bp += cur_e - cur_s
            cur_chrom, cur_s, cur_e = chrom, s, e
    merged_bp += cur_e - cur_s
    total_mb = merged_bp / 1e6
    mode = "amplicon" if total_mb < mb_threshold else "hybrid"
    return mode, {"n_regions": len(rows), "total_mb": round(total_mb, 3),
                  "mb_threshold": mb_threshold}


def load_bins(bed: Path, default_gene: str | None) -> list[tuple[str, int, int, str]]:
    """Parse a CNVkit-format BED into (chrom, start, end, gene) rows in BED order.

    For cnvkit target output, gene is in column 4. For antitarget the literal
    string 'Antitarget' is the convention; pass default_gene='Antitarget' to
    force it when the col is missing.
    """
    out: list[tuple[str, int, int, str]] = []
    with bed.open() as f:
        for line in f:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = (default_gene if default_gene is not None
                    else (parts[3] if len(parts) >= 4 else "-"))
            out.append((chrom, start, end, gene))
    return out


def write_cnn(path: Path, bins: list[tuple[str, int, int, str]],
              depths: list[float]) -> None:
    """Write a CNVkit-format coverage .cnn (chromosome/start/end/gene/depth/log2)."""
    with path.open("w") as f:
        f.write("chromosome\tstart\tend\tgene\tdepth\tlog2\n")
        for (chrom, s, e, gene), d in zip(bins, depths):
            log2 = math.log2(d) if d > 0.0 else -20.0
            f.write(f"{chrom}\t{s}\t{e}\t{gene}\t{d:.6g}\t{log2:.6g}\n")


# ----------------------------------------------------------------------------
# Per-sample synthesis: D4 -> two .cnn files
# ----------------------------------------------------------------------------


def synthesize_cnn_for_sample(
    sample: str,
    cohort_d4: Path,
    targets_bed: Path,
    antitargets_bed: Path | None,
    tgt_bins: list[tuple[str, int, int, str]],
    anti_bins: list[tuple[str, int, int, str]],
    out_dir: Path,
    skip_existing: bool,
) -> tuple[str, bool, float, str]:
    """Per-sample worker: extract per-bin mean depth from cohort.d4 and write
    target+antitarget .cnn files in the kit's coverage dir."""
    out_t = out_dir / f"{sample}.targetcoverage.cnn"
    out_a = out_dir / f"{sample}.antitargetcoverage.cnn"
    has_anti = bool(anti_bins) and antitargets_bed is not None
    if skip_existing and out_t.exists() and (not has_anti or out_a.exists()):
        return sample, True, 0.0, "skipped"

    started = time.monotonic()
    try:
        tgt_means = d4_region_means(cohort_d4, sample, targets_bed)
        tgt_depths = [tgt_means.get((c, s, e), 0.0) for (c, s, e, _) in tgt_bins]
        write_cnn(out_t, tgt_bins, tgt_depths)
        if has_anti:
            anti_means = d4_region_means(cohort_d4, sample, antitargets_bed)
            anti_depths = [anti_means.get((c, s, e), 0.0) for (c, s, e, _) in anti_bins]
            write_cnn(out_a, anti_bins, anti_depths)
    except Exception as e:
        return sample, False, round(time.monotonic() - started, 2), f"{type(e).__name__}: {e}"
    return sample, True, round(time.monotonic() - started, 2), ""


# ----------------------------------------------------------------------------
# Phases
# ----------------------------------------------------------------------------


def build_target_beds(kit_bed: Path, ref_flat: Path, access_bed: Path,
                      kit_dir: Path, mode: str, log_dir: Path) -> tuple[Path, Path]:
    """Run cnvkit target [+ antitarget]. Amplicon mode writes an empty
    antitargets.bed so the rest of the pipeline still finds the file but no
    antitarget bins are processed (matches cnv-ref-dataset's amplicon path)."""
    targets = kit_dir / "targets.bed"
    antitargets = kit_dir / "antitargets.bed"

    print(f"  cnvkit target  -> {targets.name}", flush=True)
    cnvkit(
        ["target", str(kit_bed), "--annotate", str(ref_flat), "--split",
         "-o", str(targets)],
        log=log_dir / f"cnvkit_target_{kit_dir.name}.log",
    )

    if mode == "amplicon":
        print(f"  amplicon mode: writing empty {antitargets.name}", flush=True)
        antitargets.write_text("")
    else:
        print(f"  cnvkit antitarget -> {antitargets.name}", flush=True)
        cnvkit(
            ["antitarget", str(targets), "-g", str(access_bed),
             "-o", str(antitargets)],
            log=log_dir / f"cnvkit_antitarget_{kit_dir.name}.log",
        )
    return targets, antitargets


def ensure_access_bed(ref_fa: Path, kit_dir: Path, log_dir: Path) -> Path:
    """Generate (or reuse) per-genome accessibility BED. Cached at
    <kits-dir>/.access/<ref-fa-basename>.access.bed so repeated kit builds
    against the same genome don't re-run cnvkit access (~30s)."""
    cache_dir = kit_dir.parent / ".access"
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{ref_fa.name}.access.bed"
    if out.exists() and out.stat().st_size > 0:
        return out
    print(f"  cnvkit access -> {out} (one-time per genome)", flush=True)
    cnvkit(
        ["access", str(ref_fa), "-o", str(out)],
        log=log_dir / "cnvkit_access.log",
    )
    return out


def build_pooled_reference(kit_dir: Path, ref_fa: Path, log_dir: Path) -> Path:
    """Pool the per-cohort-sample cnns into pooled_reference.cnn."""
    coverage_dir = kit_dir / "coverage"
    cnns = sorted(coverage_dir.glob("*.cnn"))
    if not cnns:
        raise RuntimeError(f"no cnn files in {coverage_dir}")
    out = kit_dir / "pooled_reference.cnn"
    print(f"  cnvkit reference ({len(cnns)} cnns) -> {out.name}", flush=True)
    cnvkit(
        ["reference", *[str(p) for p in cnns], "-f", str(ref_fa), "-o", str(out)],
        log=log_dir / f"cnvkit_reference_{kit_dir.name}.log",
    )
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cnvkit_version() -> str:
    res = subprocess.run(["cnvkit.py", "version"], capture_output=True, text=True, check=True)
    return res.stdout.strip()


def d4tools_version() -> str:
    res = subprocess.run(["d4tools", "--help"], capture_output=True, text=True)
    # `d4tools --help` prints the version on the first line; fall back to "unknown".
    line = res.stdout.splitlines()[0] if res.stdout else ""
    m = re.search(r"(\d+\.\d+\.\d+)", line)
    return m.group(1) if m else "unknown"


def write_metadata(kit_dir: Path, kit_bed: Path, cohort_d4: Path, samples: list[str],
                   ref_fa: Path, ref_flat: Path, mode: str) -> None:
    fai = ref_fa.with_suffix(ref_fa.suffix + ".fai")
    meta = {
        "kit_name": kit_dir.name,
        "kit_bed": str(kit_bed),
        "kit_bed_sha256": sha256_file(kit_bed),
        "cohort_d4": str(cohort_d4),
        "cohort_d4_sha256": sha256_file(cohort_d4),
        "n_samples": len(samples),
        "sample_ids": samples,
        "ref_fa": str(ref_fa),
        "ref_fa_fai_sha256": sha256_file(fai) if fai.exists() else "",
        "ref_flat": str(ref_flat),
        "panel_mode": mode,
        "cnvkit_version": cnvkit_version(),
        "d4tools_version": d4tools_version(),
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (kit_dir / "build.json").write_text(json.dumps(meta, indent=2))


# ----------------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------------


def cmd_ensure_kit(args: argparse.Namespace) -> int:
    kit_bed = Path(args.kit_bed).resolve()
    cohort_d4 = Path(args.cohort_d4).resolve()
    ref_fa = Path(args.ref_fa).resolve()
    ref_flat = Path(args.ref_flat).resolve()
    kits_dir = Path(args.kits_dir).resolve()

    for p, label in [(kit_bed, "kit BED"), (cohort_d4, "cohort D4"),
                     (ref_fa, "reference FASTA"), (ref_flat, "refFlat"),
                     (kits_dir, "kits dir")]:
        if not p.exists():
            print(f"{label} not found: {p}", file=sys.stderr)
            return 2

    kit_name = args.kit_name or kit_bed.stem
    kit_dir = kits_dir / kit_name
    kit_dir.mkdir(parents=True, exist_ok=True)
    coverage_dir = kit_dir / "coverage"
    coverage_dir.mkdir(parents=True, exist_ok=True)
    log_dir = kit_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    pooled_ref = kit_dir / "pooled_reference.cnn"

    # flock to serialise concurrent first-time builds for the same kit.
    lock_path = kit_dir / ".build.lock"
    print(f"Acquiring build lock {lock_path} ...", flush=True)
    with lock_path.open("w") as lock_f:
        fcntl.flock(lock_f, fcntl.LOCK_EX)
        if pooled_ref.exists() and pooled_ref.stat().st_size > 0:
            print(f"Kit reference already exists at {pooled_ref}; nothing to do.")
            return 0
        return _build_under_lock(
            kit_name, kit_dir, coverage_dir, log_dir,
            kit_bed, cohort_d4, ref_fa, ref_flat,
            mode_arg=args.mode, workers=args.workers,
            skip_existing=args.skip_existing,
        )


def _build_under_lock(
    kit_name: str, kit_dir: Path, coverage_dir: Path, log_dir: Path,
    kit_bed: Path, cohort_d4: Path, ref_fa: Path, ref_flat: Path,
    mode_arg: str, workers: int, skip_existing: bool,
) -> int:
    """Lock is held for the duration of this call."""
    if mode_arg == "auto":
        mode, stats = detect_panel_mode(kit_bed)
        print(f"Detect panel mode from BED ({kit_bed.name}):")
        print(f"  regions={stats['n_regions']}, total={stats['total_mb']} Mb "
              f"(threshold {stats['mb_threshold']} Mb)")
        print(f"  -> mode = {mode}\n")
    else:
        mode = mode_arg

    print(f"Kit:    {kit_name}")
    print(f"  BED:  {kit_bed}")
    print(f"  Mode: {mode}")
    print(f"  Out:  {kit_dir}\n")

    # Phase 1
    print(f"[1/4] cnvkit access + target{' + antitarget' if mode == 'hybrid' else ' (amplicon mode: no antitarget)'}")
    access_bed = ensure_access_bed(ref_fa, kit_dir, log_dir)
    targets_bed, antitargets_bed = build_target_beds(
        kit_bed, ref_flat, access_bed, kit_dir, mode, log_dir,
    )
    n_tgt = sum(1 for _ in targets_bed.open())
    n_anti = (sum(1 for _ in antitargets_bed.open())
              if antitargets_bed.stat().st_size > 0 else 0)
    print(f"  targets: {n_tgt} bins, antitargets: {n_anti} bins")

    # Phase 2
    samples = d4_list_tracks(cohort_d4)
    print(f"\n[2/4] D4 -> cnn synthesis ({len(samples)} samples × "
          f"{2 if mode == 'hybrid' else 1} cnn each)")
    tgt_bins = load_bins(targets_bed, default_gene=None)
    anti_bins = (load_bins(antitargets_bed, default_gene="Antitarget")
                 if mode == "hybrid" else [])
    print(f"  parsed {len(tgt_bins)} target bins, {len(anti_bins)} antitarget bins")

    started = time.monotonic()
    fails: list[str] = []
    # ThreadPoolExecutor (not Process) because each task is dominated by the
    # d4tools subprocess, which already runs in its own process.
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {
            ex.submit(
                synthesize_cnn_for_sample,
                s, cohort_d4, targets_bed,
                antitargets_bed if mode == "hybrid" else None,
                tgt_bins, anti_bins, coverage_dir, skip_existing,
            ): s for s in samples
        }
        for i, fut in enumerate(as_completed(futs), 1):
            sample, ok, elapsed, err = fut.result()
            tag = "OK " if ok else "ERR"
            print(f"  [{i:3d}/{len(samples)}] {tag} {sample:<22s} in "
                  f"{elapsed:6.2f}s" + (f"  [{err}]" if err else ""), flush=True)
            if not ok:
                fails.append(sample)
    print(f"  cnn synthesis elapsed {time.monotonic() - started:.1f}s")
    if fails:
        print(f"  Failed: {' '.join(fails)}", file=sys.stderr)
        return 1

    # Phase 3
    print("\n[3/4] cnvkit reference (pool cohort cnns)")
    started = time.monotonic()
    pooled = build_pooled_reference(kit_dir, ref_fa, log_dir)
    size_mb = pooled.stat().st_size / 1024 / 1024
    print(f"  -> {pooled} ({size_mb:.1f} MB) in {time.monotonic() - started:.1f}s")

    # Phase 4
    print("\n[4/4] build.json")
    write_metadata(kit_dir, kit_bed, cohort_d4, samples, ref_fa, ref_flat, mode)
    print(f"  -> {kit_dir / 'build.json'}")
    print(f"\nDone. Kit reference: {kit_dir}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(prog="cnv-cohort-builder")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pe = sub.add_parser("ensure-kit", help="Build the per-kit pooled reference if missing.")
    pe.add_argument("--kit-bed", required=True, help="Path to the kit's regions BED")
    pe.add_argument("--cohort-d4", required=True,
                    help="Path to the multi-track cohort D4 file")
    pe.add_argument("--ref-fa", required=True, help="Path to the reference genome FASTA")
    pe.add_argument("--ref-flat", required=True,
                    help="Path to refFlat (used by cnvkit target --annotate)")
    pe.add_argument("--kits-dir", required=True,
                    help="Directory containing per-kit subdirs; <kit-name>/ "
                         "is created here.")
    pe.add_argument("--kit-name", default=None,
                    help="Kit subdir name. Defaults to <kit-bed basename minus .bed>.")
    pe.add_argument("--mode", choices=["hybrid", "amplicon", "auto"], default="auto",
                    help="Panel chemistry. auto = pick from BED footprint "
                         "(<1 Mb merged ⇒ amplicon, else hybrid).")
    pe.add_argument("--workers", type=int, default=12,
                    help="Concurrent per-sample cnn synthesis workers.")
    pe.add_argument("--skip-existing", action="store_true",
                    help="Skip per-sample synthesis when the cnn already exists.")
    pe.set_defaults(func=cmd_ensure_kit)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
