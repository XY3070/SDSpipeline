#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKSPACE_ROOT = Path(os.environ.get("SDS_WORKSPACE_ROOT", PROJECT_ROOT.parent / "SDSworkspace")).resolve()
DEFAULT_RESULTS_ROOT = Path(os.environ.get("SDS_RESULTS_ROOT", DEFAULT_WORKSPACE_ROOT / "results")).resolve()
DEFAULT_CACHE_ROOT = Path(os.environ.get("SDS_CACHE_ROOT", DEFAULT_WORKSPACE_ROOT / "cache")).resolve()
DEFAULT_EXTERNAL_ROOT = Path(os.environ.get("SDS_EXTERNAL_ROOT", DEFAULT_WORKSPACE_ROOT / "external")).resolve()
DEFAULT_BENCH_ROOT = Path(
    os.environ.get("SDS_DEMOGRAPHY_ROOT", PROJECT_ROOT.parent / "benchmark" / "demography")
).resolve()
DEFAULT_VCF_ROOT = Path(
    os.environ.get("SDS_VCF_ROOT", DEFAULT_WORKSPACE_ROOT / "input" / "raw" / "vcf")
).resolve()
DEFAULT_SMCPP_ROOT = Path(os.environ.get("SDS_SMCPP_ROOT", DEFAULT_EXTERNAL_ROOT / "smcpp")).resolve()
DEFAULT_PROJECT_SINGULARITY = os.environ.get(
    "SDS_SMCPP_SINGULARITY_BIN",
    str(DEFAULT_EXTERNAL_ROOT / "singularity" / "bin" / "singularity"),
)
DEFAULT_PROJECT_SINGULARITY_CONFDIR = Path(
    os.environ.get("SDS_SMCPP_SINGULARITY_CONFDIR", DEFAULT_SMCPP_ROOT / "singularity_conf")
).resolve()
DEFAULT_PROJECT_SMCPP_ROOTFS = Path(
    os.environ.get("SDS_SMCPP_ROOTFS", DEFAULT_SMCPP_ROOT / "smcpp_latest.sandbox")
).resolve()
DEFAULT_PROJECT_SMCPP_MPLCONFIGDIR = Path(
    os.environ.get("SDS_SMCPP_MPLCONFIGDIR", DEFAULT_SMCPP_ROOT / "smcpp_matplotlib")
).resolve()
DEFAULT_PROJECT_SMCPP_RUNTIME_DIR = Path(
    os.environ.get("SDS_SMCPP_RUNTIME_DIR", DEFAULT_CACHE_ROOT / "smcpp_runtime")
).resolve()
DEFAULT_SMC_IMAGE = Path(os.environ.get("SDS_SMCPP_IMAGE", DEFAULT_SMCPP_ROOT / "smcpp_latest.sif")).resolve()
DEFAULT_BCFTOOLS = os.environ.get("SDS_BCFTOOLS")
if DEFAULT_BCFTOOLS is None:
    if os.environ.get("SDS_ENV_PREFIX"):
        DEFAULT_BCFTOOLS = str(Path(os.environ["SDS_ENV_PREFIX"]).expanduser() / "bin" / "bcftools")
    else:
        DEFAULT_BCFTOOLS = shutil.which("bcftools") or "bcftools"
DEFAULT_TABIX = os.environ.get("SDS_TABIX")
if DEFAULT_TABIX is None:
    if os.environ.get("SDS_ENV_PREFIX"):
        DEFAULT_TABIX = str(Path(os.environ["SDS_ENV_PREFIX"]).expanduser() / "bin" / "tabix")
    else:
        DEFAULT_TABIX = shutil.which("tabix") or "tabix"
DEFAULT_MU = 1.25e-8


def parse_chromosomes(text: str) -> list[int]:
    values: list[int] = []
    for token in text.split(","):
        item = token.strip()
        if not item:
            continue
        if "-" in item:
            start_text, end_text = item.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            if end < start:
                raise ValueError(f"Invalid chromosome range: {item}")
            values.extend(range(start, end + 1))
        else:
            values.append(int(item))
    if not values:
        raise ValueError("No chromosomes specified")
    return sorted(dict.fromkeys(values))


@dataclass(frozen=True)
class SmcppRunner:
    mode: str
    command_prefix: list[str]
    env_overrides: dict[str, str]
    metadata: dict[str, str]

    def command(self, args: list[str]) -> list[str]:
        return [*self.command_prefix, *args]


def run_checked(cmd: list[str], *, extra_env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    subprocess.run(cmd, check=True, env=env)


def read_nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def ensure_runtime_file(path: Path, content: str, *, mode: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != content:
        path.write_text(content)
    os.chmod(path, mode)
    return path


def read_distinguished_pairs(path: Path) -> list[tuple[str, str, str]]:
    pairs: list[tuple[str, str, str]] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pair_id", "sample_id_1", "sample_id_2"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Expected TSV columns {sorted(required)} in distinguished pairs file: {path}"
            )
        for row in reader:
            pairs.append((row["pair_id"], row["sample_id_1"], row["sample_id_2"]))
    if not pairs:
        raise ValueError(f"No distinguished pairs found in {path}")
    return pairs


def resolve_executable(candidate: str | Path, *, source: str) -> tuple[Path | None, str]:
    raw_value = str(candidate).strip()
    if not raw_value:
        return None, f"{source}: empty value"

    if "/" not in raw_value:
        resolved = shutil.which(raw_value)
        if resolved is None:
            return None, f"{source}: command {raw_value!r} not found on PATH"
        path = Path(resolved).resolve()
    else:
        path = Path(raw_value).expanduser().resolve()
        if not path.exists():
            return None, f"{source}: {path} does not exist"

    if not path.is_file():
        return None, f"{source}: {path} is not a file"
    if not os.access(path, os.X_OK):
        return None, f"{source}: {path} is not executable"
    return path, f"{source}: {path}"


def resolve_directory(candidate: str | Path, *, source: str) -> tuple[Path | None, str]:
    path = Path(candidate).expanduser().resolve()
    if not path.exists():
        return None, f"{source}: {path} does not exist"
    if not path.is_dir():
        return None, f"{source}: {path} is not a directory"
    return path, f"{source}: {path}"


def resolve_optional_project_directory(
    override: str | None,
    *,
    override_source: str,
    env_var: str | None,
    default_path: Path,
    default_source: str,
) -> tuple[Path | None, str | None]:
    if override:
        path, detail = resolve_directory(override, source=override_source)
        if path is None:
            raise FileNotFoundError(detail)
        return path, override_source

    if env_var and os.environ.get(env_var):
        path, detail = resolve_directory(os.environ[env_var], source=f"{env_var} environment variable")
        if path is None:
            raise FileNotFoundError(detail)
        return path, f"{env_var} environment variable"

    if default_path.exists():
        path, detail = resolve_directory(default_path, source=default_source)
        if path is None:
            raise FileNotFoundError(detail)
        return path, default_source

    return None, None


def resolve_container_runtime(override: str | None) -> tuple[Path, str]:
    if override:
        path, detail = resolve_executable(override, source="--singularity override")
        if path is None:
            raise FileNotFoundError(
                "Could not use the requested container runtime.\n"
                f"  - {detail}\n"
                "Pass a valid --singularity path or command name."
            )
        return path, "--singularity override"

    attempts: list[str] = []
    candidates: list[tuple[str, str | Path]] = []
    if os.environ.get("SINGULARITY"):
        candidates.append(("SINGULARITY environment variable", os.environ["SINGULARITY"]))
    candidates.extend(
        [
            ("PATH lookup for 'singularity'", "singularity"),
            ("PATH lookup for 'apptainer'", "apptainer"),
            ("Project fallback", DEFAULT_PROJECT_SINGULARITY),
        ]
    )

    for source, candidate in candidates:
        path, detail = resolve_executable(candidate, source=source)
        attempts.append(detail)
        if path is not None:
            return path, source

    attempted_sources = "\n".join(f"  - {detail}" for detail in attempts)
    raise FileNotFoundError(
        "Could not locate a usable container runtime for smc++.\n"
        "Tried:\n"
        f"{attempted_sources}\n"
        "Set SINGULARITY, pass --singularity /path/to/singularity, "
        "or set SMCPP_SINGULARITY before submitting the LSF wrapper."
    )


def build_rootfs_smcpp_runner(rootfs: Path, mplconfigdir: Path, rootfs_source: str) -> SmcppRunner:
    loader_candidates = sorted((rootfs / "lib" / "x86_64-linux-gnu").glob("ld-*.so"))
    if not loader_candidates:
        raise FileNotFoundError(f"No dynamic loader found under {rootfs}/lib/x86_64-linux-gnu")

    loader, detail = resolve_executable(loader_candidates[0], source="smc++ rootfs loader")
    assert loader is not None, detail
    python_bin, detail = resolve_executable(rootfs / "usr" / "bin" / "python3.8", source="smc++ rootfs python")
    assert python_bin is not None, detail
    entrypoint, detail = resolve_executable(rootfs / "usr" / "local" / "bin" / "smc++", source="smc++ rootfs entrypoint")
    assert entrypoint is not None, detail
    mpl_rc = mplconfigdir / "matplotlibrc"
    if not mpl_rc.exists():
        raise FileNotFoundError(
            f"Required matplotlibrc not found for direct smc++ runner: {mpl_rc}. "
            "Run SDSpipeline/scripts/prepare_smcpp_runtime.sh first."
        )

    library_path = ":".join(
        [
            str(rootfs / "lib" / "x86_64-linux-gnu"),
            str(rootfs / "usr" / "lib" / "x86_64-linux-gnu"),
            str(rootfs / "usr" / "lib" / "x86_64-linux-gnu" / "blas"),
            str(rootfs / "usr" / "lib" / "x86_64-linux-gnu" / "lapack"),
            str(rootfs / "usr" / "local" / "lib"),
            str(rootfs / "usr" / "lib"),
        ]
    )
    runtime_dir = DEFAULT_PROJECT_SMCPP_RUNTIME_DIR
    rootfs_pythonpath = ":".join(
        [
            str(runtime_dir),
            str(rootfs / "usr" / "local" / "lib" / "python3.8" / "dist-packages"),
            str(rootfs / "usr" / "lib" / "python3" / "dist-packages"),
        ]
    )
    env = {
        "PYTHONHOME": str(rootfs / "usr"),
        "PYTHONPATH": rootfs_pythonpath,
        "PYTHONNOUSERSITE": "1",
        "MPLCONFIGDIR": str(mplconfigdir),
        "MATPLOTLIBDATA": str(rootfs / "usr" / "share" / "matplotlib" / "mpl-data"),
        "LD_LIBRARY_PATH": library_path,
    }
    wrapper_path = ensure_runtime_file(
        DEFAULT_PROJECT_SMCPP_RUNTIME_DIR / "smcpp_rootfs_python.sh",
        "\n".join(
            [
                "#!/bin/bash",
                "set -euo pipefail",
                f"export PYTHONHOME={shlex.quote(env['PYTHONHOME'])}",
                f"export PYTHONPATH={shlex.quote(env['PYTHONPATH'])}",
                f"export PYTHONNOUSERSITE={shlex.quote(env['PYTHONNOUSERSITE'])}",
                f"export MPLCONFIGDIR={shlex.quote(env['MPLCONFIGDIR'])}",
                f"export MATPLOTLIBDATA={shlex.quote(env['MATPLOTLIBDATA'])}",
                f"export LD_LIBRARY_PATH={shlex.quote(env['LD_LIBRARY_PATH'])}",
                (
                    f"exec {shlex.quote(str(loader))} --library-path "
                    f"{shlex.quote(library_path)} {shlex.quote(str(python_bin))} \"$@\""
                ),
                "",
            ]
        ),
        mode=0o755,
    )
    sitecustomize_path = ensure_runtime_file(
        runtime_dir / "sitecustomize.py",
        "\n".join(
            [
                "import multiprocessing",
                "import os",
                "",
                "wrapper = os.environ.get('SMCPP_ROOTFS_PYTHON_WRAPPER')",
                "if wrapper:",
                "    multiprocessing.set_executable(wrapper)",
                "",
                "if os.environ.get('SMCPP_DISABLE_PROCESS_POOL', '1') == '1':",
                "    import smcpp.data_filter as _df",
                "    _df.ProcessParallelFilter.Pool = _df.ThreadParallelFilter.Pool",
                "",
            ]
        ),
        mode=0o644,
    )
    env["SMCPP_ROOTFS_PYTHON_WRAPPER"] = str(wrapper_path)
    env["SMCPP_DISABLE_PROCESS_POOL"] = "1"
    metadata = {
        "smcpp_runner_mode": "rootfs",
        "smcpp_rootfs": str(rootfs),
        "smcpp_rootfs_source": rootfs_source,
        "smcpp_loader": str(loader),
        "smcpp_python": str(python_bin),
        "smcpp_entrypoint": str(entrypoint),
        "smcpp_library_path": library_path,
        "smcpp_mplconfigdir": str(mplconfigdir),
        "smcpp_wrapper": str(wrapper_path),
        "smcpp_sitecustomize": str(sitecustomize_path),
    }
    return SmcppRunner(
        mode="rootfs",
        command_prefix=[
            str(loader),
            "--library-path",
            library_path,
            str(python_bin),
            str(entrypoint),
        ],
        env_overrides=env,
        metadata=metadata,
    )


def build_container_smcpp_runner(image: Path, singularity: Path, singularity_source: str, singularity_confdir: Path | None, singularity_confdir_source: str | None) -> SmcppRunner:
    if not image.exists():
        raise FileNotFoundError(f"Required path not found: {image}")

    env: dict[str, str] = {}
    metadata = {
        "smcpp_runner_mode": "container",
        "container_runtime": str(singularity),
        "container_runtime_source": singularity_source,
        "smcpp_image": str(image),
    }
    if singularity_confdir is not None:
        env["SINGULARITY_CONFDIR"] = str(singularity_confdir)
        metadata["singularity_confdir"] = str(singularity_confdir)
        metadata["singularity_confdir_source"] = str(singularity_confdir_source)

    return SmcppRunner(
        mode="container",
        command_prefix=[str(singularity), "exec", str(image), "smc++"],
        env_overrides=env,
        metadata=metadata,
    )


def resolve_smcpp_runner(args: argparse.Namespace, image: Path) -> SmcppRunner:
    rootfs, rootfs_source = resolve_optional_project_directory(
        args.smcpp_rootfs,
        override_source="--smcpp-rootfs override",
        env_var="SMCPP_ROOTFS",
        default_path=DEFAULT_PROJECT_SMCPP_ROOTFS,
        default_source="Project default smc++ rootfs",
    )
    if rootfs is not None:
        mplconfigdir, _ = resolve_optional_project_directory(
            args.smcpp_mplconfigdir,
            override_source="--smcpp-mplconfigdir override",
            env_var="SMCPP_MPLCONFIGDIR",
            default_path=DEFAULT_PROJECT_SMCPP_MPLCONFIGDIR,
            default_source="Project default smc++ matplotlib config",
        )
        if mplconfigdir is None:
            raise FileNotFoundError(
                "Could not locate an smc++ matplotlib config directory for the direct rootfs runner. "
                "Run SDSpipeline/scripts/prepare_smcpp_runtime.sh first or pass --smcpp-mplconfigdir."
            )
        return build_rootfs_smcpp_runner(rootfs, mplconfigdir, rootfs_source or "unknown")

    singularity, singularity_source = resolve_container_runtime(args.singularity)
    singularity_confdir, singularity_confdir_source = resolve_optional_project_directory(
        args.singularity_confdir,
        override_source="--singularity-confdir override",
        env_var="SINGULARITY_CONFDIR",
        default_path=DEFAULT_PROJECT_SINGULARITY_CONFDIR,
        default_source="Project default singularity config",
    )
    return build_container_smcpp_runner(
        image=image,
        singularity=singularity,
        singularity_source=singularity_source,
        singularity_confdir=singularity_confdir,
        singularity_confdir_source=singularity_confdir_source,
    )


def ensure_subset_vcf(
    *,
    pop: str,
    chrom: int,
    subset_samples_file: Path,
    subset_size: int,
    vcf_root: Path,
    subset_vcf_dir: Path,
    bcftools: Path,
    tabix: Path,
    threads: int,
    force: bool,
) -> Path:
    subset_vcf = subset_vcf_dir / f"UKBQC_{pop}_subset{subset_size}_chr{chrom}.vcf.gz"
    if subset_vcf.exists() and not force:
        if is_usable_bgzip_vcf(bcftools, subset_vcf):
            return subset_vcf
        print(f"[smc++] cached subset VCF is invalid, rebuilding: {subset_vcf}", flush=True)
    input_vcf = vcf_root / f"UKBQC_{pop}_chr{chrom}.vcf.gz"
    if not input_vcf.exists():
        raise FileNotFoundError(f"VCF not found: {input_vcf}")
    cmd = [
        str(bcftools),
        "view",
        "--threads",
        str(threads),
        "-S",
        str(subset_samples_file),
        "-m2",
        "-M2",
        "-v",
        "snps",
        "-O",
        "z",
        "-o",
        str(subset_vcf),
        str(input_vcf),
    ]
    run_checked(cmd)
    run_checked([str(tabix), "-f", "-p", "vcf", str(subset_vcf)])
    return subset_vcf


def infer_contig_and_length(bcftools: Path, subset_vcf: Path) -> tuple[str, int]:
    cmd = [str(bcftools), "query", "-f", "%CHROM\t%POS\n", str(subset_vcf)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    first_contig = None
    last_pos = 0
    for line in proc.stdout:
        chrom_text, pos_text = line.rstrip("\n").split("\t")
        if first_contig is None:
            first_contig = chrom_text
        last_pos = int(pos_text)
    return_code = proc.wait()
    if return_code != 0:
        raise RuntimeError(f"bcftools query failed for {subset_vcf} with exit code {return_code}")
    if first_contig is None or last_pos <= 0:
        raise ValueError(f"No variants found in subset VCF: {subset_vcf}")
    return first_contig, last_pos


def is_usable_bgzip_vcf(bcftools: Path, vcf_path: Path) -> bool:
    proc = subprocess.run(
        [str(bcftools), "view", "-h", str(vcf_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return proc.returncode == 0


def resolve_plot_csv(plot_path: Path) -> Path:
    candidates = [Path(str(plot_path) + ".csv"), plot_path.with_suffix(".csv")]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Could not find smc++ plot CSV next to {plot_path}; checked "
        + ", ".join(str(path) for path in candidates)
    )


def _write_backward_csv(plot_csv: Path, backward_csv: Path) -> None:
    """Convert smc++ plot CSV (x,y columns) to backward.py-compatible CSV (generation, mean_Ne)."""
    import csv
    rows = []
    with open(plot_csv) as fh:
        for row in csv.DictReader(fh):
            rows.append((float(row["x"]), float(row["y"])))
    rows.sort(key=lambda r: r[0])
    with open(backward_csv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["generation", "mean_Ne"])
        for x, y in rows:
            writer.writerow([x, y])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the SMC++ side of the NCN/SCN 100-sample demographic benchmark."
    )
    parser.add_argument("--pop", required=True, help="Population label, e.g. NCN or SCN.")
    parser.add_argument(
        "--subset-samples",
        default=None,
        help="Fixed subset sample file. Defaults to SDS_DEMOGRAPHY_ROOT/<POP>/subset_100.samples.txt.",
    )
    parser.add_argument(
        "--distinguished-pairs",
        default=None,
        help="TSV file with pair_id/sample_id_1/sample_id_2 columns.",
    )
    parser.add_argument(
        "--vcf-root",
        default=None,
        help="Directory containing original cohort VCFs. Defaults to SDS_VCF_ROOT/<POP>.",
    )
    parser.add_argument(
        "--subset-vcf-dir",
        default=None,
        help="Directory containing or receiving shared subset VCFs. Defaults to SDS_DEMOGRAPHY_ROOT/<POP>/subset_vcf.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for smc++ inputs and outputs. Defaults to SDS_DEMOGRAPHY_ROOT/<POP>/smcpp.",
    )
    parser.add_argument(
        "--smcpp-image",
        default=str(DEFAULT_SMC_IMAGE),
        help="Path to smcpp_latest.sif.",
    )
    parser.add_argument(
        "--singularity",
        default=None,
        help="Optional container runtime override (e.g. singularity or apptainer path/name).",
    )
    parser.add_argument(
        "--singularity-confdir",
        default=None,
        help="Optional singularity config directory override used for container-mode runs.",
    )
    parser.add_argument(
        "--smcpp-rootfs",
        default=None,
        help="Optional extracted smc++ rootfs/sandbox path for direct host execution.",
    )
    parser.add_argument(
        "--smcpp-mplconfigdir",
        default=None,
        help="Optional matplotlib config directory for direct smc++ rootfs execution.",
    )
    parser.add_argument(
        "--bcftools",
        default=str(DEFAULT_BCFTOOLS),
        help="Path to bcftools. Defaults to SDS_BCFTOOLS or the configured SDS env bcftools.",
    )
    parser.add_argument(
        "--tabix",
        default=str(DEFAULT_TABIX),
        help="Path to tabix. Defaults to SDS_TABIX or the configured SDS env tabix.",
    )
    parser.add_argument(
        "--chromosomes",
        default="1-22",
        help="Chromosome list/ranges to use, e.g. 1-22 or 1,2,22.",
    )
    parser.add_argument("--mu", type=float, default=DEFAULT_MU, help="Mutation rate per base per generation.")
    parser.add_argument("--cores", type=int, default=4, help="Cores passed to smc++ and bcftools.")
    parser.add_argument("--window-size", type=int, default=20, help="Window size passed to smc++ estimate.")
    parser.add_argument(
        "--knots",
        type=int,
        default=8,
        help="Number of spline knots for smc++ estimate (default: %(default)s).",
    )
    parser.add_argument(
        "--spline",
        choices=["cubic", "pchip", "piecewise"],
        default="piecewise",
        help="Spline type for smc++ estimate (default: %(default)s).",
    )
    parser.add_argument(
        "--timepoints",
        type=float,
        nargs=2,
        default=[10.0, 100000.0],
        metavar=("T1", "TK"),
        help="smc++ estimate --timepoints in generations (default: %(default)s).",
    )
    parser.add_argument(
        "--em-iterations",
        type=int,
        default=None,
        help="Optional smc++ estimate --em-iterations override.",
    )
    parser.add_argument(
        "--missing-cutoff",
        type=int,
        default=None,
        help="Optional smc++ vcf2smc --missing-cutoff value.",
    )
    parser.add_argument(
        "--nonseg-cutoff",
        type=int,
        default=None,
        help="Optional smc++ estimate --nonseg-cutoff value.",
    )
    parser.add_argument("--base", default=None, help="Base name for smc++ estimate outputs. Defaults to <POP>.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing subset VCFs, .smc.gz files, and rerun estimate/plot.",
    )
    parser.add_argument(
        "--probe-only",
        action="store_true",
        help="Resolve the configured smc++ runner and execute 'smc++ version', then exit.",
    )
    args = parser.parse_args()
    args.pop = args.pop.upper()
    benchmark_root = DEFAULT_BENCH_ROOT / args.pop
    if args.subset_samples is None:
        args.subset_samples = str(benchmark_root / "subset_100.samples.txt")
    if args.distinguished_pairs is None:
        args.distinguished_pairs = str(benchmark_root / "smcpp_distinguished_pairs.tsv")
    if args.vcf_root is None:
        args.vcf_root = str(DEFAULT_VCF_ROOT / args.pop)
    if args.subset_vcf_dir is None:
        args.subset_vcf_dir = str(benchmark_root / "subset_vcf")
    if args.output_dir is None:
        args.output_dir = str(benchmark_root / "smcpp")
    if args.base is None:
        args.base = args.pop
    args.chromosomes = parse_chromosomes(args.chromosomes)
    return args


def main() -> int:
    args = parse_args()
    subset_samples_file = Path(args.subset_samples).resolve()
    distinguished_pairs_file = Path(args.distinguished_pairs).resolve()
    vcf_root = Path(args.vcf_root).resolve()
    subset_vcf_dir = Path(args.subset_vcf_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    image = Path(args.smcpp_image).resolve()
    bcftools, bcftools_detail = resolve_executable(args.bcftools, source="bcftools runtime")
    if bcftools is None:
        raise FileNotFoundError(bcftools_detail)
    tabix, tabix_detail = resolve_executable(args.tabix, source="tabix runtime")
    if tabix is None:
        raise FileNotFoundError(tabix_detail)
    runner = resolve_smcpp_runner(args, image)

    print(f"[smc++] runner mode={runner.mode}", flush=True)
    for key in sorted(runner.metadata):
        print(f"[smc++] {key}={runner.metadata[key]}", flush=True)

    if args.probe_only:
        run_checked(runner.command(["version"]), extra_env=runner.env_overrides)
        return 0

    for path in [
        subset_samples_file,
        distinguished_pairs_file,
        vcf_root,
        bcftools,
        tabix,
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Required path not found: {path}")

    subset_samples = read_nonempty_lines(subset_samples_file)
    distinguished_pairs = read_distinguished_pairs(distinguished_pairs_file)
    subset_size = len(subset_samples)
    subset_sample_set = set(subset_samples)
    for pair_id, sample_a, sample_b in distinguished_pairs:
        if sample_a not in subset_sample_set or sample_b not in subset_sample_set:
            raise ValueError(
                f"Distinguished pair {pair_id} references samples outside the subset: "
                f"{sample_a}, {sample_b}"
            )

    subset_vcf_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    smc_input_dir = output_dir / "input_smc"
    analysis_dir = output_dir / "analysis"
    smc_input_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    population_arg = f"{args.pop}:{','.join(subset_samples)}"
    smc_files: list[str] = []
    subset_vcfs: list[str] = []
    contig_lengths: dict[str, dict[str, object]] = {}

    for chrom in args.chromosomes:
        subset_vcf = ensure_subset_vcf(
            pop=args.pop,
            chrom=chrom,
            subset_samples_file=subset_samples_file,
            subset_size=subset_size,
            vcf_root=vcf_root,
            subset_vcf_dir=subset_vcf_dir,
            bcftools=bcftools,
            tabix=tabix,
            threads=args.cores,
            force=args.force,
        )
        subset_vcfs.append(str(subset_vcf))
        contig_name, contig_length = infer_contig_and_length(bcftools, subset_vcf)
        contig_lengths[f"chr{chrom}"] = {
            "subset_vcf": str(subset_vcf),
            "contig_name": contig_name,
            "contig_length": contig_length,
        }
        for pair_id, sample_a, sample_b in distinguished_pairs:
            smc_path = smc_input_dir / f"{args.pop}_chr{chrom}_{pair_id}.smc.gz"
            smc_files.append(str(smc_path))
            if smc_path.exists() and not args.force:
                continue
            cmd = runner.command(
                [
                "vcf2smc",
                "--cores",
                str(args.cores),
                "-d",
                sample_a,
                sample_b,
                "--length",
                str(contig_length),
                "--ignore-missing",
                ]
            )
            if args.missing_cutoff is not None:
                cmd.extend(["--missing-cutoff", str(args.missing_cutoff)])
            cmd.extend(
                [
                    str(subset_vcf),
                    str(smc_path),
                    contig_name,
                    population_arg,
                ]
            )
            run_checked(cmd, extra_env=runner.env_overrides)

    model_json = analysis_dir / f"{args.base}.final.json"
    if args.force or not model_json.exists():
        cmd = runner.command(
            [
            "estimate",
            "--cores",
            str(args.cores),
            "-o",
            str(analysis_dir),
            "--base",
            args.base,
            "-w",
            str(args.window_size),
            "--knots",
            str(args.knots),
            "--spline",
            args.spline,
            "--timepoints",
            str(args.timepoints[0]),
            str(args.timepoints[1]),
            ]
        )
        if args.em_iterations is not None:
            cmd.extend(["--em-iterations", str(args.em_iterations)])
        if args.nonseg_cutoff is not None:
            cmd.extend(["--nonseg-cutoff", str(args.nonseg_cutoff)])
        cmd.extend([str(args.mu), *smc_files])
        run_checked(cmd, extra_env=runner.env_overrides)

    plot_path = output_dir / f"{args.base}_smcpp.png"
    if args.force or not plot_path.exists():
        cmd = runner.command(
            [
            "plot",
            "--csv",
            str(plot_path),
            str(model_json),
            ]
        )
        run_checked(cmd, extra_env=runner.env_overrides)
    plot_csv = resolve_plot_csv(plot_path)

    # Generate backward.py-compatible CSV from plot CSV
    backward_csv = output_dir / f"{args.base}_backward.csv"
    if args.force or not backward_csv.exists():
        _write_backward_csv(plot_csv, backward_csv)

    manifest = {
        "population": args.pop,
        "subset_samples_file": str(subset_samples_file),
        "distinguished_pairs_file": str(distinguished_pairs_file),
        "subset_size": subset_size,
        "chromosomes": args.chromosomes,
        "mu": args.mu,
        "cores": args.cores,
        "window_size": args.window_size,
        "knots": args.knots,
        "spline": args.spline,
        "timepoints": args.timepoints,
        "em_iterations": args.em_iterations,
        "missing_cutoff": args.missing_cutoff,
        "nonseg_cutoff": args.nonseg_cutoff,
        "subset_vcf_dir": str(subset_vcf_dir),
        "subset_vcfs": subset_vcfs,
        "smc_input_dir": str(smc_input_dir),
        "smc_files": smc_files,
        "analysis_dir": str(analysis_dir),
        "model_json": str(model_json),
        "plot_path": str(plot_path),
        "plot_csv": str(plot_csv),
        "backward_csv": str(backward_csv),
        "contigs": contig_lengths,
    }
    manifest.update(runner.metadata)
    (output_dir / f"{args.base}_smcpp_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"[Error] {exc}", file=sys.stderr)
        raise
