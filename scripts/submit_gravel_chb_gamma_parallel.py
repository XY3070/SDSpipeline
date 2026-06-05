#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

SDS_PYTHON = Path(os.environ.get("SDS_PYTHON", sys.executable)).expanduser().resolve()
if Path(sys.executable).resolve() != SDS_PYTHON:
    os.execv(str(SDS_PYTHON), [str(SDS_PYTHON), str(Path(__file__).resolve()), *sys.argv[1:]])


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKSPACE_ROOT = Path(os.environ.get("SDS_WORKSPACE_ROOT", REPO_ROOT.parent / "SDSworkspace")).expanduser().resolve()
DEFAULT_RESULTS_ROOT = Path(os.environ.get("SDS_RESULTS_ROOT", DEFAULT_WORKSPACE_ROOT / "results")).expanduser().resolve()
DEFAULT_RUNS_ROOT = Path(os.environ.get("SDS_RUNS_ROOT", DEFAULT_WORKSPACE_ROOT / "runs")).expanduser().resolve()
DEFAULT_EXTERNAL_ROOT = Path(os.environ.get("SDS_EXTERNAL_ROOT", DEFAULT_WORKSPACE_ROOT / "external")).expanduser().resolve()
DEFAULT_GAMMA_ROOT = Path(os.environ.get("SDS_GAMMA_ROOT", DEFAULT_RESULTS_ROOT / "production" / "gamma")).expanduser().resolve()
DEFAULT_DEMOGRAPHY_ROOT = Path(os.environ.get("SDS_DEMOGRAPHY_ROOT", DEFAULT_RESULTS_ROOT / "production" / "demography")).expanduser().resolve()
DEFAULT_MS_ROOT = Path(os.environ.get("SDS_MS_ROOT", DEFAULT_EXTERNAL_ROOT / "ms")).expanduser().resolve()
RUN_GAMMA_CHUNK_SCRIPT = REPO_ROOT / "scripts" / "run_single_snp_gamma_chunk.sh"
AGGREGATE_GAMMA_CHUNKS_SCRIPT = REPO_ROOT / "scripts" / "aggregate_single_snp_gamma_chunks.sh"
DEFAULT_OUTDIR = DEFAULT_GAMMA_ROOT / "gravel_chb_ne100k"
DEFAULT_DUMMY_NPZ = Path(
    os.environ.get(
        "SDS_DUMMY_NPZ",
        DEFAULT_DEMOGRAPHY_ROOT / "scenarios" / "NCN_scaled_ne0_100000.npz",
    )
).expanduser().resolve()
DEFAULT_MS_MAKE_DIR = Path(os.environ.get("SDS_MS_SCRIPTS_DIR", DEFAULT_MS_ROOT / "scripts")).expanduser().resolve()
DEFAULT_MS_BINARY = Path(os.environ.get("SDS_MS_BINARY", DEFAULT_MS_ROOT / "msdir" / "ms")).expanduser().resolve()
DEFAULT_BACKWARD_SCRIPT = Path(
    os.environ.get("SDS_BACKWARD_SCRIPT", DEFAULT_MS_MAKE_DIR / "backward.py")
).expanduser().resolve()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit chunked gamma generation for the legacy/default East-Asian model via Gravel_CHB."
    )
    parser.add_argument("--run-task-file", default=None, help="Internal worker mode: JSONL task file written by the submitter.")
    parser.add_argument("--task-index", type=int, default=None, help="1-based task index for --run-task-file.")
    parser.add_argument("--dry-run", action="store_true", help="Prepare task files and manifests without submitting jobs.")
    parser.add_argument("--queue", choices=["auto", "smp", "normal"], default="auto", help="Queue selection strategy.")
    parser.add_argument("--chunk-array-cap", type=int, default=96, help="Maximum concurrent chunk-gamma array slots.")
    parser.add_argument("--aggregate-array-cap", type=int, default=24, help="Maximum concurrent aggregate array slots.")
    parser.add_argument("--job-prefix", default="gravel_chb_ne100k_gamma", help="LSF job-name prefix.")
    parser.add_argument("--log-dir", default=str(DEFAULT_RUNS_ROOT / "gamma" / "logs"), help="Directory for LSF stdout/stderr logs.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory for chunk roots, gamma pieces, and manifests.")
    parser.add_argument("--present-ne", type=int, default=100000, help="Present-day diploid population size passed to gamma generation.")
    parser.add_argument("--sample-size", type=int, default=None, help="Optional haplotype sample size override passed into the gamma simulation Makefile.")
    parser.add_argument("--pop-model", default="Gravel_CHB", help="Backward model passed to the chunk script. Defaults to Gravel_CHB.")
    parser.add_argument("--dummy-npz", default=str(DEFAULT_DUMMY_NPZ), help="Existing NPZ path. Ignored by Gravel_CHB but required by the chunk script.")
    parser.add_argument("--sim-reps", type=int, default=1000, help="Total neutral replicates per DAF.")
    parser.add_argument("--chunk-size", type=int, default=100, help="Replicates per chunk.")
    parser.add_argument("--daf-start", type=float, default=0.05, help="Start of the DAF grid.")
    parser.add_argument("--daf-end", type=float, default=0.95, help="End of the DAF grid.")
    parser.add_argument("--daf-step", type=float, default=0.01, help="Step size of the DAF grid.")
    parser.add_argument("--ms-make-dir", default=str(DEFAULT_MS_MAKE_DIR), help="Directory containing the MS Makefile.")
    parser.add_argument("--ms-binary", default=str(DEFAULT_MS_BINARY), help="Path to the ms binary.")
    parser.add_argument("--backward-script", default=str(DEFAULT_BACKWARD_SCRIPT), help="Path to backward.py.")
    parser.add_argument("--no-reuse-existing", action="store_true", help="Force regeneration even when canonical gamma pieces already exist.")
    parser.add_argument("--chunk-worker", action="store_true", help="Generate one gamma chunk and exit.")
    parser.add_argument("--aggregate-worker", action="store_true", help="Aggregate one DAF worth of chunk outputs and exit.")
    parser.add_argument("--finalize-only", action="store_true", help="Write the final g_file from completed canonical gamma pieces and exit.")
    parser.add_argument("--daf", default=None, help="DAF string used by worker modes, e.g. 0.05.")
    parser.add_argument("--chunk-index", type=int, default=None, help="Chunk index used by --chunk-worker.")
    parser.add_argument("--chunk-start-rep", type=int, default=None, help="Inclusive replication start used by --chunk-worker.")
    parser.add_argument("--chunk-end-rep", type=int, default=None, help="Inclusive replication end used by --chunk-worker.")
    return parser.parse_args(argv)


def chunk_ranges(sim_reps: int, chunk_size: int) -> list[tuple[int, int, int]]:
    if sim_reps <= 0 or chunk_size <= 0:
        raise RuntimeError("sim_reps and chunk_size must be positive")
    ranges: list[tuple[int, int, int]] = []
    start = 1
    index = 1
    while start <= sim_reps:
        end = min(sim_reps, start + chunk_size - 1)
        ranges.append((index, start, end))
        start = end + 1
        index += 1
    return ranges


def build_daf_grid(start: float, end: float, step: float) -> list[tuple[float, str]]:
    if step <= 0.0:
        raise RuntimeError("daf_step must be positive")
    values: list[tuple[float, str]] = []
    current = start
    while current <= end + 1e-12:
        rounded = round(current + 1e-12, 2)
        values.append((float(f"{rounded:.2f}"), f"{rounded:.2f}"))
        current += step
    if not values:
        raise RuntimeError("Empty DAF grid")
    return values


def frequency_file_label(daf_label: str) -> str:
    return daf_label.replace(".", "p").replace("-", "m")


def chunk_root_for(outdir: Path, daf_label: str) -> Path:
    return outdir / "chunks" / frequency_file_label(daf_label)


def chunk_workdir_for(chunk_root: Path, chunk_index: int) -> Path:
    return chunk_root / f"chunk_{int(chunk_index):03d}"


def gamma_prefix_for(outdir: Path, present_ne: int) -> Path:
    return outdir / "final_gamma" / f"gravel_chb_present{int(present_ne)}"


def piece_path_for(outdir: Path, present_ne: int, daf_label: str) -> Path:
    return Path(f"{gamma_prefix_for(outdir, present_ne)}.{daf_label}")


def g_file_path_for(outdir: Path, present_ne: int) -> Path:
    return outdir / f"gravel_chb_present{int(present_ne)}.g_file.txt"


def chunk_result_file_count(workdir: Path, daf_label: str) -> int:
    pattern = f"res_{daf_label}_*.tab"
    return sum(1 for path in workdir.glob(pattern) if path.is_file() and path.stat().st_size > 0)


def chunk_is_complete(workdir: Path, daf_label: str, start_rep: int, end_rep: int) -> bool:
    expected = int(end_rep) - int(start_rep) + 1
    if expected <= 0:
        raise RuntimeError("chunk replicate range must be increasing")
    return chunk_result_file_count(workdir, daf_label) == expected


def load_task_record(task_file: Path, task_index: int) -> dict[str, object]:
    with task_file.open() as handle:
        for index, line in enumerate(handle, start=1):
            if index == task_index:
                return json.loads(line)
    raise RuntimeError(f"Task index {task_index} is out of range for {task_file}")


def run_task_file(task_file: Path, task_index: int | None) -> int:
    index = task_index if task_index is not None else int(os.environ["LSB_JOBINDEX"])
    payload = load_task_record(task_file, index)
    argv = [str(item) for item in payload["argv"]]
    return int(main(argv))


def probe_queue(queue_name: str) -> dict[str, object]:
    completed = subprocess.run(["bqueues", "-l", queue_name], check=True, capture_output=True, text=True)
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        parts = stripped.split()
        if len(parts) >= 10 and ":" in parts[2]:
            return {
                "name": queue_name,
                "status": parts[2],
                "open_active": parts[2] == "Open:Active",
                "njobs": int(parts[7]),
                "pend": int(parts[8]),
                "run": int(parts[9]),
            }
    raise RuntimeError(f"Could not parse queue state for {queue_name}")


def choose_queue(preference: str) -> tuple[str, dict[str, dict[str, object]]]:
    if preference in {"smp", "normal"}:
        return preference, {}
    queue_info = {name: probe_queue(name) for name in ["smp", "normal"]}
    smp = queue_info["smp"]
    normal = queue_info["normal"]
    if smp["open_active"]:
        if not normal["open_active"]:
            return "smp", queue_info
        smp_score = (int(smp["pend"]) + 1.0) / (int(smp["run"]) + 1.0)
        normal_score = (int(normal["pend"]) + 1.0) / (int(normal["run"]) + 1.0)
        if smp_score <= normal_score * 2.0:
            return "smp", queue_info
    if normal["open_active"]:
        return "normal", queue_info
    if smp["open_active"]:
        return "smp", queue_info
    raise RuntimeError("Neither smp nor normal is open and active")


def parse_job_id(text: str) -> str:
    match = re.search(r"Job <(\d+)>", text)
    if match is None:
        raise RuntimeError(f"Failed to parse LSF job id from: {text.strip()}")
    return match.group(1)


def submit_job(command: list[str]) -> tuple[str, str]:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    job_id = parse_job_id(completed.stdout)
    return job_id, completed.stdout.strip()


def dependency_expr(job_ids: list[str]) -> str:
    if not job_ids:
        return ""
    return "&&".join(f"done({job_id})" for job_id in job_ids)


def write_task_file(path: Path, tasks: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for task in tasks:
            handle.write(json.dumps(task, sort_keys=True) + "\n")


def write_plan_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "group_name",
        "task_kind",
        "action",
        "daf",
        "chunk_index",
        "chunk_start_rep",
        "chunk_end_rep",
        "piece_path",
        "sample_size",
    ]
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write(
                "\t".join(
                    [
                        str(row.get("group_name", "")),
                        str(row.get("task_kind", "")),
                        str(row.get("action", "")),
                        str(row.get("daf", "")),
                        str(row.get("chunk_index", "")),
                        str(row.get("chunk_start_rep", "")),
                        str(row.get("chunk_end_rep", "")),
                    str(row.get("piece_path", "")),
                        str(row.get("sample_size", "")),
                    ]
                )
                + "\n"
            )


def submit_array(
    *,
    group_name: str,
    tasks: list[dict[str, object]],
    queue_name: str,
    concurrency_cap: int,
    log_dir: Path,
    job_prefix: str,
    dependency: str,
) -> tuple[str, str] | None:
    if not tasks:
        return None
    task_dir = Path(str(tasks[0]["outdir"])).resolve() / "parallel_submit"
    task_file = task_dir / f"{group_name}_tasks.jsonl"
    write_task_file(task_file, tasks)
    job_name = f"{job_prefix}_{group_name}[1-{len(tasks)}]%{concurrency_cap}"
    command = [
        "bsub",
        "-q",
        queue_name,
        "-J",
        job_name,
        "-n",
        "1",
        "-R",
        "span[hosts=1]",
        "-cwd",
        str(REPO_ROOT),
        "-o",
        str(log_dir / f"{job_prefix}_{group_name}_%J_%I.out"),
        "-e",
        str(log_dir / f"{job_prefix}_{group_name}_%J_%I.err"),
    ]
    if dependency:
        command.extend(["-w", dependency])
    command.extend([sys.executable, str(Path(__file__).resolve()), "--run-task-file", str(task_file)])
    return submit_job(command)


def submit_finalize(
    *,
    args: argparse.Namespace,
    queue_name: str,
    log_dir: Path,
    job_prefix: str,
    dependency: str,
) -> tuple[str, str]:
    command = [
        "bsub",
        "-q",
        queue_name,
        "-J",
        f"{job_prefix}_finalize",
        "-n",
        "1",
        "-cwd",
        str(REPO_ROOT),
        "-o",
        str(log_dir / f"{job_prefix}_finalize_%J.out"),
        "-e",
        str(log_dir / f"{job_prefix}_finalize_%J.err"),
    ]
    if dependency:
        command.extend(["-w", dependency])
    command.extend(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--finalize-only",
            "--outdir",
            str(Path(args.outdir).resolve()),
            "--present-ne",
            str(int(args.present_ne)),
            "--pop-model",
            str(args.pop_model),
            "--dummy-npz",
            str(Path(args.dummy_npz).resolve()),
            "--sim-reps",
            str(int(args.sim_reps)),
            "--chunk-size",
            str(int(args.chunk_size)),
            "--daf-start",
            str(float(args.daf_start)),
            "--daf-end",
            str(float(args.daf_end)),
            "--daf-step",
            str(float(args.daf_step)),
            "--ms-make-dir",
            str(Path(args.ms_make_dir).resolve()),
            "--ms-binary",
            str(Path(args.ms_binary).resolve()),
            "--backward-script",
            str(Path(args.backward_script).resolve()),
            *(
                [
                    "--sample-size",
                    str(int(args.sample_size)),
                ]
                if args.sample_size is not None
                else []
            ),
        ],
    )
    return submit_job(command)


def run_checked(command: list[str], cwd: Path) -> None:
    subprocess.run(command, check=True, cwd=cwd)


def run_chunk_worker(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir).resolve()
    chunk_root = chunk_root_for(outdir, str(args.daf))
    workdir = chunk_workdir_for(chunk_root, int(args.chunk_index))
    run_checked(
        [
            "bash",
            str(RUN_GAMMA_CHUNK_SCRIPT),
            "--pop",
            str(args.pop_model),
            "--daf",
            str(args.daf),
            "--start-rep",
            str(int(args.chunk_start_rep)),
            "--end-rep",
            str(int(args.chunk_end_rep)),
            "--scenario-npz",
            str(Path(args.dummy_npz).resolve()),
            "--present-ne",
            str(int(args.present_ne)),
            "--workdir",
            str(workdir),
            "--ms-make-dir",
            str(Path(args.ms_make_dir).resolve()),
            "--ms-binary",
            str(Path(args.ms_binary).resolve()),
            "--backward-script",
            str(Path(args.backward_script).resolve()),
            *(
                [
                    "--sample-size",
                    str(int(args.sample_size)),
                ]
                if args.sample_size is not None
                else []
            ),
        ],
        cwd=REPO_ROOT,
    )
    print(workdir)


def run_aggregate_worker(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir).resolve()
    chunk_root = chunk_root_for(outdir, str(args.daf))
    gamma_prefix = gamma_prefix_for(outdir, int(args.present_ne))
    piece_path = piece_path_for(outdir, int(args.present_ne), str(args.daf))
    piece_path.parent.mkdir(parents=True, exist_ok=True)
    expected_ranges = chunk_ranges(int(args.sim_reps), int(args.chunk_size))
    missing = [
        f"chunk_{chunk_index:03d}"
        for chunk_index, start_rep, end_rep in expected_ranges
        if not chunk_is_complete(chunk_workdir_for(chunk_root, chunk_index), str(args.daf), start_rep, end_rep)
    ]
    if missing:
        raise RuntimeError(f"Missing chunks for DAF {args.daf}: {', '.join(missing)}")
    run_checked(
        [
            "bash",
            str(AGGREGATE_GAMMA_CHUNKS_SCRIPT),
            "--daf",
            str(args.daf),
            "--chunk-root",
            str(chunk_root),
            "--gamma-prefix",
            str(gamma_prefix),
        ],
        cwd=REPO_ROOT,
    )
    if not piece_path.exists() or piece_path.stat().st_size == 0:
        raise RuntimeError(f"Aggregate did not create a non-empty piece file: {piece_path}")
    print(piece_path)


def run_finalize(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    grid = build_daf_grid(float(args.daf_start), float(args.daf_end), float(args.daf_step))
    lines: list[str] = []
    missing: list[str] = []
    for _, daf_label in grid:
        piece_path = piece_path_for(outdir, int(args.present_ne), daf_label)
        if not piece_path.exists() or piece_path.stat().st_size == 0:
            missing.append(daf_label)
            continue
        content = piece_path.read_text().strip()
        if not content:
            missing.append(daf_label)
            continue
        lines.append(content)
    if missing:
        raise RuntimeError(f"Cannot finalize g_file; missing DAF pieces: {', '.join(missing)}")
    g_file_path = g_file_path_for(outdir, int(args.present_ne))
    g_file_path.write_text("\n".join(lines) + "\n")
    summary = {
        "pop_model": str(args.pop_model),
        "present_ne": int(args.present_ne),
        "sample_size": int(args.sample_size) if args.sample_size is not None else None,
        "sim_reps": int(args.sim_reps),
        "chunk_size": int(args.chunk_size),
        "piece_count": len(lines),
        "g_file": str(g_file_path),
    }
    summary_path = outdir / "finalize_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(g_file_path)


def orchestrate(args: argparse.Namespace) -> int:
    if args.run_task_file is not None:
        return run_task_file(Path(args.run_task_file), args.task_index)
    if args.chunk_worker:
        run_chunk_worker(args)
        return 0
    if args.aggregate_worker:
        run_aggregate_worker(args)
        return 0
    if args.finalize_only:
        run_finalize(args)
        return 0

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if not Path(args.dummy_npz).exists():
        raise RuntimeError(f"Required path does not exist: {args.dummy_npz}")
    for required in [RUN_GAMMA_CHUNK_SCRIPT, AGGREGATE_GAMMA_CHUNKS_SCRIPT, Path(args.ms_make_dir), Path(args.ms_binary), Path(args.backward_script)]:
        if not required.exists():
            raise RuntimeError(f"Required path does not exist: {required}")
    grid = build_daf_grid(float(args.daf_start), float(args.daf_end), float(args.daf_step))
    ranges = chunk_ranges(int(args.sim_reps), int(args.chunk_size))

    chunk_tasks: list[dict[str, object]] = []
    aggregate_tasks: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []

    for _, daf_label in grid:
        piece_path = piece_path_for(outdir, int(args.present_ne), daf_label)
        if piece_path.exists() and piece_path.stat().st_size > 0 and not args.no_reuse_existing:
            manifest_rows.append(
                {
                    "group_name": "pieces",
                    "task_kind": "piece",
                    "action": "skip",
                    "daf": daf_label,
                    "chunk_index": "",
                    "chunk_start_rep": "",
                    "chunk_end_rep": "",
                    "piece_path": str(piece_path),
                    "sample_size": str(int(args.sample_size)) if args.sample_size is not None else "",
                }
            )
            continue

        chunk_root = chunk_root_for(outdir, daf_label)
        for chunk_index, start_rep, end_rep in ranges:
            workdir = chunk_workdir_for(chunk_root, chunk_index)
            reusable_chunk = chunk_is_complete(workdir, daf_label, start_rep, end_rep) and not args.no_reuse_existing
            manifest_rows.append(
                {
                    "group_name": "chunks",
                    "task_kind": "chunk",
                    "action": "skip" if reusable_chunk else "submit",
                    "daf": daf_label,
                    "chunk_index": chunk_index,
                    "chunk_start_rep": start_rep,
                    "chunk_end_rep": end_rep,
                    "piece_path": str(piece_path),
                    "sample_size": str(int(args.sample_size)) if args.sample_size is not None else "",
                }
            )
            if reusable_chunk:
                continue
            chunk_tasks.append(
                {
                    "outdir": str(outdir),
                    "argv": [
                        "--chunk-worker",
                        "--outdir",
                        str(outdir),
                        "--present-ne",
                        str(int(args.present_ne)),
                        "--pop-model",
                        str(args.pop_model),
                        "--dummy-npz",
                        str(Path(args.dummy_npz).resolve()),
                        "--sim-reps",
                        str(int(args.sim_reps)),
                        "--chunk-size",
                        str(int(args.chunk_size)),
                        "--ms-make-dir",
                        str(Path(args.ms_make_dir).resolve()),
                        "--ms-binary",
                        str(Path(args.ms_binary).resolve()),
                        "--backward-script",
                        str(Path(args.backward_script).resolve()),
                        *(
                            [
                                "--sample-size",
                                str(int(args.sample_size)),
                            ]
                            if args.sample_size is not None
                            else []
                        ),
                        "--daf",
                        daf_label,
                        "--chunk-index",
                        str(chunk_index),
                        "--chunk-start-rep",
                        str(start_rep),
                        "--chunk-end-rep",
                        str(end_rep),
                    ],
                }
            )

        manifest_rows.append(
            {
                "group_name": "aggregates",
                "task_kind": "aggregate",
                "action": "submit",
                "daf": daf_label,
                "chunk_index": "",
                "chunk_start_rep": "",
                "chunk_end_rep": "",
                "piece_path": str(piece_path),
                "sample_size": str(int(args.sample_size)) if args.sample_size is not None else "",
            }
        )
        aggregate_tasks.append(
            {
                "outdir": str(outdir),
                "argv": [
                    "--aggregate-worker",
                    "--outdir",
                    str(outdir),
                    "--present-ne",
                    str(int(args.present_ne)),
                    "--pop-model",
                    str(args.pop_model),
                    "--dummy-npz",
                    str(Path(args.dummy_npz).resolve()),
                    "--sim-reps",
                    str(int(args.sim_reps)),
                    "--chunk-size",
                    str(int(args.chunk_size)),
                    "--ms-make-dir",
                    str(Path(args.ms_make_dir).resolve()),
                    "--ms-binary",
                    str(Path(args.ms_binary).resolve()),
                    "--backward-script",
                    str(Path(args.backward_script).resolve()),
                    *(
                        [
                            "--sample-size",
                            str(int(args.sample_size)),
                        ]
                        if args.sample_size is not None
                        else []
                    ),
                    "--daf",
                    daf_label,
                ],
            }
        )

    submit_dir = outdir / "parallel_submit"
    submit_dir.mkdir(parents=True, exist_ok=True)
    write_plan_manifest(submit_dir / "task_plan.tsv", manifest_rows)
    write_task_file(submit_dir / "chunks_tasks.jsonl", chunk_tasks)
    write_task_file(submit_dir / "aggregates_tasks.jsonl", aggregate_tasks)

    queue_name, queue_info = choose_queue(args.queue)
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)
    plan_payload = {
        "outdir": str(outdir),
        "queue_name": queue_name,
        "queue_info": queue_info,
        "present_ne": int(args.present_ne),
        "sample_size": int(args.sample_size) if args.sample_size is not None else None,
        "pop_model": str(args.pop_model),
        "chunk_task_count": len(chunk_tasks),
        "aggregate_task_count": len(aggregate_tasks),
        "expected_g_file": str(g_file_path_for(outdir, int(args.present_ne))),
    }
    with (submit_dir / "submission_plan.json").open("w") as handle:
        json.dump(plan_payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if args.dry_run:
        print(json.dumps(plan_payload, indent=2, sort_keys=True))
        return 0

    submitted_jobs: dict[str, dict[str, str]] = {}
    chunk_submission = submit_array(
        group_name="chunks",
        tasks=chunk_tasks,
        queue_name=queue_name,
        concurrency_cap=int(args.chunk_array_cap),
        log_dir=log_dir,
        job_prefix=args.job_prefix,
        dependency="",
    )
    if chunk_submission is not None:
        submitted_jobs["chunks"] = {"job_id": chunk_submission[0], "raw": chunk_submission[1]}

    aggregate_dependency = dependency_expr([submitted_jobs["chunks"]["job_id"]]) if "chunks" in submitted_jobs else ""
    aggregate_submission = submit_array(
        group_name="aggregates",
        tasks=aggregate_tasks,
        queue_name=queue_name,
        concurrency_cap=int(args.aggregate_array_cap),
        log_dir=log_dir,
        job_prefix=args.job_prefix,
        dependency=aggregate_dependency,
    )
    if aggregate_submission is not None:
        submitted_jobs["aggregates"] = {"job_id": aggregate_submission[0], "raw": aggregate_submission[1]}

    finalize_dependency = dependency_expr([payload["job_id"] for payload in submitted_jobs.values()])
    finalize_submission = submit_finalize(
        args=args,
        queue_name=queue_name,
        log_dir=log_dir,
        job_prefix=args.job_prefix,
        dependency=finalize_dependency,
    )
    submitted_jobs["finalize"] = {"job_id": finalize_submission[0], "raw": finalize_submission[1]}

    with (submit_dir / "submitted_jobs.json").open("w") as handle:
        json.dump({"queue_name": queue_name, "jobs": submitted_jobs}, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"outdir": str(outdir), "queue_name": queue_name, "jobs": submitted_jobs}, indent=2, sort_keys=True))
    return 0


def validate_args(args: argparse.Namespace) -> None:
    mode_count = sum(bool(flag) for flag in [args.run_task_file is not None, args.chunk_worker, args.aggregate_worker, args.finalize_only])
    if mode_count > 1:
        raise RuntimeError("Use only one of --run-task-file, --chunk-worker, --aggregate-worker, or --finalize-only")
    if int(args.sim_reps) <= 0 or int(args.chunk_size) <= 0:
        raise RuntimeError("--sim-reps and --chunk-size must be positive")
    if args.chunk_worker and any(value is None for value in [args.daf, args.chunk_index, args.chunk_start_rep, args.chunk_end_rep]):
        raise RuntimeError("--chunk-worker requires --daf, --chunk-index, --chunk-start-rep, and --chunk-end-rep")
    if args.aggregate_worker and args.daf is None:
        raise RuntimeError("--aggregate-worker requires --daf")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    validate_args(args)
    return orchestrate(args)


if __name__ == "__main__":
    raise SystemExit(main())
