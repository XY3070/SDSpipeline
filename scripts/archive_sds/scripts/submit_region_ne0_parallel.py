#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

SDS_PYTHON = Path("/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python")
if Path(sys.executable).resolve() != SDS_PYTHON.resolve():
    os.execv(str(SDS_PYTHON), [str(SDS_PYTHON), str(Path(__file__).resolve()), *sys.argv[1:]])

from benchmark_adh_ne0 import (
    REPO_ROOT,
    chunk_is_complete,
    chunk_ranges,
    chunk_root_for,
    chunk_workdir_for,
    discover_piece_fragment,
    frequency_key,
    main as benchmark_main,
    parse_args as parse_benchmark_args,
    prepare_benchmark,
    purpose_sim_reps,
    scenario_requires_smoke_prepass,
)


DEFAULT_RESCUE_OUTDIR = REPO_ROOT / "tmp" / "region_ne0_positive_20260426"
BENCHMARK_SCRIPT = Path(__file__).with_name("benchmark_adh_ne0.py").resolve()


def parse_args(argv: list[str] | None = None) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(
        description="Submit the region Ne(0) positive-control rescue as parallel LSF arrays plus a dependent finalize job."
    )
    parser.add_argument("--run-task-file", default=None, help="Internal worker mode: JSONL task file written by the submitter.")
    parser.add_argument("--task-index", type=int, default=None, help="1-based task index for --run-task-file.")
    parser.add_argument("--dry-run", action="store_true", help="Prepare the outdir and task plan without submitting jobs.")
    parser.add_argument("--queue", choices=["auto", "smp", "normal"], default="auto", help="Queue selection strategy.")
    parser.add_argument("--full-array-cap", type=int, default=32, help="Maximum concurrent full-gamma array slots.")
    parser.add_argument("--smoke-array-cap", type=int, default=16, help="Maximum concurrent smoke-gamma array slots.")
    parser.add_argument("--chunk-array-cap", type=int, default=64, help="Maximum concurrent chunk-gamma array slots.")
    parser.add_argument("--job-prefix", default="region_ne0_parallel", help="LSF job-name prefix.")
    parser.add_argument("--log-dir", default=str(REPO_ROOT / "logs"), help="Directory for LSF stdout/stderr logs.")
    args, benchmark_argv = parser.parse_known_args(argv)
    return args, benchmark_argv


def ensure_default_resume_outdir(argv: list[str]) -> list[str]:
    has_outdir = any(
        token == "--resume-outdir"
        or token.startswith("--resume-outdir=")
        or token == "--outdir"
        or token.startswith("--outdir=")
        for token in argv
    )
    if has_outdir:
        return list(argv)
    return list(argv) + ["--resume-outdir", str(DEFAULT_RESCUE_OUTDIR)]


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
    return int(benchmark_main(argv))


def probe_queue(queue_name: str) -> dict[str, object]:
    completed = subprocess.run(
        ["bqueues", "-l", queue_name],
        check=True,
        capture_output=True,
        text=True,
    )
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
        "scenario_label",
        "purpose",
        "frequency",
        "chunk_index",
        "chunk_start_rep",
        "chunk_end_rep",
        "piece_status",
        "piece_path",
        "source_label",
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
                        str(row.get("scenario_label", "")),
                        str(row.get("purpose", "")),
                        str(row.get("frequency", "")),
                        str(row.get("chunk_index", "")),
                        str(row.get("chunk_start_rep", "")),
                        str(row.get("chunk_end_rep", "")),
                        str(row.get("piece_status", "")),
                        str(row.get("piece_path", "")),
                        str(row.get("source_label", "")),
                    ]
                )
                + "\n"
            )


def build_piece_worker_argv(common_benchmark_argv: list[str], scenario_label: str, purpose: str, frequency: float) -> list[str]:
    return list(common_benchmark_argv) + [
        "--piece-worker",
        "--scenario-label",
        scenario_label,
        "--gamma-purpose",
        purpose,
        "--frequency",
        frequency_key(frequency),
    ]


def build_chunk_worker_argv(
    common_benchmark_argv: list[str],
    scenario_label: str,
    purpose: str,
    frequency: float,
    chunk_index: int,
    start_rep: int,
    end_rep: int,
) -> list[str]:
    return list(common_benchmark_argv) + [
        "--chunk-worker",
        "--scenario-label",
        scenario_label,
        "--gamma-purpose",
        purpose,
        "--frequency",
        frequency_key(frequency),
        "--chunk-index",
        str(int(chunk_index)),
        "--chunk-start-rep",
        str(int(start_rep)),
        "--chunk-end-rep",
        str(int(end_rep)),
    ]


def build_aggregate_worker_argv(common_benchmark_argv: list[str], scenario_label: str, purpose: str, frequency: float) -> list[str]:
    return list(common_benchmark_argv) + [
        "--aggregate-worker",
        "--scenario-label",
        scenario_label,
        "--gamma-purpose",
        purpose,
        "--frequency",
        frequency_key(frequency),
    ]


def build_finalize_argv(common_benchmark_argv: list[str]) -> list[str]:
    return list(common_benchmark_argv) + ["--finalize-only"]


def collect_piece_tasks(
    prepared,
    benchmark_args,
    common_benchmark_argv: list[str],
) -> tuple[dict[str, list[dict[str, object]]], list[dict[str, object]]]:
    chunked_mode = benchmark_args.gamma_generation_mode == "chunked"
    tasks_by_group: dict[str, list[dict[str, object]]] = (
        {
            "smoke_chunks": [],
            "smoke_aggregates": [],
            "full_after_smoke_chunks": [],
            "full_after_smoke_aggregates": [],
            "full_direct_chunks": [],
            "full_direct_aggregates": [],
        }
        if chunked_mode
        else {
            "smoke": [],
            "full_after_smoke": [],
            "full_direct": [],
        }
    )
    manifest_rows: list[dict[str, object]] = []
    reuse_existing = not benchmark_args.no_reuse_existing_gamma

    for scenario in prepared.scenarios:
        if scenario.smoke_only:
            scheduled_groups = [("smoke", "smoke")]
        elif scenario_requires_smoke_prepass(scenario, benchmark_args):
            scheduled_groups = [("smoke", "smoke"), ("full_after_smoke", "full")]
        else:
            scheduled_groups = [("full_direct", "full")]

        for group_name, purpose in scheduled_groups:
            for frequency in prepared.unique_frequencies:
                fragment = discover_piece_fragment(
                    prepared=prepared,
                    scenario=scenario,
                    purpose=purpose,
                    frequency=frequency,
                    args=benchmark_args,
                    allow_legacy_backfill=True,
                )
                should_skip = reuse_existing and fragment.gamma_value is not None
                action = "skip" if should_skip else ("submit_forced" if fragment.gamma_value is not None else "submit")
                manifest_rows.append(
                    {
                        "group_name": group_name,
                        "task_kind": "piece",
                        "action": action,
                        "scenario_label": scenario.label,
                        "purpose": purpose,
                        "frequency": frequency_key(frequency),
                        "chunk_index": "",
                        "chunk_start_rep": "",
                        "chunk_end_rep": "",
                        "piece_status": fragment.status,
                        "piece_path": str(fragment.piece_path),
                        "source_label": fragment.source_label or "",
                    }
                )
                if should_skip:
                    continue
                if not chunked_mode:
                    tasks_by_group[group_name].append(
                        {
                            "scenario_label": scenario.label,
                            "purpose": purpose,
                            "frequency": float(frequency),
                            "piece_path": str(fragment.piece_path),
                            "argv": build_piece_worker_argv(common_benchmark_argv, scenario.label, purpose, frequency),
                            "outdir": str(prepared.outdir),
                        }
                    )
                    continue

                chunk_group_name = f"{group_name}_chunks"
                aggregate_group_name = f"{group_name}_aggregates"
                sim_reps = purpose_sim_reps(purpose, benchmark_args)
                ranges = chunk_ranges(sim_reps, int(benchmark_args.gamma_chunk_size))
                chunk_root = chunk_root_for(scenario.scenario_dir, purpose, frequency)
                for chunk_index, start_rep, end_rep in ranges:
                    workdir = chunk_workdir_for(chunk_root, chunk_index)
                    chunk_complete = chunk_is_complete(workdir, frequency, start_rep, end_rep)
                    manifest_rows.append(
                        {
                            "group_name": chunk_group_name,
                            "task_kind": "chunk",
                            "action": "skip" if chunk_complete else "submit",
                            "scenario_label": scenario.label,
                            "purpose": purpose,
                            "frequency": frequency_key(frequency),
                            "chunk_index": chunk_index,
                            "chunk_start_rep": start_rep,
                            "chunk_end_rep": end_rep,
                            "piece_status": "reused_chunk" if chunk_complete else "missing_chunk",
                            "piece_path": str(fragment.piece_path),
                            "source_label": fragment.source_label or "",
                        }
                    )
                    if chunk_complete:
                        continue
                    tasks_by_group[chunk_group_name].append(
                        {
                            "scenario_label": scenario.label,
                            "purpose": purpose,
                            "frequency": float(frequency),
                            "piece_path": str(fragment.piece_path),
                            "argv": build_chunk_worker_argv(
                                common_benchmark_argv,
                                scenario.label,
                                purpose,
                                frequency,
                                chunk_index,
                                start_rep,
                                end_rep,
                            ),
                            "outdir": str(prepared.outdir),
                        }
                    )

                manifest_rows.append(
                    {
                        "group_name": aggregate_group_name,
                        "task_kind": "aggregate",
                        "action": "submit",
                        "scenario_label": scenario.label,
                        "purpose": purpose,
                        "frequency": frequency_key(frequency),
                        "chunk_index": "",
                        "chunk_start_rep": "",
                        "chunk_end_rep": "",
                        "piece_status": fragment.status,
                        "piece_path": str(fragment.piece_path),
                        "source_label": fragment.source_label or "",
                    }
                )
                tasks_by_group[aggregate_group_name].append(
                    {
                        "scenario_label": scenario.label,
                        "purpose": purpose,
                        "frequency": float(frequency),
                        "piece_path": str(fragment.piece_path),
                        "argv": build_aggregate_worker_argv(common_benchmark_argv, scenario.label, purpose, frequency),
                        "outdir": str(prepared.outdir),
                    }
                )

    return tasks_by_group, manifest_rows


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
    finalize_argv: list[str],
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
    command.extend([sys.executable, str(BENCHMARK_SCRIPT), *finalize_argv])
    return submit_job(command)


def orchestrate(argv: list[str] | None = None) -> int:
    args, benchmark_argv = parse_args(argv)
    if args.run_task_file is not None:
        return run_task_file(Path(args.run_task_file), args.task_index)

    common_benchmark_argv = ensure_default_resume_outdir(benchmark_argv)
    benchmark_args = parse_benchmark_args(common_benchmark_argv)
    if (
        benchmark_args.prepare_only
        or benchmark_args.piece_worker
        or benchmark_args.chunk_worker
        or benchmark_args.aggregate_worker
        or benchmark_args.finalize_only
    ):
        raise RuntimeError("Submitter input should not include benchmark mode flags")
    prepared = prepare_benchmark(benchmark_args)

    tasks_by_group, manifest_rows = collect_piece_tasks(prepared, benchmark_args, common_benchmark_argv)
    submit_dir = prepared.outdir / "parallel_submit"
    submit_dir.mkdir(parents=True, exist_ok=True)
    write_plan_manifest(submit_dir / "task_plan.tsv", manifest_rows)

    queue_name, queue_info = choose_queue(args.queue)
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)
    finalize_argv = build_finalize_argv(common_benchmark_argv)

    plan_payload = {
        "outdir": str(prepared.outdir),
        "queue_name": queue_name,
        "queue_info": queue_info,
        "task_counts": {group: len(tasks) for group, tasks in tasks_by_group.items()},
        "finalize_argv": finalize_argv,
    }
    with (submit_dir / "submission_plan.json").open("w") as handle:
        json.dump(plan_payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if args.dry_run:
        print(json.dumps(plan_payload, indent=2, sort_keys=True))
        return 0

    submitted_jobs: dict[str, dict[str, str]] = {}
    if benchmark_args.gamma_generation_mode == "chunked":
        smoke_chunks_submission = submit_array(
            group_name="smoke_chunks",
            tasks=tasks_by_group["smoke_chunks"],
            queue_name=queue_name,
            concurrency_cap=args.chunk_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency="",
        )
        if smoke_chunks_submission is not None:
            submitted_jobs["smoke_chunks"] = {"job_id": smoke_chunks_submission[0], "raw": smoke_chunks_submission[1]}

        smoke_aggregate_dependency = dependency_expr([submitted_jobs["smoke_chunks"]["job_id"]]) if "smoke_chunks" in submitted_jobs else ""
        smoke_aggregates_submission = submit_array(
            group_name="smoke_aggregates",
            tasks=tasks_by_group["smoke_aggregates"],
            queue_name=queue_name,
            concurrency_cap=args.full_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency=smoke_aggregate_dependency,
        )
        if smoke_aggregates_submission is not None:
            submitted_jobs["smoke_aggregates"] = {"job_id": smoke_aggregates_submission[0], "raw": smoke_aggregates_submission[1]}

        full_after_smoke_dependency = dependency_expr([submitted_jobs["smoke_aggregates"]["job_id"]]) if "smoke_aggregates" in submitted_jobs else ""
        full_after_smoke_chunks_submission = submit_array(
            group_name="full_after_smoke_chunks",
            tasks=tasks_by_group["full_after_smoke_chunks"],
            queue_name=queue_name,
            concurrency_cap=args.chunk_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency=full_after_smoke_dependency,
        )
        if full_after_smoke_chunks_submission is not None:
            submitted_jobs["full_after_smoke_chunks"] = {
                "job_id": full_after_smoke_chunks_submission[0],
                "raw": full_after_smoke_chunks_submission[1],
            }

        full_after_smoke_aggregate_dependency = (
            dependency_expr([submitted_jobs["full_after_smoke_chunks"]["job_id"]])
            if "full_after_smoke_chunks" in submitted_jobs
            else ""
        )
        full_after_smoke_aggregates_submission = submit_array(
            group_name="full_after_smoke_aggregates",
            tasks=tasks_by_group["full_after_smoke_aggregates"],
            queue_name=queue_name,
            concurrency_cap=args.full_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency=full_after_smoke_aggregate_dependency,
        )
        if full_after_smoke_aggregates_submission is not None:
            submitted_jobs["full_after_smoke_aggregates"] = {
                "job_id": full_after_smoke_aggregates_submission[0],
                "raw": full_after_smoke_aggregates_submission[1],
            }

        full_direct_chunks_submission = submit_array(
            group_name="full_direct_chunks",
            tasks=tasks_by_group["full_direct_chunks"],
            queue_name=queue_name,
            concurrency_cap=args.chunk_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency="",
        )
        if full_direct_chunks_submission is not None:
            submitted_jobs["full_direct_chunks"] = {"job_id": full_direct_chunks_submission[0], "raw": full_direct_chunks_submission[1]}

        full_direct_aggregate_dependency = (
            dependency_expr([submitted_jobs["full_direct_chunks"]["job_id"]])
            if "full_direct_chunks" in submitted_jobs
            else ""
        )
        full_direct_aggregates_submission = submit_array(
            group_name="full_direct_aggregates",
            tasks=tasks_by_group["full_direct_aggregates"],
            queue_name=queue_name,
            concurrency_cap=args.full_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency=full_direct_aggregate_dependency,
        )
        if full_direct_aggregates_submission is not None:
            submitted_jobs["full_direct_aggregates"] = {
                "job_id": full_direct_aggregates_submission[0],
                "raw": full_direct_aggregates_submission[1],
            }

        finalize_dependency = dependency_expr(
            [
                payload["job_id"]
                for group_name, payload in submitted_jobs.items()
                if group_name.endswith("_aggregates")
            ]
        )
    else:
        smoke_submission = submit_array(
            group_name="smoke",
            tasks=tasks_by_group["smoke"],
            queue_name=queue_name,
            concurrency_cap=args.smoke_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency="",
        )
        if smoke_submission is not None:
            submitted_jobs["smoke"] = {"job_id": smoke_submission[0], "raw": smoke_submission[1]}

        smoke_dependency = dependency_expr([submitted_jobs["smoke"]["job_id"]]) if "smoke" in submitted_jobs else ""
        full_after_smoke_submission = submit_array(
            group_name="full_after_smoke",
            tasks=tasks_by_group["full_after_smoke"],
            queue_name=queue_name,
            concurrency_cap=args.full_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency=smoke_dependency,
        )
        if full_after_smoke_submission is not None:
            submitted_jobs["full_after_smoke"] = {
                "job_id": full_after_smoke_submission[0],
                "raw": full_after_smoke_submission[1],
            }

        full_direct_submission = submit_array(
            group_name="full_direct",
            tasks=tasks_by_group["full_direct"],
            queue_name=queue_name,
            concurrency_cap=args.full_array_cap,
            log_dir=log_dir,
            job_prefix=args.job_prefix,
            dependency="",
        )
        if full_direct_submission is not None:
            submitted_jobs["full_direct"] = {"job_id": full_direct_submission[0], "raw": full_direct_submission[1]}

        finalize_dependency = dependency_expr([payload["job_id"] for payload in submitted_jobs.values()])

    finalize_submission = submit_finalize(
        finalize_argv=finalize_argv,
        queue_name=queue_name,
        log_dir=log_dir,
        job_prefix=args.job_prefix,
        dependency=finalize_dependency,
    )
    submitted_jobs["finalize"] = {"job_id": finalize_submission[0], "raw": finalize_submission[1]}

    with (submit_dir / "submitted_jobs.json").open("w") as handle:
        json.dump({"queue_name": queue_name, "jobs": submitted_jobs}, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps({"outdir": str(prepared.outdir), "queue_name": queue_name, "jobs": submitted_jobs}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(orchestrate())
