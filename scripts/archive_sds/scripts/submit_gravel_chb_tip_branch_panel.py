#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
SDS_PYTHON = Path("/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python")
RUN_PIECE = REPO_ROOT / "scripts" / "run_single_snp_gamma_piece.sh"
SUMMARIZE = REPO_ROOT / "scripts" / "summarize_tip_branch_panel.py"
DEFAULT_DUMMY_NPZ = REPO_ROOT / "tmp" / "region_ne0_positive_20260426" / "scaled_ne0_100000" / "NCN_scaled_ne0_100000.npz"
DEFAULT_OUTDIR = REPO_ROOT / "tmp" / "gravel_chb_tip_branch_ne_panel_20260525"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Submit a small Gravel_CHB Ne(0) tip-branch diagnostic panel.")
    parser.add_argument("--run-task-file", default=None, help="Internal worker mode: JSONL task file.")
    parser.add_argument("--task-index", type=int, default=None, help="1-based task index for worker mode.")
    parser.add_argument("--queue", choices=["auto", "smp", "normal"], default="auto")
    parser.add_argument("--array-cap", type=int, default=6)
    parser.add_argument("--present-ne-grid", default="60000,80000,100000,150000,200000,270000")
    parser.add_argument("--daf-grid", default="0.05,0.07,0.10")
    parser.add_argument("--sim-reps", type=int, default=100)
    parser.add_argument("--pop-model", default="Gravel_CHB")
    parser.add_argument("--dummy-npz", default=str(DEFAULT_DUMMY_NPZ))
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    parser.add_argument("--log-dir", default=str(REPO_ROOT / "logs"))
    parser.add_argument("--job-prefix", default="gravel_tip_panel")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def parse_int_grid(text: str) -> list[int]:
    values = []
    for item in text.split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    if not values:
        raise RuntimeError("Empty present-ne grid")
    return values


def parse_float_grid(text: str) -> list[float]:
    values = []
    for item in text.split(","):
        item = item.strip()
        if item:
            values.append(float(item))
    if not values:
        raise RuntimeError("Empty daf grid")
    return values


def load_task_record(task_file: Path, task_index: int) -> dict[str, object]:
    with task_file.open() as handle:
        for index, line in enumerate(handle, start=1):
            if index == task_index:
                return json.loads(line)
    raise RuntimeError(f"Task index {task_index} out of range for {task_file}")


def run_task_file(task_file: Path, task_index: int | None) -> int:
    index = task_index if task_index is not None else int(os.environ["LSB_JOBINDEX"])
    payload = load_task_record(task_file, index)
    argv = [str(item) for item in payload["argv"]]
    completed = subprocess.run(argv, check=False)
    return completed.returncode


def probe_queue(queue_name: str) -> dict[str, object]:
    completed = subprocess.run(["bqueues", "-l", queue_name], check=True, capture_output=True, text=True)
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        parts = stripped.split()
        if len(parts) >= 10 and ":" in parts[2]:
            return {
                "name": queue_name,
                "open_active": parts[2] == "Open:Active",
                "pend": int(parts[8]),
                "run": int(parts[9]),
            }
    raise RuntimeError(f"Could not parse queue state for {queue_name}")


def choose_queue(preference: str) -> str:
    if preference in {"smp", "normal"}:
        return preference
    try:
        smp = probe_queue("smp")
        normal = probe_queue("normal")
    except subprocess.CalledProcessError:
        return "smp"
    if smp["open_active"] and smp["pend"] == 0:
        return "smp"
    if smp["open_active"] and not normal["open_active"]:
        return "smp"
    if normal["open_active"]:
        return "normal"
    if smp["open_active"]:
        return "smp"
    raise RuntimeError("Neither smp nor normal is open and active")


def parse_job_id(text: str) -> str:
    match = re.search(r"Job <(\d+)>", text)
    if match is None:
        raise RuntimeError(f"Failed to parse job id from {text!r}")
    return match.group(1)


def submit_job(command: list[str]) -> tuple[str, str]:
    completed = subprocess.run(
        ["/bin/bash", "-lc", shlex.join(command)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"LSF submission failed (code={completed.returncode})\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return parse_job_id(completed.stdout), completed.stdout.strip()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.run_task_file is not None:
        return run_task_file(Path(args.run_task_file), args.task_index)

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    ne_grid = parse_int_grid(args.present_ne_grid)
    daf_grid = parse_float_grid(args.daf_grid)
    queue = choose_queue(args.queue)

    tasks: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    for present_ne in ne_grid:
        for daf in daf_grid:
            scenario_label = f"gravel_chb_ne{present_ne}_daf{str(daf).replace('.', 'p')}"
            workdir = outdir / scenario_label / "workdir"
            gamma_prefix = outdir / scenario_label / "gamma" / scenario_label
            gamma_output = Path(f"{gamma_prefix}.{daf:.2f}")
            argv_task = [
                "bash",
                str(RUN_PIECE),
                "--pop",
                args.pop_model,
                "--daf",
                f"{daf:.2f}",
                "--scenario-npz",
                str(Path(args.dummy_npz).resolve()),
                "--present-ne",
                str(present_ne),
                "--gamma-prefix",
                str(gamma_prefix),
                "--workdir",
                str(workdir),
                "--sim-reps",
                str(args.sim_reps),
                "--dump-tip-branches",
            ]
            tasks.append({"argv": argv_task})
            manifest_rows.append(
                {
                    "scenario_label": scenario_label,
                    "present_ne": present_ne,
                    "daf": f"{daf:.2f}",
                    "sim_reps": args.sim_reps,
                    "workdir": str(workdir),
                    "gamma_prefix": str(gamma_prefix),
                    "gamma_output": str(gamma_output),
                }
            )

    task_file = outdir / "parallel_submit" / "tasks.jsonl"
    task_file.parent.mkdir(parents=True, exist_ok=True)
    with task_file.open("w") as handle:
        for task in tasks:
            handle.write(json.dumps(task, sort_keys=True) + "\n")

    manifest = outdir / "tip_branch_panel_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    summary_prefix = outdir / "tip_branch_panel"

    queue_info = {"selected_queue": queue, "n_tasks": len(tasks)}
    (outdir / "submission_context.json").write_text(json.dumps(queue_info, indent=2) + "\n")

    if args.dry_run:
        print(json.dumps({"manifest": str(manifest), "task_file": str(task_file), **queue_info}, indent=2))
        return 0

    worker_cmd = (
        f'"{SDS_PYTHON}" "{SCRIPT_PATH}" --run-task-file "{task_file}" --task-index "$LSB_JOBINDEX"'
    )
    array_cmd = [
        "bsub",
        "-q",
        queue,
        "-n",
        "2",
        "-R",
        "span[hosts=1]",
        "-J",
        f"{args.job_prefix}[1-{len(tasks)}]%{args.array_cap}",
        "-o",
        str(log_dir / f"{args.job_prefix}_%I.out"),
        "-e",
        str(log_dir / f"{args.job_prefix}_%I.err"),
        "/bin/bash",
        "-lc",
        worker_cmd,
    ]
    array_job, array_text = submit_job(array_cmd)

    summary_cmd = [
        "bsub",
        "-q",
        queue,
        "-w",
        f"done({array_job})",
        "-n",
        "1",
        "-R",
        "span[hosts=1]",
        "-J",
        f"{args.job_prefix}_summary",
        "-o",
        str(log_dir / f"{args.job_prefix}_summary.out"),
        "-e",
        str(log_dir / f"{args.job_prefix}_summary.err"),
        "/bin/bash",
        "-lc",
        f'"{SDS_PYTHON}" "{SUMMARIZE}" --manifest "{manifest}" --output-prefix "{summary_prefix}"',
    ]
    summary_job, summary_text = submit_job(summary_cmd)

    print(
        json.dumps(
            {
                "manifest": str(manifest),
                "task_file": str(task_file),
                "array_job": array_job,
                "array_submit": array_text,
                "summary_job": summary_job,
                "summary_submit": summary_text,
                "queue": queue,
                "n_tasks": len(tasks),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
