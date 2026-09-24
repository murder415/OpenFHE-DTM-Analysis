"""Run four sequential Windows benchmark processes and preserve their evidence.

Cold/warm means first/subsequent application execution. The OS file cache is not
flushed. Packing reuses saved keys AND saved ciphertext packs on runs 2--4.
Analysis uses a fresh OpenFHE context/key pair in every process.
"""
from __future__ import annotations

import argparse
import ctypes
from ctypes import wintypes
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parent


def now() -> str:
    return dt.datetime.now().astimezone().isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, object]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def write_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def cache_inventory(root: Path) -> dict[str, object]:
    # Metadata only. Do not read/hash multi-GB key files between timed runs.
    files = []
    if root.exists():
        for path in sorted(root.rglob("*")):
            if path.is_file():
                stat = path.stat()
                files.append({"name": str(path.relative_to(root)), "bytes": stat.st_size,
                              "modified_ns": stat.st_mtime_ns})
    return {"file_count": len(files), "bytes": sum(f["bytes"] for f in files), "files": files}


class ProcessMemoryCountersEx(ctypes.Structure):
    _fields_ = [
        ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t),
        ("PrivateUsage", ctypes.c_size_t),
    ]


class WindowsMemoryReader:
    def __init__(self) -> None:
        self.psapi = ctypes.WinDLL("psapi", use_last_error=True)
        self.query = self.psapi.GetProcessMemoryInfo
        self.query.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessMemoryCountersEx), wintypes.DWORD]
        self.query.restype = wintypes.BOOL

    def read(self, process: subprocess.Popen) -> tuple[dict[str, int] | None, int | None]:
        counters = ProcessMemoryCountersEx()
        counters.cb = ctypes.sizeof(counters)
        # CPython retains this native process handle through wait(), enabling a
        # final OS high-water-counter query after termination when supported.
        handle = wintypes.HANDLE(int(process._handle))
        if not self.query(handle, ctypes.byref(counters), counters.cb):
            return None, ctypes.get_last_error()
        return {"peak_working_set_bytes": int(counters.PeakWorkingSetSize),
                "working_set_bytes": int(counters.WorkingSetSize),
                "private_commit_bytes": int(counters.PrivateUsage),
                "peak_pagefile_bytes": int(counters.PeakPagefileUsage)}, None


def measure(command: list[str], output: Path, run: int, timeout: float,
            poll_seconds: float, environment: dict[str, str]) -> dict[str, object]:
    stdout_path = output / f"run_{run}.stdout.txt"
    stderr_path = output / f"run_{run}.stderr.txt"
    memory = WindowsMemoryReader()
    observations = 0
    peak_working_set = peak_private_commit_observed = peak_pagefile = 0
    query_errors: dict[str, int] = {}
    final_query = None
    process = None
    timed_out = cancelled = False

    def observe() -> dict[str, int] | None:
        nonlocal observations, peak_working_set, peak_private_commit_observed, peak_pagefile
        values, error = memory.read(process)
        if values is None:
            key = str(error)
            query_errors[key] = query_errors.get(key, 0) + 1
            return None
        observations += 1
        peak_working_set = max(peak_working_set, values["peak_working_set_bytes"])
        peak_private_commit_observed = max(peak_private_commit_observed, values["private_commit_bytes"])
        peak_pagefile = max(peak_pagefile, values["peak_pagefile_bytes"])
        return values

    started_at = now()
    start = time.perf_counter()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, cwd=str(Path(command[0]).parent), env=environment,
                                   stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
                                   creationflags=subprocess.CREATE_NO_WINDOW)
        try:
            while True:
                observe()
                elapsed = time.perf_counter() - start
                if elapsed >= timeout:
                    timed_out = True
                    process.kill()
                    process.wait()
                    break
                try:
                    process.wait(timeout=min(poll_seconds, timeout - elapsed))
                    break
                except subprocess.TimeoutExpired:
                    pass
        except KeyboardInterrupt:
            cancelled = True
            process.kill()
            process.wait()
        except BaseException:
            process.kill()
            process.wait()
            raise
        finally:
            wall_seconds = time.perf_counter() - start
            final_query = observe()
    return {
        "run": run, "label": "first_execution" if run == 1 else "subsequent_execution",
        "started_at": started_at, "completed_at": now(), "pid": process.pid,
        "command": command, "working_directory": str(Path(command[0]).parent),
        "exit_code": process.returncode, "timed_out": timed_out, "cancelled": cancelled,
        "wall_seconds": wall_seconds,
        "peak_working_set_bytes": peak_working_set,
        "peak_working_set_MB_decimal": peak_working_set / 1_000_000,
        "peak_working_set_MiB": peak_working_set / (1024 * 1024),
        "peak_private_commit_observed_bytes": peak_private_commit_observed,
        "peak_pagefile_bytes": peak_pagefile,
        "memory_observations": observations, "memory_query_errors": query_errors,
        "final_memory_query": final_query,
        "final_peak_counter_available": final_query is not None and final_query["peak_working_set_bytes"] > 0,
        "stdout": str(stdout_path), "stderr": str(stderr_path),
    }


def parse_program_output(mode: str, stdout: str) -> dict[str, object]:
    if mode == "analysis":
        records = [json.loads(line[len("RESULT_JSON "):]) for line in stdout.splitlines()
                   if line.startswith("RESULT_JSON ")]
        if len(records) != 1:
            raise RuntimeError("Expected one RESULT_JSON record from final_analysis_demo")
        return records[0]
    parsed = {}
    for name in ("pack", "pack16"):
        match = re.search(rf"^{name}: total=(\d+) encrypted=(\d+) skipped=(\d+) (?:mask/list|split-stitches)=(\d+)$",
                          stdout, flags=re.MULTILINE)
        if not match:
            raise RuntimeError(f"Expected {name} summary from packed_store_validation")
        parsed[name] = dict(zip(("total", "encrypted", "skipped", "stitches"), map(int, match.groups())))
    if "real/imag split verified" not in stdout:
        raise RuntimeError("Packing program did not report successful split verification")
    return parsed


def system_metadata() -> dict[str, object]:
    data = {"platform": platform.platform(), "processor": platform.processor(),
            "machine": platform.machine(), "logical_cpus": os.cpu_count(),
            "python_version": sys.version, "python_executable": sys.executable}
    try:
        import psutil
        data.update({"physical_cpus": psutil.cpu_count(logical=False),
                     "physical_memory_bytes": psutil.virtual_memory().total,
                     "psutil_version": psutil.__version__})
    except ImportError:
        data["psutil_version"] = None
    return data


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("analysis", "packing"), required=True)
    parser.add_argument("--tdb", type=Path, default=ROOT / "data/taiwan.tdb")
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--output", type=Path, help="New or empty output directory; defaults to a unique log directory")
    parser.add_argument("--runs", type=int, default=4)
    parser.add_argument("--height", type=float, default=14.76)
    parser.add_argument("--angle", type=float, default=45.0)
    parser.add_argument("--grid-interval", type=float)
    parser.add_argument("--threads", type=int, help="Set OMP_NUM_THREADS; default preserves inherited setting")
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--poll-ms", type=float, default=10.0)
    args = parser.parse_args()
    if os.name != "nt":
        parser.error("This benchmark uses Windows process memory counters")
    if not 2 <= args.runs <= 100 or args.timeout <= 0 or args.poll_ms <= 0:
        parser.error("Require 2..100 runs and positive timeout/poll interval")
    if args.threads is not None and args.threads <= 0:
        parser.error("--threads must be positive")
    executable = (args.exe or ROOT / ("build/final_analysis_demo.exe" if args.mode == "analysis"
                                      else "build/upstream/packed_store_validation.exe")).resolve()
    dataset = args.tdb.resolve()
    if not executable.is_file() or not dataset.is_file():
        parser.error(f"Missing executable or dataset: {executable}, {dataset}")
    benchmark_root = ROOT / "logs/benchmarks"
    benchmark_root.mkdir(parents=True, exist_ok=True)
    output = (args.output or benchmark_root / (
        args.mode + "_" + dt.datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8])).resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f"Output directory must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    # A common lock prevents this script's analysis and packing modes from
    # overlapping. Unrelated builds/programs must also be stopped by the caller.
    lock_path = benchmark_root / ".benchmark_process.lock"
    try:
        lock_descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        parser.error(f"Another benchmark may be running; inspect {lock_path}")
    with os.fdopen(lock_descriptor, "w", encoding="utf-8") as lock_file:
        json.dump({"pid": os.getpid(), "started_at": now(), "mode": args.mode, "output": str(output)}, lock_file)
    environment = os.environ.copy()
    if args.threads is not None:
        environment["OMP_NUM_THREADS"] = str(args.threads)
    cache = output / "packing_cache"
    metadata: dict[str, object] = {}
    runs: list[dict[str, object]] = []
    try:
        metadata = {
            "mode": args.mode, "created_at": now(), "requested_runs": args.runs,
            "executable": file_identity(executable), "dataset": file_identity(dataset),
            "benchmark_script": file_identity(Path(__file__).resolve()),
            "system": system_metadata(),
            "thread_environment": {key: environment.get(key) for key in (
                "OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PROC_BIND", "OMP_PLACES",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
            "timer": "perf_counter around fresh Popen through process termination, including startup/context/keys, excluding log parsing",
            "memory": "Windows GetProcessMemoryInfo PeakWorkingSetSize for the single executable process, including its threads; bytes/decimal MB/MiB are separately labeled",
            "memory_poll_ms": args.poll_ms,
            "warm_definition": "Runs 2 onward; every run is a new process. OS cache is not flushed and binary/dataset hashing occurs before run 1.",
            "cache_policy": ("No saved context/key directory; all runs generate fresh keys" if args.mode == "analysis" else
                             "Run 1 uses a unique absent application cache; runs 2 onward reuse saved keys and existing ciphertext packs. No existing user cache is modified."),
            "workload": ("Final-decrypt sightline plus boundary-slope earthwork over the demo rectangle; desired height and excavation angle are direct arguments" if args.mode == "analysis" else
                         "Original packed_store_validation: Taiwan DEM 2-stitch packing plus 256x256 synthetic-fixture 4-stitch complex packing; includes complex probe, save/load/split, cache-clear reload, and right-half slope parity check"),
            "packing_cache": str(cache) if args.mode == "packing" else None,
            "processes_are_sequential": True,
        }
        write_json(output / "metadata.json", metadata)
        for run in range(1, args.runs + 1):
            command = [str(executable), str(dataset)]
            before = None
            if args.mode == "analysis":
                command += [str(args.height), str(args.angle)]
                if args.grid_interval is not None:
                    command.append(str(args.grid_interval))
            else:
                before = cache_inventory(cache)
                if run == 1 and (cache.exists() or before["file_count"]):
                    raise RuntimeError("Cold packing cache must not preexist")
                command.append(str(cache))
            result = measure(command, output, run, args.timeout, args.poll_ms / 1000.0, environment)
            if args.mode == "packing":
                result["cache_before"] = before
                result["cache_after"] = cache_inventory(cache)
            runs.append(result)
            write_json(output / "results.json", runs)
            if result["exit_code"] != 0 or result["timed_out"] or result["cancelled"]:
                raise RuntimeError(f"Run {run} failed; see {result['stderr']}")
            if not result["peak_working_set_bytes"]:
                raise RuntimeError(f"Run {run}: no valid peak working set was captured")
            text = Path(result["stdout"]).read_text(encoding="utf-8", errors="replace")
            result["program_result"] = parse_program_output(args.mode, text)
            if args.mode == "packing":
                for name in ("pack", "pack16"):
                    counts = result["program_result"][name]
                    if run == 1 and (counts["encrypted"] != counts["total"] or counts["skipped"] != 0):
                        raise RuntimeError(f"Cold run unexpectedly reused {name} packs")
                    if run > 1 and (counts["skipped"] != counts["total"] or counts["encrypted"] != 0):
                        raise RuntimeError(f"Warm run unexpectedly regenerated {name} packs")
            write_json(output / "results.json", runs)
            print(json.dumps({key: result[key] for key in (
                "run", "exit_code", "wall_seconds", "peak_working_set_MB_decimal",
                "peak_working_set_MiB", "final_peak_counter_available")}), flush=True)
        warm = runs[1:]
        summary = {
            "passed": True, "mode": args.mode, "completed_at": now(), "runs": len(runs),
            "first_wall_seconds": runs[0]["wall_seconds"],
            "subsequent_mean_wall_seconds": statistics.mean(r["wall_seconds"] for r in warm),
            "subsequent_range_wall_seconds": [min(r["wall_seconds"] for r in warm), max(r["wall_seconds"] for r in warm)],
            "all_runs_max_peak_working_set_bytes": max(r["peak_working_set_bytes"] for r in runs),
            "all_runs_max_peak_working_set_MB_decimal": max(r["peak_working_set_MB_decimal"] for r in runs),
            "all_runs_max_peak_working_set_MiB": max(r["peak_working_set_MiB"] for r in runs),
            "all_final_peak_counters_available": all(r["final_peak_counter_available"] for r in runs),
            "interpretation": metadata["cache_policy"],
            "workload": metadata["workload"],
            "warning": "Observed whole-process costs are specific to this executable, machine, input, and cache policy; they do not isolate a single cryptographic operation.",
        }
        write_json(output / "summary.json", summary)
        print("SUMMARY " + json.dumps(summary), flush=True)
        print(str(output), flush=True)
        return 0
    except BaseException as error:
        write_json(output / "summary.json", {"passed": False, "mode": args.mode,
                   "completed_at": now(), "completed_run_records": len(runs),
                   "error": f"{type(error).__name__}: {error}"})
        raise
    finally:
        lock_path.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
