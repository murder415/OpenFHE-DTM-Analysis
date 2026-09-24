"""Independent post-run checks of sightline operation DAGs and numeric CSVs.

This does not edit the C++ implementation, rerun HE, or turn logged values into
claims about survey accuracy. Run again after regenerating the executable logs.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def fingerprint(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}


def verify_trace(path: Path, samples: int, slots: int) -> dict[str, object]:
    events = read_csv(path)
    batches = (samples + slots - 1) // slots
    nodes: dict[int, dict[str, object]] = {}
    roots: list[int] = []
    counts: collections.Counter[str] = collections.Counter()
    seen_decrypt = False
    for event_number, row in enumerate(events):
        event, output, a, b = (int(row[name]) for name in ("event", "output_cipher", "input_a", "input_b"))
        op = row["operation"]
        check(event == event_number, f"{path}: event sequence gap")
        check(not seen_decrypt or op == "decrypt", f"{path}: HE work after first decrypt")
        for input_id in (a, b):
            check(not input_id or input_id in nodes, f"{path}: missing/forward ciphertext input {input_id}")
        counts[op] += 1
        if op == "decrypt":
            check(output == 0 and a > 0 and b == 0, f"{path}: malformed decrypt event")
            roots.append(a)
            seen_decrypt = True
        elif op == "encode":
            check(output == a == b == 0, f"{path}: malformed encode event")
        else:
            check(output > 0 and output not in nodes, f"{path}: repeated/absent output identity")
            nodes[output] = {"op": op, "a": a, "b": b}
    expected_counts = collections.Counter({
        "encode": 17 * batches,
        "encrypt": 8 * batches,
        "mulPlain": 8 * batches,
        "add": 6 * batches,
        "subPlain": batches,
        "sub": batches,
        "decrypt": batches,
    })
    check(counts == expected_counts, f"{path}: operation counts {counts} != {expected_counts}")
    check(len(roots) == batches and len(set(roots)) == batches, f"{path}: repeated/missing terminal batch")
    check(len(nodes) == 24 * batches, f"{path}: unexpected ciphertext node count")

    def interpolation_tree(root: int) -> set[int]:
        visited: set[int] = set()
        tree_counts: collections.Counter[str] = collections.Counter()

        def visit(node_id: int) -> None:
            check(node_id not in visited, f"{path}: reused node within interpolation sum")
            visited.add(node_id)
            node = nodes[node_id]
            op, a, b = node["op"], int(node["a"]), int(node["b"])
            tree_counts[str(op)] += 1
            if op == "add":
                check(a > 0 and b > 0, f"{path}: malformed interpolation add")
                visit(a)
                visit(b)
            elif op == "mulPlain":
                check(a > 0 and b == 0, f"{path}: malformed plaintext multiplication")
                check(nodes[a]["op"] == "encrypt", f"{path}: interpolation leaf is not encrypt->mulPlain")
                visit(a)
            elif op == "encrypt":
                check(a == b == 0, f"{path}: malformed source encryption")
            else:
                raise AssertionError(f"{path}: unexpected op {op} in interpolation tree")

        check(nodes[root]["op"] == "add", f"{path}: interpolation root is not an addition")
        visit(root)
        check(tree_counts == {"add": 3, "mulPlain": 4, "encrypt": 4},
              f"{path}: interpolation is not a sum of four encrypted weighted corners: {tree_counts}")
        return visited

    all_used: set[int] = set()
    batch_records = []
    for batch_number, root in enumerate(roots):
        final = nodes[root]
        check(final["op"] == "sub", f"{path}: decrypted value is not final height subtraction")
        plane_id, ground_id = int(final["a"]), int(final["b"])
        plane = nodes[plane_id]
        check(plane["op"] == "subPlain" and plane["b"] == 0,
              f"{path}: final left branch is not base minus negative plaintext rise")
        base_id = int(plane["a"])
        ground_nodes, base_nodes = interpolation_tree(ground_id), interpolation_tree(base_id)
        check(not (ground_nodes & base_nodes), f"{path}: ground/base interpolation trees overlap")
        batch_nodes = ground_nodes | base_nodes | {root, plane_id}
        check(len(batch_nodes) == 24, f"{path}: batch topology has extra/shared nodes")
        check(not (all_used & batch_nodes), f"{path}: batches share or reuse ciphertext nodes")
        all_used |= batch_nodes
        batch_records.append({
            "batch": batch_number,
            "active_samples": min(slots, samples - batch_number * slots),
            "decrypted_root": root,
            "ground_interpolation_root": ground_id,
            "base_interpolation_root": base_id,
            "ciphertext_nodes": 24,
            "valid": True,
        })
    check(all_used == set(nodes), f"{path}: unused ciphertext branch was hidden from terminal topology")
    return {
        "passed": True,
        "counts": dict(counts),
        "batch_details": batch_records,
        "all_he_before_first_decrypt": True,
        "unused_or_unverified_ciphertext_nodes": 0,
    }


def metrics(errors: list[float]) -> dict[str, float]:
    return {
        "max_abs_error_m": max(errors),
        "mae_m": math.fsum(errors) / len(errors),
        "rmse_m": math.sqrt(math.fsum(e * e for e in errors) / len(errors)),
    }


def verify_numeric(path: Path, case: dict[str, object]) -> dict[str, object]:
    rows = read_csv(path)
    check(len(rows) == int(case["samples"]), f"{path}: sample count mismatch")
    tangent = math.tan(float(case["angle_degrees"]) * math.pi / 180.0)
    same_raw_errors, same_height_errors, independent_height_errors = [], [], []
    distance_differences, baseline_differences, original_csv_reconstruction_differences = [], [], []
    near_zero_sign_mismatches = 0
    for i, row in enumerate(rows):
        check(int(row["sample"]) == i, f"{path}: sample ordering mismatch")
        values = {key: float(value) for key, value in row.items()}
        check(all(math.isfinite(v) for v in values.values()), f"{path}: nonfinite input")
        base, ground = values["plain_base_m"], values["plain_ground_m"]
        same_diff = (base + values["returned_distance_m"] * tangent) - ground
        same_height = max(0.0, same_diff)
        independent_diff = (base + values["independent_distance_m"] * tangent) - ground
        raw, actual = values["decoded_difference_m"], values["he_height_m"]
        check(actual == max(0.0, raw), f"{path}: output was not exactly max(0,decoded)")
        same_raw_errors.append(abs(raw - same_diff))
        same_height_errors.append(abs(actual - same_height))
        independent_height_errors.append(abs(actual - values["plain_height_m"]))
        distance_differences.append(abs(values["returned_distance_m"] - values["independent_distance_m"]))
        baseline_differences.append(abs(same_diff - values["plain_difference_m"]))
        original_csv_reconstruction_differences.append(abs(independent_diff - values["plain_difference_m"]))
        if abs(same_diff) <= 1e-7 and ((same_diff > 0) != (raw > 0)):
            near_zero_sign_mismatches += 1
    check(max(same_raw_errors) < 1e-6 and max(same_height_errors) < 1e-6,
          f"{path}: same-distance validation tolerance exceeded")
    check(math.isclose(max(independent_height_errors), float(case["max_abs_error_m"]),
                       rel_tol=1e-12, abs_tol=1e-24), f"{path}: original reported error does not match CSV")
    return {
        "passed": True,
        "same_returned_distance_raw_difference": metrics(same_raw_errors),
        "same_returned_distance_final_height": metrics(same_height_errors),
        "independent_distance_final_height": metrics(independent_height_errors),
        "max_distance_oracle_difference_m": max(distance_differences),
        "max_plaintext_baseline_difference_m": max(baseline_differences),
        "max_python_vs_recorded_plaintext_reconstruction_difference_m": max(original_csv_reconstruction_differences),
        "same_distance_near_zero_sign_mismatches": near_zero_sign_mismatches,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log_root", nargs="?", type=Path,
                        default=Path(__file__).parent / "logs" / "sightline")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.log_root / "independent_trace_audit.json"
    records, sources = [], []
    report_files = sorted(args.log_root.glob("run_*/results.json"))
    check(bool(report_files), "No per-run results were found")
    for report_path in report_files:
        report = json.loads(report_path.read_text(encoding="utf-8-sig"))
        sources.append(fingerprint(report_path))
        case_names = {case["name"] for case in report["cases"]}
        traces = list(report_path.parent.glob("*_trace.csv"))
        check({p.name.removesuffix("_trace.csv") for p in traces} == case_names,
              f"{report_path}: trace/case inventory differs")
        for case in report["cases"]:
            trace = report_path.parent / (str(case["name"]) + "_trace.csv")
            samples = report_path.parent / (str(case["name"]) + "_samples.csv")
            sources.extend([fingerprint(trace), fingerprint(samples)])
            records.append({
                "run": report["run"], "name": case["name"], "samples": case["samples"],
                "topology": verify_trace(trace, int(case["samples"]), int(report["slots"])),
                "numeric": verify_numeric(samples, case),
            })
    taiwan = [record for record in records if record["name"] == "taiwan_demo_10"]
    result = {
        "passed": True,
        "runs": len(report_files), "cases": len(records),
        "terminal_decrypts_verified": sum(r["topology"]["counts"]["decrypt"] for r in records),
        "scope": "Read-only postprocessing of captured backend interface traces and numeric CSVs. No C++ changes or HE reruns.",
        "topology_contract": "For every decrypted root: sub(subPlain(base weighted-four-corner interpolation, negative rise), ground weighted-four-corner interpolation). Every interpolation has exactly 4 encrypt, 4 mulPlain, 3 add; all ciphertext nodes belong to exactly one terminal batch and all HE calls precede the first decrypt.",
        "limitations": [
            "Trace does not record plaintext encode arrays/weights or primitive operations internal to OpenFHE. It verifies the exposed IHeBackend call graph; correctness of slot values is additionally supported by the numerical tests and source review.",
            "The same-returned-distance baseline removes the independent geodesic distance difference, but still includes ordinary floating-point rounding and possible Python/C++ libm differences. Neither metric isolates pure CKKS error.",
            "Stored intermediate plaintext baseline elevations are external validation data, not output from the production sightline method.",
            "These traces may precede input-guard-only changes; rerun this script after regenerating final executable logs.",
        ],
        "taiwan_same_distance_final_height_max_error_range_m": [
            min(r["numeric"]["same_returned_distance_final_height"]["max_abs_error_m"] for r in taiwan),
            max(r["numeric"]["same_returned_distance_final_height"]["max_abs_error_m"] for r in taiwan),
        ],
        "taiwan_same_distance_raw_difference_max_error_range_m": [
            min(r["numeric"]["same_returned_distance_raw_difference"]["max_abs_error_m"] for r in taiwan),
            max(r["numeric"]["same_returned_distance_raw_difference"]["max_abs_error_m"] for r in taiwan),
        ],
        "records": records, "inputs": sources,
    }
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "passed", "runs", "cases", "terminal_decrypts_verified",
        "taiwan_same_distance_final_height_max_error_range_m",
        "taiwan_same_distance_raw_difference_max_error_range_m")}, ensure_ascii=False))
    print(str(output.resolve()))


if __name__ == "__main__":
    main()
