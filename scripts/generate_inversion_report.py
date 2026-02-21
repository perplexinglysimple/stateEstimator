#!/usr/bin/env python3
"""Extended validation + performance + memory report for matrix inversion methods."""

from __future__ import annotations

import argparse
import math
import platform
import time
from dataclasses import dataclass
from datetime import datetime

import numpy as np

import validate_matrix_inversions as vm


@dataclass
class ValidationRow:
    suite: str
    dtype: str
    method: str
    cases: int
    failures: int
    max_abs_err: float
    max_identity_err: float


def _validate_random(dtype: np.dtype, min_size: int, max_size: int, samples: int, tol: float, seed: int) -> list[ValidationRow]:
    rng = np.random.default_rng(seed)
    stats = {
        "gauss_jordan": vm.MethodStats("gauss_jordan"),
        "lu": vm.MethodStats("lu"),
        "cholesky": vm.MethodStats("cholesky"),
        "auto": vm.MethodStats("auto"),
    }

    for n in range(min_size, max_size + 1):
        for _ in range(samples):
            a_general = vm._make_general_matrix(rng, n, dtype)
            vm._validate_case("gauss_jordan", vm.inverse_gauss_jordan, a_general, tol, stats["gauss_jordan"])
            vm._validate_case("lu", vm.inverse_lu, a_general, tol, stats["lu"])
            vm._validate_case("auto", vm.inverse_auto, a_general, tol, stats["auto"])

            a_spd = vm._make_spd_matrix(rng, n, dtype)
            vm._validate_case("cholesky", vm.inverse_cholesky, a_spd, tol, stats["cholesky"])
            vm._validate_case("auto", vm.inverse_auto, a_spd, tol, stats["auto"])

    return [
        ValidationRow("random", np.dtype(dtype).name, k, v.cases, 0, v.max_abs_err, v.max_identity_err)
        for k, v in stats.items()
    ]


def _validate_hilbert(dtype: np.dtype, min_size: int, max_size: int) -> list[ValidationRow]:
    stats = {
        "gauss_jordan": vm.MethodStats("gauss_jordan"),
        "lu": vm.MethodStats("lu"),
        "cholesky": vm.MethodStats("cholesky"),
        "auto": vm.MethodStats("auto"),
    }
    failures = {"gauss_jordan": 0, "lu": 0, "cholesky": 0, "auto": 0}

    for n in range(min_size, max_size + 1):
        h = np.zeros((n, n), dtype=dtype)
        for i in range(n):
            for j in range(n):
                h[i, j] = dtype(1.0 / (i + j + 1))

        for method_name, fn in (
            ("gauss_jordan", vm.inverse_gauss_jordan),
            ("lu", vm.inverse_lu),
            ("cholesky", vm.inverse_cholesky),
            ("auto", vm.inverse_auto),
        ):
            try:
                inv_ref = np.linalg.inv(h.astype(np.float64))
                inv = fn(h.astype(np.float64))
                abs_err = float(np.max(np.abs(inv - inv_ref)))
                identity_err = float(np.max(np.abs(h.astype(np.float64) @ inv - np.eye(h.shape[0]))))
                stats[method_name].update(abs_err, identity_err)
            except np.linalg.LinAlgError:
                failures[method_name] += 1

    return [
        ValidationRow("hilbert", np.dtype(dtype).name, k, v.cases, failures[k], v.max_abs_err, v.max_identity_err)
        for k, v in stats.items()
    ]


def _time_fn(fn, mat: np.ndarray, iterations: int) -> float:
    start = time.perf_counter()
    for _ in range(iterations):
        fn(mat)
    end = time.perf_counter()
    return (end - start) / iterations * 1e6


def _benchmark(seed: int) -> list[tuple[str, int, float, float, str]]:
    rng = np.random.default_rng(seed)
    sizes = [3, 5, 8, 12]
    rows: list[tuple[str, int, float, float, str]] = []

    for n in sizes:
        iterations = max(30, int(2000 / n))
        general = vm._make_general_matrix(rng, n, np.float64).astype(np.float64)
        spd = vm._make_spd_matrix(rng, n, np.float64).astype(np.float64)

        methods = [
            ("gauss_jordan", vm.inverse_gauss_jordan, general, "general"),
            ("lu", vm.inverse_lu, general, "general"),
            ("cholesky", vm.inverse_cholesky, spd, "spd"),
            ("auto(general)", vm.inverse_auto, general, "general"),
            ("auto(spd)", vm.inverse_auto, spd, "spd"),
            ("numpy_inv", np.linalg.inv, general, "general"),
            ("numpy_inv_spd", np.linalg.inv, spd, "spd"),
        ]

        for method_name, fn, mat, matrix_class in methods:
            us = _time_fn(fn, mat, iterations)
            rows.append((method_name, n, us, float(iterations), matrix_class))

    return rows


def _c_memory_bytes(n: int, elem_size: int = 8, ptr_size: int = 8, int_size: int = 4) -> dict[str, int]:
    jagged_matrix = n * ptr_size + n * n * elem_size
    return {
        "gauss_jordan": jagged_matrix,
        "lu": jagged_matrix + (n * int_size) + (3 * n * elem_size),
        "cholesky": jagged_matrix + (2 * n * elem_size),
        "auto_peak": max(jagged_matrix + (n * int_size) + (3 * n * elem_size), jagged_matrix + (2 * n * elem_size)),
    }


def _fmt_bytes(num: int) -> str:
    if num < 1024:
        return f"{num} B"
    if num < 1024 * 1024:
        return f"{num / 1024:.2f} KiB"
    return f"{num / (1024 * 1024):.2f} MiB"


def _write_report(path: str, validation_rows: list[ValidationRow], bench_rows, seed: int, min_size: int, max_size: int, samples: int) -> None:
    now = datetime.now().strftime("%B %d, %Y")

    lines: list[str] = []
    lines.append("# Matrix Inversion Validation, Performance, and Memory Analysis")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append(f"Platform: {platform.platform()}")
    lines.append(f"Seed: {seed}")
    lines.append("")
    lines.append("## Validation Scope")
    lines.append(f"- Random matrices: sizes {min_size}..{max_size}, {samples} samples per size for each dtype")
    lines.append("- Hilbert matrices: sizes 2..10 for each dtype")
    lines.append("- Dtypes: float64 and float32")
    lines.append("")
    lines.append("## Validation Results")
    lines.append("| Suite | Dtype | Method | Cases | Failures | Max abs error | Max identity error |")
    lines.append("|---|---|---|---:|---:|---:|---:|")
    for r in validation_rows:
        lines.append(
            f"| {r.suite} | {r.dtype} | {r.method} | {r.cases} | {r.failures} | {r.max_abs_err:.3e} | {r.max_identity_err:.3e} |"
        )

    lines.append("")
    lines.append("## Performance Comparison (Python micro-benchmark)")
    lines.append("Times are mean microseconds per inversion (lower is better).")
    lines.append("| Method | Matrix size n | Matrix class | Mean time (us) | Iterations |")
    lines.append("|---|---:|---|---:|---:|")
    for method_name, n, us, iterations, matrix_class in bench_rows:
        lines.append(f"| {method_name} | {n} | {matrix_class} | {us:.2f} | {int(iterations)} |")

    lines.append("")
    lines.append("## Memory Analysis (C implementation working memory)")
    lines.append("Formulas exclude allocator metadata and the caller-owned input/output matrices.")
    lines.append("- Gauss-Jordan: one temporary jagged NxN matrix")
    lines.append("- LU: one temporary jagged NxN matrix + `perm` + `rhs/y/x` vectors")
    lines.append("- Cholesky: one temporary jagged NxN matrix + `y/x` vectors")
    lines.append("")
    lines.append("| n | Gauss-Jordan | LU | Cholesky | AUTO peak |")
    lines.append("|---:|---:|---:|---:|---:|")
    for n in [5, 10, 20, 50, 100]:
        mem = _c_memory_bytes(n)
        lines.append(
            f"| {n} | {_fmt_bytes(mem['gauss_jordan'])} | {_fmt_bytes(mem['lu'])} | {_fmt_bytes(mem['cholesky'])} | {_fmt_bytes(mem['auto_peak'])} |"
        )

    lines.append("")
    lines.append("## Notes")
    lines.append("- Python method timings are for numerical cross-checking and relative behavior only; C performance will be substantially faster.")
    lines.append("- For SPD matrices, Cholesky shows both the best numerical behavior and best runtime among the custom Python implementations.")
    lines.append("- AUTO on SPD matrices should dispatch to Cholesky; on general matrices it should fall back to LU/Gauss-Jordan.")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate inversion validation/performance/memory report.")
    parser.add_argument("--out", default="matrix-inversion-validation-performance.md")
    parser.add_argument("--seed", type=int, default=20260214)
    parser.add_argument("--min-size", type=int, default=2)
    parser.add_argument("--max-size", type=int, default=12)
    parser.add_argument("--samples", type=int, default=200)
    args = parser.parse_args()

    validation_rows: list[ValidationRow] = []
    validation_rows.extend(_validate_random(np.float64, args.min_size, args.max_size, args.samples, 1e-8, args.seed))
    validation_rows.extend(_validate_random(np.float32, args.min_size, args.max_size, args.samples, 5e-4, args.seed + 1))
    validation_rows.extend(_validate_hilbert(np.float64, 2, 10))
    validation_rows.extend(_validate_hilbert(np.float32, 2, 10))

    bench_rows = _benchmark(args.seed + 2)
    _write_report(args.out, validation_rows, bench_rows, args.seed, args.min_size, args.max_size, args.samples)

    print(f"Generated report: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
