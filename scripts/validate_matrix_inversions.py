#!/usr/bin/env python3
"""Validate matrix inversion algorithms against NumPy references.

This script mirrors the methods introduced in src/matrixInversion.c:
- Gauss-Jordan with partial pivoting
- LU with partial pivoting
- Cholesky (SPD)
- AUTO dispatch (Cholesky -> LU -> Gauss-Jordan)
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass
class MethodStats:
    name: str
    cases: int = 0
    max_abs_err: float = 0.0
    max_identity_err: float = 0.0

    def update(self, abs_err: float, identity_err: float) -> None:
        self.cases += 1
        self.max_abs_err = max(self.max_abs_err, abs_err)
        self.max_identity_err = max(self.max_identity_err, identity_err)


def _pivot_epsilon(dtype: np.dtype) -> float:
    if dtype == np.float32:
        return 1e-6
    return 1e-12


def _is_symmetric(a: np.ndarray, tol: float) -> bool:
    return np.max(np.abs(a - a.T)) <= tol


def inverse_gauss_jordan(a: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    eps = _pivot_epsilon(a.dtype)
    work = a.astype(a.dtype, copy=True)
    inv = np.eye(n, dtype=a.dtype)

    for i in range(n):
        pivot_row = i + int(np.argmax(np.abs(work[i:, i])))
        pivot = work[pivot_row, i]
        if abs(float(pivot)) <= eps:
            raise np.linalg.LinAlgError("singular in Gauss-Jordan")

        if pivot_row != i:
            work[[i, pivot_row], :] = work[[pivot_row, i], :]
            inv[[i, pivot_row], :] = inv[[pivot_row, i], :]

        pivot = work[i, i]
        work[i, :] /= pivot
        inv[i, :] /= pivot

        for r in range(n):
            if r == i:
                continue
            factor = work[r, i]
            if abs(float(factor)) <= eps:
                work[r, i] = 0.0
                continue
            work[r, :] -= factor * work[i, :]
            inv[r, :] -= factor * inv[i, :]
            work[r, i] = 0.0

    return inv


def _lu_decompose(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = a.shape[0]
    eps = _pivot_epsilon(a.dtype)
    lu = a.astype(a.dtype, copy=True)
    perm = np.arange(n)

    for k in range(n):
        pivot_row = k + int(np.argmax(np.abs(lu[k:, k])))
        pivot = lu[pivot_row, k]
        if abs(float(pivot)) <= eps:
            raise np.linalg.LinAlgError("singular in LU")

        if pivot_row != k:
            lu[[k, pivot_row], :] = lu[[pivot_row, k], :]
            perm[[k, pivot_row]] = perm[[pivot_row, k]]

        for i in range(k + 1, n):
            lu[i, k] /= lu[k, k]
            lu[i, k + 1 :] -= lu[i, k] * lu[k, k + 1 :]

    return lu, perm


def inverse_lu(a: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    eps = _pivot_epsilon(a.dtype)
    lu, perm = _lu_decompose(a)
    inv = np.zeros_like(a)

    for col in range(n):
        rhs = (perm == col).astype(a.dtype)

        y = np.zeros(n, dtype=a.dtype)
        for i in range(n):
            y[i] = rhs[i] - np.dot(lu[i, :i], y[:i])

        x = np.zeros(n, dtype=a.dtype)
        for i in range(n - 1, -1, -1):
            pivot = lu[i, i]
            if abs(float(pivot)) <= eps:
                raise np.linalg.LinAlgError("singular in LU back solve")
            x[i] = (y[i] - np.dot(lu[i, i + 1 :], x[i + 1 :])) / pivot

        inv[:, col] = x

    return inv


def inverse_cholesky(a: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    eps = _pivot_epsilon(a.dtype)
    sym_tol = 10.0 * eps
    if not _is_symmetric(a, sym_tol):
        raise np.linalg.LinAlgError("not symmetric for Cholesky")

    l = np.zeros_like(a)
    for i in range(n):
        for j in range(i + 1):
            s = a[i, j] - np.dot(l[i, :j], l[j, :j])
            if i == j:
                if s <= eps:
                    raise np.linalg.LinAlgError("not SPD for Cholesky")
                l[i, j] = np.sqrt(s)
            else:
                if abs(float(l[j, j])) <= eps:
                    raise np.linalg.LinAlgError("zero diagonal in Cholesky")
                l[i, j] = s / l[j, j]

    inv = np.zeros_like(a)
    for col in range(n):
        rhs = np.zeros(n, dtype=a.dtype)
        rhs[col] = 1.0

        y = np.zeros(n, dtype=a.dtype)
        for i in range(n):
            y[i] = (rhs[i] - np.dot(l[i, :i], y[:i])) / l[i, i]

        x = np.zeros(n, dtype=a.dtype)
        for i in range(n - 1, -1, -1):
            x[i] = (y[i] - np.dot(l[i + 1 :, i], x[i + 1 :])) / l[i, i]

        inv[:, col] = x

    return 0.5 * (inv + inv.T)


def inverse_auto(a: np.ndarray) -> np.ndarray:
    for fn in (inverse_cholesky, inverse_lu, inverse_gauss_jordan):
        try:
            return fn(a)
        except np.linalg.LinAlgError:
            continue
    raise np.linalg.LinAlgError("AUTO failed for matrix")


def _make_general_matrix(rng: np.random.Generator, n: int, dtype: np.dtype) -> np.ndarray:
    for _ in range(200):
        a = rng.normal(0.0, 1.0, size=(n, n)).astype(dtype)
        a += np.eye(n, dtype=dtype) * dtype(0.2)
        if abs(np.linalg.det(a.astype(np.float64))) > 1e-10:
            return a
    raise RuntimeError(f"could not generate invertible general matrix (n={n})")


def _make_spd_matrix(rng: np.random.Generator, n: int, dtype: np.dtype) -> np.ndarray:
    m = rng.normal(0.0, 1.0, size=(n, n)).astype(dtype)
    a = (m @ m.T).astype(dtype)
    a += np.eye(n, dtype=dtype) * dtype(0.5)
    return a


def _validate_case(
    method_name: str,
    fn: Callable[[np.ndarray], np.ndarray],
    a: np.ndarray,
    tol: float,
    stats: MethodStats,
) -> None:
    inv_ref = np.linalg.inv(a.astype(np.float64))
    inv = fn(a.astype(np.float64))

    abs_err = float(np.max(np.abs(inv - inv_ref)))
    identity_err = float(np.max(np.abs(a.astype(np.float64) @ inv - np.eye(a.shape[0]))))
    stats.update(abs_err, identity_err)

    if abs_err > tol or identity_err > tol:
        raise AssertionError(
            f"{method_name} failed: abs_err={abs_err:.3e}, identity_err={identity_err:.3e}, tol={tol:.3e}"
        )


def run_validation(
    min_size: int,
    max_size: int,
    samples_per_size: int,
    seed: int,
    tol: float,
) -> list[MethodStats]:
    rng = np.random.default_rng(seed)
    stats = {
        "gauss_jordan": MethodStats("gauss_jordan"),
        "lu": MethodStats("lu"),
        "cholesky": MethodStats("cholesky"),
        "auto": MethodStats("auto"),
    }

    for n in range(min_size, max_size + 1):
        for _ in range(samples_per_size):
            a_general = _make_general_matrix(rng, n, np.float64)
            _validate_case("gauss_jordan", inverse_gauss_jordan, a_general, tol, stats["gauss_jordan"])
            _validate_case("lu", inverse_lu, a_general, tol, stats["lu"])
            _validate_case("auto", inverse_auto, a_general, tol, stats["auto"])

            a_spd = _make_spd_matrix(rng, n, np.float64)
            _validate_case("cholesky", inverse_cholesky, a_spd, tol, stats["cholesky"])
            _validate_case("auto", inverse_auto, a_spd, tol, stats["auto"])

    return [stats["gauss_jordan"], stats["lu"], stats["cholesky"], stats["auto"]]


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate new matrix inversion methods.")
    parser.add_argument("--min-size", type=int, default=2)
    parser.add_argument("--max-size", type=int, default=8)
    parser.add_argument("--samples", type=int, default=100, help="samples per matrix size")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tol", type=float, default=1e-8)
    args = parser.parse_args()

    if args.min_size < 2 or args.max_size < args.min_size:
        raise ValueError("invalid size bounds")

    stats = run_validation(
        min_size=args.min_size,
        max_size=args.max_size,
        samples_per_size=args.samples,
        seed=args.seed,
        tol=args.tol,
    )

    print("Matrix inversion validation passed.")
    for s in stats:
        print(
            f"- {s.name}: cases={s.cases}, max_abs_err={s.max_abs_err:.3e}, "
            f"max_identity_err={s.max_identity_err:.3e}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
