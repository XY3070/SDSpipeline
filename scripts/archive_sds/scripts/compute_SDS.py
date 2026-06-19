#!/usr/bin/env python3
"""
Optimized Singleton Density Score (SDS) computation.

- Feather for test SNP intermediate exchange
- Parquet for archived final output
- CSV for small summaries/validation
- Pickle cache for expensive Python-side parsed inputs
"""

import argparse
import csv
import fcntl
import gzip
import hashlib
import io
import logging
import os
import pickle
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, Optional

import numpy as np
import pyarrow as pa
import pyarrow.feather as feather
import pyarrow.parquet as pq
from scipy.optimize import minimize

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
repo_root_str = str(REPO_ROOT)
if repo_root_str not in sys.path:
    sys.path.insert(0, repo_root_str)

try:
    import numba
    from numba import prange
except ImportError:
    numba = None
    prange = range  # fallback: plain range when numba is not available

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

RESULT_COLUMNS = ['ID', 'AA', 'DA', 'POS', 'DAF', 'nG0', 'nG1', 'nG2', 'rSDS', 'SuggestedInitPoint']


def _jit_or_noop(*jit_args, **jit_kwargs):
    if numba is None:
        def _decorator(func):
            return func
        return _decorator
    return numba.jit(*jit_args, **jit_kwargs)


@_jit_or_noop(nopython=True, cache=True)
def _logsum_vec(log_a, log_b):
    """Element-wise log(exp(a) + exp(b)) for vector + scalar."""
    n = len(log_a)
    out = np.empty(n)
    for k in range(n):
        hi = max(log_a[k], log_b)
        lo = min(log_a[k], log_b)
        out[k] = hi + np.log1p(np.exp(lo - hi))
    return out


@_jit_or_noop(nopython=True, cache=True)
def _logsum_scalar(a, b):
    hi = max(a, b)
    lo = min(a, b)
    return hi + np.log1p(np.exp(lo - hi))


@_jit_or_noop(nopython=True, cache=True)
def _logsum_pairwise(log_a, log_b):
    """Element-wise log(exp(a) + exp(b)) for vector + vector."""
    n = len(log_a)
    out = np.empty(n)
    for k in range(n):
        hi = max(log_a[k], log_b[k])
        lo = min(log_a[k], log_b[k])
        out[k] = hi + np.log1p(np.exp(lo - hi))
    return out


@_jit_or_noop(nopython=True, cache=True)
def _nll_kernel(params, dat0, dat1, dat2, A1, A2, eps, ls_const1, ls_const2):
    logE1, logE2 = params[0], params[1]
    logA1, logA2 = np.log(A1), np.log(A2)
    logB1, logB2 = logA1 - logE1, logA2 - logE2
    ll = 0.0
    n0 = len(dat0)
    n1 = len(dat1)
    n2 = len(dat2)

    if n0 > 0:
        log_d0 = np.log(dat0 + eps)
        ls0 = _logsum_vec(log_d0, logB1)
        ll += n0 * (
            2.0 * A1 * (logB1 - np.mean(ls0))
            + np.mean(log_d0)
            + np.log(2.0)
            + logA1
            - 2.0 * np.mean(ls0)
            + ls_const1
        )
    if n2 > 0:
        log_d2 = np.log(dat2 + eps)
        ls2 = _logsum_vec(log_d2, logB2)
        ll += n2 * (
            2.0 * A2 * (logB2 - np.mean(ls2))
            + np.mean(log_d2)
            + np.log(2.0)
            + logA2
            - 2.0 * np.mean(ls2)
            + ls_const2
        )
    if n1 > 0:
        log_d1 = np.log(dat1 + eps)
        ls1B1 = _logsum_vec(log_d1, logB1)
        ls1B2 = _logsum_vec(log_d1, logB2)
        sub1 = -2.0 * ls1B1 + logA1 + _logsum_scalar(logA1, 0.0)
        sub2 = -2.0 * ls1B2 + logA2 + _logsum_scalar(logA2, 0.0)
        sub3 = np.log(2.0) + logA1 + logA2 - ls1B1 - ls1B2
        ll += n1 * (
            A1 * (logB1 - np.mean(ls1B1))
            + A2 * (logB2 - np.mean(ls1B2))
            + np.mean(log_d1)
            + np.mean(_logsum_pairwise(_logsum_pairwise(sub1, sub2), sub3))
        )

    if np.isnan(ll) or np.isinf(ll):
        return 1e15
    return -ll


@_jit_or_noop(nopython=True, cache=True, parallel=True)
def _nll_kernel_batch(param_matrix, dat0, dat1, dat2, A1, A2, eps, ls_const1, ls_const2):
    """[OPT] parallel prange over grid points."""
    k = param_matrix.shape[0]
    out = np.empty(k)
    for i in prange(k):
        out[i] = _nll_kernel(param_matrix[i], dat0, dat1, dat2, A1, A2, eps, ls_const1, ls_const2)
    return out


@_jit_or_noop(nopython=True, cache=True)
def _nll_grad_kernel(params, dat0, dat1, dat2, A1, A2, eps, ls_const1, ls_const2):
    logE1, logE2 = params[0], params[1]
    logA1 = np.log(A1)
    logA2 = np.log(A2)
    logB1 = logA1 - logE1
    logB2 = logA2 - logE2
    n0 = len(dat0)
    n1 = len(dat1)
    n2 = len(dat2)
    ll = 0.0
    dll_dlogE1 = 0.0
    dll_dlogE2 = 0.0

    if n0 > 0:
        log_d0 = np.log(dat0 + eps)
        ls0 = _logsum_vec(log_d0, logB1)
        mean_ls0 = np.mean(ls0)
        mean_log_d0 = np.mean(log_d0)
        ll += n0 * (
            2.0 * A1 * (logB1 - mean_ls0)
            + mean_log_d0
            + np.log(2.0)
            + logA1
            - 2.0 * mean_ls0
            + ls_const1
        )
        w0 = np.mean(np.exp(logB1 - ls0))
        dll_dlogE1 += n0 * (-2.0 * A1 * (1.0 - w0) + 2.0 * w0)

    if n2 > 0:
        log_d2 = np.log(dat2 + eps)
        ls2 = _logsum_vec(log_d2, logB2)
        mean_ls2 = np.mean(ls2)
        mean_log_d2 = np.mean(log_d2)
        ll += n2 * (
            2.0 * A2 * (logB2 - mean_ls2)
            + mean_log_d2
            + np.log(2.0)
            + logA2
            - 2.0 * mean_ls2
            + ls_const2
        )
        w2 = np.mean(np.exp(logB2 - ls2))
        dll_dlogE2 += n2 * (-2.0 * A2 * (1.0 - w2) + 2.0 * w2)

    if n1 > 0:
        log_d1 = np.log(dat1 + eps)
        ls1B1 = _logsum_vec(log_d1, logB1)
        ls1B2 = _logsum_vec(log_d1, logB2)
        w1_1 = np.exp(logB1 - ls1B1)
        w1_2 = np.exp(logB2 - ls1B2)
        mean_ls1B1 = np.mean(ls1B1)
        mean_ls1B2 = np.mean(ls1B2)
        mean_log_d1 = np.mean(log_d1)

        sub1 = -2.0 * ls1B1 + logA1 + _logsum_scalar(logA1, 0.0)
        sub2 = -2.0 * ls1B2 + logA2 + _logsum_scalar(logA2, 0.0)
        sub3 = np.log(2.0) + logA1 + logA2 - ls1B1 - ls1B2
        lse12 = _logsum_pairwise(sub1, sub2)
        lse123 = _logsum_pairwise(lse12, sub3)
        mean_lse123 = np.mean(lse123)

        ll += n1 * (
            A1 * (logB1 - mean_ls1B1)
            + A2 * (logB2 - mean_ls1B2)
            + mean_log_d1
            + mean_lse123
        )

        p12 = np.exp(sub1 - lse12)
        p123_12 = np.exp(lse12 - lse123)
        p123_3 = np.exp(sub3 - lse123)
        d_mean_lse123_dlogE1 = np.mean(p123_12 * p12 * 2.0 * w1_1 + p123_3 * w1_1)
        dll_dlogE1 += n1 * (-A1 * (1.0 - np.mean(w1_1)) + d_mean_lse123_dlogE1)

        q12 = np.exp(sub2 - lse12)
        d_mean_lse123_dlogE2 = np.mean(p123_12 * q12 * 2.0 * w1_2 + p123_3 * w1_2)
        dll_dlogE2 += n1 * (-A2 * (1.0 - np.mean(w1_2)) + d_mean_lse123_dlogE2)

    if np.isnan(ll) or np.isinf(ll):
        grad = np.zeros(2, dtype=np.float64)
        return 1e15, grad

    grad = np.empty(2, dtype=np.float64)
    grad[0] = -dll_dlogE1
    grad[1] = -dll_dlogE2
    return -ll, grad


@dataclass
class SDSConfig:
    debug_mode: bool = False
    precision: int = 4
    e_grid_num_points: int = 50  # Retained for CLI compatibility; fine grid has been removed.
    e_grid_scale_factor: float = 20.0
    optim_num_iterations: int = 5  # Retained for CLI compatibility; optimization now uses fixed+random starts.
    optimizer_method: str = "bfgs"  # one of: nelder-mead, bfgs, lbfgsb
    optimizer_maxiter: int = 1000
    optimizer_gtol: float = 1e-5
    lbfgsb_bounds_scale: float = 1.0
    skip_boundary_missing_fraction: float = 0.10
    boundary_missing_mode: str = "skip"
    max_singletons_per_indv: int = 10000
    progress_every: int = 50000
    parquet_chunk_size: int = 50000  # [OPT] larger chunks reduce flush overhead


class PickleCache:
    def __init__(self, cache_dir: Optional[str] = None):
        root = cache_dir or os.environ.get("SDS_PICKLE_CACHE_DIR") or ".cache/sds_pickle"
        self.cache_dir = Path(root)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _cache_path(self, label: str, filepath: str, extra: str = "") -> Path:
        src = Path(filepath)
        stat = src.stat()
        digest = hashlib.sha256(
            f"{src.resolve()}|{stat.st_mtime_ns}|{stat.st_size}|{extra}".encode("utf-8")
        ).hexdigest()
        return self.cache_dir / f"{label}_{digest}.pkl"

    @staticmethod
    def _safe_load(cache_path: Path):
        with cache_path.open("rb") as handle:
            return pickle.load(handle)

    def load_or_create(self, label: str, filepath: str, builder: Callable[[], object], extra: str = ""):
        cache_path = self._cache_path(label, filepath, extra)
        lock_path = cache_path.with_suffix(f"{cache_path.suffix}.lock")

        try:
            if cache_path.exists():
                return self._safe_load(cache_path)
        except (pickle.UnpicklingError, EOFError, OSError, AttributeError) as exc:
            logger.warning("Discarding corrupt cache %s before rebuild: %s", cache_path, exc)
            cache_path.unlink(missing_ok=True)

        with lock_path.open("w") as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)

            try:
                if cache_path.exists():
                    try:
                        return self._safe_load(cache_path)
                    except (pickle.UnpicklingError, EOFError, OSError, AttributeError) as exc:
                        logger.warning("Discarding corrupt cache %s inside lock: %s", cache_path, exc)
                        cache_path.unlink(missing_ok=True)

                value = builder()
                tmp_path = cache_path.with_suffix(f"{cache_path.suffix}.tmp.{os.getpid()}")
                try:
                    with tmp_path.open("wb") as handle:
                        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
                        handle.flush()
                        os.fsync(handle.fileno())
                    os.replace(tmp_path, cache_path)
                finally:
                    if tmp_path.exists():
                        tmp_path.unlink(missing_ok=True)
                return value
            finally:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


class FileReader:
    @staticmethod
    def open_file(filepath: str, mode: str = 'r'):
        if filepath == '-':
            return sys.stdin if 'r' in mode else sys.stdout
        path = Path(filepath)
        if path.suffix == '.gz':
            return gzip.open(filepath, mode + 't')
        return open(filepath, mode)

    @staticmethod
    def _safe_float(v: str) -> float:
        if not v or v == 'NA':
            return np.nan
        try:
            return float(v)
        except ValueError:
            return np.nan

    @staticmethod
    def read_matrix(filepath: str, max_cols: Optional[int] = None) -> np.ndarray:
        data = []
        with FileReader.open_file(filepath) as f:
            for line in f:
                values = line.strip().split()
                if not values:
                    continue
                row = [FileReader._safe_float(v) for v in values]
                if max_cols:
                    row.extend([np.nan] * (max_cols - len(row)))
                    row = row[:max_cols]
                data.append(row)
        if not data:
            raise ValueError(f"No data found in file: {filepath}")
        return np.array(data)

    @staticmethod
    def read_vector(filepath: str) -> np.ndarray:
        with FileReader.open_file(filepath) as f:
            line = f.readline().strip()
            if not line:
                raise ValueError(f"Empty file: {filepath}")
            return np.array([FileReader._safe_float(v) for v in line.split()])


@_jit_or_noop(nopython=True, cache=True, parallel=True)
def _searchsorted_2d(finite, value):
    """[OPT] Vectorized searchsorted over rows using numba prange, replacing Python for loop."""
    n = finite.shape[0]
    out = np.empty(n, dtype=np.intp)
    for i in prange(n):
        out[i] = np.searchsorted(finite[i], value)
    return out


class SingletonData:
    def __init__(self, data: np.ndarray):
        self.data = data
        self.n_individuals = self.data.shape[0]
        # [OPT] Precompute finite array once, reused for every SNP's get_intervals
        self._finite = np.where(np.isnan(data), np.inf, data)

    @classmethod
    def from_file(cls, filepath: str, max_cols: int, cache: PickleCache):
        logger.info(f"Loading singletons from {filepath}")

        def _build():
            matrix = FileReader.read_matrix(filepath, max_cols)
            non_nan_mask = ~np.all(np.isnan(matrix), axis=0)
            return matrix[:, non_nan_mask]

        data = cache.load_or_create("singletons", filepath, _build, extra=str(max_cols))
        logger.info(f"Loaded {data.shape[0]} individuals")
        return cls(data)

    def get_intervals(self, test_position, boundary_upstream, boundary_downstream):
        data = self.data
        n, m = data.shape
        # [OPT] Use precomputed finite array
        idx = _searchsorted_2d(self._finite, test_position)

        upstream = np.full(n, np.nan)
        downstream = np.full(n, np.nan)

        has_prev = idx > 0
        prev_idx = np.where(has_prev, idx - 1, 0)
        prev_pos = data[np.arange(n), prev_idx]
        valid_up = has_prev & ~np.isnan(prev_pos) & (prev_pos >= boundary_upstream)
        upstream[valid_up] = test_position - prev_pos[valid_up]

        has_next = idx < m
        next_idx = np.where(has_next, idx, m - 1)
        next_pos = data[np.arange(n), next_idx]
        valid_dn = has_next & ~np.isnan(next_pos) & (next_pos <= boundary_downstream)
        downstream[valid_dn] = next_pos[valid_dn] - test_position
        return upstream, downstream


class GammaShapeInterpolator:
    def __init__(self, frequencies: np.ndarray, shapes: np.ndarray):
        self.frequencies = frequencies
        self.shapes = shapes

    @classmethod
    def from_file(cls, filepath: str, cache: PickleCache):
        logger.info(f"Loading gamma shape parameters from {filepath}")

        def _build():
            data = FileReader.read_matrix(filepath)
            sort_idx = np.argsort(data[:, 0])
            return data[sort_idx, 0], data[sort_idx, 1]

        frequencies, shapes = cache.load_or_create("gamma_shape", filepath, _build)
        return cls(frequencies, shapes)

    def get_shape(self, frequency: float) -> float:
        idx = np.searchsorted(self.frequencies, frequency)
        if idx == 0:
            return self.shapes[0]
        if idx >= len(self.frequencies):
            return self.shapes[-1]
        x1, x2 = self.frequencies[idx - 1], self.frequencies[idx]
        y1, y2 = self.shapes[idx - 1], self.shapes[idx]
        return y1 + (y2 - y1) * (frequency - x1) / (x2 - x1)


class SDSComputer:
    def __init__(self, singletons, observability, gamma_shape, boundaries, init_guess, config):
        self.singletons = singletons
        self.observability = observability / np.mean(observability)
        self.gamma_shape = gamma_shape
        self.boundaries = boundaries
        self.init_guess = init_guess
        self.config = config
        self.boundary_idx = 0
        self.eps = 1e-9

    def logsum(self, log_a, log_b):
        return np.maximum(log_a, log_b) + np.log1p(np.exp(-np.abs(log_a - log_b)))

    def compute_neg_log_likelihood(self, params, dat0, dat1, dat2, A1, A2, ls_const1, ls_const2):
        # Keep signature for compatibility with existing callsites.
        _ = ls_const1, ls_const2
        logE1, logE2 = params
        logA1, logA2 = np.log(A1), np.log(A2)
        logB1, logB2 = logA1 - logE1, logA2 - logE2
        n0, n1, n2 = len(dat0), len(dat1), len(dat2)
        ll = 0.0

        if n0 > 0:
            log_dat0 = np.log(dat0 + self.eps)
            ll += n0 * (
                2.0 * A1 * (logB1 - np.mean(self.logsum(log_dat0, logB1)))
                + np.mean(log_dat0)
                + np.log(2.0)
                + logA1
                - 2.0 * np.mean(self.logsum(log_dat0, logB1))
                + self.logsum(np.log(2.0) + logA1, 0.0)
            )
        if n2 > 0:
            log_dat2 = np.log(dat2 + self.eps)
            ll += n2 * (
                2.0 * A2 * (logB2 - np.mean(self.logsum(log_dat2, logB2)))
                + np.mean(log_dat2)
                + np.log(2.0)
                + logA2
                - 2.0 * np.mean(self.logsum(log_dat2, logB2))
                + self.logsum(np.log(2.0) + logA2, 0.0)
            )
        if n1 > 0:
            log_dat1 = np.log(dat1 + self.eps)
            sub1 = -2.0 * self.logsum(log_dat1, logB1) + logA1 + self.logsum(logA1, 0.0)
            sub2 = -2.0 * self.logsum(log_dat1, logB2) + logA2 + self.logsum(logA2, 0.0)
            sub3 = np.log(2.0) + logA1 + logA2 - self.logsum(log_dat1, logB1) - self.logsum(log_dat1, logB2)
            ll += n1 * (
                A1 * (logB1 - np.mean(self.logsum(log_dat1, logB1)))
                + A2 * (logB2 - np.mean(self.logsum(log_dat1, logB2)))
                + np.mean(log_dat1)
                + np.mean(self.logsum(self.logsum(sub1, sub2), sub3))
            )

        if np.isnan(ll) or np.isinf(ll):
            return 1e15
        return -ll

    def optimize_likelihood(self, dat0, dat1, dat2, A1, A2, ls_const1, ls_const2):
        log_center = np.log(self.init_guess)
        log_range = np.log(self.config.e_grid_scale_factor)

        # Diagonal parameterization keeps broad rSDS coverage with compact initialization cost.
        n_mu = 3
        n_delta = 34
        mu_grid = np.linspace(log_center - 1.0, log_center + 1.0, n_mu)
        delta_grid = np.linspace(-log_range, log_range, n_delta)
        mu_arr, delta_arr = np.meshgrid(mu_grid, delta_grid, indexing='ij')
        logE1_arr = mu_arr + delta_arr
        logE2_arr = mu_arr - delta_arr
        pm_coarse = np.ascontiguousarray(
            np.column_stack([logE1_arr.ravel(), logE2_arr.ravel()]),
            dtype=np.float64,
        )
        fixed_init = np.array([log_center, log_center], dtype=np.float64)
        random_init = pm_coarse[np.random.randint(pm_coarse.shape[0])]
        candidate_inits = (fixed_init, random_init)
        method = self.config.optimizer_method.lower()

        def _nll_and_grad(params):
            nll, grad = _nll_grad_kernel(
                np.ascontiguousarray(params, dtype=np.float64),
                dat0,
                dat1,
                dat2,
                A1,
                A2,
                self.eps,
                ls_const1,
                ls_const2,
            )
            return float(nll), np.ascontiguousarray(grad, dtype=np.float64)

        def _run_solver(init):
            if method == "nelder-mead":
                return minimize(
                    self.compute_neg_log_likelihood,
                    init,
                    args=(dat0, dat1, dat2, A1, A2, ls_const1, ls_const2),
                    method='Nelder-Mead',
                    options={'maxiter': int(self.config.optimizer_maxiter)},
                )
            if method == "bfgs":
                return minimize(
                    _nll_and_grad,
                    init,
                    method='BFGS',
                    jac=True,
                    options={'maxiter': int(self.config.optimizer_maxiter), 'gtol': float(self.config.optimizer_gtol)},
                )
            if method == "lbfgsb":
                bound_range = float(log_range) * float(self.config.lbfgsb_bounds_scale)
                bounds = [
                    (log_center - bound_range, log_center + bound_range),
                    (log_center - bound_range, log_center + bound_range),
                ]
                return minimize(
                    _nll_and_grad,
                    init,
                    method='L-BFGS-B',
                    jac=True,
                    bounds=bounds,
                    options={'maxiter': int(self.config.optimizer_maxiter), 'gtol': float(self.config.optimizer_gtol)},
                )
            raise ValueError(f"Unsupported optimizer method: {self.config.optimizer_method}")

        # Always keep a robust fallback that matches legacy behavior.
        def _run_nelder_fallback(init):
            return minimize(
                self.compute_neg_log_likelihood,
                init,
                args=(dat0, dat1, dat2, A1, A2, ls_const1, ls_const2),
                method='Nelder-Mead',
                options={'maxiter': int(self.config.optimizer_maxiter)},
            )

        best_res = None
        best_fun = np.inf
        for init in candidate_inits:
            res = _run_solver(init)
            if (not np.isfinite(res.fun)) or (not np.all(np.isfinite(res.x))):
                res = _run_nelder_fallback(init)
            if np.isfinite(res.fun) and res.fun < best_fun:
                best_fun = float(res.fun)
                best_res = res

        # Fallback: if both starts fail, use center start.
        if best_res is None:
            return fixed_init[0], fixed_init[1], np.exp(np.mean(fixed_init))
        logE1, logE2 = best_res.x
        return logE1, logE2, np.exp(np.mean(best_res.x))

    def process_test_snp(self, snp_id, allele_anc, allele_der, position, genotypes):
        seed_val = int(position) % (2**32 - 1)
        np.random.seed(seed_val)

        while self.boundary_idx < len(self.boundaries) and self.boundaries[self.boundary_idx, 1] < position:
            self.boundary_idx += 1

        if self.boundary_idx >= len(self.boundaries) or self.boundaries[self.boundary_idx, 0] > position:
            return None

        valid_mask = genotypes != -1
        if not np.any(valid_mask):
            return None

        daf = np.mean(genotypes[valid_mask]) / 2.0
        if daf <= 0.0 or daf >= 1.0:
            return None

        boundary_up, boundary_down = self.boundaries[self.boundary_idx]
        upstream, downstream = self.singletons.get_intervals(position, boundary_up, boundary_down)
        up_nan_frac = np.isnan(upstream).mean()
        dn_nan_frac = np.isnan(downstream).mean()

        if self.config.boundary_missing_mode == "skip":
            if (
                up_nan_frac > self.config.skip_boundary_missing_fraction
                or dn_nan_frac > self.config.skip_boundary_missing_fraction
            ):
                return None
        elif self.config.boundary_missing_mode != "cap_to_boundary":
            raise ValueError(f"Unsupported boundary_missing_mode: {self.config.boundary_missing_mode}")

        # When a side has no observed singleton at all, fall back to the arm boundary
        # distance instead of propagating all-NaN intervals into the likelihood step.
        up_boundary_cap = max(float(position - boundary_up), self.eps)
        dn_boundary_cap = max(float(boundary_down - position), self.eps)
        if self.config.boundary_missing_mode == "cap_to_boundary":
            upstream[np.isnan(upstream)] = up_boundary_cap
            downstream[np.isnan(downstream)] = dn_boundary_cap
        else:
            up_fill = up_boundary_cap if np.all(np.isnan(upstream)) else np.nanmax(upstream)
            dn_fill = dn_boundary_cap if np.all(np.isnan(downstream)) else np.nanmax(downstream)
            upstream[np.isnan(upstream)] = up_fill
            downstream[np.isnan(downstream)] = dn_fill
        intervals = (upstream + downstream) * self.observability

        dat0 = intervals[genotypes == 0]
        dat1 = intervals[genotypes == 1]
        dat2 = intervals[genotypes == 2]

        if len(dat0) + len(dat1) + len(dat2) < 10:
            return None

        A1 = self.gamma_shape.get_shape(1 - daf)
        A2 = self.gamma_shape.get_shape(daf)
        ls_const1 = float(_logsum_scalar(np.log(2.0) + np.log(A1), 0.0))
        ls_const2 = float(_logsum_scalar(np.log(2.0) + np.log(A2), 0.0))
        logE1, logE2, suggested_init = self.optimize_likelihood(
            dat0,
            dat1,
            dat2,
            A1,
            A2,
            ls_const1,
            ls_const2,
        )

        return {
            'ID': snp_id,
            'AA': allele_anc,
            'DA': allele_der,
            'POS': int(position),
            'DAF': round(daf, self.config.precision),
            'nG0': len(dat0),
            'nG1': len(dat1),
            'nG2': len(dat2),
            'rSDS': round(logE1 - logE2, self.config.precision),
            'SuggestedInitPoint': f"1e{int(np.log10(suggested_init))}",
        }


def _iter_text_t_file(filepath: str) -> Iterator[tuple[str, str, str, float, np.ndarray]]:
    """[OPT] Batch-read lines and use numpy for fast genotype parsing."""
    with FileReader.open_file(filepath) as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(None, 4)  # split only first 4 whitespace
            if len(parts) < 5:
                continue
            geno_str = parts[4]
            # [OPT] Use numpy fromstring for bulk int conversion, handle NA separately
            has_na = 'N' in geno_str  # quick check for NA
            if has_na:
                tokens = geno_str.split()
                # Use numpy array with bytestring comparison for speed
                arr = np.empty(len(tokens), dtype=np.int64)
                for j, t in enumerate(tokens):
                    arr[j] = -1 if t == "NA" else int(t)
                genotypes = arr
            else:
                # Fast path: no NA, pure numeric - use np.fromstring
                genotypes = np.fromstring(geno_str, dtype=np.int64, sep=' ')
            yield parts[0], parts[1], parts[2], float(parts[3]), genotypes


def _iter_feather_t_file(filepath: str) -> Iterator[tuple[str, str, str, float, np.ndarray]]:
    table = feather.read_table(filepath)
    for batch in table.to_batches():
        id_col = batch.column(0)
        aa_col = batch.column(1)
        da_col = batch.column(2)
        pos_col = batch.column(3)
        genotype_col = batch.column(4)

        for row_index in range(batch.num_rows):
            genotype_str = genotype_col[row_index].as_py() or ""
            # [OPT] Same fast-path parsing as text file
            if not genotype_str:
                genotypes = np.array([], dtype=np.int64)
            elif 'N' in genotype_str:
                genotype_tokens = genotype_str.split("\t")
                arr = np.empty(len(genotype_tokens), dtype=np.int64)
                for j, t in enumerate(genotype_tokens):
                    arr[j] = -1 if t == "NA" else int(t)
                genotypes = arr
            else:
                genotypes = np.fromstring(genotype_str.replace('\t', ' '), dtype=np.int64, sep=' ')
            yield (
                id_col[row_index].as_py(),
                aa_col[row_index].as_py(),
                da_col[row_index].as_py(),
                float(pos_col[row_index].as_py()),
                genotypes,
            )


def iter_test_snps(filepath: str) -> Iterator[tuple[str, str, str, float, np.ndarray]]:
    suffix = Path(filepath).suffix.lower()
    if suffix == ".feather":
        return _iter_feather_t_file(filepath)
    return _iter_text_t_file(filepath)


def write_tsv(results: list[dict], output_path: str):
    with open(output_path, "w") as handle:
        print("\t".join(RESULT_COLUMNS), file=handle)
        for row in results:
            print("\t".join(str(row[col]) for col in RESULT_COLUMNS), file=handle)


def write_parquet(results: list[dict], output_path: str):
    schema = pa.schema(
        [
            ("ID", pa.string()),
            ("AA", pa.string()),
            ("DA", pa.string()),
            ("POS", pa.int64()),
            ("DAF", pa.float64()),
            ("nG0", pa.int64()),
            ("nG1", pa.int64()),
            ("nG2", pa.int64()),
            ("rSDS", pa.float64()),
            ("SuggestedInitPoint", pa.string()),
        ]
    )
    pq.write_table(pa.Table.from_pylist(results, schema=schema), output_path)


def write_summary_csv(summary_path: str, rows_in: int, rows_out: int, t_file: str, parquet_path: Optional[str]):
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["t_file", "rows_in", "rows_out", "parquet_output"])
        writer.writeheader()
        writer.writerow(
            {
                "t_file": t_file,
                "rows_in": rows_in,
                "rows_out": rows_out,
                "parquet_output": parquet_path or "",
            }
        )


class StreamingTSVWriter:
    """[OPT] Buffered TSV writer - flushes every N rows to reduce I/O syscalls."""
    _FLUSH_EVERY = 500

    def __init__(self, output_path: str):
        self.handle = open(output_path, "w", buffering=1 << 16)  # 64KB buffer
        self.handle.write("\t".join(RESULT_COLUMNS) + "\n")
        self._buf: list[str] = []

    def write_row(self, row: dict):
        self._buf.append("\t".join(str(row[col]) for col in RESULT_COLUMNS))
        if len(self._buf) >= self._FLUSH_EVERY:
            self._flush_buf()

    def _flush_buf(self):
        if self._buf:
            self.handle.write("\n".join(self._buf) + "\n")
            self._buf.clear()

    def close(self):
        self._flush_buf()
        if not self.handle.closed:
            self.handle.close()


class StreamingParquetWriter:
    schema = pa.schema(
        [
            ("ID", pa.string()),
            ("AA", pa.string()),
            ("DA", pa.string()),
            ("POS", pa.int64()),
            ("DAF", pa.float64()),
            ("nG0", pa.int64()),
            ("nG1", pa.int64()),
            ("nG2", pa.int64()),
            ("rSDS", pa.float64()),
            ("SuggestedInitPoint", pa.string()),
        ]
    )

    def __init__(self, output_path: str, chunk_size: int):
        self.output_path = output_path
        self.chunk_size = chunk_size
        self.buffer: list[dict] = []
        self.writer: Optional[pq.ParquetWriter] = None

    def write_row(self, row: dict):
        self.buffer.append(row)
        if len(self.buffer) >= self.chunk_size:
            self.flush()

    def flush(self):
        if not self.buffer:
            return
        table = pa.Table.from_pylist(self.buffer, schema=self.schema)
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)
        self.writer.write_table(table)
        self.buffer.clear()

    def close(self):
        self.flush()
        if self.writer is None:
            pq.write_table(pa.Table.from_pylist([], schema=self.schema), self.output_path)
            return
        self.writer.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('s_file')
    parser.add_argument('t_file')
    parser.add_argument('o_file')
    parser.add_argument('b_file')
    parser.add_argument('g_file')
    parser.add_argument('init', type=float)
    parser.add_argument('s_file_ncol', nargs='?', type=int, default=10000)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--output', '-o', required=True)
    parser.add_argument('--output-parquet')
    parser.add_argument('--summary-csv')
    parser.add_argument('--pickle-cache-dir')
    parser.add_argument('--progress-every', type=int, default=50000)
    parser.add_argument('--parquet-chunk-size', type=int, default=10000)
    parser.add_argument('--optimizer-method', choices=['nelder-mead', 'bfgs', 'lbfgsb'], default='bfgs')
    parser.add_argument('--optimizer-maxiter', type=int, default=1000)
    parser.add_argument('--optimizer-gtol', type=float, default=1e-5)
    parser.add_argument('--lbfgsb-bounds-scale', type=float, default=1.0)
    parser.add_argument('--skip-boundary-missing-fraction', type=float, default=0.10)
    parser.add_argument('--boundary-missing-mode', choices=['skip', 'cap_to_boundary'], default='skip')
    args = parser.parse_args()

    config = SDSConfig(
        debug_mode=args.debug,
        max_singletons_per_indv=args.s_file_ncol,
        optimizer_method=args.optimizer_method,
        optimizer_maxiter=max(1, args.optimizer_maxiter),
        optimizer_gtol=max(1e-12, args.optimizer_gtol),
        lbfgsb_bounds_scale=max(0.01, args.lbfgsb_bounds_scale),
        skip_boundary_missing_fraction=min(max(0.0, args.skip_boundary_missing_fraction), 1.0),
        boundary_missing_mode=args.boundary_missing_mode,
        progress_every=max(1, args.progress_every),
        parquet_chunk_size=max(1, args.parquet_chunk_size),
    )
    cache = PickleCache(args.pickle_cache_dir)

    try:
        if numba is None:
            logger.warning(
                "Numba is not available; using pure-Python fallback kernels. "
                "Performance will be significantly lower."
            )
        singletons = SingletonData.from_file(args.s_file, config.max_singletons_per_indv, cache)
        observability = cache.load_or_create("observability", args.o_file, lambda: FileReader.read_vector(args.o_file))
        gamma_shape = GammaShapeInterpolator.from_file(args.g_file, cache)
        boundaries = cache.load_or_create("boundaries", args.b_file, lambda: FileReader.read_matrix(args.b_file))
        computer = SDSComputer(singletons, observability, gamma_shape, boundaries, args.init, config)
        dummy = np.array([1.0, 2.0, 3.0])
        _nll_kernel(np.array([0.0, 0.0]), dummy, dummy, dummy, 1.0, 1.0, 1e-9, np.log(3.0), np.log(3.0))
        _nll_kernel_batch(
            np.array([[0.0, 0.0], [0.1, 0.1]]),
            dummy,
            dummy,
            dummy,
            1.0,
            1.0,
            1e-9,
            np.log(3.0),
            np.log(3.0),
        )
        if config.optimizer_method in {"bfgs", "lbfgsb"}:
            _nll_grad_kernel(
                np.array([0.0, 0.0]),
                dummy,
                dummy,
                dummy,
                1.0,
                1.0,
                1e-9,
                np.log(3.0),
                np.log(3.0),
            )
        logger.info("Numba JIT warm-up complete (method=%s)", config.optimizer_method)

        rows_in = 0
        rows_out = 0
        start_time = time.time()
        last_log_time = start_time
        tsv_writer = StreamingTSVWriter(args.output)
        parquet_writer = None
        if args.output_parquet:
            parquet_writer = StreamingParquetWriter(args.output_parquet, config.parquet_chunk_size)

        for snp_id, allele_anc, allele_der, position, genotypes in iter_test_snps(args.t_file):
            rows_in += 1
            res = computer.process_test_snp(snp_id, allele_anc, allele_der, position, genotypes)
            if res:
                rows_out += 1
                tsv_writer.write_row(res)
                if parquet_writer is not None:
                    parquet_writer.write_row(res)

            if rows_in % config.progress_every == 0:
                now = time.time()
                elapsed = max(now - start_time, 1e-9)
                delta = max(now - last_log_time, 1e-9)
                logger.info(
                    "Progress: processed=%d kept=%d elapsed=%.1fs avg_rate=%.1f rows/s window_rate=%.1f rows/s",
                    rows_in,
                    rows_out,
                    elapsed,
                    rows_in / elapsed,
                    config.progress_every / delta,
                )
                last_log_time = now

        tsv_writer.close()
        if parquet_writer is not None:
            parquet_writer.close()

        elapsed = time.time() - start_time
        logger.info(
            "Finished: processed=%d kept=%d elapsed=%.1fs avg_rate=%.1f rows/s output=%s",
            rows_in,
            rows_out,
            elapsed,
            rows_in / max(elapsed, 1e-9),
            args.output,
        )
        if args.summary_csv:
            write_summary_csv(args.summary_csv, rows_in, rows_out, args.t_file, args.output_parquet)
    except Exception as exc:
        logger.error(f"Error: {exc}", exc_info=True)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
