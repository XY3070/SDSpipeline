#!/usr/bin/env python3
"""
Singleton Density Score (SDS) Computation
Fixed: 'NA' parsing error in FileReader
"""

import sys
import gzip
import argparse
from pathlib import Path
from typing import Optional, Tuple, List
import numpy as np
from scipy.optimize import minimize
from dataclasses import dataclass
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class SDSConfig:
    debug_mode: bool = False
    precision: int = 4
    e_grid_num_points: int = 50
    e_grid_scale_factor: float = 20.0
    optim_num_iterations: int = 5
    skip_boundary_missing_fraction: float = 0.10
    max_singletons_per_indv: int = 10000

class FileReader:
    @staticmethod
    def open_file(filepath: str, mode: str = 'r'):
        if filepath == '-': return sys.stdin if 'r' in mode else sys.stdout
        path = Path(filepath)
        if path.suffix == '.gz': return gzip.open(filepath, mode + 't')
        return open(filepath, mode)
    
    @staticmethod
    def _safe_float(v: str) -> float:
        """Helper to convert string to float, handling 'NA' as NaN"""
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
                if not values: continue
                # FIX: Handle 'NA' string explicitly
                row = [FileReader._safe_float(v) for v in values]
                if max_cols:
                    row.extend([np.nan] * (max_cols - len(row)))
                    row = row[:max_cols]
                data.append(row)
        if not data: raise ValueError(f"No data found in file: {filepath}")
        return np.array(data)
    
    @staticmethod
    def read_vector(filepath: str) -> np.ndarray:
        with FileReader.open_file(filepath) as f:
            line = f.readline().strip()
            if not line: raise ValueError(f"Empty file: {filepath}")
            # FIX: Handle 'NA' string explicitly
            return np.array([FileReader._safe_float(v) for v in line.split()])

class SingletonData:
    def __init__(self, filepath: str, max_cols: int):
        logger.info(f"Loading singletons from {filepath}")
        self.data = FileReader.read_matrix(filepath, max_cols)
        self.n_individuals = self.data.shape[0]
        self.current_indices = np.zeros(self.n_individuals, dtype=int)
        
        non_nan_mask = ~np.all(np.isnan(self.data), axis=0)
        self.data = self.data[:, non_nan_mask]
        logger.info(f"Loaded {self.n_individuals} individuals")
    
    def get_intervals(self, test_position, boundary_upstream, boundary_downstream):
        upstream = np.full(self.n_individuals, np.nan)
        downstream = np.full(self.n_individuals, np.nan)
        n_cols = self.data.shape[1]
        
        for i in range(self.n_individuals):
            # If inputs are not sorted/reset logic
            if self.current_indices[i] < n_cols and not np.isnan(self.data[i, self.current_indices[i]]) and self.data[i, self.current_indices[i]] > test_position:
                 self.current_indices[i] = 0

            while (self.current_indices[i] < n_cols and
                   not np.isnan(self.data[i, self.current_indices[i]]) and
                   self.data[i, self.current_indices[i]] < test_position):
                self.current_indices[i] += 1
            
            curr_idx = self.current_indices[i]
            if curr_idx > 0:
                prev_pos = self.data[i, curr_idx - 1]
                if not np.isnan(prev_pos) and prev_pos >= boundary_upstream:
                    upstream[i] = test_position - prev_pos
            if curr_idx < n_cols:
                next_pos = self.data[i, curr_idx]
                if not np.isnan(next_pos) and next_pos <= boundary_downstream:
                    downstream[i] = next_pos - test_position
        return upstream, downstream

class GammaShapeInterpolator:
    def __init__(self, filepath: str):
        logger.info(f"Loading gamma shape parameters from {filepath}")
        data = FileReader.read_matrix(filepath)
        self.frequencies = data[:, 0]
        self.shapes = data[:, 1]
        sort_idx = np.argsort(self.frequencies)
        self.frequencies, self.shapes = self.frequencies[sort_idx], self.shapes[sort_idx]
    
    def get_shape(self, frequency: float) -> float:
        idx = np.searchsorted(self.frequencies, frequency)
        if idx == 0: return self.shapes[0]
        if idx >= len(self.frequencies): return self.shapes[-1]
        x1, x2 = self.frequencies[idx-1], self.frequencies[idx]
        y1, y2 = self.shapes[idx-1], self.shapes[idx]
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
    
    def compute_neg_log_likelihood(self, params, dat0, dat1, dat2, A1, A2):
        logE1, logE2 = params
        logA1, logA2 = np.log(A1), np.log(A2)
        logB1, logB2 = logA1 - logE1, logA2 - logE2
        n0, n1, n2 = len(dat0), len(dat1), len(dat2)
        LL = 0.0
        
        if n0 > 0:
            log_dat0 = np.log(dat0 + self.eps)
            LL += n0 * (2.0*A1*(logB1 - np.mean(self.logsum(log_dat0, logB1))) + np.mean(log_dat0) + np.log(2)+logA1 - 2.0*np.mean(self.logsum(log_dat0, logB1)) + np.mean(self.logsum(np.log(2)+logA1, 0)))
        if n2 > 0:
            log_dat2 = np.log(dat2 + self.eps)
            LL += n2 * (2.0*A2*(logB2 - np.mean(self.logsum(log_dat2, logB2))) + np.mean(log_dat2) + np.log(2)+logA2 - 2.0*np.mean(self.logsum(log_dat2, logB2)) + np.mean(self.logsum(np.log(2)+logA2, 0)))
        if n1 > 0:
            log_dat1 = np.log(dat1 + self.eps)
            sub1 = -2.0*self.logsum(log_dat1, logB1) + logA1 + self.logsum(logA1, 0)
            sub2 = -2.0*self.logsum(log_dat1, logB2) + logA2 + self.logsum(logA2, 0)
            sub3 = (np.log(2)+logA1+logA2 - self.logsum(log_dat1, logB1) - self.logsum(log_dat1, logB2))
            LL += n1 * (A1*(logB1 - np.mean(self.logsum(log_dat1, logB1))) + A2*(logB2 - np.mean(self.logsum(log_dat1, logB2))) + np.mean(log_dat1) + np.mean(self.logsum(self.logsum(sub1, sub2), sub3)))
        
        if np.isnan(LL) or np.isinf(LL): return 1e15
        return -LL
    
    def optimize_likelihood(self, dat0, dat1, dat2, A1, A2):
        log_center = np.log(self.init_guess)
        log_range = np.log(self.config.e_grid_scale_factor)
        best_ll, best_params = -np.inf, None
        
        for _ in range(self.config.optim_num_iterations):
            init = np.random.uniform(log_center - log_range, log_center + log_range, 2)
            res = minimize(self.compute_neg_log_likelihood, init, args=(dat0, dat1, dat2, A1, A2), method='Nelder-Mead', options={'maxiter': 1000})
            if -res.fun > best_ll: best_ll, best_params = -res.fun, res.x
            
        res = minimize(self.compute_neg_log_likelihood, np.array([log_center]*2), args=(dat0, dat1, dat2, A1, A2), method='Nelder-Mead', options={'maxiter': 1000})
        if -res.fun > best_ll: best_params = res.x
        return best_params[0], best_params[1], np.exp(np.mean(best_params))
    
    def process_test_snp(self, snp_id, allele_anc, allele_der, position, genotypes):
        seed_val = int(position) % (2**32 - 1)
        np.random.seed(seed_val)

        while (self.boundary_idx < len(self.boundaries) and self.boundaries[self.boundary_idx, 1] < position):
            self.boundary_idx += 1
        
        if (self.boundary_idx >= len(self.boundaries) or self.boundaries[self.boundary_idx, 0] > position):
            return None
        
        valid_mask = genotypes != -1
        if not np.any(valid_mask): return None
            
        daf = np.mean(genotypes[valid_mask]) / 2.0
        if daf <= 0.0 or daf >= 1.0: return None

        boundary_up, boundary_down = self.boundaries[self.boundary_idx]
        upstream, downstream = self.singletons.get_intervals(position, boundary_up, boundary_down)
        
        if (np.isnan(upstream).mean() > self.config.skip_boundary_missing_fraction or 
            np.isnan(downstream).mean() > self.config.skip_boundary_missing_fraction):
            return None
        
        upstream[np.isnan(upstream)] = np.nanmax(upstream)
        downstream[np.isnan(downstream)] = np.nanmax(downstream)
        intervals = (upstream + downstream) * self.observability
        
        dat0 = intervals[genotypes == 0]
        dat1 = intervals[genotypes == 1]
        dat2 = intervals[genotypes == 2]
        
        if len(dat0) + len(dat1) + len(dat2) < 10: return None

        A1 = self.gamma_shape.get_shape(1 - daf)
        A2 = self.gamma_shape.get_shape(daf)
        logE1, logE2, suggested_init = self.optimize_likelihood(dat0, dat1, dat2, A1, A2)
        
        return {
            'ID': snp_id, 'AA': allele_anc, 'DA': allele_der, 'POS': int(position),
            'DAF': round(daf, self.config.precision), 'nG0': len(dat0), 'nG1': len(dat1), 'nG2': len(dat2),
            'rSDS': round(logE1 - logE2, self.config.precision),
            'SuggestedInitPoint': f"1e{int(np.log10(suggested_init))}"
        }

def main():
    args = argparse.ArgumentParser()
    args.add_argument('s_file'); args.add_argument('t_file'); args.add_argument('o_file')
    args.add_argument('b_file'); args.add_argument('g_file'); args.add_argument('init', type=float)
    args.add_argument('s_file_ncol', nargs='?', type=int, default=10000)
    args.add_argument('--debug', action='store_true'); args.add_argument('--output', '-o')
    args = args.parse_args()
    
    config = SDSConfig(debug_mode=args.debug, max_singletons_per_indv=args.s_file_ncol)
    try:
        singletons = SingletonData(args.s_file, config.max_singletons_per_indv)
        observability = FileReader.read_vector(args.o_file)
        gamma_shape = GammaShapeInterpolator(args.g_file)
        boundaries = FileReader.read_matrix(args.b_file)
        computer = SDSComputer(singletons, observability, gamma_shape, boundaries, args.init, config)
        
        output = open(args.output, 'w') if args.output else sys.stdout
        print('\t'.join(['ID', 'AA', 'DA', 'POS', 'DAF', 'nG0', 'nG1', 'nG2', 'rSDS', 'SuggestedInitPoint']), file=output)
        
        with FileReader.open_file(args.t_file) as f:
            for line in f:
                if not line.strip(): continue
                parts = line.split()
                if len(parts) < 4: continue
                genotypes = []
                for g in parts[4:]:
                    try: genotypes.append(int(g))
                    except: genotypes.append(-1)
                res = computer.process_test_snp(parts[0], parts[1], parts[2], float(parts[3]), np.array(genotypes))
                if res:
                    print('\t'.join(str(res[k]) for k in ['ID', 'AA', 'DA', 'POS', 'DAF', 'nG0', 'nG1', 'nG2', 'rSDS', 'SuggestedInitPoint']), file=output)

        if output != sys.stdout: output.close()
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True); return 1
    return 0

if __name__ == '__main__': sys.exit(main())