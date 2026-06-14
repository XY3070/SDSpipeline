#!/data/home/grp-wangyf/intern/miniforge3/envs/ms/bin/python

import os, sys, re, math
import argparse
import numpy as np
import csv

import simuOpt
simuOpt.setOptions(quiet=True)

import simuPOP as sim
import math
from simuPOP.utils import *


#-----------------#
# Parse arguments #
#-----------------#

parser = argparse.ArgumentParser(description="Backward simulation of selection using the simuPOP package. The input includes the demographic model, the current derived allele frequency, the selection coefficient, and the backward time window at which selection was active.")

parser.add_argument("-g", "--number_of_generations", 
                    type=int,
                    default=10000,
                    help="The number of generations simulated \
                          (default = %(default)s).")

parser.add_argument("-te", "--selection_end_time", 
                    type=int,
                    default=0,
                    help="The generation (backward in time) at which selection ended \
                          (default = %(default)s, i.e., current time).")

parser.add_argument("-ts", "--selection_start_time", 
                    type=int,
                    default=100,
                    help="The generation (backward in time) at which selection started \
                          (default = %(default)s).")

parser.add_argument("-ps", "--selection_sigma_paramter",
                    type=float,
                    default=0.0,
                    help="The selection \"sigma\" paramter \
                          (default = %(default)s).")

parser.add_argument("-ph", "--selection_dominance_paramter",
                    type=float,
                    default=0.5,
                    help="The selection h paramter \
                          (default = %(default)s).")

parser.add_argument("-f", "--current_snp_frequency",
                    type=float,
                    default=0.6,
                    help="The frequency of the snp at current time \
                          (default = %(default)s).")

parser.add_argument("-sNe", "--scale_Ne", 
                    type=float,
                    default=1.0,
                    help="A scale factor for Ne, the effective population_size, \
                          (default = %(default)s).")

parser.add_argument("-m", "--population_size_model",
                    choices=["Tennessen_CEU","Gravel_AFR","Gravel_CEU","Gravel_CHB","Const_1E4", "NCN", "SCN", "NCN_SMC", "SCN_SMC"],
                    default="Tennessen_CEU",
                    help="Define the population size history model. \
                          (default = %(default)s).")

parser.add_argument("--npz_path", type=str, default="china_pop_history.npz",
                    help="Path to the pre-extracted China population history NPZ file")

parser.add_argument("--ne_stat", choices=["median", "lower", "upper"], default="median",
                    help="Which Ne statistic to use (NPZ path only)")

parser.add_argument("--smcpp_csv", type=str, default=None,
                    help="Path to SMC++ averaged_Ne.csv (for NCN_SMC/SCN_SMC models)")

parser.add_argument("--ne0_override", type=float, default=None,
                    help="Override Ne(0) to this value (e.g. 100000). "
                         "Scales entire SMC++ trajectory proportionally.")

parser.add_argument("--smcpp_stat", type=str, default="mean_Ne",
                    help="Column in SMC++ CSV to use (mean_Ne, shuffle_00, ...)")

args = parser.parse_args()

#-------------------------------
#
#-------------------------------

OUTPUT = sys.stdout

num_of_generations = args.number_of_generations-1

#

def genBackwards(gen_forward):
    return int(abs(gen_forward-(num_of_generations+1)))

#

def Nt(gen):
    if gen==0:
        return 1
    tmp_baseNt = baseNt(gen)
    return int(args.scale_Ne * tmp_baseNt)


POP_DATA = None
SMC_DATA = None


def load_npz_models(path):
    global POP_DATA
    try:
        POP_DATA = np.load(path)
    except Exception as e:
        sys.exit(f"Failed to load NPZ file: {e}")


def load_smcpp_csv(path, stat):
    global SMC_DATA
    try:
        with open(path) as fh:
            # Sniff delimiter: tab or comma
            sample = fh.read(8192)
            fh.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
            reader = csv.DictReader(fh, dialect=dialect)
            rows = list(reader)
    except Exception as e:
        sys.exit(f"Failed to load SMC++ CSV: {e}")
    if not rows:
        sys.exit(f"SMC++ CSV is empty: {path}")
    header = list(rows[0].keys())
    # Auto-detect format: raw smc++ plot CSV (x,y) vs resampled CSV (generation, mean_Ne, ...)
    if "x" in header and "y" in header:
        t = np.array([float(r["x"]) for r in rows])
        ne = np.array([float(r["y"]) for r in rows])
    elif "generation" in header and stat in header:
        t = np.array([float(r["generation"]) for r in rows])
        ne = np.array([float(r[stat]) for r in rows])
    elif "generation" in header:
        sys.exit(f"SMC++ CSV has 'generation' but missing '{stat}' column. Available: {header}")
    else:
        sys.exit(f"SMC++ CSV must have 'x'+'y' or 'generation'+'{stat}' columns. Found: {header}")
    order = np.argsort(t, kind="stable")
    SMC_DATA = (t[order], ne[order])


def baseNt(gen):
    model = args.population_size_model
    if model in ["NCN", "SCN"]:
        return baseNt_China_Custom(gen, model)
    elif model in ["NCN_SMC", "SCN_SMC"]:
        return baseNt_SMC(gen)
    elif model == "Tennessen_CEU":
        return baseNt_Tennessen_CEU(gen)
    elif model == "Gravel_AFR":
        return baseNt_Gravel_AFR(gen)
    elif model == "Gravel_CEU":
        return baseNt_Gravel_CEU(gen)
    elif model == "Gravel_CHB":
        return baseNt_Gravel_CHB(gen)
    elif model == "Const_1E4":
        return baseNt_Const_1E4(gen)


def baseNt_China_Custom(gen, pop_name):
    my_gen_backwards = genBackwards(gen)

    # get the array of cooresponding pop and stats
    t_array = POP_DATA[f'{pop_name}_t']
    ne_array = POP_DATA[f'{pop_name}_{args.ne_stat}']

    # efficient search: because t_array is ordered, np.searchsorted is O(logN)
    idx = np.searchsorted(t_array, my_gen_backwards, side='right') - 1
    idx = max(0, min(idx, len(ne_array) - 1))

    return int(ne_array[idx])


def baseNt_SMC(gen):
    my_gen_backwards = genBackwards(gen)
    t_grid, ne_grid = SMC_DATA
    if args.ne0_override is not None:
        scale = args.ne0_override / ne_grid[0]
        ne_grid = ne_grid * scale
    # piecewise-constant: SMC++ step function
    idx = np.searchsorted(t_grid, my_gen_backwards, side='right') - 1
    idx = max(0, min(idx, len(ne_grid) - 1))
    return int(ne_grid[idx])


def baseNt_Tennessen_CEU(gen):
    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards > 2040:
        return 14474
    if my_gen_backwards > 920:
        return 1861
    if my_gen_backwards == 920:
        return 1032
    if my_gen_backwards >= 205:
        return int(1032*(1.0*9300/1032)**(1.0*(920-my_gen_backwards)/(920-205)))
    if my_gen_backwards >= 1:
        return int(9300*(1.0*512000/9300)**(1.0*(205-my_gen_backwards)/(205-1)))


def baseNt_Gravel_AFR(gen):
    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards >= 1:
        return 14474


def baseNt_Gravel_CEU(gen):
    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards > 2040:
        return 14474
    if my_gen_backwards > 920:
        return 1861
    if my_gen_backwards == 920:
        return 1032
    if my_gen_backwards >= 1:
        return int(1032*(1.0038)**(1.0*(920-my_gen_backwards)))


def baseNt_Gravel_CHB(gen):
    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards > 5920:
        return 7310
    if my_gen_backwards > 2040:
        return 14474
    if my_gen_backwards > 920:
        return 1861
    if my_gen_backwards == 920:
        return 550
    if my_gen_backwards >= 1:
        return int(550*(1.0048)**(1.0*(920-my_gen_backwards)))


def baseNt_Const_1E4(gen):
    return 10000


#

def Ft(gen, subPop):
    my_gen_backwards = genBackwards(gen)

    if my_gen_backwards < args.selection_end_time:
        return (1.0, 1.0, 1.0)
    if my_gen_backwards > args.selection_start_time:
        return (1.0, 1.0, 1.0)
    else:
        #older version had a bug... fixed: Sep 4, 2015
        #return (1.0, 1.0+(2.0*args.selection_sigma_paramter*args.selection_dominance_paramter), 1.0+(2.0*args.selection_sigma_paramter))
        return (1.0-args.selection_sigma_paramter, 1.0-(args.selection_sigma_paramter*args.selection_dominance_paramter), 1.0)


if args.population_size_model in ["NCN", "SCN"]:
    load_npz_models(args.npz_path)
elif args.population_size_model in ["NCN_SMC", "SCN_SMC"]:
    if args.smcpp_csv is None:
        sys.exit("--smcpp_csv is required for NCN_SMC/SCN_SMC models")
    load_smcpp_csv(args.smcpp_csv, args.smcpp_stat)

traj = simulateBackwardTrajectory(N=Nt, fitness=Ft, endGen=num_of_generations, endFreq=args.current_snp_frequency)
#traj = simulateBackwardTrajectory(N=Nt, fitness=[1,1,1], endGen=num_of_generations, endFreq=0.6)

#for i in range(0,1+num_of_generations):
for i in range(num_of_generations,0,-1):
    my_gen = str(genBackwards(i))
    my_popsize = str(Nt(i))
    my_fitness = Ft(i,0)
    my_fitness0 = str(my_fitness[0])
    my_fitness1 = str(my_fitness[1])
    my_fitness2 = str(my_fitness[2])
    my_freq = str(traj.freq(i,0)[0])
    OUTPUT.write( my_gen + "\t" + my_popsize + "\t" + my_fitness0 + "\t" + my_fitness1 + "\t" + my_fitness2 + "\t" + my_freq  +"\n")
