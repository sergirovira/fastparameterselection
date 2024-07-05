import math
import csv
import sys

sys.path.append('../lattice-estimator')


try:
    from estimator import *
except ImportError:
    print("Warning: Failed to import lattice_estimator, some options will not work")

#Auxiliary functions needed in estimate.py


def load_all_from_csv(file_path):
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        entries = [row for row in reader]
    return entries

def closest_power_of_2(n):
    if n <= 0:
        raise ValueError("Input must be a positive number.")
    
    # Calculate the power of 2 just below and above the number
    lower_pow = 2 ** math.floor(math.log2(n))
    upper_pow = 2 ** math.ceil(math.log2(n))
    
    # Determine which is closer
    if abs(n - lower_pow) < abs(n - upper_pow):
        return lower_pow
    else:
        return upper_pow

def print_table(headers, rows):
    # Calculate the maximum width for each column
    col_widths = [max(len(str(cell)) for cell in col) for col in zip(headers, *rows)]
    
    # Create a format string for each row
    row_format = " | ".join(["{:<" + str(width) + "}" for width in col_widths])
    
    # Print the header
    print(row_format.format(*headers))
    
    # Print the separator
    print("-+-".join(['-' * width for width in col_widths]))
    
    # Print the rows
    for row in rows:
        print(row_format.format(*row))

def parse_logq(logq_str):
    logq = []
    parts = logq_str.split(';')
    for part in parts:
        if '-' in part:
            start, end = map(int, part.split('-'))
            logq.extend(range(start, end + 1))
        else:
            logq.append(int(part))
    return logq

def run_verification(lq,secret,est_usvp,est_bdd,est_usvp_pow,est_bdd_pow):
    lwe_parameters_usvp = []
    lwe_parameters_bdd = []
    if(secret == 'binary'):
        lwe_parameters_usvp = LWE.Parameters(est_usvp, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd = LWE.Parameters(est_bdd, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
        lwe_parameters_usvp_pow = LWE.Parameters(est_usvp_pow, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd_pow = LWE.Parameters(est_bdd_pow, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
    else:
        lwe_parameters_usvp = LWE.Parameters(est_usvp, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd = LWE.Parameters(est_bdd, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
        lwe_parameters_usvp_pow = LWE.Parameters(est_usvp_pow, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd_pow = LWE.Parameters(est_bdd_pow, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))

    lwe_usvp = math.floor(math.log2(LWE.primal_usvp(lwe_parameters_usvp)["rop"]))
    lwe_bdd = math.floor(math.log2(LWE.primal_bdd(lwe_parameters_bdd)["rop"]))
    lwe_usvp_pow = math.floor(math.log2(LWE.primal_usvp(lwe_parameters_usvp_pow)["rop"]))
    lwe_bdd_pow = math.floor(math.log2(LWE.primal_bdd(lwe_parameters_bdd_pow)["rop"]))

    return lwe_usvp, lwe_bdd, lwe_usvp_pow, lwe_bdd_pow