import math
import csv
import sys

sys.path.append('../lattice-estimator')


try:
    from estimator import *
except ImportError:
    print("Warning: Failed to import lattice_estimator, some options will not work")

#Auxiliary functions needed in estimate.py

#Exctracted from the Lattice Estimator
def UniformModStd(q):
    a = -(q // 2)
    b = -a -1 if q % 2 == 0 else -a

    if b < a:
        raise ValueError(f"upper limit must be larger than lower limit but got: {b} < {a}")
    m = b - a + 1
    mean = (a + b) / float(2)
    stddev = math.sqrt((m**2 - 1) / float(12))

    return stddev

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

def helper():
    #print('python3 estimate.py --param "lambda" --file "example_lambda_binary.csv"')
    # print('python3 estimate.py --param "n" --file "example_n_binary.csv"')
    # print('python3 estimate.py --param "logq" --file "example_logq.csv"')
    # print('python3 estimate.py --param "error" --file "example_err.csv"')
    print('python3 estimate.py --param "lambda" --n "1024" --logq "20-30;35;40-60" --secret "binary" --error "3.19"')
    print('python3 estimate.py --param "n" --lambda "80" --logq "20-30" --secret "binary" --error "3.19"')
    print('python3 estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19"')
    print('python3 estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --secret "binary"')
    print('You can add  --verify 1 to any of the above commands to check the results against the Lattice Estimator')
    sys.exit()

paper = 'https://eprint.iacr.org/2024/1001'

def create_explanation_dict(headers):
    explanations = {
        "Secret dist.": "The distribution of the secret (can be either binary or ternary)",
        "LWE dim.": "The Learning With Errors (LWE) dimension",
        "lambda": "The security level",
        "log q": "The size of the modulus q in bits",
        "usvp_s (Eq. 21)": "The output of Eq. 21 of " + paper,
        "lwe est": "The output of running the Lattice Estimator using the output of our formulas and the rest of the LWE parameters",
        "usvp_s pow2": "Closest power of 2 to the output of Eq. 21",
        "bdd_s (Eq. 22)": "The output of Eq. 22 of " + paper,
        "bdd_s pow2": "Closest power of 2 to the output of Eq. 22",
        "bdd": "The output of Eq. XX of " + paper, #TODO update this reference
        "bdd pow2": "Closest power of 2 to the output of Eq. XX", #TODO update this reference
        "usvp (Eq. 14)": "The output of Eq. 14 of " + paper,
        "usvp_s (Eq. 16)": "The output of Eq. 16 of " + paper,
        "bdd (Eq. 17)": "The output of Eq. 17 of " + paper,
        "bdd_s (Eq. 20)": "The output of Eq. 20 of " + paper,
        "logq usvp": "The result of numerically approximating log q using usvp",
        "logq bdd": "The result of numerically approximating log q using bdd",
        "std_e usvp": "The result of numerically approximating the standard deviation of the error using usvp",
        "std_e bdd": "The result of numerically approximating the standard deviation of the error using bdd",
        "bdd 3.19": "The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_bdd",
        "usvp 3.19": "The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_usvp",
        "diff": "The difference between the output of the previous column and the output of the Lattice Estimator"
    }

    # Create a dictionary using the headers and explanations
    explanation_dict = {}
    for header in headers:
        # Add the explanation if it exists in the explanations dictionary, otherwise use a default message
        explanation_dict[header] = explanations.get(header, "No explanation available for this header.")

    return explanation_dict

def helper_headers(header):
    explanation_dict = create_explanation_dict(header)

    max_length = max(len(header) for header in explanation_dict.keys())
    max_length_exp = max(len(explanation) for explanation in explanation_dict.values())

    # Print each header and its explanation with proper formatting
    for header, explanation in explanation_dict.items():
        print(f"{header:<{max_length}}: {explanation}")

    print("." * max_length_exp)
    print('\n')
