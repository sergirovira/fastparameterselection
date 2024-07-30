import math
import sqlite3
import csv
from os.path import exists
import ast
import sys, getopt
import pickle
import json
import os
from pathlib import Path

from formula_params import *
from formulas import *
from aux import *
from numerical_solver import *

sys.path.append('../lattice-estimator')
estimator_installed = 1

try:
    from estimator import *
except ImportError:
    print("Warning: Failed to import lattice_estimator, some options will not work")
    estimator_installed = 0

def main(argv):
    secret = "binary"
    param = None
    lwe_d = None
    logq = 64
    file_path = None
    verify = 0
    lwe_parameters = []
    std_s = 0.5
    std_e = 3.19
    secret_q = 2

    headers = []
    data = []

    try:
        opts, args = getopt.getopt(argv, "h", ["secret=", "error=", "param=", "n=", "lambda=", "logq=", "verify=", "file="])
    except getopt.GetoptError:
        print('python3 estimate.py --param "lambda" --file "example_lambda_binary.csv"')
        print('python3 estimate.py --param "lambda" --n "1024" --logq "20-30;35;40-60" --dist "binary"')
        print('python3 estimate.py --param "n" --lambda "80" --logq "20-30" --dist "binary" --verify 1')
        print('python3 estimate.py --param "n" --file "example_n_ternary.csv"')
        print('python3 estimate.py --param "logq" --lambda "80" --n "1024" --dist "binary"')
        print('python3 estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --dist "binary"')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('python3 estimate.py --param "lambda" --file "example_lambda_binary.csv"')
            print('python3 estimate.py --param "lambda" --n "1024" --logq "20-30;35;40-60" --dist "binary"')
            print('python3 estimate.py --param "n" --lambda "80" --logq "20-30" --dist "binary" --verify 1')
            print('python3 estimate.py --param "n" --file "example_n_ternary.csv"')
            print('python3 estimate.py --param "logq" --lambda "80" --n "1024" --dist "binary"')
            print('python3 estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --dist "binary"')
            sys.exit()
        elif opt == '--secret':
            secret = arg
            if secret == 'binary': 
                std_s = UniformModStd(2)
                secret_q = 2
            elif secret == 'ternary': 
                std_s = UniformModStd(3)
                secret_q = 3
            else: 
                print("Secret distribution not supported")
                sys.exit() 
        elif opt == '--error':
            std_e = float(arg)
        elif opt == '--param':
            param = arg
        elif opt == '--n':
            lwe_d = int(arg)
        elif opt == '--lambda':
            l = int(arg)
        elif opt == '--logq':
            logq = parse_logq(arg)
        elif opt == '--verify':
            verify = int(arg)
        elif opt == '--file':
            file_path = arg

    if secret == "binary":
        lambda_usvp = lambda_usvp_bin
        lambda_usvp_s = lambda_usvp_s_bin
        lambda_bdd = lambda_bdd_bin
        lambda_bdd_s = lambda_bdd_s_bin
        n_usvp_s = n_usvp_s_bin
        n_bdd_s = n_bdd_s_bin
    else:
        lambda_usvp = lambda_usvp_ter
        lambda_usvp_s = lambda_usvp_s_ter
        lambda_bdd = lambda_bdd_ter
        lambda_bdd_s = lambda_bdd_s_ter
        n_usvp_s = n_usvp_s_ter
        n_bdd_s = n_bdd_s_ter


    # If we select to run the formulas for the LWE dimension, we get an output of the following form:
    #
    #Secret dist. | lambda | log q | usvp_s (Eq. 21) | lwe est | usvp_s pow2 | lwe est | bdd_s (Eq. 22) | lwe est | bdd_s pow2 | lwe est
    #-------------+--------+-------+-----------------+---------+-------------+---------+----------------+---------+------------+--------
    #binary       | 80     | 20    | 501             | 84      | 512         | 85      | 531            | 86      | 512        | 84     
    #binary       | 80     | 21    | 525             | 83      | 512         | 81      | 554            | 86      | 512        | 80     
    #binary       | 80     | 22    | 550             | 83      | 512         | 78      | 577            | 85      | 512        | 76  
    #
    # The column 'lwe est' is the output of the Lattice Estimator run with the parameters in columns 'Secret dist.', 'log q' and the previous 
    # column. For example, The value 84 in the fourth column is obtained by running the Lattice Estimator with n = 501. 

    if param == 'n':
        
        if(verify):
            headers = ["Secret dist.", "lambda", "log q", "usvp_s (Eq. 21)", "lwe est", "usvp_s pow2", "lwe est", "bdd_s (Eq. 22)", "lwe est", "bdd_s pow2", "lwe est"]
        else:
            headers = ["Secret dist.", "lambda", "log q", "usvp_s (Eq. 21)", "usvp_s pow2", "usvp_s num", "bdd_s (Eq. 22)", "bdd_s pow2", "bdd_s num"]

        data = []
        if file_path:
            entries = load_all_from_csv(file_path)
            for entry in entries:
                secret = entry['secret']
                l = int(entry['lambda'])
                logq = int(entry['logq'])
                est_usvp = int(math.ceil(model_n_usvp(l, logq, n_usvp_s)))
                est_bdd = int(math.ceil(model_n_bdd(l, logq, std_s, std_e, n_bdd_s)))
                est_usvp_pow = closest_power_of_2(est_usvp)
                est_bdd_pow = closest_power_of_2(est_bdd)
                if(verify and estimator_installed):
                    lwe_usvp, lwe_bdd, lwe_usvp_pow, lwe_bdd_pow = run_verification(logq,secret,est_usvp,est_bdd,est_usvp_pow,est_bdd_pow)
                    data_point = [secret, l, logq, est_usvp, lwe_usvp, est_usvp_pow, lwe_usvp_pow, est_bdd, lwe_bdd, est_bdd_pow, lwe_bdd_pow]
                else:
                    data_point = [secret, l, logq, est_usvp, est_usvp_pow, est_bdd, est_bdd_pow]
                data.append(data_point)
        else:
            for lq in logq:
                est_usvp = int(math.ceil(model_n_usvp(l, lq, n_usvp_s)))
                est_bdd = int(math.ceil(model_n_bdd(l, lq,  std_s, std_e, n_bdd_s)))
                est_usvp_pow = closest_power_of_2(est_usvp)
                est_bdd_pow = closest_power_of_2(est_bdd)
                std_s = 0.816496580927726
                std_e = 3.19
                est_usvp_numerical = int(math.ceil(numerical_n_usvp(est_usvp, lq,  std_s, std_e)))
                est_bdd_numerical = int(math.ceil(numerical_n_bdd(est_bdd, lq,  std_s, std_e)))

                if(verify and estimator_installed):
                    lwe_usvp, lwe_bdd, lwe_usvp_pow, lwe_bdd_pow = run_verification(lq,secret,est_usvp,est_bdd,est_usvp_pow,est_bdd_pow)
                    data_point = [secret, l, lq, est_usvp, lwe_usvp, est_usvp_pow, lwe_usvp_pow, est_bdd, lwe_bdd, est_bdd_pow, lwe_bdd_pow]

                else:
                    data_point = [secret, l, lq, est_usvp, est_usvp_pow, est_usvp_numerical, est_bdd, est_bdd_pow, est_bdd_numerical]

                data.append(data_point)

    if param == 'logq':
        
        if(verify):
            headers = ["Secret dist.", "lambda", "n", "logq usvp", "logq bdd", "lwe est"]
        else:
            headers = ["Secret dist.", "lambda", "n", "logq usvp", "logq bdd"]

        if file_path:
            entries = load_all_from_csv(file_path)
            for entry in entries:
                secret = entry['secret']
                l = int(entry['lambda'])
                logq = int(entry['logq'])
                est_usvp = int(math.ceil(model_n_usvp(l, logq, n_usvp_s_bin)))
                est_bdd = int(math.ceil(model_n_bdd(l, logq, secret, n_bdd_s_bin)))
                est_usvp_pow = closest_power_of_2(est_usvp)
                est_bdd_pow = closest_power_of_2(est_bdd)
                if(verify and estimator_installed):
                    lwe_usvp, lwe_bdd, lwe_usvp_pow, lwe_bdd_pow = run_verification(logq,secret,est_usvp,est_bdd,est_usvp_pow,est_bdd_pow)
                    data_point = [secret, l, logq, est_usvp, lwe_usvp, est_usvp_pow, lwe_usvp_pow, est_bdd, lwe_bdd, est_bdd_pow, lwe_bdd_pow]
                else:
                    data_point = [secret, l, logq, est_usvp, est_usvp_pow, est_bdd, est_bdd_pow]
                data.append(data_point)
        else:
                est_usvp_numerical = int(math.ceil(numerical_logq_usvp(l, lwe_d, std_s, std_e)))
                est_bdd_numerical = int(math.ceil(numerical_logq_bdd(l, lwe_d, std_s, std_e)))

                if(verify and estimator_installed):
                    lwe_parameters_bdd = LWE.Parameters(lwe_d, 2 ** est_bdd_numerical, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
                    lwe_bdd = math.floor(math.log2(LWE.primal_bdd(lwe_parameters_bdd)["rop"]))
                    data_point = [secret, l, lwe_d, est_usvp_numerical, est_bdd_numerical,lwe_bdd]
                else:
                    data_point = [secret, l, lwe_d, est_usvp_numerical, est_bdd_numerical]

                data.append(data_point)

    if param == 'std_e':
        
        if(verify):
            headers = ["Secret dist.", "lambda", "n", "logq", "std_e bdd", "lwe est"]
        else:
            headers = ["Secret dist.", "lambda", "n", "logq", "std_e usvp", "std_e bdd"]

        if file_path:
            entries = load_all_from_csv(file_path)
            for entry in entries:
                secret = entry['secret']
                l = int(entry['lambda'])
                logq = int(entry['logq'])
                est_usvp = int(math.ceil(model_n_usvp(l, logq, n_usvp_s_bin)))
                est_bdd = int(math.ceil(model_n_bdd(l, logq, secret, n_bdd_s_bin)))
                est_usvp_pow = closest_power_of_2(est_usvp)
                est_bdd_pow = closest_power_of_2(est_bdd)
                if(verify and estimator_installed):
                    lwe_usvp, lwe_bdd, lwe_usvp_pow, lwe_bdd_pow = run_verification(logq,secret,est_usvp,est_bdd,est_usvp_pow,est_bdd_pow)
                    data_point = [secret, l, logq, est_usvp, lwe_usvp, est_usvp_pow, lwe_usvp_pow, est_bdd, lwe_bdd, est_bdd_pow, lwe_bdd_pow]
                else:
                    data_point = [secret, l, logq, est_usvp, est_usvp_pow, est_bdd, est_bdd_pow]
                data.append(data_point)
        else:

            for lq in logq:

                std_s = 0.816496580927726
                est_usvp_numerical = numerical_std_e_usvp(lwe_d, lq, std_s)
                est_bdd_numerical = numerical_std_e_bdd(lwe_d, lq, std_s)

                if(verify and estimator_installed):
                    lwe_parameters_bdd = LWE.Parameters(lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(est_bdd_numerical))
                    lwe_bdd = math.floor(math.log2(LWE.primal_bdd(lwe_parameters_bdd)["rop"]))
                    data_point = [secret, l, lwe_d, lq, est_bdd_numerical,lwe_bdd]
                else:
                    data_point = [secret, l, lwe_d, lq, est_usvp_numerical, est_bdd_numerical]

                data.append(data_point)

    # If we select to run the formulas for the security level, we get an output of the following form:
    #
    #Secret dist. | LWE dim. | log q | usvp (Eq. 14) | diff | usvp_s (Eq. 16) | diff | bdd (Eq. 17) | diff | bdd_s (Eq. 20) | diff
    #-------------+----------+-------+---------------+------+-----------------+------+--------------+------+----------------+-----
    #binary       | 1024     | 20    | 170           | 1    | 172             | 3    | 173          | 6    | 173            | 6   
    #binary       | 1024     | 24    | 139           | -2   | 142             | 1    | 142          | 3    | 142            | 3   
    #binary       | 1024     | 25    | 133           | -2   | 136             | 1    | 136          | 3    | 136            | 3  
    #
    # The column 'diff' is the output of the difference between the previous column and the output of the Lattice Estimator run with 
    # the parameters in columns 'Secret dist.', 'log q' and 'LWE dim.'For example, The value 1 in the fourth column is obtained by 
    # running the Lattice Estimator and substructing the output to 170. 


    if param == 'lambda':

        if(verify):
            headers = ["Secret dist.", "LWE dim.", "log q", "usvp (Eq. 14)", "diff", "usvp_s (Eq. 16)", "diff", "bdd (Eq. 17)", "diff", "bdd_s (Eq. 20)", "diff", "Estimator"]
        else:
            headers = ["Secret dist.", "LWE dim.", "log q", "usvp (Eq. 14)", "usvp_s (Eq. 16)", "bdd (Eq. 17)", "bdd_s (Eq. 20)"]

        if file_path:
            entries = load_all_from_csv(file_path)
            for entry in entries:
                secret = entry['secret']
                lwe_d = int(entry['lwe_d'])
                logq = int(entry['logq'])
                est_usvp = int(round(model_lambda_usvp(lwe_d, logq, secret, lambda_usvp_bin)))
                est_usvp_s = int(round(model_lambda_usvp_s(lwe_d, logq, lambda_usvp_s_bin)))
                est_bdd = 0
                try: 
                    est_bdd = int(round(model_lambda_bdd(lwe_d, logq, secret, lambda_bdd_bin)[0].real))
                except Exception as e:
                    pass
                est_bdd_s = int(round(model_lambda_bdd_s(lwe_d, logq, lambda_bdd_s_bin)))
                if(verify and estimator_installed):
                    lwe_usvp, lwe_bdd, lwe_usvp_s, lwe_bdd_s = run_verification(logq,secret,lwe_d,lwe_d,lwe_d,lwe_d)
                    data_point = [secret, lwe_d, logq, est_usvp, est_usvp-lwe_usvp, est_usvp_s, est_usvp_s-lwe_usvp_s, est_bdd, est_bdd-lwe_bdd, est_bdd_s, est_bdd_s-lwe_bdd_s]
                else:
                    data_point = [secret, lwe_d, logq, est_usvp, est_usvp_s, est_bdd, est_bdd_s]
                data.append(data_point)
        
        else:
            for lq in logq:

                est_usvp = int(round(model_lambda_usvp(lwe_d, lq, std_s, std_e, lambda_usvp)))
                est_usvp_s = int(round(model_lambda_usvp_s(lwe_d, lq, lambda_usvp_s)))
                est_bdd = 0
                try: 
                    est_bdd = int(round(model_lambda_bdd(lwe_d, lq, std_s, std_e, lambda_bdd)[0].real))
                except Exception as e:
                    pass
                est_bdd_s = int(round(model_lambda_bdd_s(lwe_d, lq, lambda_bdd_s)))

                if(verify and estimator_installed):
                    lwe_usvp, lwe_bdd, lwe_usvp_s, lwe_bdd_s = run_verification(lq,secret,lwe_d,lwe_d,lwe_d,lwe_d)
                    data_point = [secret, lwe_d, lq, est_usvp, est_usvp-lwe_usvp, est_usvp_s, est_usvp_s-lwe_usvp_s, est_bdd, est_bdd-lwe_bdd, est_bdd_s, est_bdd_s-lwe_bdd_s]
                else:
                    data_point = [secret, lwe_d, lq, est_usvp, est_usvp_s, est_bdd, est_bdd_s]

                data.append(data_point)

    print_table(headers,data)
    
    if param == "est":
        for lq in logq:
            parameters = LWE.Parameters(lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
            LWE.estimate(parameters)

    
    print("\n")
    if(verify and not estimator_installed): 
        print("Warning: Verification not possible, Lattice Estimator not installed")
    if(param == "lambda" and not estimator_installed):
        print("Warning: Bdd set to 0")
    print("\n")

if __name__ == "__main__":
    main(sys.argv[1:])
