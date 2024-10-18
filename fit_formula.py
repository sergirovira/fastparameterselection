#!/usr/bin/env python3
import sys
sys.path.append('../lattice-estimator')

from estimator import *

import lmfit
import math
import numpy as np
import sqlite3
import csv
from os.path import exists
import ast
import sys, getopt
import time

import pickle
import json
import matplotlib.lines as mlines
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

from multiprocessing import Pool, Manager
import itertools
import matplotlib.colors as mcolors
import random

from config import *
from formulas import *
from aux import *

secret = "binary"
param = None
lwe_d = None
logQ = 64
name_file = None
verify = 0
lwe_parameters = []
std_s = 0.5
std_e = 3.19
secret_q = 2
ntru_flag = False
simpl = 0

headers = []
data = []
params_list = [('P0', 0.07,0), ('P1', 0.34,0), ('P2', 1,0), ('P3', 1,0), ('P4', 1,0)] #, ('P5', 1,0)]

verbose = 0

def write_tuples_to_txt(tuples_list, filename):
    with open(filename, 'w') as file:
        for tup in tuples_list:
            file.write(','.join(map(str, tup)) + '\n')

def merge_txt_files(input_files, output_file):
    with open(output_file, 'w') as out_file:
        for file_name in input_files:
            with open(file_name, 'r') as in_file:
                out_file.write(in_file.read())

#In the function 'model_n_bdd' these are called by their position.
#E.g, alpha = params[0], beta = params[1], gamma = params[2], delta = params[3].
#If you add a new parameter and don't use it in 'model_n_bdd', a math domain error will triger. 


def random_color_generator():
    color = random.choice(list(mcolors.XKCD_COLORS.keys()))
    return color

#following eq 6 of the paper
def compute_lambda(d, logq, secret="u3", error="dg3.19"):
    sigma = 3.19
    eta = 2*sigma

    if(secret == "u3"): eta = math.sqrt(float(3)/2)*sigma

    if(logq >= d): return {}

    numerator = 1.75 * d * math.log(float(d) / logq)
    A = logq + math.sqrt(numerator / logq) - math.log(2 * math.pi * math.e * sigma)
    #print("Numerator: ", numerator)
    #print("A: ", A)
    lowerbound = (2 * d * (logq - math.log(eta)) * math.log(numerator / (2 * math.pi * math.e * logq))) / A**2 #This is the right-hand side of eq. 6
    #print("lowerbound: ", numerator / (2 * math.pi * math.e * logq))
    lmbda = 0.292 * lowerbound #this should give an approximation of lambda
    dic = {}
    dic['bdd'] = lmbda
    return dic

#"binary" "u3"
def estimate(d, logq, attacks=None, secret="u3", error="dg3.19", force=False, verbose=True):
    '''
    Database-caching wrapper for the Lattice Estimator

    Parameters:
      d:       polynomial degree
      logq:    bit size of the ciphertext modulus
      attacks: list of attacks to estimate, supported values are:
               "arora-gb", "coded-bkw", "usvp", "bdd", "bdd_hybrid", "dual", "dual_hybrid"
      secret:  the secret distribution, supported values are:
               "cbd1", "cbd21", "dg3.19", "u3"

    Return Value:
      A dictonary with estimates for each attack.
    '''
    db = sqlite3.connect(database, timeout=30.0)

    # input casting (for sage)
    # convert sage.rings.integer.Integer to int
    d = int(d)
    #d = math.log(2,int(d))
    logq = int(logq)

    # input validation
    allAttacks = ['arora-gb', 'coded-bkw', 'usvp', 'bdd', 'bdd_hybrid', 'dual', 'dual_hybrid']
    skipAttacks = ['arora-gb', 'coded-bkw','dual_hybrid']


    versionID = db.execute("SELECT versionID FROM Versions WHERE hash = ?", (version,)).fetchone()[0]
    if versionID == None:
        raise Exception("unsupported version", version)
    versionID = 3

    attackIDs = []
    if attacks == None:
        attacks = allAttacks
    for attack in attacks:
        if attack not in skipAttacks:
            attackID = db.execute("SELECT attackID FROM attacks WHERE name = ?", (attack,)).fetchone()
            if attackID == None:
                raise Exception("unsupported attack", attack)
            attackIDs.append(attackID[0])


    secretID = db.execute("SELECT distributionID FROM distributions WHERE name = ?", (secret,)).fetchone()
    if secret == "umod" or secretID == None:
        raise Exception("unsupported secret distribution", secret)
    secretID = 4

    errorID = db.execute("SELECT distributionID FROM distributions WHERE name = ?", (error,)).fetchone()
    if error == "umod" or errorID == None:
        raise Exception("unsupported error distribution", error)
    errorID = errorID[0]

    # setup estimator
    dists = [
        ND.CenteredBinomial(1),
        ND.CenteredBinomial(21),
        ND.DiscreteGaussian(3.19),
        ND.UniformMod(3),
        ND.UniformMod(2),
    ]
    params = LWE.Parameters(d, 2**logq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))

    estimates = {}
    estimates_list = []

    for attackID in attackIDs:
        # check for cached result
        #print("versionID {}, attackID {}, secretID {}, errorID {}, d {}, logq {}".format(versionID, attackID, secretID, errorID, d, logq))
        estimate = db.execute("SELECT estimate FROM Estimates WHERE (versionID, attackID, secretID, errorID, d, logq) = (?,?,?,?,?,?)", (versionID, attackID, secretID, errorID, d, logq)).fetchone()
        #print("estimate", estimate)
        if estimate != None:
            estimates[allAttacks[attackID - 1]] = estimate[0]
        continue 
        # run the lattice estimator
        if verbose:
            print(f"Running lattice estimator for {allAttacks[attackID - 1]} with degree {d} and modulus size {logq}")

        # insert null for the parameter crashing
        #db.execute("INSERT INTO Estimates (versionID, attackID, secretID, errorID, d, logq, estimate) VALUES (?,?,?,?,?,?,?)", (versionID, attackID, secretID, errorID, d, logq, None))
        #db.commit()

        #red_cost_model=RC.BDGL16
        rop = 0

        r_error = 0;

        try: 
            if attackID == 1:
                rop = LWE.arora_gb(params).rop
            elif attackID == 2:
                rop = LWE.coded_bkw(params).rop
            elif attackID == 3:
                rop = LWE.primal_usvp(params,red_cost_model=RC.BDGL16).rop
            elif attackID == 4:
                rop = LWE.primal_bdd(params,red_cost_model=RC.BDGL16).rop
            elif attackID == 5:
                rop = LWE.primal_hybrid(params,red_cost_model=RC.BDGL16).rop
            elif attackID == 6:
                rop = LWE.dual(params,red_cost_model=RC.BDGL16).rop
            elif attackID == 7:
                rop = LWE.dual_hybrid(params).rop
            estimate = math.log(rop, 2)
            estimates[allAttacks[attackID - 1]] = estimate
            r_error = 0
        except:
            print("runtime error")
            r_error = 1
            pass

        #print("estimate, versionID, attackID, secretID, errorID, d, logq", estimate, versionID, attackID, secretID, errorID, d, logq)
        if(r_error == 0):
            estimates_list.append((estimate, versionID, attackID, secretID, errorID, d, logq))

    db.close() 

    return estimates, estimates_list


def degree_loss(model, d):
    diff = d-model
    # for i, x in enumerate(diff):
    #     if x > 0:
    #         diff[i] = 0

    return diff

def error(model, d):
    mape = np.nanmean(np.abs((d - model)/d))*100
    return mape

def merge(list1, list2):
    
    merged_list = []
    for i in range(0, len(list1)):
        #if not math.isnan(list1[i]): merged_list.append((list1[i],list2[i]))
        merged_list.append((list1[i],list2[i]))
    return merged_list


    residues = ()
    if d60 != None:
        residues = np.concatenate((residues, degree_loss(model60, d60)))
    if d70 != None:
        residues = np.concatenate((residues, degree_loss(model70, d70)))
    if d80 != None:
        residues = np.concatenate((residues, degree_loss(model80, d80)))
    if d100 != None:
        residues = np.concatenate((residues, degree_loss(model100, d100)))
    if d128 != None:
        residues = np.concatenate((residues, degree_loss(model128, d128)))
    if d150 != None:
        residues = np.concatenate((residues, degree_loss(model150, d150)))
    if d192 != None:
        residues = np.concatenate((residues, degree_loss(model192, d192)))

    return residues

def degree_fit(params, logq, levels, e_std, s_std, std_s_num):

    residues = ()
    model = None
    new_params = []

    for p in params_list:
        new_params.append(params[p[0]])

    for i, level in enumerate(levels):
        #print("level fit: ", level)
        #print("attack: ", attack)

        #print("degrees", degrees)

        #print("i: ", i)

        if param == 'n' and simpl == '0':
            if attack == 'bdd' :
                model  = model_n_bdd(security_levels[i], logq, std_s, std_e, new_params) #, params['delta'])
            if attack == 'usvp':
                model  = model_n_usvp(security_levels[i], logq, new_params) #, params['delta'])

        if param == 'n' and simpl == '1':
            if attack == 'bdd' :
                model  = model_n_bdd_s(security_levels[i], logq, std_s, std_e, new_params) #, params['delta'])
            if attack == 'usvp':
                model  = model_n_usvp_s(security_levels[i], logq, new_params) #, params['delta'])

        if param == 'lambda' and simpl == '0':
            if attack == 'usvp':
                model  = model_lambda_usvp(degrees[i], logq, s_std, e_std, new_params) #, params['delta'])
            if attack == 'bdd':
                model = model_lambda_bdd(degrees[i],  logq,  s_std, e_std, std_s_num, new_params) #,delta)
        if param == 'lambda' and simpl == '1':
            if attack == 'usvp':
                model  = model_lambda_usvp_s(degrees[i], logq, new_params) #, params['delta'])
            if attack == 'bdd':
                model  = model_lambda_bdd_s(degrees[i], logq, new_params) #, params['delta'])

        if level != None:
            residues = np.concatenate((residues, degree_loss(model, level)))

    return residues

def convert(array):
    li = [float(m) if m != 'nan' else np.nan for m in array]
    return li

def convert_int(array):
    li = [int(m) if m != 'nan' else np.nan for m in array]
    return li

def get_position_in_bound(degree):
    if(degree >= 2**10 and degree < 2**11): return 0
    if(degree >= 2**11 and degree < 2**12): return 1
    if(degree >= 2**12 and degree < 2**13): return 2
    if(degree >= 2**13 and degree < 2**(13.5)): return 3
    if(degree >= 2**(13.5) and degree < 2**14): return 4
    if(degree >= 2**14 and degree < 2**(14.5)): return 5
    if(degree >= 2**(14.5) and degree <=2**15): return 6
    return -1



def check_bound(degree, q, bounds):
    pos = get_position_in_bound(degree)

    if(pos == -1): 
        print("No bound given for d >= ", degree)
        exit(0)

    min_q = bounds[pos][0]
    max_q = bounds[pos][1]

    if(q > min_q and q <= max_q): return 1
    else: return 0

def plot_values(data):
    
    for v in data:
        x_value = v[1]
        y_value = v[0]
        plt.scatter(x_value, y_value, label=f"{v[2]}")

    plt.xlabel('a')
    plt.ylabel('t1')
    plt.title('Scatter Plot')
    plt.legend()
    plt.show()

def save_dict_to_file(dictionary, filename):
    with open(filename, 'wb') as file:
        pickle.dump(dictionary, file)

def read_dict_from_file(filename):
    with open(filename, 'rb') as file:
        data = pickle.load(file)
        return data

def fit_formula(points_est, e_std, s_std, std_s_num, params):

    args = {'logq': logQ}

    args['levels'] = points_est
    args['e_std'] = e_std
    args['s_std'] = s_std
    args['std_s_num'] = std_s_num 

    if(verbose): print(points_est)

    fit = []

    fit = lmfit.minimize(degree_fit, params, nan_policy='omit', kws=args)
    #print(lmfit.fit_report(fit))

    fit_results = []

    for p in params_list:
        res = fit.params[p[0]].value
        fit_results.append(res)

    return fit_results

def plot_points_est(param, points_atk, points_est, points_secret_dist, fit_results, std_s, std_e):
    fig, ax = plt.subplots(figsize=(11,8))

    fig.gca().set_ylabel(r'$d$')
    fig.gca().set_xlabel(r'$\log q$')
    
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)

    j = 0
    offsetx = 3
    offsety = 0.5
    parity = 0

    colours = []

    markers = ['x', 'o', '^', 'v', '1', 'P', '*']

    if(attack == 'bdd'):
        markers = ['x','o']
    if(attack == 'usvp'):
        markers = ['o','x']

    plotted_points = [[] for _ in range(len(points_atk))]
    
    num_points = 0

    for j, d_attack_level in enumerate(points_atk):

        color = random_color_generator()
        colours.append(color)

        for i, a in enumerate(d_attack_level):
        
            estimate = points_est[j][i]

            if not np.isnan(estimate):
                if(verbose): print("point: ", "logq", logQ[i], "estimate", estimate, "color", color, "n: ", degrees[j])
                num_points += 1
            
            if a == 0:
                ax.plot(logQ[i], estimate, marker='x', linestyle='', color=color)
                plotted_points[j].append((logQ[i],estimate,color))
            elif a == 1:
                ax.plot(logQ[i], estimate, marker='o', linestyle='', color=color)
                plotted_points[j].append((logQ[i],estimate,color)) 
            elif a == 2:
                ax.plot(logQ[i], estimate, marker='^', linestyle='', color=color)
            elif a == 3:
                ax.plot(logQ[i], estimate, marker='v', linestyle='', color=color)
            elif a == 4:
                ax.plot(logQ[i], estimate, marker='1', linestyle='', color=color)
            elif a == 5:
                ax.plot(logQ[i], estimate, marker='P', linestyle='', color=color)
            elif a == 6:
                ax.plot(logQ[i], estimate, marker='*', linestyle='', color=color)

    # for j, d_attack_level in enumerate(points_atk):

    #     color = random_color_generator()
    #     colours.append(color)

    #     for i, a in enumerate(d_attack_level):
        
    #         estimate = points_est[j][i]
    #         sigma_e = points_secret_dist[j][i]

    #         if not np.isnan(estimate):
    #             print("point: ", "sigma_e", sigma_e, "estimate", estimate, "color", color, "n: ", degrees[j])
    #             num_points += 1
            
    #         if a == 0:
    #             ax.plot(sigma_e, estimate, marker='x', linestyle='', color=color)
    #             plotted_points[j].append((sigma_e,estimate,color))
    #         elif a == 1:
    #             ax.plot(sigma_e, estimate, marker='o', linestyle='', color=color)
    #             plotted_points[j].append((sigma_e,estimate,color)) 
    #         elif a == 2:
    #             ax.plot(logQ[i], estimate, marker='^', linestyle='', color=color)
    #         elif a == 3:
    #             ax.plot(logQ[i], estimate, marker='v', linestyle='', color=color)
    #         elif a == 4:
    #             ax.plot(logQ[i], estimate, marker='1', linestyle='', color=color)
    #         elif a == 5:
    #             ax.plot(logQ[i], estimate, marker='P', linestyle='', color=color)
    #         elif a == 6:
    #             ax.plot(logQ[i], estimate, marker='*', linestyle='', color=color)

    
    if(verbose): print("plotted points: ", plotted_points)

    modelQ = list(range(xmin,logQ[0])) + logQ + list(range(logQ[-1] + 1,xmax))

    error_file = {}

    if os.path.isfile('dict_errors.pkl'):
        os.remove('dict_errors.pkl')

    if param == 'n':
        for i, level in enumerate(security_levels):
            if attack == 'bdd' and simpl == '0':
                model = model_n_bdd(level,  modelQ, std_s, std_e, fit_results) #,delta)
                model_error = model_n_bdd(level,  logQ, std_s, std_e, fit_results) #,delta)
            if attack == 'bdd' and simpl == '1':
                model = model_n_bdd_s(level,  modelQ, fit_results) #,delta)
                model_error = model_n_bdd_s(level,  logQ, fit_results) #,delta)
            if attack == 'usvp' and simpl == '0': 
                model = model_n_usvp(level,  modelQ, fit_results) #,delta)
                model_error = model_n_usvp(level,  logQ, fit_results) #,delta)
            if attack == 'usvp' and simpl == '1': 
                model = model_n_usvp_s(level,  modelQ, fit_results) #,delta)
                model_error = model_n_usvp_s(level,  logQ, fit_results) #,delta)


            # error_file[level] = {}

            # for j, q in enumerate(logQ):
            #     if not math.isnan(points_est[i][j]):
            #         error_file[level][q] = (points_est[i][j], model_error[j])

            #print(error_file)
            ax.plot(modelQ, model,  linewidth=2, color=colours[i], linestyle="dotted", label='$lambda=' + str(level) + '$')

    if param == 'lambda':
        for i, level in enumerate(degrees):
            if attack == 'usvp' and simpl == '0':
                model = model_lambda_usvp(level,  modelQ, secret, fit_results) #,delta)
                model_error = model_lambda_usvp(level,  logQ, secret, fit_results) #,delta)
            if attack == 'usvp' and simpl == '1':
                model = model_lambda_usvp_s(level,  modelQ, fit_results) #,delta)
                model_error = model_lambda_usvp_s(level,  logQ, fit_results) #,delta)
            if attack == 'bdd' and simpl == '0':
                model = model_lambda_bdd(level,  modelQ, points_secret_dist[i], fit_results) #,delta)
                model_error = model_lambda_bdd(level,  logQ, points_secret_dist[i], fit_results) #,delta)
                print("model: ", model)
            if attack == 'bdd' and simpl == '1':
                model = model_lambda_bdd_s(level,  modelQ, fit_results) #,delta)
                model_error = model_lambda_bdd_s(level,  logQ, fit_results) #,delta)

            error_file[level] = {}

            for j, q in enumerate(logQ):
                if not math.isnan(points_est[i][j]):
                    error_file[level][q] = (points_est[i][j], model_error[j])

            #print(error_file)
            ax.plot(modelQ, model,  linewidth=2, color=colours[i], linestyle="dotted", label='$lambda=' + str(level) + '$')

    with open('dict_errors.pkl', 'wb') as fp:
        pickle.dump(error_file, fp)

    if(attack == 'bdd'):
        markers = ['o']
    if(attack == 'usvp'):
        markers = ['x']


    handles = []

    handles.append(mlines.Line2D([], [], color='purple', marker=markers[0], linestyle='None',
                          markersize=10, label=attack))

    legend1 = plt.legend(handles=handles, loc=1)

    plt.gca().add_artist(legend1)

    ax.legend(loc=2)

    ax.grid(color="grey", linestyle='dashed', linewidth = 1)
    #plt.axhline(y = 2296, color = 'r', linestyle = '-')
    #plt.axvline(x = 95, color = 'r', linestyle = '-')
    plt.savefig('degree_ter_bdd.png', dpi=1000)

    plt.show()

def filter_points(d_estimates):

    d_estimates_filtered = []

    for d in d_estimates:

        d = sorted(d, key=lambda tup: (tup[0], -tup[1]))

        filtered_list = [tup for i, tup in enumerate(d) if i == 0 or (tup[0] != d[i-1][0] or tup[1] > d[i-1][1])]
        filtered_list2 = [tup for i, tup in enumerate(filtered_list) if i == 0 or tup[1] != filtered_list[i-1][1]]

        d_estimates_filtered.append(filtered_list2)

    return d_estimates_filtered

# Function to load dictionary from file
def load_dictionary(filename):
    with open(filename, 'rb') as file:
        dictionary = pickle.load(file)
        return dictionary

# Function to save dictionary to file
def save_dict(dictionary, filename):
    try:
        with open(filename, 'wb') as file:
            pickle.dump(dictionary, file)
            print("Dictionary saved successfully to", filename)
    except Exception as e:
        print("Error saving dictionary:", e)


def main(argv):

    global simpl
    global param
    global attack
    global logQ
    global degrees
    
    # argv = sys.argv[1:]
   
    # opts, args = getopt.getopt(argv,"hb",["xmin=","xmax=","ymin=","ymax=", "yminbound=", "bound=","levels=", "attack=", "dist=", "param=", "simpl="])

    verbose = 0

    # for opt in opts:
    #     if opt[0] == '-h':
    #         print('python3 find_d_HOPE_linear_sL.py --xmin <min x coordinate> --xmax <max x coordinate> --ymin <min y coordinate> --ymax <max x coordinate> --bound <bound> --levels "level1 level2 ... leveln" ')
    #         sys.exit()
    #     elif opt[0] in ("--xmin"):
    #         xmin = int(opt[1])
    #     elif opt[0] in ("--xmax"):
    #         xmax = int(opt[1])
    #     elif opt[0] in ("--ymin"):
    #         ymin = int(opt[1])
    #     elif opt[0] in ("--ymax"):
    #         ymax = int(opt[1])
    #     elif opt[0] in ("-b,--bound"):
    #         bound = float(opt[1])
    #     elif opt[0] in ("--yminbound"):
    #         yminbound = int(opt[1])
    #     elif opt[0] in ("--attack"):
    #         attack = opt[1]
    #     elif opt[0] in ("--dist"):
    #         secret = opt[1]
    #     elif opt[0] in ("--param"):
    #         param = opt[1]
    #     elif opt[0] in ("--simpl"):
    #         simpl = int(opt[1])
    #     elif opt[0] in ("--levels"):
    #         security_levels = opt[1].split()
    #         security_levels = [eval(i) for i in security_levels]
    #         security_levels.sort(reverse=True)

    #try:
    opts, args = getopt.getopt(argv, "h", ["attack=", "dist=", "secret=", "error=", "param=", "n=", "lambda=", "file=", "simpl="])
    #except:
    #helper_fit()

    #if (len(opts) == 0): helper_fit()

    output_dict = {}

    for opt, arg in opts:
        if opt == '-h':
            helper()
        elif opt == '--attack':
             attack = arg
        elif opt == '--secret':
            secret = arg
            if secret == 'binary': 
                std_s = UniformModStd(2)
                secret_q = 2
                output_dict['std_s'] = 0.5
            elif secret == 'ternary': 
                std_s = UniformModStd(3)
                secret_q = 3
                output_dict['std_s'] = math.sqrt(2./3)
            else: 
                print("Secret distribution not supported")
                sys.exit() 
        elif opt == '--error':
            std_e = float(arg)
            output_dict['std_e'] = std_e
        elif opt == "--param":
            param = arg
        elif opt == '--file':
            name_file = arg
        elif opt == '--simpl':
            simpl = arg
        else:
            helper_fit()

    if param == 'n':
        if secret == "binary":
            degrees = range(2**10, 2**12, 32)
            logQ = list(range(20,65))
            ymin = 1000
            ymax = 2050
            xmin = 20
            xmax = 65
            if attack == 'bdd':
                name_file = 'data_binary_bdd.pkl'
            elif attack == 'usvp':
                name_file = 'data_binary_usvp.pkl'
            else:
                print("Attack " + attack + " not considered")
                exit(1) 
        else:
            degrees = list(range(2**10,2**15,2**2))
            logQ = list(range(10,200)) + list(range(200,500, 10)) + list(range(500,1000,10)) + list(range(1000,1600,50))
            ymin = 15000
            ymax = 32390
            xmin = 400
            xmax = 1400
            name_file = 'data_ternary.pkl'

    if param == "lambda":
        if secret == "binary":
            degrees = [2**10, 2**11]
            logQ = list(range(20,65))
            ymin = 80
            ymax = 256
            xmin = 10
            xmax = 90
            if attack == 'bdd':
                name_file = 'data_binary_bdd.pkl'
            elif attack == 'usvp':
                name_file = 'data_binary_usvp.pkl'
            else:
                print("Attack " + attack + " not considered")
                exit(1) 
        if secret == "ternary":
            degrees = [2**10, 1696, 2**11, 2032, 2**13,9216,10240,11264, 12609, 13633,  14657,15681, 16384,17408, 19456, 22528, 24194, 25218, 28290, 32386, 32768]
            #degrees = [10240,11264, 12609, 13633, 14657, 1568]
            logQ = list(range(10,200)) + list(range(200,500, 10)) + list(range(500,1000,10)) + list(range(1000,1600,50))
            ymin = 80
            ymax = 256
            xmin = 20
            xmax = 1400
            name_file = 'data_ternary.pkl'
        if secret == "tfhe":
            #degrees = [1517, 1567, 1493, 1476, 1542, 1443, 1418, 1369, 1394, 1336, 1319, 1262, 1287, 1237, 1212, 1163, 1517, 1567, 1493, 1476, 1542, 1443, 1418, 1369, 1394, 1336, 1319, 1262, 1287, 1237, 1212, 1163, 1188, 1138, 1113, 1089, 1015, 1064, 990, 1039, 965, 916, 940, 891, 866, 841, 767, 817, 742, 792, 718, 668, 643, 693, 619, 594, 520, 569, 495, 545, 429, 470, 404, 446, 379, 355, 313, 330, 288, 264, 239, 190, 214]
            degrees = [1517, 1567, 1493, 1476, 1542, 2087, 2163, 2061, 2121, 2020]
            logQ = [64]
            ymin = 75
            ymax = 195
            xmin = 1
            xmax = 64
            name_file = 'combined.pkl'
    
    points_est = []
    points_atk = []
    points_secret_dist = []
    number_of_levels = len(security_levels)

    max_length = 0

    estimates_list = []

    items = []
    
    degrees_dict = {}
    
    d_estimates = []

# Read the data back using pickle
    with open(name_file, 'rb') as file:
        d_estimates = pickle.load(file)

    if(verbose): 
        print("file name", name_file)
        print("d_estimates", d_estimates)

    if secret == 'tfhe':
        d_estimates_filtered = [d_estimates]
        degrees = []
        for d in d_estimates:
            degrees.append(d[0])
    else:
        d_estimates_filtered = filter_points(d_estimates)

    #d_estimates_filtered.pop()

    params = lmfit.Parameters()

    for p in params_list:
        params.add(p[0],value = p[1])

    est_dict = {}

    if(param == 'n'):

        for level in security_levels:
            points_est.append([np.nan] * len(logQ))
            points_atk.append([np.nan] * len(logQ))
            points_secret_dist.append([np.nan] * len(logQ))

        for i,d in enumerate(security_levels):
            for tup in d_estimates_filtered[i]:
                if(verbose): print(tup)
                if(attack == 'bdd' and tup[3] == 1):
                    points_est[i][logQ.index(tup[1])] = tup[0]
                    points_atk[i][logQ.index(tup[1])] = tup[3]
                    # if secret == 'tfhe':
                    #     points_secret_dist[i][logQ.index(tup[1])] = 2**tup[4]
                    # else:
                    #     points_secret_dist[i][logQ.index(tup[1])] = 3.19
                if(attack == 'usvp' and tup[3] == 0):
                    points_est[i][logQ.index(tup[1])] = tup[0]
                    points_atk[i][logQ.index(tup[1])] = tup[3]
                    # if secret == 'tfhe':
                    #     points_secret_dist[i][logQ.index(tup[1])] = 2**tup[4]
                    # else:
                    #     points_secret_dist[i][logQ.index(tup[1])] = 3.19

    if(param == 'lambda'):

        for level in degrees:
            points_est.append([np.nan] * len(logQ))
            points_atk.append([np.nan] * len(logQ))
            points_secret_dist.append([np.nan] * len(logQ))

        for filtered in d_estimates_filtered:
            for tup in filtered:
                if(verbose): print("tup: ", tup)

                for i,d in enumerate(degrees):
                    if tup[0] == d and attack == 'usvp' and tup[3] == 0: 
                        points_est[i][logQ.index(tup[1])] = tup[2]
                        points_atk[i][logQ.index(tup[1])] = tup[3]
                        # if secret == 'tfhe':
                        #     points_secret_dist[i][logQ.index(tup[1])] = 2**tup[4]
                        # else:
                        #     points_secret_dist[i][logQ.index(tup[1])] = 3.19
                    if tup[0] == d and attack == 'bdd' and tup[3] == 1: 
                        points_est[i][logQ.index(tup[1])] = tup[2]
                        points_atk[i][logQ.index(tup[1])] = tup[3]
                        # if secret == 'tfhe':
                        #     points_secret_dist[i][logQ.index(tup[1])] = 2**tup[4]
                        # else:
                        #     points_secret_dist[i][logQ.index(tup[1])] = 3.19

    if(verbose): print("points_est",points_est)

    results = fit_formula(points_est, std_e, std_s, secret_q, params)

    #results = [0.28891458, 0.878668965, 19.1069565,1,1]
    #results = [0.26497, 3.25511,-13.69437,1,1]

    if(secret == 'tfhe'):
        for p in d_estimates:
            print(p)
            print("res: ", model_lambda_bdd(p[0],  [p[1]], [2**p[4]], results))
            params_est = LWE.Parameters(p[0], 2 ** 64, ND.UniformMod(2), ND.DiscreteGaussian(2**p[4]))
            print("lattice estimator: ")
            print(LWE.estimate(params_est))
        exit(0)
        #model_lambda_bdd(level,  modelQ, points_secret_dist[i], fit_results)

    #plot_points_est(param, points_atk, points_est, points_secret_dist, results, std_s, std_e)

    #s = param + attack + str(simpl) + secret

    #name = "params_" + s + ".txt"
    #np.savetxt(name, np.array(results))

    print("\nFormula: ", opts)
    print("Params: ", results, "\n")


if __name__ == "__main__":
    main(sys.argv[1:])

