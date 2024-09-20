#!/usr/bin/env python3
from sage.all import *

import lmfit
import math
import numpy as np
import sqlite3
import csv
from os.path import exists
import ast
import sys, getopt
#sys.path.append('lattice-estimator')
#from estimator import *

import istarmap

from tqdm import tqdm

import pickle
import json
import matplotlib.lines as mlines
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

from multiprocessing import Pool
import itertools
import matplotlib.colors as mcolors
import random

from scipy.optimize import fsolve

#from mplcursors import cursor

database = "estimates.db"
#version  = "93bc553f0661c9651147d5b8930395d6a13184e6" # "b53aff022614373a8f629961ab91813a7ff09750"  #general
version = "77de68daecd823babbb58edb1c8e14d7106e83bb" #using BDGL16
rubblue  = "#003560"
rubgreen = "#8dae10"
rubgray  = "#e7e7e7"

ymin = 1024
ymax = 2048
xmin = 30
xmax = 65

bound = 0.5
yminbound = 0

sigma = 3.19

security_levels = [150, 145, 140, 135, 128, 125, 120, 115, 110, 105, 100, 95, 90, 85, 80]
degrees = [1024, 2048]

attacks = ['usvp','bdd']
secret = "binary"
error = "standard_gaussian"

bounds = []

bounds.append((0,100))     #2**10 and degree < 2**11
bounds.append((0,200))     #2**11 and degree < 2**12
bounds.append((110,190))   #>= 2**12 and degree < 2**13
bounds.append((150,430,40))   #>= 2**13 and degree < 2**(13.5)
bounds.append((320,770,40))   #>= 2**(13.5) and degree < 2**14)
bounds.append((440,780,50))  #>= 2**14 and degree < 2**(14.5)
bounds.append((640,1600,50))  #>= 2**(14.5) and degree <= 2**15