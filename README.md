Guidance for Efficient Selection of Secure
Parameters for FHE
=======================================

This repository implements the formulas of the paper [Guidance for Efficient Selection of Secure Parameters for Fully Homomorphic Encryption](https://eprint.iacr.org/2024/1001).

Usage
-----
You can use the formulas to find an estimation of the security level by running:
   ````
  python3 estimate.py --param "lambda" --n "1024" --logq "20-30;35;40-60" --dist "ternary"
   ````
  or
  ````
  python3 estimate.py --param "lambda" --file "example_lambda_ternary.csv"
 ````
You can use the formulas to find an estimation of the LWE dimension required to obtain a given security level by running:
 ````
 python3 estimate.py --param "n" --lambda "80" --logq "20-30" --dist "binary" 
 ````
or 
 ````
 python3 estimate.py --param "n" --file "example_n_binary.csv"
 ````
Note: you can add the option ````--verify 1```` to any of the commands to compare the output of the formulas against the Lattice Estimator (see Dependencies).

Dependencies
---------
We have added the functionality to compare the output of our formulas against the [Lattice Estimator](https://github.com/malb/lattice-estimator). 
Please download the Estimator if you want to use such functionality. 
Note: At present, the Estimator is also needed to run one of the formulas, this will be fixed shortly. 

[Numpy](https://numpy.org/) is required.

TODOs
---------

We will address the following tasks in the upcoming weeks:

- ``[ ]`` Allow the user to choose the distribution of the error. At the moment is set to 3.19 as default. 
         
Bugs
----

Please report bugs through the [GitHub issue tracker](https://github.com/sergirovira/fastparameterselection/issues)

Citing
------

The paper will be presented at Africacrypt 2024, before the proceedings are published, a pre-print is available as

    | Cryptology ePrint Archive, Report 2024/1001, 2024. https://eprint.iacr.org/2024/1001
