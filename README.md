Guidance for Efficient Selection of Secure
Parameters for FHE
=======================================

This repository implements the formulas of the paper [Guidance for Efficient Selection of Secure Parameters for Fully Homomorphic Encryption](https://eprint.iacr.org/2024/1001).

Usage (Parameter estimation)
-----
Find an estimation of the security level by running:
   ````
   python3 estimate.py --param "lambda" --n "1024" --logq "20-30;35;40-60" --secret "binary" --error "3.19"
   ````
Find an estimation of the LWE dimension required to obtain a given security level:
 ````
 python3 estimate.py --param "n" --lambda "80" --logq "20-30" --secret "binary" --error "3.19"
 ````
Find an estimation of the size of the modulus q:
 ````
 python3 estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19"
 ````
Find an estimation of the standard deviation of the error distribution:
````
python3 estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --secret "binary" --error "3.19"
 ````

Note: you can add the option ````--verify 1```` to any of the commands to compare the output of the formulas against the Lattice Estimator (see Dependencies).

Usage (Fitting function)
-----
We include the file fit_formula.py which computes the best coefficients to fit our formulas with the output of the lattice estimator. 
To use it, you do the following:

   ````
   python3 fit_formula.py --param A --attack B --dist C --simpl D
   ````
where

A is in {'lambda','n'}, B is in {'bdd', 'usvp'}, C is in {'binary', 'ternary'} and D is in {0,1}.

For example, 

   ````
   python3 fit_formula.py --param 'n' --attack 'usvp' --dist 'ternary' --simpl 0
   ````
will find the best parameters for the LWE dimension n, considering the usvp attack, ternary distribution and the formula containing beta. If we set --simpl 1
we consider the formula where the dependency on beta has been removed. 

Dependencies
---------
We have added the functionality to compare the output of our formulas against the [Lattice Estimator](https://github.com/malb/lattice-estimator). 
Please download the Estimator if you want to use such functionality. 
Note: At present, the Estimator is also needed to run one of the formulas, this will be fixed shortly. 

[Numpy](https://numpy.org/) is required.
     
Bugs
----

Please report bugs through the [GitHub issue tracker](https://github.com/sergirovira/fastparameterselection/issues)

Citing
------

The paper was presented at Africacrypt 2024:

 ````
 @inproceedings{kirshanova2024guidance,
  			title={Guidance for efficient selection of secure parameters for fully homomorphic encryption},
  			author={Kirshanova, Elena and Marcolla, Chiara and Rovira, Sergi},
  			booktitle={International Conference on Cryptology in Africa},
  			pages={376--400},
  			year={2024},
  			organization={Springer}
			}
 ````
  
  			
A pre-print is available as
	 
	 	
 | Cryptology ePrint Archive, Report [2024/1001](https://eprint.iacr.org/2024/1001), 2024.

