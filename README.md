A Tool for Fast and Secure LWE Parameter
Selection
=======================================
In this repository, we offer a tool to select secure parameters for LWE-based applications in a fast and flexible way. The tool can provide you with any of the following parameters: security level, size of the ciphertext modulus, LWE dimension and standard deviation of the error distribution. 

Our tool is constructed by studying the uSVP and BDD attacks against LWE. From this study, we derive formulas which describe each of the aforementioned parameters as a function of the others. You can find all the details in this paper [A Tool for Fast and Secure LWE Parameter Selection: the FHE case](https://eprint.iacr.org/2024/1895). Section 6 of the paper comprehensively explains the tool's usage, which complements this readme.  

Usage (Parameter estimation)
-----
Find an estimation of the security level by running:
   ````
   python3 src/estimate.py --param "lambda" --n "1024" --logq "20-30;35;40-60" --secret "binary" --error "3.19"
   ````
Find an estimation of the LWE dimension required to obtain a given security level:
 ````
 python3 src/estimate.py --param "n" --lambda "80" --logq "20-30" --secret "binary" --error "3.19"
 ````
Find an estimation of the size of the modulus q:
 ````
 python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19"
 ````
Find an estimation of the standard deviation of the error distribution:
````
python3 src/estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --secret "binary" --error "3.19"
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
[SciPy](https://scipy.org/) is required. 

You can run the following command to install:
````
pip install -r requirements.txt
````

Use with Docker
-------
You can build and run the repository with [Docker](...) with the command
````
docker-compose -f ./docker/docker-compose.yaml up --build
````
To only run, use
````
docker-compose -f ./docker/docker-compose.yaml up
````

Currently, it runs estimate.py to obtain the parameter lambda, given n = 1024, logq = 35, binary secret distribution and standard deviation of the error distribution 3.19. To run the estimation with your parameters, you can modify the command line in docker-compose.yaml as follows.
Find an estimation of the security level:
````
command: [ "sage", "--python3", "src/estimate.py", "--param", "lambda",  "--n", "1024", "--logq", "20-30\\;35\\;40-60", "--secret", "binary", "--error", "3.19"]
````
Find an estimation of the security level and verify it against the Lattice Estimator:
````
command: [ "sage", "--python3", "src/estimate.py", "--param", "lambda",  "--n", "1024", "--logq", "20-30\\;35\\;40-60", "--secret", "binary", "--error", "3.19", "--verify", "1" ]
````
Find an estimation of the LWE dimension:
````
command: ["sage", "--python3", "src/estimate.py", "--param", "n",  "--lambda", "80", "--logq", "20", "--secret", "binary", "--error", "3.19"]
````
Find an estimation of the size of the modulus q:
````
command: ["sage", "--python3", "src/estimate.py", "--param", "logq",  "--lambda", "80", "--n", "1024", "--secret", "binary", "--error", "3.19"]
````
Find an estimation of the standard deviation of the error distribution:
````
command: ["sage", "--python3", "src/estimate.py", "--param", "std_e",  "--lambda", "80", "--n", "1024", "--logq", "20", "--secret", "binary"]
````
Find an estimation of the security level, given the example parameters in example_lambda_binary.csv: 
````
command: ["sage", "--python3", "src/estimate.py", "--param", "lambda", "--file", "./examples/example_lambda_binary.csv", "--verify", "1"]
````
Find an estimation of the LWE dimension, given the example parameters in example_n_ternary.csv:
````
command: ["sage", "--python3", "src/estimate.py", "--param", "n", "--file", "./examples/example_n_ternary.csv"]
````
Find an estimation of the size of the modulus q, given the example parameters in example_logq_binary.csv:
````
command: ["sage", "--python3", "src/estimate.py", "--param", "logq", "--file", "./examples/example_logq_binary.csv"]
````
Find an estimation of the standard deviation of the error distribution, given the example parameters in example_error_binary.csv:
````
command: ["sage", "--python3", "src/estimate.py", "--param", "error", "--file", "./examples/example_error_binary.csv"]
````

Note: in the Docker version, we applied a change to the [Lattice Estimator](...), to address an imprecision in the case where the standard deviation of the error distribution is much larger than 3.19.
     
Bugs
----

Please report bugs through the [GitHub issue tracker](https://github.com/sergirovira/fastparameterselection/issues)

Citing
------
 ````
@misc{cryptoeprint:2024/1895,
      author = {Beatrice Biasioli and Elena Kirshanova and Chiara Marcolla and Sergi Rovira},
      title = {A Tool for Fast and Secure {LWE} Parameter Selection: the {FHE} case},
      howpublished = {Cryptology {ePrint} Archive, Paper 2024/1895},
      year = {2024},
      url = {https://eprint.iacr.org/2024/1895}
}
 ````
  
The paper associated with our tool is a follow-up and extension of the following paper presented at Africacrypt 2024.   			
The pre-prints of both papers are available at
	 
 | Cryptology ePrint Archive, Report [2024/1001](https://eprint.iacr.org/2024/1001), 2024. 	
 | Cryptology ePrint Archive, Report [2024/1895](https://eprint.iacr.org/2024/1895), 2024.

