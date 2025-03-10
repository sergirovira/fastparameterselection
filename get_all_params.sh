#!/bin/bash

python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "binary" --error "3.19" --simpl 0
python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "binary" --error "3.19" --simpl 0
python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "binary" --error "3.19" --simpl 1
python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "binary" --error "3.19" --simpl 1

python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 0
#python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 0 #breaks approximating beta since n < lnq at some point
python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 1
python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 1 #it works if we fix the value of P1 to always be positive

python3 fit_formula.py --param "n" --attack 'usvp' --secret "binary" --error "3.19" --simpl 0
python3 fit_formula.py --param "n" --attack 'bdd' --secret "binary" --error "3.19" --simpl 0 #Need to verify where to place the parameters
#python3 fit_formula.py --param "n" --attack 'usvp' --secret "binary" --error "3.19" --simpl 1 #Not defined
python3 fit_formula.py --param "n" --attack 'bdd' --secret "binary" --error "3.19" --simpl 1

python3 fit_formula.py --param "n" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 0
python3 fit_formula.py --param "n" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 0 #Need to verify where to place the parameters
python3 fit_formula.py --param "n" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 1 #not defined
python3 fit_formula.py --param "n" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 1