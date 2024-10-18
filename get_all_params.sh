#!/bin/bash

python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "binary" --error "3.19" --simpl 0
python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "binary" --error "3.19" --simpl 0
python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "binary" --error "3.19" --simpl 1
python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "binary" --error "3.19" --simpl 1

#python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 0
#python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 0
#python3 fit_formula.py --param "lambda" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 1
#python3 fit_formula.py --param "lambda" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 1

python3 fit_formula.py --param "n" --attack 'usvp' --secret "binary" --error "3.19" --simpl 0
#python3 fit_formula.py --param "n" --attack 'bdd' --secret "binary" --error "3.19" --simpl 0
#python3 fit_formula.py --param "n" --attack 'usvp' --secret "binary" --error "3.19" --simpl 1
python3 fit_formula.py --param "n" --attack 'bdd' --secret "binary" --error "3.19" --simpl 1

python3 fit_formula.py --param "n" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 0
#python3 fit_formula.py --param "n" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 0
#python3 fit_formula.py --param "n" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 1
python3 fit_formula.py --param "n" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 1