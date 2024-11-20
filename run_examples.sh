#!/bin/bash

python3 estimate.py --param "lambda" --file "example_lambda_binary.csv"

python3 estimate.py --param "lambda" --file "example_lambda_binary.csv" --verify 1

python3 estimate.py --param "lambda" --n "1024" --logq "20-30;35;40-60" --secret "binary" --error "3.19"

python3 estimate.py --param "n" --lambda "80" --logq "20-30" --secret "binary" --error "3.19" --verify 1

python3 estimate.py --param "n" --file "example_n_ternary.csv"

python3 estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19"

python3 estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --secret "binary"

