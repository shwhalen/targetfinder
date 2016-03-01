#!/bin/sh

./generate_enhancers.py $1
./generate_promoters.py $1
./generate_pairs.py $1
./generate_training.py $1
