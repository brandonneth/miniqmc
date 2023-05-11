#!/bin/bash

# To be run from build directory

for x in 1 2 3; do
for y in 1 2 3; do
for z in 1 2 3; do
./bin/miniqmc -g "$x $y $z"
done
done
done

