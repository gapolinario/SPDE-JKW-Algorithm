#!/bin/bash

run="python3 statistics.py"

numcores=8
prefix=(0)

N=( 7 8 9 10 )
dtconst=( 1. )
nu=( 0.01 )

parallel -P $numcores $run '{1} {2} {3} {4}' ::: "${prefix[@]}" ::: "${N[@]}" ::: "${dtconst[@]}" ::: "${nu[@]}"
