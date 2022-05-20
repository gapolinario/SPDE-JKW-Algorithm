#!/bin/bash

run="python3 algorithm.py"

numcores=8
prefix=0

size=16
ensemble_init=$(($prefix*1000))
ensemble_finl=$(($prefix*1000 + $size -1))

R=`seq $ensemble_init $ensemble_finl`
N=( 7 8 9 10 )
dtconst=( 0.02 )
nu=( 0.01 )

parallel -P $numcores $run '{1} {2} {3} {4}' ::: "${R[@]}" ::: "${N[@]}" ::: "${dtconst[@]}" ::: "${nu[@]}"
