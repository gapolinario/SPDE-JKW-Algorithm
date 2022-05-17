numcores=8
prefix=0

size=16
ensemble_init=$(($prefix*1000))
ensemble_finl=$(($prefix*1000 + $size -1))

N=8
dtconst=0.1
nu=0.01

seq $ensemble_init $ensemble_finl | xargs -t -I {} -P $numcores python3 algorithm.py {}
