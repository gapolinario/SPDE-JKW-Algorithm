numcores=8
prefix=1

size=16
ensemble_init=$(($prefix*1000))
ensemble_finl=$(($prefix*1000+$size))

seq $ensemble_init $ensemble_finl | xargs -t -I {} -P $numcores python3 algorithm.py {}

#seq $ensemble_init $ensemble_finl | xargs -t -I {} -P $numcores python3 algorithm_sde.py {}
