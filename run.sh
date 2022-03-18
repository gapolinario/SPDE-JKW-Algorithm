numcores=8
prefix=8

ensemble_init=$(($prefix*1000))
ensemble_finl=$(($prefix*1000+95))

seq $ensemble_init $ensemble_finl | xargs -t -I {} -P $numcores python3 algorithm.py {}
