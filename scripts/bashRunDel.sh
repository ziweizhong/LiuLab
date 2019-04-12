#!bin/bash

#$ -N CAcount-ALL
#$ -q pub64
#$ -pe openmp 28
#$ -R y

module load anaconda/3.7-5.3.0
cp /pub/ziweiz2/LiuLab/scripts/$1 . 
python $@
rm ./$1
