set -x
sh revision.sh $1 1 0 1
sh revision.sh $1 1 0.1 1
sh revision.sh $1 1 0.05 1
sh revision.sh $1 1 0.01 1
sh revision.sh $1 1 0.005 1
sh revision.sh $1 1 0.001 1
sh revision.sh $1 1 0.0005 1
sh revision.sh $1 1 0.0001 1

sh revision.sh $1 2 0 0.5
sh revision.sh $1 5 0 0.2
sh revision.sh $1 10 0 0.1
sh revision.sh $1 20 0 0.05
sh revision.sh $1 33 0 0.03
sh revision.sh $1 50 0 0.02
sh revision.sh $1 100 0 0.01
