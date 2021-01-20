# 0 - Relaxation run of 10MLN  steps to try obliterating the initial rosettes;
# 1 - Relaxation run of 100MLN   steps saving every 1000 steps.

for replica in $(seq 2 1 2);
do
    replicadir=replica_${replica}
    mkdir -p ${replicadir}
    cd ${replicadir}
    sed -e "s/XXXreplicaXXX/${replica}/g" ../relaxation.py > relaxation_replica_${replica}.py
    mpirun -np 16 python relaxation_replica_${replica}.py > output_replica_${replica}.log &
    cd ..
done # Close cycle over replica

