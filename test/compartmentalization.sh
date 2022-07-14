for epsilon in $(seq 0.0 0.1 1.0);
do

    mkdir -p epsilon_${epsilon}kBT
    cd epsilon_${epsilon}kBT
    pwd

    sed "s/XXXepsilonXXX/${epsilon}/g" ../compartmentalization.py | python > output_${epsilon}kBT.log &
    
    cd ..
    
done
