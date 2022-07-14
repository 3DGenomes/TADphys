chrlength=20000

#for nextruders in $(seq 1 10 100);
#do
for nbarriers in 1 ; # 10 20 30 40 50 60 70 80 90 100;
do
    #dir=nextruders_${nextruders}_nbarriers_${nbarriers}
    dir=nbarriers_${nbarriers}

    if [[ -d $dir ]];
    then
	continue
    fi

    echo $nbarriers
    
    barriers=$(awk -v cl=$chrlength -v nb=${nbarriers} 'BEGIN{printf("["); db=int(cl/(nb)); for(i=db/2;i<cl-db/2; i+=db) if(i-1>0) printf("%d,", int(i));  printf("%d]",cl-db/2)}')
    echo $barriers
    
    mkdir -p ${dir}
    cd $dir
    pwd
    
    sed -e "s/XXXnextrudersXXX/${nextruders}/g" -e "s/XXXbarriersXXX/${barriers}/g" ../loop_extrusion.py > loop_extrusion_${nbarriers}.py
    python loop_extrusion_${nbarriers}.py > output_nbarriers_${nbarriers}.log &
    cd ..
    
done # Close cycle over barriers
#done # Close cycle over extruders
