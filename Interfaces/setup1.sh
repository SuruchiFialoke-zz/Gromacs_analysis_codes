Nrange=$(seq 2400 80 4800)

for i in $Nrange
do
    mkdir N$i
    
    sed "s/NNN/N$i/g" zcut_fixed_pwII.cpp > N$i/zcut_fixed_pwII_N$i.cpp
    sed "s/NNN/N$i/g" submit.sh > N$i/submit1.sh
    cp miscfun.h N$i
    cp Marching_Cube.h N$i
 
    cd N$i
    wdir=`pwd`
    sed "s|WORKDIR|$wdir|g" < submit1.sh >submit.sh
    rm submit1.sh
    qsub submit.sh
    cd ..

done
