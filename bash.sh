#!/bin/bash
#Runs simulations for GED model

runDir="/mnt/e/runs/"
baseName="squareLattice/R1/N4900/moreExchanges"

N=4900
f=0.1
mu=0.1

sq=1
er=0
sf=0

ring=0

sfc=0

baseDir=$runDir$baseName
mainDir=`pwd`
echo $mainDir
echo $baseDir

echo -n `date` >> runLog.txt
echo -e -n "\t" >> runLog.txt
echo $baseName >> runLog.txt

for log10oneminuslambda in 0 -0.5 -1 -1.5 -2 -2.5 -3
#for clone in {1..1000}
do
    #log10oneminuslambda=-3
    lambda=`echo "scale=10; 1-e(l(10)*$log10oneminuslambda)" | bc -l`
    ##lambda=$log10oneminuslambda
    #echo $lambda
    path=$baseDir/l$log10oneminuslambda
    #path=$baseDir/$clone/l$log10oneminuslambda

    nano=`date +%N`
    randChoice=$((${nano#0} % 5001))
    echo $randChoice
    loadWealths="squareLattice/R1/N4900/l-3/data/wealthDistribution-"
    txt=".txt"
    loadWealths=$runDir$loadWealths$randChoice$txt


    mkdir -p $path
    cp exchange.h $path
    cp exchange.cpp $path
    cp main.cpp $path
    cp read.in $path
    cp bash.sh $path
    cd $path
    mkdir data

    #Now we need to insert all of our variable values into our copy of read.in:
    sed -i -e 's/?p?/0.0/g' read.in
    sed -i -e 's/?lambda?/'$lambda'/g' read.in
    sed -i -e 's/?T?/1000000/g' read.in
    sed -i -e 's/?sf?/'$sf'/g' read.in
    sed -i -e 's/?a?/1.5/g' read.in
    sed -i -e 's/?er?/'$er'/g' read.in
    sed -i -e 's/?r?/'`echo "scale=10; 4/$N" | bc`'/g' read.in
    sed -i -e 's/?N?/'$N'/g' read.in
    sed -i -e 's/?mu?/'$mu'/g' read.in
    sed -i -e 's/?breakRate?/0.0/g' read.in
    sed -i -e 's/?minDegree?/2.0/g' read.in
    sed -i -e 's/?intensity?/0.0/g' read.in
    sed -i -e 's/?noise?/0.0/g' read.in
    sed -i -e 's/?an?/0.0/g' read.in
    sed -i -e 's/?doFit?/0/g' read.in
    sed -i -e 's/?doEG?/0/g' read.in
    sed -i -e 's/?hubify?/0/g' read.in
    sed -i -e 's/?cml?/0/g' read.in
    sed -i -e 's/?sq?/'$sq'/g' read.in
    sed -i -e 's/?rsq?/1.0/g' read.in
    sed -i -e 's/?f?/'$f'/g' read.in
    sed -i -e 's/?nzo?/0/g' read.in
    sed -i -e 's/?ring?/'$ring'/g' read.in
    sed -i -e 's/?rring?/1.0/g' read.in
    sed -i -e 's/?sfc?/'$sfc'/g' read.in
    sed -i -e 's/?tde?/0/g' read.in
    sed -i -e 's/?halfk?/2/g' read.in
    sed -i -e 's/?doLearn?/0/g' read.in
    sed -i -e 's/?lS?/3.0/g' read.in
    sed -i -e 's/?nuke?/'$N'/g' read.in
    sed -i -e 's|?lW?|none|g' read.in
    sed -i -e 's/?seed?/0/g' read.in
    sed -i -e 's/?nData?/5000/g' read.in
    sed -i -e 's/?mX?/1/g' read.in


    g++ -O3 exchange.h exchange.cpp main.cpp -o a.out

    #echo we compiled! Time to run! 
    ./a.out >> log.txt
    cd $mainDir
done
