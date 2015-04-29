#! /bin/bash
source ~/.bashrc
gmx455
lpy
openmm601

function runSim {
simdir=$1
syst=$2
datadir=$3

    if [ ! -d $simdir ]; then
        mkdir $simdir
    fi
    
    cd $simdir
    grompp -f ${datadir}/sim0.mdp -c ${datadir}/${syst}.gro -p ${datadir}/${syst}.top -o sim.tpr
    mdrun -deffnm sim -v
    echo 0 | trjconv -f sim.trr -s sim.tpr -o bondedPout.pdb
    printf "Potential\n0\n" | g_energy -f sim.edr -o potential.xvg
    cd ..
}


function runOpenmm {
    ./simOpenmm.py



}

function main {
    datadir="${PWD}/data"
    rundir="sim"
    syst="xqCG"
    runSim $rundir $syst $datadir
    runOpenmm
}

main

