#!/usr/bin/env python2
#SSchott
import os, sys, glob, re, math
import argparse, subprocess

explanation = "Script to setup lateral pressure profile calculations with NAMD from an AMBER equilibrium simulation."

parser = argparse.ArgumentParser(prog="setup_latpress", description = explanation)
parser.add_argument("-d","--dir",type=str, default="md/", help="path to directory with replicas ('01','02', and so on...).")
parser.add_argument("-c","--cpus",type=int, default=8, help="number of cpus for each individual NAMD run. Lateral pressure profiles run only on CPUs!")
parser.add_argument("-l","--length",type=int, default=1, help="number of ns on which the profile will be calculated")
parser.add_argument("-f","--freq",type=int, default=1000, help="frequency in number of steps on which the profile will be recorded. The same will be used for trajectory saving")
parser.add_argument("-s","--slabsfac",type=int, default=1, help="number of slabs per angstrom")
parser.add_argument("-n","--name",type=str, default="latpress", help="name of the slurm job. Defaults to 'latpress'. Number of the replica will be appended at the end")


def is_number(num):
    try:
        float(num)
        return True
    except:
        return False

def restrt_to_pdb(restart, topology, output="latpress.pdb"):
    cpptraj_cmd = "center ^1\nimage\ncenter !:WAT,K+,Cl-,Na+\nimage\ntrajout %s\nexit" % output
#    cpptraj_cmd = "trajout %s\nexit" % output
    with open("pdb.cpptraj","w") as f: f.write(cpptraj_cmd)
#    os.system(("ambpdb -p %s -c %s > %s") % (topology, restart, output))
    os.system(("cpptraj -p %s -y %s -i %s > cpp.log") % (topology, restart, "pdb.cpptraj"))
    return output

args = parser.parse_args()

if "AMBERHOME" not in os.environ:
    print "Load the AMBER module before proceeding"
    raise RuntimeError


dirs = glob.glob(args.dir+"/*/")
replica_paths = [d for d in dirs if is_number(d.split("/")[-2])]
cwd = os.getcwd()
namd_input1 = "npt_latpress.namd"
namd_input2 = "npt_latpress_ewald.namd"
slabs = 0

for replica in replica_paths:
    r_num = replica.split("/")[-2]
    os.chdir(replica)
    if not os.path.isdir("latpress"):
        os.mkdir("latpress")
    os.chdir("latpress")
    full_path = os.getcwd()
    slurm_script = """#!/bin/bash -x
#SBATCH --job-name=%s
#SBATCH --ntasks=%s
#SBATCH --output=%s.%%j
#SBATCH --error=%s.%%j
#SBATCH --time=24:00:00
#SBATCH --partition=rest-of-gpu1,rest-of-gpu2,rest-of-gpu3,rest-of-gpu4

source /usr/lib/Modules/init/bash
module purge
module load NAMD/git20190206
EXE=namd2
DO_PARALLEL="srun -n $SLURM_NTASKS "

WORKDIR=%s

cd $WORKDIR

INPUT=%s
$DO_PARALLEL $EXE $WORKDIR/$INPUT > ${INPUT}.out
INPUT=%s
$DO_PARALLEL $EXE $WORKDIR/$INPUT > ${INPUT}.out


exit
""" % (args.name+r_num, args.cpus, args.name, args.name, full_path, namd_input1, namd_input2)
    slurm_file = open("run_latpress.pbs","w")
    slurm_file.write(slurm_script)
    slurm_file.close()

    if os.path.isfile("../prod/run_prod.pbs"):
        with open("../prod/run_prod.pbs") as f:
            for l in f:
                if re.match("PRMTOP", l):
                    topology = l.split("=")[-1].strip()
                    break
    else:
        raise IOError("No run_prod.pbs file in prod folder")

    restart = sorted(glob.glob("../prod/*.restrt"))[-1]

    latpress_pdb = restrt_to_pdb(restart, topology, output="latpress_%s.pdb" % r_num)

    with open(latpress_pdb) as f:
        for l in f:
            if l.startswith("CRYST1"):
                v = l.strip().split()
                x,y,z = map(float,v[1:4])
                break

    if slabs == 0:
        slabs = math.ceil(z)/args.slabsfac

    namd_script1 = """amber              on
parmfile           %s
coordinates        %s
readexclusions     on	                            ;# Read exclusions just in case
exclude        	   scaled1-4                        ;# AMBER style
1-4scaling     	   0.833333                         ;# SCEE value. 1.2 default in AMBER 1-4scaling=1/SCEE
switching          off                              ;# OFF in AMBER

cutoff             10                               ;# as used before
pairListDist       12                               ;# cutoff+skinnb (2 by default in AMBER)
LJcorrection       on                               ;# correct tailing values after cutoff
zeroMomentum       on                               ;# removes PME drift. Experimental if restraints?

rigidBonds         all                              ;# SHAKE it!
useSettle          on                               ;# SETTLE for waters
rigidTolerance     0.00001                          ;# AMBER tolerance

set temperature    300                              ;# target temperature used several times below

# starting from scratch
temperature        $temperature                     ;# initialize velocities randomly

outputName         npt_latpress   ;# base name for output from this run

restartfreq        500000     ;# 500 steps = every 1ps
dcdfreq                %i
xstFreq             10000

outputEnergies     10000     ;# 10000 steps = every 20 ps
outputTiming       1000      ;# shows time per step and time to completion every 2 ps

stepspercycle       10   ;# redo pairlists every ten steps

# Integrator Parameters
timestep            2.0  ;# 2fs/step
nonbondedFreq       1    ;# nonbonded forces every step
fullElectFrequency  2    ;# PME only every other step

# Constant Temperature Control
langevin            on            ;# langevin dynamics
langevinDamping     1.            ;# damping coefficient of 1/ps
langevinTemp        $temperature  ;# random noise at this level
langevinHydrogen    no            ;# don't couple bath to hydrogens

# Periodic Boundary conditions  
cellBasisVector1    %f  0.    0.  ;# vector to the next image
cellBasisVector2     0.   %f  0.
cellBasisVector3     0.     0.   %f
cellOrigin           0.     0.    0.  ;# the *center* of the cell

#PME (for full-system periodic electrostatics)
PME                 yes

# let NAMD determine grid
PMEGridSpacing      1.0
PMETolerance        1E-5           ;# AMBER tolerance

# Constant Pressure Control (variable volume)
useGroupPressure      yes                           ;# needed for rigid bonds
useFlexibleCell       yes                           ;# allow changes in every axis no for water box, yes for membrane
useConstantRatio       no                           ;# keep x/y // needs useflexiblecell
useConstantArea        no                           ;# keeps xy fixed. As in amber, not working with berendsen

BerendsenPressure                   on
BerendsenPressureTarget             1.              ;# 1 bar
BerendsenPressureRelaxationTime     1000            ;# taup=1
BerendsenPressureCompressibility    4.46E-5         ;# AMBER comp default
BerendsenPressureFreq               10

pressureProfile                     on              ;#
pressureProfileSlabs                %i
pressureProfileFreq                 %d              ;# Every 10th step?

reinitvels          $temperature   ;# restarting from AMBER // assume zeros velocities
run %d                          """ % (topology, latpress_pdb, args.freq, x, y, z, slabs, args.freq, args.length*500000)


    namd_script2 = """amber              on
parmfile           %s
coordinates        %s
readexclusions     on	                            ;# Read exclusions just in case
exclude        	   scaled1-4                        ;# AMBER style
1-4scaling     	   0.833333                         ;# SCEE value. 1.2 default in AMBER 1-4scaling=1/SCEE
switching          off                              ;# OFF in AMBER

cutoff             10                               ;# as used before
pairListDist       12                               ;# cutoff+skinnb (2 by default in AMBER)
LJcorrection       on                               ;# correct tailing values after cutoff
zeroMomentum       on                               ;# removes PME drift. Experimental if restraints?

rigidBonds         all                              ;# SHAKE it!
useSettle          on                               ;# SETTLE for waters
rigidTolerance     0.00001                          ;# AMBER tolerance

set temperature    300                              ;# target temperature used several times below

# starting from scratch
temperature        $temperature                     ;# initialize velocities randomly

outputName         npt_latpress_ewald   ;# base name for output from this run

restartfreq        500000     ;# 500 steps = every 1ps
dcdfreq                %i
xstFreq             10000

outputEnergies     10000     ;# 10000 steps = every 20 ps
outputTiming       1000      ;# shows time per step and time to completion every 2 ps

stepspercycle       10   ;# redo pairlists every ten steps

# Integrator Parameters
timestep            2.0  ;# 2fs/step
nonbondedFreq       1    ;# nonbonded forces every step
fullElectFrequency  2    ;# PME only every other step

# Constant Temperature Control
langevin            on            ;# langevin dynamics
langevinDamping     1.            ;# damping coefficient of 1/ps
langevinTemp        $temperature  ;# random noise at this level
langevinHydrogen    no            ;# don't couple bath to hydrogens

# Periodic Boundary conditions  
cellBasisVector1    %f  0.    0.  ;# vector to the next image
cellBasisVector2     0.   %f  0.
cellBasisVector3     0.     0.   %f
cellOrigin           0.     0.    0.  ;# the *center* of the cell

#PME (for full-system periodic electrostatics)
PME                 yes

# let NAMD determine grid
PMEGridSpacing      1.0
PMETolerance        1E-5           ;# AMBER tolerance

# Constant Pressure Control (variable volume)
useGroupPressure      yes                           ;# needed for rigid bonds
useFlexibleCell       yes                           ;# allow changes in every axis no for water box, yes for membrane
useConstantRatio       no                           ;# keep x/y // needs useflexiblecell
useConstantArea        no                           ;# keeps xy fixed. As in amber, not working with berendsen

BerendsenPressure                   on
BerendsenPressureTarget             1.              ;# 1 bar
BerendsenPressureRelaxationTime     1000            ;# taup=1
BerendsenPressureCompressibility    4.46E-5         ;# AMBER comp default
BerendsenPressureFreq               10

pressureProfile                     on              ;#
pressureProfileSlabs                %i
pressureProfileFreq                 %d              ;# Every 10th step?

pressureProfileEwald  on
pressureProfileEwaldX  16
pressureProfileEwaldY  16
pressureProfileEwaldZ  16

set ts 0
firstTimestep $ts

coorfile open dcd npt_latpress.dcd
while { [coorfile read] != -1 } {
  incr ts %i
  firstTimestep $ts
  run 0
}
coorfile close
                  """ % (topology, latpress_pdb, args.freq, x, y, z, slabs, args.freq, args.freq)

    namd_file = open(namd_input1, "w")
    namd_file.write(namd_script1)
    namd_file.close()

    namd_file = open(namd_input2, "w")
    namd_file.write(namd_script2)
    namd_file.close()

    subprocess.check_output('/bin/bash -c "sbatch run_latpress.pbs"', shell=True)

    os.chdir(cwd)
