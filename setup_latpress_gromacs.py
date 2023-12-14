#!/usr/bin/env python2
#SSchott
import os, sys, glob, re, math
import argparse, subprocess

explanation = "Script to setup lateral pressure profile calculations with AMBER and postprocessing with GROMACS."

parser = argparse.ArgumentParser(prog="setup_latpress", description = explanation)
parser.add_argument("-d","--dir",type=str, default="md/", help="path to directory with replicas ('01','02', and so on...).")
parser.add_argument("-l","--length",type=int, default=100, help="number of ns to run. Devs recommend 100ns")
parser.add_argument("-c","--chunk",type=int, default=10, help="chunk size in ns")
parser.add_argument("-f","--freq",type=int, default=2500, help="frequency in number of steps on which the trajectory will be recorded. Devs recommend 5 to 10 ps")
parser.add_argument("-n","--name",type=str, default="G-LS", help="name of the slurm job. Defaults to 'G-LS'. Number of the replica will be appended at the end")


def is_number(num):
    try:
        float(num)
        return True
    except:
        return False

args = parser.parse_args()

if "AMBERHOME" not in os.environ:
    print "Load the AMBER module before proceeding"
    raise RuntimeError


dirs = glob.glob(args.dir+"/*/")
replica_paths = [d for d in dirs if is_number(d.split("/")[-2])]
cwd = os.getcwd()
queue_input = "run_vel.pbs"
md_input = "md_npt_vel_0001.in"

for replica in replica_paths:
    r_num = replica.split("/")[-2]
    os.chdir(replica)
    if not os.path.isdir("g_lstress"):
        os.mkdir("g_lstress")
    os.chdir("g_lstress")
    full_path = os.getcwd()

    if os.path.isfile("../prod/run_prod.pbs"):
        with open("../prod/run_prod.pbs") as f:
            for l in f:
                if re.match("PRMTOP", l):
                    topology = l.split("=")[-1].strip()
                    break
    else:
        raise IOError("No run_prod.pbs file in prod folder")

    restart = sorted(glob.glob("../prod/*.restrt"))[-1]
    slurm_script = """#!/bin/bash -x
#SBATCH --job-name=%s
#SBATCH --output=%s.%%j
#SBATCH --error=%s.%%j
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu-ser
#SBATCH --gres=gpu:1

#
# --- Prepare calculation
#

echo "------------------------------------------------------------"
echo "Job-ID               = $SLURM_JOB_ID"
echo "Host                 = $SLURMD_NODENAME"
echo "Nodelist             = $SLURM_JOB_NODELIST"
echo "Nodes                = $SLURM_JOB_NUM_NODES"
echo "GPUs                 = $SLURM_GPUS"
echo "Submission Directory = $SLURM_SUBMIT_DIR"
echo "------------------------------------------------------------" 

module load AMBER/20.1
EXE=pmemd.cuda
DO_PARALLEL="srun -n 1 "

#################################
# --- DEFINE YOUR VARIABLES --- # 
#################################
#
WORKDIR=%s
cd $WORKDIR

RESOLD=%s 	#OLD RESTRT
PRMTOP=%s 	#TOPOLOGY
BASE=%s 	#BASENAME
SCRIPT=%s  	#PATH TO THIS SCRIPT

START=1
echo "Start is "$START
END=$(($START + 1))
echo "End is "$END
MAXREP=%s #NUMBER OF CYCLES

PRMTOP=$(readlink -f $PRMTOP)

##########################################################
# --- HERE SHOULD NOT BE ANY NEED FOR CHANGES BELOW --- #
##########################################################

# --- Init for loop
#

COUNT=$START
COUNTEXP=$(printf %%04d $COUNT)

# --- Save copy of this script
#

cp $SCRIPT $SCRIPT"_"$COUNTEXP

#
# --- Maximum count
#
if [ $START -gt $MAXREP ]
then
  exit
fi

#
# --- Processing
#

TMPFOLDER=$(mktemp -d)
cd $TMPFOLDER

$DO_PARALLEL $EXE -O -i ${WORKDIR}/$BASE"_"$COUNTEXP.in -o $BASE"_"$COUNTEXP.out -c ${WORKDIR}/${RESOLD} -p $PRMTOP -r $BASE"_"$COUNTEXP.restrt -x $BASE"_"$COUNTEXP.mdcrd

if [ ! -f $BASE"_"$COUNTEXP.restrt ]
then
echo "Restart file $BASE"_"$COUNTEXP.restrt not created!"
exit
fi

cd $WORKDIR
cp ${TMPFOLDER}/* .
rm -rf $TMPFOLDER

#UNCOMMENT IF YOU WANT GZIP
#echo "Zipping mdcrd file"
#gzip -9 -f $BASE"_"$COUNTEXP.mdcrd

###########################
### Prepare for new run ###
###########################

# --- Modify input file
#
OLD=$BASE"_"$COUNTEXP".in"

COM1="s/^RESOLD=.*/RESOLD="$BASE"_"$COUNTEXP".restrt""/"

COUNT=$END
COUNTEXP=$(printf %%04d $COUNT)

cp $OLD $BASE"_"$COUNTEXP".in"
NEW_INPUT=$BASE"_"$COUNTEXP".in"

echo Modifying $SCRIPT
COM2="s/START="$START"/START="$COUNT"/"

cat $SCRIPT | sed -e "$COM1" | sed -e "$COM2" > tmp
mv tmp $SCRIPT

#
# --- Re-queue
#
sbatch $SCRIPT
# --- Exit this script
#
exit
""" % (args.name+r_num, args.name, args.name, full_path, restart, topology, md_input[:-8], full_path+"/"+queue_input, args.length/args.chunk)
    slurm_file = open(queue_input,"w")
    slurm_file.write(slurm_script)
    slurm_file.close()


    mdin_script = """NTP production with velocities using pmemd.cuda
 &cntrl
  ntx = 5,
  irest = 1,
  ntpr = %s,
  ntwr = %s,
  iwrap = 1,
  ntwx = %s,
  ntwv = -1,

  nstlim = %s,
  nscm = 0,
  dt = 0.002,

  ntt = 3,
  tempi = 300.0,
  ig = -1,
  gamma_ln = 1,

  ntp = 3,
  csurften = 3,

  ntc = 2,

  ntf = 2,
  cut = 10.0,
 

  nmropt = 0,

 &end

""" % (args.freq, 500000*args.chunk, args.freq, 500000*args.chunk)
    mdin_file = open(md_input, "w")
    mdin_file.write(mdin_script)
    mdin_file.close()

    subprocess.check_output('/bin/bash -c "sbatch '+queue_input+'"', shell=True)

    os.chdir(cwd)
