#!/usr/bin/env python2
#SSchott
import os, sys, glob, re, math
import argparse, subprocess
import parmed as pmd

explanation = "Script to setup lateral pressure profile postprocessing with GROMACS."

parser = argparse.ArgumentParser(prog="setup_latpress", description = explanation)
parser.add_argument("-d","--dir",type=str, default="md/", help="path to directory with replicas ('01','02', and so on...).")
parser.add_argument("-c","--chunk",type=int, default=250, help="chunk size for GROMACS stress in ps. Number of frames depends on MD run.")
parser.add_argument("-f","--freq",type=int, default=2500, help="frequency in number of steps on which the trajectory will be recorded. Devs recommend 5 to 10 ps")
parser.add_argument("-n","--name",type=str, default="G-LS", help="name of the slurm job. Defaults to 'G-LS'. Number of the replica will be appended at the end")
parser.add_argument("--center", default="!r WAT&!r K+", help="selection for centering atoms with trjconv. Defaults to not WAT not K+")


def is_number(num):
    try:
        float(num)
        return True
    except:
        return False


def amb2gro(topo, crd):
    amber = pmd.load_file(topo, crd)
    new_topo = (".").join(topo.split(".")[:-1])+".gromacs.top"
    new_gro  = (".").join(topo.split(".")[:-1])+".gromacs.gro"
    if not os.path.exists(new_topo):
        amber.save(new_topo)
    if not os.path.exists(new_gro):
        amber.save(new_gro)
    return new_topo, new_gro

args = parser.parse_args()

if "AMBERHOME" not in os.environ:
    print("Load the AMBER module before proceeding")
    raise RuntimeError

if "GMXBIN" not in os.environ:
    print("Load the GROMACS-LS module before proceeding")
    raise RuntimeError


dirs = glob.glob(args.dir+"/*/")
replica_paths = [d for d in dirs if is_number(d.split("/")[-2])]
cwd = os.getcwd()
mdp_input = "stress.mdp"

for replica in replica_paths:
    r_num = replica.split("/")[-2]
    print("Working on replica: %s" % (r_num))
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

    restart = sorted(glob.glob("*.restrt"))[-1]

    topo, gro = amb2gro(topology, restart)

    trr_out = "md_npt.trr"
    if not os.path.isfile(trr_out) and not glob.glob(trr_out.replace(".trr","*?.trr")):
        print("Concatenating traj into TRR file")
        cmd1="""srun --mem=5G -p rest-of-gpu1,rest-of-gpu2,rest-of-gpu3,rest-of-gpu4 cpptraj -p %s -y %s << EOF > cpp1.log
center ^1
image
center :OL
image
center :OL
image
trajout %s
go
EOF
""" % (topology, "*.mdcrd", trr_out)
#        cmd2="""srun --mem=10G -p rest-of-gpu1,rest-of-gpu2,rest-of-gpu3,rest-of-gpu4 cpptraj -p %s -y temp.nc << EOF > cpp2.log
#center :OL,PA origin
#image origin
#trajout %s
#go
#EOF
#""" % (topology, trr_out)
#
        subprocess.check_output(cmd1, shell=True)
#        subprocess.check_output(cmd2, shell=True)
#        os.remove("temp.nc")

    mdp_script = """integrator = sd
bd-fric = 1
dt = 0.002
nsteps = 500000
nstlog = 0
nstenergy = 0
continuation = no
constraints = h-bonds
constraint-algorithm = SHAKE
coulombtype = Cut-off
rcoulomb = 2.0
vdwtype = Cut-off
DispCorr = EnerPres
comm-mode = Linear
nstcomm = 0
comm_grps = System
tinit = 0.000
nstxout = 0
nstvout = 0
pbc = xyz
rlist = 2.0
fourierspacing = 0.1
pme_order = 4
ewald_geometry = 3d
ewald-rtol = 1e-5
rvdw-switch = 0.0
rvdw = 2.0
tcoupl = Berendsen
tc_grps = System
tau_t = 1.0
ref_t = 300
Pcoupl = Berendsen
pcoupltype = semiisotropic
tau_p = 1.0 1.0
compressibility = 4.46e-5 4.46e-5
ref_p = 1.0 1.0
gen_vel = no
gen-temp = 300"""

    mdp_file = open(mdp_input,"w")
    mdp_file.write(mdp_script)
    mdp_file.close()

    tpr_input = "md_npt.tpr"

    if not os.path.exists(tpr_input):
        print("Making TPR file")
        os.system("grompp_LS -f %s -c %s -p %s -o %s" % (mdp_input, gro, topo, tpr_input))

    cmd="""make_ndx_LS -f %s << EOF &> make_ndx.log
%s
quit
EOF
""" % (tpr_input, args.center)

    if not os.path.exists("index.ndx"):
        print("Making index file")
        subprocess.check_output(cmd, shell=True)

    trr_centered = "md_npt_centered.trr"
    if not glob.glob(trr_centered.replace(".trr","*?.trr")):
        print("Centering and splitting TRR")
        cmd="""srun --mem=5G -p rest-of-gpu1,rest-of-gpu2,rest-of-gpu3,rest-of-gpu4 trjconv_LS -split %s -f %s -o %s -center -s %s -n index.ndx << EOF &> trjconv.log
%s
0
EOF
""" % (args.chunk, trr_out, trr_centered, tpr_input, args.center.replace("r ","").replace("&","_&_"))
        subprocess.check_output(cmd, shell=True)
        print("Deleting original TRR")
        os.remove(trr_out)

    if not os.path.isdir("stress"):
        os.mkdir("stress")
    os.chdir("stress")
    stress_path = os.getcwd()

    if not os.path.exists("RUNNING"):
        print("Generating queue files")
        for trj in glob.glob("../"+trr_centered.replace(".trr","*?.trr")):
            trr_num = ''.join(c for c in os.path.basename(trj) if c.isdigit())
            if os.path.exists("ls_%s.dat0" % (trr_num)):
                continue
            queue_input = "run_gls_"+trr_num+".pbs"
            slurm_script = """#!/bin/bash -x
#SBATCH --job-name=%s
#SBATCH --output=%s.%%j
#SBATCH --error=%s.%%j
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=8:00:00
#SBATCH --partition=rest-of-gpu1,rest-of-gpu2,rest-of-gpu3,rest-of-gpu4

#
# --- Prepare calculation
#

echo "------------------------------------------------------------"
echo "Job-ID               = $SLURM_JOB_ID"
echo "Host                 = $SLURMD_NODENAME"
echo "Nodelist             = $SLURM_JOB_NODELIST"
echo "Nodes                = $SLURM_JOB_NUM_NODES"
echo "Submission Directory = $SLURM_SUBMIT_DIR"
echo "------------------------------------------------------------" 

module load GROMACS/4.5.5-ls
EXE=mdrun_LS
DO_PARALLEL="srun -n 1 "

#################################
# --- DEFINE YOUR VARIABLES --- # 
#################################
#
WORKDIR=%s
cd $WORKDIR

TPR=$(readlink -f ../%s)
TRR=$(readlink -f %s)
BASE=$(basename ${TRR//[a-z,.,_]/})


TMPFOLDER=$(mktemp -d)
cd $TMPFOLDER

$DO_PARALLEL $EXE -s $TPR -rerun $TRR -lsgridx 1 -lsgridy 1 -e ener${BASE}.edr -ols ls_${BASE}.dat -g md_${BASE}.log &> log${BASE}

cd $WORKDIR
cp ${TMPFOLDER}/* .
rm -rf $TMPFOLDER

if [ ! -f ls_${BASE}.dat ]
then
echo "Local stress file ls_${BASE}.dat not created!"
exit
fi

exit
""" % (args.name+r_num+trr_num, 
    args.name, args.name, 
    stress_path, 
    tpr_input, 
    trj)

            slurm_file = open(queue_input,"w")
            slurm_file.write(slurm_script)
            slurm_file.close()


            subprocess.check_output('/bin/bash -c "sbatch '+queue_input+'"', shell=True)
    
        os.system("touch RUNNING")
    os.chdir(cwd)


