#PBS -S /bin/bash
#PBS -q debug
#PBS -l mppwidth=32
#PBS -l walltime=00:30:00

# This script does not work with hyper-threading.
# This script also assumes that mppwidth and OMP_NUM_THREADS are set to reasonable numbers.
# For example, mppwidth=31 will probably break the script.

MY_EXE="./main.Linux.Intel.prof.mpi.omp.exe"
MY_INPUTSFILE="inputs_SMC"

#MY_NUM_CORES=32  # Usually you don't need to set this.
                  # However, this can be set to a number smaller than mppwidth,
                  # if you wish to have more memory per core.

export OMP_NUM_THREADS=8    # set to 1 got pure MPI jobs
#export OMP_STACKSIZE=100M

if [ "$NERSC_HOST" != "edison" ]; then
    echo "This script is for NERSC Edison only"
    exit 1
fi

if [ x"$PBS_NP" == "x" ]; then
    echo "PBS_NP is not set."
    exit 1
fi 

if [ x"$MY_NUM_CORES" == "x" ]; then
    MY_NUM_CORES="$PBS_NP"
elif [ "$MY_NUM_CORES" -gt "$PBS_NP" ]; then
    echo "Please change mppwidth to at least $MY_NUM_CORES"
    exit 1
fi

let "MY_PES = MY_NUM_CORES / OMP_NUM_THREADS"
if [ "$MY_PES" -lt "1" ]; then
    echo "How did this happen?  MY_NUM_CORES=$MY_NUM_CORES  OMP_NUM_THREADS=$OMP_NUM_THREADS"
fi

let "MY_NUM_NODES = PBS_NP / 16"  # 16 cores per node
let "MY_PES_PER_NODE = MY_PES / MY_NUM_NODES"
let "MY_PES_PER_NUMA_NODE = MY_PES_PER_NODE / 2"  # 2 numa nodes per node
if [ "$MY_PES_PER_NUMA_NODE" -lt "1" ]; then
    MY_PES_PER_NUMA_NODE=1
fi

if [[ "$MY_EXE" == *.Intel.* ]] 
then
    if [ "$MY_PES" != "1" ]; then  # MPI/OMP
	export KMP_AFFINITY=compact
	MY_CC_KEYWORD=numa_node
	MY_NUMA_MEMORY=-ss
    else  # pure OMP job 
	export KMP_AFFINITY=scatter
	MY_CC_KEYWORD=none
    fi
else
    MY_CC_KEYWORD=cpu    
    if [ "$MY_PES" != "1" ]; then  # MPI/OMP
	MY_NUMA_MEMORY=-ss
    fi
fi


cd $PBS_O_WORKDIR

MY_APRUN="-n $MY_PES -N $MY_PES_PER_NODE -S $MY_PES_PER_NUMA_NODE -d $OMP_NUM_THREADS -cc $MY_CC_KEYWORD $MY_NUMA_MEMORY $MY_EXE $MY_INPUTSFILE"
echo "aprun $MY_APRUN"
echo ""

aprun $MY_APRUN
