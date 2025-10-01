#!/bin/csh
#

# cd $TMPDIR

# # module load molpro/2020.1/openmp

# module load mpi/openmpi/3.1.4/intel-2018
# export OMP_NUM_THREADS=14
# export MKL_NUM_THREADS=14
module load cfour_mpi

# 1 12640
set point = $1
set endPoint = $2
echo $point 
while ($point < $endPoint)
    set line = `awk 'NR=='"$point" /scratch/vol1/asmola/HOCl/molpro/HOCL-GRID-LinearFix.txt`

    if ("$line" != "") then
        csh -f /scratch/vol1/asmola/HOCl/DBOCs/Cl35/RecomputeDBOCs/AUG-PCVTZ/runGenerateCfourScript.csh $line
    endif
@ point++
end
xclean