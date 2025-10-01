#!/bin/csh
#
# Generation of the input file for MOLPRO
#


setenv PATH /home/asmola/.local/bin:$PATH

set pwd = `pwd`

set point = $1        
set directory = /scratch/vol1/asmola/HOCl/DBOCs/Cl35/RecomputeDBOCs/AUG-PCVTZ/Computed
set fname = cfour_HOCL37_DBOC_AUG-PCVTZ_${point}

cat<<endb> ZMAT
Hypochlorous acid DBOC calculation                                      
O                                                                               
H 1 ROH                                                                         
CL 1 ROCL 2 A                                                                   
                                                                                
ROH=$2                                                                        
ROCL=$3                                                                      
A=$4                                                                    
                                                                                
                                                                                
*CFOUR(CALC=CCSD,BASIS=AUG-PCVTZ,SCF_CONV=10
CC_CONV=10,DBOC=ON,MEMORY=500000000)

%isotopes
16
1
37
endb

# ls /home/zcaposm/bin/
# ls /home/zcaposm
xcfour > ${fname}.out
rm ZMAT
cp ${fname}.out ${directory}
rm ${fname}.out
xclean
