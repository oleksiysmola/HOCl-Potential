#!/bin/bash -l 
#
# Generation of the input file for MOLPRO 
#

cd $2
export pwd=`pwd`
echo $pwd

export irun=$1
export ijob=$3

export grid="1d"


export outdir=qz.41.${ijob}


if [ ! -e $pwd/${outdir}-dat ]; then
   mkdir $pwd/${outdir}-dat
fi

if [ ! -e $pwd/${outdir}-out ]; then
   mkdir $pwd/${outdir}-out
fi



echo $outdir


export fname=hclo.${irun}.${outdir}
export datname=${outdir}.${irun}

echo $fname
#
cat<<endb> $fname.com

emp20 = -535.418740621004

gthresh,energy=1.d-10,zero=1.d-10,thrint=1.d-10,oneint=1.d-14,twoint=1.d-10,prefac=1.d-16
memory,512,m

PROC calc-mp2
 {hf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,17;wf,orbital}
 {mp2}
ENDPROC


PROC CCT_opt
 {hf;orbprint,-1;maxit,100;accu,18}
 {ccsd(t); thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 {STATUS,nocheck}
 {OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 {readvar,'xxx.opt'}
!{frequencies;}
ENDPROC



PROC CCT-f12
 !{rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;}
 !{ccsd(t)-f12b,thrden=1.0d-12,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.2;
 !thresh,energy=1.d-10,coeff=1.d-10,thrint=1.d-10,zero=1.d-11;maxit,100}
 {ccsd(t)-f12b,thrden=1.0d-10,ri_basis=optri,df_basis=awcv5z/mp2fit,df_basis_exch=v5z/jkfit,gem_beta=1.0;
  thresh,energy=1.d-10,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,800}
 !{STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC


PROC CCT
 !{rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,18;}
 {ccsd(t);thresh,energy=1.d-10,coeff=1.d-11,thrint=1.d-10,zero=1.d-10;maxit,100}
 !{STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC



symmetry,nosym
geometry={angstrom
H 
O,   1, rho,
Cl,  2, rclo, 1, alpha}


!basis = vqz
basis= vqz-f12
ansatz=3c,fix=1,canonical=1


i1 = 1
field=[-0.005, 0.005]

rhoe = 0.96
rcloe = 1.683
alphae = 105 


grida   = [107,105,110,100,115,95,120,90,130,80,140,70,150,60,160,50,170,40,180,40,30,20,15,10,8,6,4,2,0]

gridd1=  [0.0,0.01,-0.01,0.02,-0.02,0.05,-0.05,0.10,-0.10,-0.15,0.15,0.20,-0.20,0.30,-0.25,0.40,-0.30,0.50,0.60,0.8,1.0,1.2,1.5,1.6,1.8] Ang
gridd2=  [0.0,0.01,-0.01,0.02,-0.02,0.05,-0.05,0.10,-0.10,-0.15,0.15,0.20,-0.20,0.25,-0.30,0.40,-0.35,0.50,0.60,0.8,1.0,1.2,1.5,1.6,1.8] Ang

i1 = 1
field=[-0.005, 0.005]

N1max= #gridd1
N2max= #gridd2
N3max= #grida
!
!
imin =     ${irun}
imax =     ${irun} 
!
i  = 0
i1 = 0
n1 = 0

NOGPRINT,VARIABLE
Nmax = 60
do qN=0,Nmax
 do k3=1,N3max
   do q1=1,qN+1
       !
       q2 = qN-(q1-1)+1
       !
       show,k3,q1,q2
       !
       i=i+1
       !
       if (q2.le.N2max.and.q2.ge.0.and.q1.le.N1max.and.q1.ge.0) then
        irunjob = 6;
       else
        irunjob = 0;
       end if;
       !
       SHOW,irunjob;
       !
       if (irunjob.ne.0) then
         !
         i1 = i1+1
         !
         show,i1,i
         !
         if (imin.le.i1.and.i1.le.imax) then
            !
            n1 = n1 + 1 
            !
            SHOW,n1;
            !
            rho  = rhoe+gridd1(q1)
            rclo  = rcloe+gridd2(q2)
            alpha  = grida(k3)
            !
            calc-mp2
            !
            !MP2
            !
            emp2cm = (energy - emp20)*tocm
            !
            if (emp2cm.lt.80000.0) then
             !
             CCT-f12
             eeee(n1) = energy
             !
             SHOW,energy;
             !
             xh(n1)   =  rho
             xcl(n1)  =  rclo 
             xalpha(n1) =  alpha 
             !
             xxx(n1) = 'opqr'
             yyy(n1) = 'abcd'
             !
             text ### HClO QZ
             table,xxx,xh,xcl,xalpha,eeee
             DIGITS, 3,  8, 8,   5,  10
             !save,$fname.dat,new
             !
             goto,END:
             !
            endif
            !
         end if
       end if
  end do   
 end do    
end do 

END: TEXT ### HClO


text ### HClO
table,yyy,xh,xcl,xalpha,eeee
DIGITS, 3,  8,  8,   5,  10
save,$fname.dat,new

---

endb

export OMP_NUM_THREADS=1

#limit
#echo $TMP

export TMP=$TMPDIR
export wdir=$TMPDIR


# Run Molpro
echo "System TMPDIR = $TMPDIR"
echo "wdir = $wdir"

cd $wdir

#module unload mpi/intel/2017/update1/intel
#module unload compilers/intel/2017/update1
#module unload default-modules/2017
#module load gcc-libs/4.9.2
##module load compilers/gnu/4.9.2
#module load mpi/openmpi/1.8.4/gnu-4.9.2
#module load molpro/2012.1.25/gnu-4.9.2

module load molpro/2015.1.5/intel-2015-update2


#df

echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out

#$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $wdir/${fname}.out

#ls *.dat

/bin/cp *.dat $pwd/${outdir}-dat

#ls $datname.dat

if [ -e $pwd/$datname.dat ]; then
   echo "pwd-dat:" $datname.dat
   /bin/cp  $pwd/$datname.dat  $pwd/${outdir}-dat
fi

if [ -e $datname.dat ]; then
   echo "dat:" $datname.dat
   /bin/mv  $datname.dat  $pwd/${outdir}-dat
fi

if [ -e $datname.xyz.dat ]; then
   /bin/cp  $datname.xyz.dat  $pwd/${outdir}-dat
fi

if [ -e $pwd/$fname.out ]; then
    gzip $pwd/$fname.out
   /bin/mv  $pwd/$fname.out.gz  $pwd/${outdir}-out
fi

if [ -e $pwd/$fname.com ]; then
     /bin/rm $pwd/$fname.com
fi



