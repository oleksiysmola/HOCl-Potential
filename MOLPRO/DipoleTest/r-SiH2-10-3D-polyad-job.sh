#!/bin/bash -l
#
# Generation of the input file for MOLPRO 
# run as ./scriptname igeometry  DIR OPTION
#

#Ntot = 1/2*(Pol+2)*(Pol+1)


cd $2
export pwd=`pwd`
echo $pwd

export irun=$1
export job=$3

export grid="3d"

export outdir=QZ-3d-${job}.10

echo $outdir


if [ ! -e $pwd/${outdir}-dat ]; then
   mkdir $pwd/${outdir}-dat
fi

if [ ! -e $pwd/${outdir}-out ]; then
   mkdir $pwd/${outdir}-out
fi

export fname=sih2.${irun}.${outdir}
export datname=${irun}.${outdir}

echo $fname
#




cat<<endb> $fname.com
***,PO

gthresh,energy=1.d-10,zero=1.d-14,thrint=1.d-14,oneint=1.d-14,twoint=1.d-14,prefac=1.d-20
memory,4112,m

!gexpec,dm,rel,darwin,massv



JOB=$3


PROC CCT-f12
 {rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,17;wf,orbital}
 !{ccsd(t)-f12b,thrden=1.0d-10,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.2;maxit,100}
 {ccsd(t)-f12c,thrden=1.0d-10,ri_basis=optri,df_basis=awcv5z/mp2fit,df_basis_exch=v5z/jkfit,gem_beta=1.0;maxit,100}
 !{mp2}
 !{STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}



ENDPROC
!
!ANSATZ=3C,fix=1,hybrid=1,canonical=1
!

basis=vqz-f12
ansatz=3c,fix=1,canonical=1


field=[-0.005, 0.005]

j1 = 0

rh1 = 1.52 Ang
rh1 = 1.52 Ang
ah1h2 = 94.5
!


!stretching
gR=[1.52,1.48, 1.55, 1.45, 1.57, 1.40, 1.60,1.35,1.65,1.30, 1.70, 1.25, 1.75, 1.20, 1.30, 1.10, 1.40,1.00,1.50,1.60,1.80,2.0,2.2,2.5,3.0] Ang
!bending
gAA=[94.5,96.,90.0, 100., 85., 105, 80.0, 110., 70., 120.0,60.,130.0,50.0,140.0, 40.0,150.0,30.0,160.0 ]


field=[-0.005, 0.005]

Vsmax1= #gR
Vsmax2= #gR
Vbmax=  #gAA
Pmin = 0 ! $job
Pmax = $job

imin =     ${irun}
imax =     ${irun} 

i=0
ic = 0
i1 = 0
do pol = 0,Pmax
 do j1=0,pol
    jt1 = j1+1;
  do j2=0,pol-j1
      jt2 = j2+1
      !
      k1 = pol-j1-j2+1
      !
      i=i+1
      !
      irunjob = 1
      !
      if (jt1.gt.Vsmax1.or.jt2.gt.Vsmax2.or.k1.gt.Vbmax) then
       irunjob = 0;
      else
       irunjob = 1;
      end if;
      !
      SHOW,irunjob;
      !
      if (irunjob.ne.0) then
          !
          ic = ic+1
          !
          if (imin.le.ic.and.ic.le.imax) then
            !
            i1 = i1 + 1
            !
            rh1  = gR(jt1)
            rh2   = gR(jt2)
            ah1h2 = gAA(k1)
            !
            r1(i1) =  rh1
            r2(i1) =  rh2
            r3(i1) =  ah1h2
            !
            symmetry,y
            ORIENT, NOORIENT
            geometry={angstrom
            Si,                  0.d0   ,    0.d0,         0.d0           ,
            H,     rh1*cos(0.5d0*ah1h2) ,    0.d0,   rh1*sin(0.5d0*ah1h2) ,
            H,     rh2*cos(0.5d0*ah1h2) ,    0.d0,  -rh2*sin(0.5d0*ah1h2) ,
            }
            !
            gexpec, rel,darwin,massv,lx,ly,lz,dmx,dmy,dmz;
            !
            CCT-f12
            !
            ener(i1) = energy(1)
            xxyy = 'xxyyE'
            !
            table, xxyy,r1,r2,r3,ener
            DIGITS,   0, 4, 4, 4,  12
            save,$fname.dat
            !
            !put, xyz, $datname.x.dat, old
            !
            goto,END:
            !
          endif
      endif
   enddo
 enddo
enddo

!
END: text ### MRCI
          
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

module load molpro/2015.1.5/intel-2015-update2

module list

#df

echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

time molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $pwd/${fname}.out

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



