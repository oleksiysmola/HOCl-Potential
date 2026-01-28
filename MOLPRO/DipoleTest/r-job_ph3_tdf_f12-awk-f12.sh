#!/bin/csh 
#
# Generation of the input file for MOLPRO 
#

cd $2
set pwd = `pwd`
echo $pwd

set irun = $1


set fname = so3.f12.dip.${irun}.tdf.sg.th
#
cat<<endb> $fname.com

gthresh,energy=1.d-10,zero=1.d-14,thrint=1.d-14,oneint=1.d-14,twoint=1.d-14,prefac=1.d-20



memory,1600,m
PROC CCT_opt
 {hf;orbprint,-1;maxit,100;accu,18;wf,charge=1}
 {ccsd(t); thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 {STATUS,nocheck}
 {OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 {readvar,'xxx.opt'}
!{frequencies;}
ENDPROC

PROC CCT
 {hf;orbprint,-1;maxit,100;accu,18;wf,charge=1}
 {ccsd(t); thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 {STATUS,nocheck}
ENDPROC


PROC HF-proc
 {rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,17;wf,40,1,orbital,IGNORE_ERROR}
ENDPROC


PROC CCT-f12
 {rhf,ipnit=1,maxdis=30;orbprint,-1;maxit,150;accu,17;wf,40,1,orbital,IGNORE_ERROR}
 {ccsd(t)-f12b,thrden=1.0d-12,ri_basis=ri_basis_set,df_basis=AWCV5Z/MP2FIT,df_basis_exch=V5Z/JKFIT,gem_beta=1.2;
 thresh,energy=1.d-9,coeff=1.d-9,thrint=1.d-10,zero=1.d-10;maxit,100}
 !{STATUS,nocheck}
 !{OPTG,grad=5,energy=8,saveact='xxx.opt';print,history;inactive,alpha}
 !{readvar,'xxx.opt'}
ENDPROC


ANSATZ=3C,fix=1,hybrid=1,canonical=1



BASIS
!
H=aug-cc-pVQZ

P=aug-cc-pV(Q+d)Z
!
set,ri_basis_set
!
! K.E. Yousaf and K.A. Peterson, Chem. Phys. Lett. 476, 303 (2009).
! Auxiliary RI (OptRI) matched to the aug-cc-pV(T+d)Z Basis Set for P
!
!This should be used within the CABS framework with the above orbital set
!
! aug-cc-pVQZ
!
!
s,P,1.063375E+01,1.709283E+00,7.460733E-01,3.328648E-01,1.398463E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
!
p,P,1.557718E+01,6.900878E+00,1.503643E+00,8.394805E-01,3.147842E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
c,5.5,1.000000E+00;
!
d,P,2.105815E+01,5.928729E+00,1.525675E+00,5.426733E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
!
f,P,4.171852E+00,4.698425E-01,1.849512E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
!
g,P,8.909170E-01,3.980391E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
!
h,P,7.072552E-01;
c,1.1,1.000000E+00;
!
s,H,7.744478E+00,2.171468E+00,4.626660E-01,1.722015E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
!
p,H,1.109436E+01,7.317025E+00,1.259771E+00,5.442124E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
c,4.4,1.000000E+00;
!
d,H,3.086886E+00,1.146690E+00,4.119373E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
c,3.3,1.000000E+00;
!
f,H,2.250697E+00,7.980389E-01;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
!
g,H,2.069345E+00,1.030651E+00;
c,1.1,1.000000E+00;
c,2.2,1.000000E+00;
!
END



symmetry,nosym;
orient,noorient;

geometry={angstrom;
 S;
 O, 1, R1SO;
 O, 1, R2SO, 2, AA3;
 O, 1, R3SO, 2, AA2, 3, AA1, 1;
}

!BASIS=VDZ-F12
!ANSATZ=3c,fix=1,canonical=1
!orbital,2102.2,IGNORE_ERROR


hartree = 219474.63
emp20 = -623.083175267015

!stretching
jt0 = 1
gRR=[1.417, 1.43, 1.41, 1.44, 1.40, 1.46, 1.39, 1.50, 1.37, 1.55, 1.60]
!bending
gAA=[120., 119., 121.0, 118.0, 116., 110.0, 105.0, 100.0, 95.0, 90.0 ]

 geometry={
 nosym;
 noorient;
 P;
 H, 1, R1PH(i1);
 H, 1, R2PH(i1), 2, AA3(i1);
 H, 1, R3PH(i1), 2, AA2(i1), 3, AA1(i1), 1
 }

!stretching
RR=[1.20, 1.30, 1.40, 1.50, 1.60, 1.70] ANG
!bending
AA=[119., 110., 100., 90., 80., 70., 60.]




field=[-0.001d0, 0.001d0]

Vsmax= #gRR
Vbmax= #gAA


imin =     ${irun}
imax =     ${irun} 

i=0
i1 = 0
n1 = 0
do jt1=1,Vsmax
 do jt2=1,Vsmax
 if (gRR(jt1).le.gRR(jt2)) then 
 do jt3=1,Vsmax
   if (gRR(jt2).le.gRR(jt3)) then 
   do k1=1,Vbmax 
   do k2=1,Vbmax
   if (gAA(k1).le.gAA(k2)) then 
   do k3=1,Vbmax
   if (gAA(k2).le.gAA(k3)) then 
      !
      i=i+1
      !
      if (jt1.eq.jt2.and.jt1.eq.jt3.and.k1+k2+k3.le.11) then
       irunjob = 1;
      else if (k1.eq.k2.and.k1.eq.k3.and.jt1+jt2+jt3.le.12) then
       irunjob = 2;
      else if ((jt1+jt2+jt3).le.12) then
       irunjob = 3;
      else if ((k1+k2+k3).le.11) then
       irunjob = 4;
      else if (abs(jt1+jt2+jt3+k1+k2+k3).le.12) then
       irunjob = 5;
      else if (mod(i,8).eq.0) then
       irunjob = 10;
      else
       irunjob = 0;
      end if;
      !
      SHOW,irunjob;
      !
      if (irunjob.ne.0) then
         !
         i1 = i1+1
         if (imin.le.i1.and.i1.le.imax) then
            !
            !n1 = n1 + 1
            !
            n1 = 1 
            !
            SHOW,n1;
            !
            R1SO = gRR(jt1)
            R2SO = gRR(jt2)
            R3SO = gRR(jt3)
            AA1  = gAA(k1)
            AA2  = gAA(k2)
            AA3  = gAA(k3)
            !
            R1SOt(n1) =gRR(jt1)
            R2SOt(n1) =gRR(jt2)
            R3SOt(n1) =gRR(jt3)
            AAt1(n1)  = gAA(k1)
            AAt2(n1)  = gAA(k2)
            AAt3(n1)  = gAA(k3)
            !
            !gexpec,rel,darwin,massv;
            !
            HF-proc
            !
            MP2
            !
            emp2cm = (energy - emp20)*tocm
            !
            if (emp2cm.lt.50000.0) then
                !
                dip,0.0d0
                CCT-f12
                eeee(n1) = energy
                !
                dip,field(1)
                CCT-f12
                ex6d1=energy
                dip,field(2)
                CCT-f12
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dx(n1)=dipau*TODEBYE
                !             
                dip,,field(1)
                CCT-f12
                ex6d1=energy
                dip,,field(2)
                CCT-f12
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dy(n1)=dipau*TODEBYE
                !
                dip,,,field(1)
                CCT-f12
                ex6d1=energy
                dip,,,field(2)
                CCT-f12
                ex6d2=energy
                dipau=(ex6d2-ex6d1)/(field(2)-field(1))
                dip6dz(n1)=dipau*TODEBYE
                !
                put, xyz, $fname.xyz.dat, old
                !
                SHOW,energy;
                !
                xxx(n1) = 'ttt'
                yyy(n1) = 'ppp'
                !
                text ### SO3 energies 
                table,xxx,R1SOt,R2SOt,R3SOt,AAt1,AAt2,AAt3,eeee,dip6dx,dip6dy,dip6dz
                DIGITS, 3,    5,    5,    5,   5,   5,   5,  12,    12,    12,    12 
                !
                goto,END:
                !
            endif
            !
         end if
      end if
   end if
   end do
   end if
   end do
   end do
 end if
 end do
 end if
 end do
end do
          
END: text ### SO3 CCSD(T)-f12
          
text ###  SO3 nergies 
table,yyy,R1SOt,R2SOt,R3SOt,AAt1,AAt2,AAt3,eeee,dip6dx,dip6dy,dip6dz
DIGITS, 3,    5,    5,    5,   5,   5,   5,  12,    12,    12,    12 
          
save,$fname.dat,new
          
---       
endb


setenv OMP_NUM_THREADS 1

limit
#echo $TMP

setenv wdir $TMPDIR


# Run Molpro
echo "System TMPDIR = $TMPDIR"
echo "wdir = $wdir"

cd $wdir

module load mkl/10.2.5/035
module load compilers/intel/11.1/072
module load mpi/openmpi/1.4.1/intel

setenv molproexe /shared/ucl/apps/Molpro/2012.1-bindist/bin/

#setenv molproexe /shared/ucl/apps/Molpro/2009.1/bin/
#setenv PATH /shared/ucl/apps/Molpro/2009.1/bin/:"${PATH}"


echo "Running molpro -n 1 -d $TMPDIR -W $wdir < $pwd/$fname.com > $pwd/$fname.out"

#$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $pwd/${fname}.out

$molproexe/molpro -n 1 -I $TMPDIR -d $TMPDIR -W $wdir < $pwd/${fname}.com > $wdir/${fname}.out

if (-e $fname.xyz.dat) then
   /bin/cp  $fname.xyz.dat  $pwd
endif

if (-e $fname.dat) then
   /bin/cp  $fname.dat  $pwd
endif

if (-e $fname.out) then
   /bin/cp  $fname.out  $pwd
endif


