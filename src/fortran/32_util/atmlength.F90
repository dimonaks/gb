!{\src2tex{textfont=tt}}
!!****f* ABINIT/atmlength
!! NAME
!! atmlength
!!
!! FUNCTION
!! Return atomic decay length for one given type of atom.
!! This length is used to generate an approximate atomic gaussian density
!! in reciprocal space:   n^AT(G)=exp[-(2pi.length.G)^2]
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (DCA,XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! densty=parameter for initialisation of the density of this atom type
!!        if densty>0, returned decay length if densty !
!! zion=charge on current type of atom (real number)
!! znucl=atomic number, for current type of atom
!!
!! OUTPUT
!! length=decay lenth
!!
!! PARENTS
!!      initro
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine atmlength(densty,length,zion,znucl)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: densty,zion,znucl
 real(dp),intent(out) :: length

!Local variables-------------------------------
!scalars
 integer :: nval
 real(dp) :: coreel
!arrays
 real(dp) :: data_length(16)

! *************************************************************************

!Either use the input value, or the default value, tabulated now.
 if(abs(densty)>tol10)then
   length=densty
 else

!  Count the number of core electrons.
   coreel=znucl-zion
!  Round the number of valence electrons
   nval=nint(zion)

!  For each set of core electron numbers, there are different decay lengths,
!  they start from nval=1, and proceed by group of 5, until a default is used

   if (nval==0) then
     length=zero

!    Bare ions : adjusted on 1h and 2he only
   else if(coreel<0.5)then
     data_length(1:4)=(/ .6_dp,.4_dp,.3_dp,.25_dp /)
     length=.2_dp
     if(nval<=4)length=data_length(nval)

!    1s2 core : adjusted on 3li, 6c, 7n, and 8o
   else if(coreel<2.5)then
     data_length(1:8)=(/ 1.8_dp,1.4_dp,1.0_dp ,.7_dp,.6_dp,&
&     .5_dp, .4_dp, .35_dp /)
     length=.3_dp
     if(nval<=8)length=data_length(nval)

!    Ne core (1s2 2s2 2p6) : adjusted on 11na, 13al, 14si and 17cl
   else if(coreel<10.5)then
     data_length(1:10)=(/ 2.0_dp,1.6_dp,1.25_dp,1.1_dp,1.0_dp,&
&     .9_dp, .8_dp, .7_dp , .7_dp, .7_dp  /)
     length=.6_dp
     if(nval<=10)length=data_length(nval)

!    Mg core (1s2 2s2 2p6 3s2) : adjusted on 19k, and on coreel==10
   else if(coreel<12.5)then
     data_length(1:10)=(/ 1.9_dp,1.5_dp,1.15_dp,1.0_dp,0.9_dp,&
&     .8_dp, .7_dp, .6_dp , .6_dp, .6_dp  /)
     length=.5_dp
     if(nval<=10)length=data_length(nval)

!    Ar core (Ne + 3s2 3p6) : adjusted on 20ca, 25mn and 30zn
   else if(coreel<18.5)then
     data_length(1:12)=(/ 2.0_dp ,1.8_dp ,1.5_dp,1.2_dp ,1.0_dp,&
&     .9_dp , .85_dp, .8_dp, .75_dp, .7_dp,&
&     .65_dp, .65_dp /)
     length=.6_dp
     if(nval<=12)length=data_length(nval)

!    Full 3rd shell core (Ar + 3d10) : adjusted on 31ga, 34se and 38sr
   else if(coreel<28.5)then
     data_length(1:14)=(/ 1.5_dp ,1.25_dp,1.15_dp,1.05_dp,1.00_dp,&
&     .95_dp, .95_dp, .9_dp , .9_dp , .85_dp,&
&     .85_dp, .80_dp, .8_dp , .75_dp         /)
     length=.7_dp
     if(nval<=14)length=data_length(nval)

!    Krypton core (Ar + 3d10 4s2 4p6) : adjusted on 39y, 42mo and 48cd
   else if(coreel<36.5)then
     data_length(1:12)=(/ 2.0_dp ,2.00_dp,1.60_dp,1.40_dp,1.25_dp,&
&     1.10_dp,1.00_dp, .95_dp, .90_dp, .85_dp,&
&     .80_dp, .75_dp /)
     length=.7_dp
     if(nval<=12)length=data_length(nval)

!    For the remaining elements, consider a function of nval only
   else
     data_length(1:12)=(/ 2.0_dp ,2.00_dp,1.55_dp,1.25_dp,1.15_dp,&
&     1.10_dp,1.05_dp,1.0_dp , .95_dp , .9_dp,&
&     .85_dp, .85_dp /)
     length=.8_dp
     if(nval<=12)length=data_length(nval)

   end if

!  End the choice between default and no-default
 end if

!DEBUG
!Here, use the previous default
!length=1.2_dp
!ENDDEBUG

end subroutine atmlength
!!***
