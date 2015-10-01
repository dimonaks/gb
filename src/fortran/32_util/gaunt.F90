!{\src2tex{textfont=tt}}
!!****f* ABINIT/gaunt
!! NAME
!! gaunt
!!
!! FUNCTION
!! Returns gaunt coefficient, i.e.
!!   the integral of Sqrt[4 \pi] Y*(l_i,m_i) Y*(ll,mm) Y(l_j,m_j)
!!   See the 3-j and 6-j symbols by Rotenberg, etc., (Technology Press, 1959), pg.5.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!   ll,mm,l1,l2,m1,m2= six quantum numbers defining the Gaunt coef.
!!
!! OUTPUT
!!   gaunt(ll,mm,l1,l2,m1,m2)=the value of the integral
!!
!! PARENTS
!!      realgaunt
!!
!! CHILDREN
!!      factorial,permutations
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function gaunt(ll,mm,l1,m1,l2,m2)

 use defs_basis

 use m_special_funcs,  only : factorial

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util, except_this_one => gaunt
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: l1,l2,ll,m1,m2,mm
 real(dp) :: gaunt

!Local variables ------------------------------
!scalars
 integer :: i1,i2,j1,j1half,j2,j2half,j3,j3half,j_half,jj,k1,k2,n1,n2
 real(dp) :: argument,sign,sum,xx,yy
 logical :: ok
!************************************************************************
 gaunt=zero;sum=zero;ok =.true.

 if((-m1-mm+m2) /= 0) ok = .false.
 if(abs(m1) > l1) ok = .false.
 if(abs(mm) > ll) ok = .false.
 if(abs(m2) > l2) ok = .false.

 jj = l1 + ll + l2

 if (mod(jj,2)/=0) ok = .false.
 j1 = jj-2*l2
 j2 = jj-2*ll
 j3 = jj-2*l1

 if (j1<0 .or. j2<0 .or. j3<0) ok = .false.

 if (ok) then

   xx = (2 * l1 + 1) * (2 * ll + 1) * (2 * l2 + 1)

   j1half = j1/2
   j2half = j2/2
   j3half = j3/2
   j_half = jj/2

   gaunt = (-1)**j1half * sqrt(xx)
   gaunt = gaunt * factorial(j2) * factorial(j3) / factorial(jj+1)
   gaunt = gaunt * factorial(j_half)/(factorial(j1half)&
&   * factorial(j2half)*factorial(j3half))

   yy = factorial(l2 + m2) * factorial(l2 - m2)

   if (mm>=0) then
     yy = yy * permutations(ll+mm,2*mm)
   else
     yy = yy / permutations(ll-mm,-2*mm)
   end if

   if (m1>=0) then
     yy = yy / permutations(l1+m1,2*m1)
   else
     yy = yy * permutations(l1-m1,-2*m1)
   end if

   gaunt = gaunt * sqrt(yy)

   i1 = l2 - ll - m1
   i2 = l2 - l1 + mm
   k1 = -min(0, i1, i2)

   n1 = l1 + m1
   n2 = ll - mm
   k2 = min(j1, n1, n2)

   sign = 1._dp
   if(k1>0) sign = (-1._dp)**k1

   argument = sign     * permutations(n1,k1)/factorial(k1)
   argument = argument * permutations(n2,k1)/factorial(i1 + k1)
   argument = argument * permutations(j1,k1)/factorial(i2 + k1)
   sum = sum + argument

   sign = -sign
   k1 = k1 + 1

   do while(k1 <= k2)
     argument = sign     * permutations(n1, k1)/factorial(k1)
     argument = argument * permutations(n2, k1)/factorial(i1 + k1)
     argument = argument * permutations(j1, k1)/factorial(i2 + k1)
     sum = sum + argument

     sign = -sign
     k1 = k1 + 1

   end do

 end if

 gaunt = gaunt * sum

 end function gaunt
!!***
