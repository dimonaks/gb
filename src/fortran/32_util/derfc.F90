!{\src2tex{textfont=tt}}
!!****f* ABINIT/derfc
!! NAME 
!! derfc
!!
!! FUNCTION
!! Evaluates the complementary error function in real(dp).
!! Same implementation as imsl.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! yy
!!
!! OUTPUT
!! derfc_yy=complementary error function of yy
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!      ewald,ewald2,ewald3,ewald4,ewald9,make_efg_ion,psp2lo
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine derfc(derfc_yy,yy)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: yy
 real(dp),intent(out) :: derfc_yy

!Local variables-------------------------------
  !no_abirules
  integer          ::  done,ii,isw
  real(dp), parameter :: &
       ! coefficients for 0.0 <= yy < .477
       &  pp(5)=(/ 113.8641541510502e0_dp, 377.4852376853020e0_dp,  &
       &           3209.377589138469e0_dp, .1857777061846032e0_dp,  &
       &           3.161123743870566e0_dp /)
  real(dp), parameter :: &
       &  qq(4)=(/ 244.0246379344442e0_dp, 1282.616526077372e0_dp,  &
       &           2844.236833439171e0_dp, 23.60129095234412e0_dp/)
  ! coefficients for .477 <= yy <= 4.0
  real(dp), parameter :: &
       &  p1(9)=(/ 8.883149794388376e0_dp, 66.11919063714163e0_dp,  &
       &           298.6351381974001e0_dp, 881.9522212417691e0_dp,  &
       &           1712.047612634071e0_dp, 2051.078377826071e0_dp,  &
       &           1230.339354797997e0_dp, 2.153115354744038e-8_dp, &
       &           .5641884969886701e0_dp /)
  real(dp), parameter :: &
       &  q1(8)=(/ 117.6939508913125e0_dp, 537.1811018620099e0_dp,  &
       &           1621.389574566690e0_dp, 3290.799235733460e0_dp,  &
       &           4362.619090143247e0_dp, 3439.367674143722e0_dp,  &
       &           1230.339354803749e0_dp, 15.74492611070983e0_dp/)
  ! coefficients for 4.0 < y,
  real(dp), parameter :: &
       &  p2(6)=(/ -3.603448999498044e-01_dp, -1.257817261112292e-01_dp,   &
       &           -1.608378514874228e-02_dp, -6.587491615298378e-04_dp,   &
       &           -1.631538713730210e-02_dp, -3.053266349612323e-01_dp/)
  real(dp), parameter :: &
       &  q2(5)=(/ 1.872952849923460e0_dp   , 5.279051029514284e-01_dp,    &
       &           6.051834131244132e-02_dp , 2.335204976268692e-03_dp,    &
       &           2.568520192289822e0_dp /)
  real(dp), parameter :: &
       &  sqrpi=.5641895835477563e0_dp, xbig=13.3e0_dp, xlarge=6.375e0_dp, xmin=1.0e-10_dp
  real(dp) ::  res,xden,xi,xnum,xsq,xx

!******************************************************************
!BEGIN EXECUTABLE SECTION

!DEBUG
!write(6,*)' imsl derfc routine '
!stop
!ENDDEBUG

 xx = yy
 isw = 1
!Here change the sign of xx, and keep track of it thanks to isw
 if (xx<0.0e0_dp) then
   isw = -1
   xx = -xx
 end if

 done=0

!Residual value, if yy < -6.375e0_dp
 res=2.0e0_dp

!abs(yy) < .477, evaluate approximation for erfc
 if (xx<0.477e0_dp) then
!  xmin is a very small number
   if (xx<xmin) then
     res = xx*pp(3)/qq(3)
   else
     xsq = xx*xx
     xnum = pp(4)*xsq+pp(5)
     xden = xsq+qq(4)
     do ii = 1,3
       xnum = xnum*xsq+pp(ii)
       xden = xden*xsq+qq(ii)
     end do
     res = xx*xnum/xden
   end if
   if (isw==-1) res = -res
   res = 1.0e0_dp-res
   done=1
 end if

!.477 < abs(yy) < 4.0 , evaluate approximation for erfc
 if (xx<=4.0e0_dp .and. done==0 ) then
   xsq = xx*xx
   xnum = p1(8)*xx+p1(9)
   xden = xx+q1(8)
   do ii=1,7
     xnum = xnum*xx+p1(ii)
     xden = xden*xx+q1(ii)
   end do
   res = xnum/xden
   res = res* exp(-xsq)
   if (isw.eq.-1) res = 2.0e0_dp-res
   done=1
 end if

!y > 13.3e0_dp
 if (isw > 0 .and. xx > xbig .and. done==0 ) then
   res = 0.0e0_dp
   done=1
 end if

!4.0 < yy < 13.3e0_dp  .or. -6.375e0_dp < yy < -4.0
!evaluate minimax approximation for erfc
 if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
   xsq = xx*xx
   xi = 1.0e0_dp/xsq
   xnum= p2(5)*xi+p2(6)
   xden = xi+q2(5)
   do ii = 1,4
     xnum = xnum*xi+p2(ii)
     xden = xden*xi+q2(ii)
   end do
   res = (sqrpi+xi*xnum/xden)/xx
   res = res* exp(-xsq)
   if (isw.eq.-1) res = 2.0e0_dp-res
 end if

!All cases have been investigated
 derfc_yy = res

end subroutine derfc
!!***
