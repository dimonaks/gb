!{\src2tex{textfont=tt}}
!!****f* ABINIT/plm_coeff
!! NAME
!! plm_coeff
!!
!! FUNCTION
!! Compute coefficients depending on Plm and its derivatives where P_lm is a legendre polynome.
!! They are used to compute the second derivatives of spherical harmonics
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mpsang=1+ maximum l quantum number
!!  xx= input value
!!
!! OUTPUT
!!  blm(5,mpsang*mpsang)=coefficients depending on Plm and its derivatives where P_lm is a legendre polynome
!!
!! NOTES
!!
!!
!! PARENTS
!!      initylmg
!!
!! CHILDREN
!!      leave_new,pl_deriv,plm_d2theta,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine plm_coeff(blm,mpsang,xx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
 use interfaces_32_util, except_this_one => plm_coeff
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpsang
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: blm(5,mpsang*mpsang)

!Local variables ---------------------------------------
!scalars
 integer :: il,ilm,ilm0,ilm1,im
 real(dp) :: dplm_dt,d2plm_dt2,llp1,onemx2,plm,sqrx,xsqrx,xx2,yy
 logical :: is_one
 character(len=500) :: message
!arrays
 real(dp) :: pl_d2(mpsang),plm_d2t(mpsang*mpsang)

!************************************************************************

 if (abs(xx).gt.1.d0) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' plm_d2theta : ERROR -',ch10,&
&   '   xx > 1 !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 blm=zero
 is_one=(abs(abs(xx)-one)<=tol12)
 xx2=xx**2
 onemx2=abs(one-xx2)
 sqrx=sqrt(onemx2)
 xsqrx=xx*sqrt(onemx2)

 call plm_d2theta(mpsang,plm_d2t,xx)
 if (is_one) then
   yy=sign(one,xx)
   call pl_deriv(mpsang,pl_d2,yy)
 end if


 do il=0,mpsang-1
   llp1=dble(il*(il+1))
   ilm0=il*il+il+1
   do im=0,il
     ilm=ilm0+im;ilm1=ilm0-im

     plm      =(-1)**im*ass_leg_pol(il,im,xx)
     dplm_dt  =(-1)**im*plm_dtheta(il,im,xx)
     d2plm_dt2=         plm_d2t(ilm)

     blm(1,ilm)=         two*xsqrx    *dplm_dt+onemx2*d2plm_dt2
     blm(2,ilm)=         (one-two*xx2)*dplm_dt-xsqrx *d2plm_dt2
     blm(3,ilm)=llp1*plm+                             d2plm_dt2
     blm(4,ilm)=        -two*xsqrx    *dplm_dt+xx2   *d2plm_dt2


     if (is_one) then
       if (im==1) then
         blm(5,ilm)=llp1*plm+d2plm_dt2
       end if
       if (im==2) then
         blm(5,ilm)=d2plm_dt2-three*pl_d2(il+1)
       end if
     else
       if(im>0) then
         blm(5,ilm)=plm/onemx2-dplm_dt*xx/sqrx
       end if
     end if

     if (im>0) then
       blm(1,ilm1)=blm(1,ilm)
       blm(2,ilm1)=blm(2,ilm)
       blm(3,ilm1)=blm(3,ilm)
       blm(4,ilm1)=blm(4,ilm)
       blm(5,ilm1)=blm(5,ilm)
     end if

   end do
 end do

 end subroutine plm_coeff
!!***
