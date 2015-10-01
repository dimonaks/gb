!{\src2tex{textfont=tt}}
!!****f* ABINIT/plm_d2theta
!! NAME
!! plm_d2theta
!!
!! FUNCTION
!! Compute d2(Plm (cos(theta)))/d(theta)2  where P_lm is a legendre polynome
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
!!  plm_d2t(mpsang*mpsang)
!!
!! NOTES
!!
!!
!! PARENTS
!!      plm_coeff
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine plm_d2theta(mpsang,plm_d2t,xx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
 use interfaces_32_util, except_this_one => plm_d2theta
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpsang
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: plm_d2t(mpsang*mpsang)

!Local variables ---------------------------------------
!scalars
 integer :: il,ilm,ilmm1,ilmm2,im
 real(dp) :: sqrx
 character(len=500) :: message

!************************************************************************
 if (abs(xx).gt.1.d0) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' plm_d2theta : ERROR -',ch10,&
&   '   xx > 1 !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if


 plm_d2t=zero
 if (mpsang>1) then
   sqrx=sqrt(abs((1.d0-xx)*(1.d0+xx)))

   do il=1,mpsang-1
     ilm=il*il+2*il+1
     ilmm1=(il-1)*(il-1)+2*(il-1)+1
!    terme d2(Pll)/dtet2
     plm_d2t(ilm)=(2*il-1)*(sqrx*(plm_d2t(ilmm1)-(-1)**(il-1)*ass_leg_pol(il-1,il-1,xx))+&
&     2.d0*xx*(-1)**(il-1)*plm_dtheta(il-1,il-1,xx))
     plm_d2t(ilm-2*il)=plm_d2t(ilm)
!    terme d2(Pl(l-1))/dtet2
     plm_d2t(ilm-1)=(2*il-1)*(xx*(plm_d2t(ilmm1)-(-1)**(il-1)*ass_leg_pol(il-1,il-1,xx))-&
&     2.d0*sqrx*(-1)**(il-1)*plm_dtheta(il-1,il-1,xx))
     if(il>1) plm_d2t(il*il+2)=plm_d2t(ilm-1)
   end do
!  terme d2(Plm)/dtet2
   if(mpsang>2) then
     do il=2,mpsang-1
       do im=0,il-2
         ilm=il*il+il+1+im
         ilmm1=(il-1)*(il-1)+il+im
         ilmm2=(il-2)*(il-2)+il-1+im
         plm_d2t(ilm)=dfloat(2*il-1)/dfloat(il-im)*(xx*(plm_d2t(ilmm1)-(-1)**im*ass_leg_pol(il-1,im,xx))-&
&         2.d0*sqrx*(-1)**im*plm_dtheta(il-1,im,xx))-&
&         dfloat(il+im-1)/dfloat(il-im)*plm_d2t(ilmm2)
         plm_d2t(ilm-2*im)=plm_d2t(ilm)
       end do
     end do
   end if
 end if

 end subroutine plm_d2theta
!!***
