!{\src2tex{textfont=tt}}
!!****f* ABINIT/plm_dtheta
!! NAME
!! plm_dtheta
!!
!! FUNCTION
!! Compute -(1-x^2)^1/2*d/dx{P_lm(x)} where P_lm is a legendre polynome
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (NH, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ll= l quantum number
!!  mm= m quantum number
!!  xx= input value
!!
!! OUTPUT
!!  plm_dtheta(xx)
!!
!! NOTES
!!  This routine comes from Function Der_Theta_P(L,m,x) (pwpaw code from N. Holzwarth,
!!                                                       implemented by Y. Abraham))
!!
!! PARENTS
!!      initylmg,initylmr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function plm_dtheta(ll,mm,xx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,mm
 real(dp) :: plm_dtheta
 real(dp),intent(in) :: xx

!Local variables ---------------------------------------
!scalars
 integer :: il,im
 real(dp) :: dosomx2,dpll,dpmm,dpmmp1,fact,pll,pmm,pmmp1,somx2
 character(len=500) :: message

! *********************************************************************

 if (mm.lt.0.or.mm.gt.ll.or.abs(xx).gt.1.d0) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' plm_dtheta : ERROR -',ch10,&
&   '   mm < 0 or mm > ll or xx > 1 !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 plm_dtheta=zero
 pmm=one
 dpmm=one
 dosomx2=one
 somx2=sqrt((1-xx)*(1+xx))
 if(mm==0)then
   dpmm=zero
 elseif (mm > 0) then
   fact=one
   do im=1,mm
     pmm=-pmm*fact*somx2
     dpmm=-dpmm*fact
     fact=fact+2
   end do
   if(mm>1)then
     do im=2,mm
       dosomx2=dosomx2*somx2
     end do
   end if
   dpmm= dpmm*mm*xx*dosomx2
 end if
 if(ll==mm)then
   plm_dtheta=dpmm
 else
   pmmp1=xx*(2*mm+1)*pmm
   dpmmp1=-(2*mm+1)*somx2*pmm+xx*(2*mm+1)*dpmm
   if(ll==mm+1) then
     plm_dtheta=dpmmp1
   else if(ll>=mm+2)then
     do il=mm+2,ll
       pll=(xx*(2*il-1)*pmmp1-(il+mm-1)*pmm)/(il-mm)
       dpll=(-somx2*(2*il-1)*pmmp1+(xx*(2*il-1)*dpmmp1-(il+mm-1)*dpmm))/(il-mm)
       pmm=pmmp1
       pmmp1=pll
       dpmm=dpmmp1
       dpmmp1=dpll
     end do
     plm_dtheta=dpll
   end if
 end if

 end function plm_dtheta
!!***
