!{\src2tex{textfont=tt}}
!!****f* ABINIT/plm_dphi
!! NAME
!! plm_dphi
!!
!! FUNCTION
!! Compute  m*P_lm(x)/sqrt((1-x^2)where P_lm is a legendre polynome
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
!!  plm_dphi(xx)
!!
!! NOTES
!!  This routine comes from Function Der_Phi_P(L,m,x) (pwpaw code from N. Holzwarth,
!!                                                          implemented by Y. Abraham))
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

function plm_dphi(ll,mm,xx)

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
 real(dp) :: plm_dphi
 real(dp),intent(in) :: xx

!Local variables ---------------------------------------
!scalars
 integer :: il,im
 real(dp) :: dosomx2,fact,pll,pmm,pmmp1,somx2
 character(len=500) :: message

! *********************************************************************

 if (mm.lt.0.or.mm.gt.ll.or.abs(xx).gt.1.d0) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' plm_dphi : ERROR -',ch10,&
&   '   mm < 0 or mm > ll or xx > 1 !'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 plm_dphi=zero
 if (mm==0) return

 pmm=one
 dosomx2=one
 if (mm > 0) then
   somx2=sqrt((1-xx)*(1+xx))
   fact=one
   do im=1,mm
     pmm=-pmm*fact
     fact=fact+2
   end do
   if (mm > 1) then
     do im=2,mm
       dosomx2=somx2*dosomx2
     end do
   end if
   pmm=pmm*dosomx2 !due to one more term (-1^M)
 end if
 if(ll==mm) then
   plm_dphi=pmm*mm
 else
   pmmp1=xx*(2*mm+1)*pmm
   if(ll==mm+1) then
     plm_dphi=pmmp1*mm
   else if(ll>=mm+2) then
     do il=mm+2,ll
       pll=(xx*(2*il-1)*pmmp1-(il+mm-1)*pmm)/(il-mm)
       pmm=pmmp1
       pmmp1=pll
     end do
     plm_dphi=pll*mm
   end if
 end if

 end function plm_dphi
!!***
