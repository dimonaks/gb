!{\src2tex{textfont=tt}}
!!****f* ABINIT/solvbes
!! NAME
!! solvbes
!!
!! FUNCTION
!!    Find nq first roots of instrinsic equation:
!!               alpha.jl(Q) + beta.Q.djl/dr(Q) = 0
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  alpha,beta= factors in intrinsic equation
!!  ll= l quantum number
!!  nq= number of roots to find
!!
!! OUTPUT
!!  root(nq)= roots of instrinsic equation
!!
!! PARENTS
!!      shapebes
!!
!! CHILDREN
!!      jbessel
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine solvbes(root,alpha,beta,ll,nq)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util, except_this_one => solvbes
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: ll,nq
 real(dp) :: alpha,beta
!arrays
 real(dp) :: root(nq)

!Local variables-------------------------------
!scalars
 integer :: nroot
 real(dp),parameter :: dh=0.1_dp,tol=tol12
 real(dp) :: dum,hh,jbes,jbesp,qq,qx,y1,y2

! *************************************************************************

 qq=dh;nroot=0

 do while (nroot<nq)
   call jbessel(jbes,jbesp,dum,ll,1,qq)
   y1=alpha*jbes+beta*qq*jbesp
   qq=qq+dh
   call jbessel(jbes,jbesp,dum,ll,1,qq)
   y2=alpha*jbes+beta*qq*jbesp

   do while (y1*y2>=zero)
     qq=qq+dh
     call jbessel(jbes,jbesp,dum,ll,1,qq)
     y2=alpha*jbes+beta*qq*jbesp
   end do

   hh=dh;qx=qq
   do while (hh>tol)
     hh=half*hh
     if (y1*y2<zero) then
       qx=qx-hh
     else
       qx=qx+hh
     end if
     call jbessel(jbes,jbesp,dum,ll,1,qx)
     y2=alpha*jbes+beta*qx*jbesp
   end do
   nroot=nroot+1
   root(nroot)=qx

 end do

end subroutine solvbes
!!***
