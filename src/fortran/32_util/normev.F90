!{\src2tex{textfont=tt}}
!!****f* ABINIT/normev
!! NAME
!! normev
!!
!! FUNCTION
!! Normalize a set of num eigenvectors of complex length ndim
!! (real length 2*ndim) and set phases to make evec(i,i) real and positive.
!! Near convergence, evec(i,j) approaches delta(i,j).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  evec(2*ndim,num)=num unnormalized eigenvectors
!!  ndim=dimension of evec as shown
!!  num=number of eigenvectors and complex length thereof.
!!
!! OUTPUT
!!  evec(2*ndim,num)=num normalized eigenvectors
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      subdiago
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine normev(evec,ndim,num)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim,num
!arrays
 real(dp),intent(inout) :: evec(2*ndim,num)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: den,evim,evre,phim,phre,xnorm

! *************************************************************************
!
!Loop over vectors
 do ii=1,num
!  find norm
   xnorm=0.0d0
   do jj=1,2*ndim
     xnorm=xnorm+evec(jj,ii)**2
   end do
   xnorm=1.0d0/sqrt(xnorm)
!  Set up phase
   phre=evec(2*ii-1,ii)
   phim=evec(2*ii,ii)
   if (phim/=0.0d0) then
     den=1.0d0/sqrt(phre**2+phim**2)
     phre=phre*xnorm*den
     phim=phim*xnorm*den
   else
!    give xnorm the same sign as phre (negate if negative)
     phre=sign(xnorm,phre)
     phim=0.0d0
   end if
!  normalize with phase change
   do jj=1,2*ndim,2
     evre=evec(jj,ii)
     evim=evec(jj+1,ii)
     evec(jj,ii)=phre*evre+phim*evim
     evec(jj+1,ii)=phre*evim-phim*evre
   end do
 end do

end subroutine normev
!!***
