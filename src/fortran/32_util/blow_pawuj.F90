!{\src2tex{textfont=tt}}
!!****f* ABINIT/blow_pawuj
!!
!! NAME
!! blow_pawuj
!!
!! FUNCTION
!! This subroutine reads a real nxn matrice and appends lines n+1 and clumn n+1 containing 
!! the sum of the lines 
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DJA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt . 
!!
!! INPUTS
!!  mat(nj,nj) matrix to be completed 
!!
!! OUTPUT
!!  matt(nj+1,nj+1) completed matrix 
!!
!! PARENTS
!!      pawuj_utils
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine blow_pawuj(mat,nj,matt)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)        :: nj
!arrays
 real(dp),intent(in)       :: mat(nj,nj)
 real(dp),intent(out)      :: matt(nj+1,nj+1)

!Local variables-------------------------------
!scalars
 integer                   :: ii
!arrays

! *************************************************************************

 matt(1:nj,1:nj)=mat
 do  ii = 1,nj
   matt(ii,nj+1)=-sum(matt(ii,1:nj))
 end do

 do  ii = 1,nj+1
   matt(nj+1,ii)=-sum(matt(1:nj,ii))
 end do
 
end subroutine blow_pawuj 
!!***
