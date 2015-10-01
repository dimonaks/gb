!{\src2tex{textfont=tt}}
!!****f* ABINIT/permutations
!! NAME
!! permutations
!!
!! FUNCTION
!! Returns N!/(N-k)!  if N>=0 and N-k>0
!!                    otherwise 0 is returned
!! Output is real
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!   kk=number k to use
!!   nn=number N to use
!!
!! OUTPUT
!!   permutations= n!/(n-k)! (real)
!!
!! PARENTS
!!      gaunt
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function permutations(nn,kk)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kk,nn
 real(dp) :: permutations

!Local variables ---------------------------------------
!scalars
 integer :: ii
 real(dp) :: pp

! *********************************************************************

 if ((nn>=0).and.((nn-kk)>=0)) then
   pp=one
   do ii=nn-kk+1,nn
     pp=pp*ii
   end do
 else
   pp=zero
 end if

 permutations=pp

end function permutations
!!***
