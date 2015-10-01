!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_kptrank
!!
!! NAME
!! m_kptrank
!!
!! FUNCTION
!! This module deals with rank objects for hashing k-point vector lists
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_kptrank

use defs_basis

implicit none

! structure to contain a rank/inverse rank pair of arrays, with dimensions
 type kptrank_type
   integer :: max_linear_density
   integer :: max_rank
   integer :: npoints
   integer, pointer :: invrank(:),rank(:)
 end type kptrank_type

contains
!!***

!!****f* ABINIT/mkkptrank
!!
!! NAME
!! mkkptrank
!!
!! FUNCTION
!! This routine sets up the kpt ranks for comparing kpts
!!
!! COPYRIGHT
!! Copyright (C) 2005-2010 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt = number of kpoints
!!  kpt = coordinates of kpoints
!!
!! OUTPUT
!!  kptrank_t = object containing ranking and inverse ranking
!!
!! NOTES
!!
!! PARENTS
!!      bfactor,get_full_kgrid,get_tetra,mkfskgrid,mkqptequiv,order_fs_kpts
!!      read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkkptrank (kpt,nkpt,kptrank_t)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
!arrays
 type(kptrank_type), intent(out) :: kptrank_t
 real(dp),intent(in) :: kpt(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: ikpt
 real(dp) :: smallestlen
!arrays

! *********************************************************************

! find smallest linear length
 smallestlen = one
 do ikpt=1, nkpt
   if (abs(kpt(1,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(1,ikpt)))
   if (abs(kpt(2,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(2,ikpt)))
   if (abs(kpt(3,ikpt)) > tol10) &
&     smallestlen = min(smallestlen, abs(kpt(3,ikpt)))
 end do

 kptrank_t%max_linear_density = int(one/smallestlen)+1
 kptrank_t%max_rank = 2*kptrank_t%max_linear_density**3
 kptrank_t%npoints = nkpt

 allocate (kptrank_t%rank(nkpt))
 allocate (kptrank_t%invrank(kptrank_t%max_rank))
 kptrank_t%invrank(:) = -1

!Ensure kpt(i)+one is positive, and the smallest
!difference between kpts should be larger than 1/100
!ie ngkpt < 100.
 do ikpt=1,nkpt
   call get_rank_1kpt (kpt(:,ikpt), kptrank_t%rank(ikpt), kptrank_t)

   if (kptrank_t%rank(ikpt) > kptrank_t%max_rank .or. kptrank_t%rank(ikpt) < 1) then
     write (*,*) 'ikpt, rank ', ikpt, kptrank_t%rank(ikpt)
     stop 'Error in mkkptrank : max rank exceeded or < 1'
   end if
   kptrank_t%invrank(kptrank_t%rank(ikpt)) = ikpt
 end do

end subroutine mkkptrank
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_rank_1kpt
!!
!! NAME
!! get_rank_1kpt
!!
!! FUNCTION
!! This routine calculates the rank for one kpt
!!
!! COPYRIGHT
!! Copyright (C) 2005-2010 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kpt = coordinates of kpoints
!!  kptrank_t = rank object for the k-grid we are using
!!
!! OUTPUT
!!  rank = rank of the kpoint
!!
!! NOTES
!!
!! PARENTS
!!      bfactor,elphon,get_full_kgrid,get_tetra,integrate_gamma
!!      integrate_gamma_alt,k_neighbors,m_kptrank,mkfskgrid,mkqptequiv
!!      read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE
subroutine get_rank_1kpt (kpt,rank,kptrank_t)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: rank
 type(kptrank_type), intent(in) :: kptrank_t
!arrays
 real(dp),intent(in) :: kpt(3)
! local
 real(dp) :: redkpt(3), res

 call wrap2_zero_one(kpt(1),redkpt(1),res)
 call wrap2_zero_one(kpt(2),redkpt(2),res)
 call wrap2_zero_one(kpt(3),redkpt(3),res)

 rank = int(real(kptrank_t%max_linear_density)*(redkpt(3)+half+tol8 +&
&           real(kptrank_t%max_linear_density)*(redkpt(2)+half+tol8 +&
&           real(kptrank_t%max_linear_density)*(redkpt(1)+half+tol8))))

 if (rank > kptrank_t%max_rank) then
   write (*,*) ' mkkptrank : error : rank should be inferior to ', kptrank_t%max_rank
   stop
 end if

end subroutine get_rank_1kpt
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/copy_kptrank
!!
!! NAME
!! copy_kptrank
!!
!! FUNCTION
!! This routine deallocates the arrays in a kptrank_type structure
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  kptrank_t = object containing ranking and inverse ranking, to be deallocated
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine copy_kptrank (kptrank_t_in, kptrank_t_out)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(kptrank_type), intent(in) :: kptrank_t_in
 type(kptrank_type), intent(out) :: kptrank_t_out

!Local variables -------------------------

! *********************************************************************
 kptrank_t_out%max_linear_density = kptrank_t_in%max_linear_density
 kptrank_t_out%max_rank = kptrank_t_in%max_rank
 kptrank_t_out%npoints = kptrank_t_in%npoints
 
 allocate (kptrank_t_out%rank(kptrank_t_out%npoints))
 kptrank_t_out%rank = kptrank_t_in%rank
 
 allocate (kptrank_t_out%invrank(kptrank_t_out%max_rank))
 kptrank_t_out%invrank = kptrank_t_in%invrank

end subroutine copy_kptrank
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/destroy_kptrank
!!
!! NAME
!! destroy_kptrank
!!
!! FUNCTION
!! This routine deallocates the arrays in a kptrank_type structure
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  kptrank_t = object containing ranking and inverse ranking, to be deallocated
!!
!! NOTES
!!
!! PARENTS
!!      bfactor,elphon,get_full_kgrid,get_tetra,mkfskgrid,mkqptequiv
!!      order_fs_kpts,read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_kptrank (kptrank_t)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(kptrank_type), intent(out) :: kptrank_t

!Local variables -------------------------

! *********************************************************************

 if (associated(kptrank_t%rank)) deallocate(kptrank_t%rank)
 if (associated(kptrank_t%invrank)) deallocate(kptrank_t%invrank)

end subroutine destroy_kptrank
!!***

end module m_kptrank
!!***
