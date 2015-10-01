!{\src2tex{textfont=tt}}
!!****f* ABINIT/slxyzs
!! NAME
!! slxyzs
!!
!! FUNCTION
!! computes the matrix element <Sl'm'|L_idir|Slm>
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   integer :: lp,mp,idir,ll,mm
!!
!! OUTPUT
!!   complex(dpc) :: sls_val
!! 
!! NOTES
!! Slm is the real spherical harmonic used througout abinit, L_idir is a component of the angular
!! momentum operator. The subroutine computes <S_l'm'|L_idir|S_lm>
!!
!! PARENTS
!!      gipaw_aug_fields
!!
!! CHILDREN
!!      lxyz,ys
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine slxyzs(lp,mp,idir,ll,mm,sls_val)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util, except_this_one => slxyzs
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: idir,ll,lp,mm,mp
 complex(dpc),intent(out) :: sls_val

!Local variables ---------------------------------------
!scalars
 integer :: lpp,lppp,mpp,mppp
 complex(dpc) :: lidir,sy_val,ys_val

! *********************************************************************

 sls_val = czero
 
 if (lp == ll) then
   lpp  = ll
   lppp = ll
   do mpp = -lpp, lpp
     call ys(lpp,mpp,lp,mp,sy_val)
     do mppp = -lppp, lppp
       call lxyz(lpp,mpp,idir,lppp,mppp,lidir)
       call ys(lppp,mppp,ll,mm,ys_val)
       sls_val = sls_val + conjg(sy_val)*lidir*ys_val
     end do
   end do
 end if
 
end subroutine slxyzs

!!***
