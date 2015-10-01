!{\src2tex{textfont=tt}}
!!****f* ABINIT/ys
!! NAME
!! ys
!!
!! FUNCTION
!! computes the matrix element <Yl'm'|Slm>
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   integer :: l',m',l,m
!!
!! OUTPUT
!!   complex(dpc) :: ys_val 
!! 
!! NOTES
!! Ylm is the standard complex-valued spherical harmonic, Slm is the real spherical harmonic
!! used througout abinit. <Yl'm'|Slm> is their overlap.
!!
!! PARENTS
!!      slxyzs
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ys(lp,mp,ll,mm,ys_val)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,lp,mm,mp
 complex(dpc),intent(out) :: ys_val

!Local variables ---------------------------------------
!scalars
 real(dp) :: mone

! *********************************************************************

 mone = -one
 ys_val = czero
 
 if (lp == ll) then 
   select case (mm)
     case (0) ! case for S_l0
       if (mp == mm) ys_val = cmplx(1.0,0.0)
     case (:-1) ! case for S_lm with m < 0
       if (mp == -mm) ys_val = cmplx(0.0,-mone**mm*sqrthalf)
       if (mp == mm) ys_val = cmplx(0.0,sqrthalf)
     case (1:) ! case for S_lm with m > 0
       if (mp == mm) ys_val = cmplx(mone**mm*sqrthalf,0.0)
       if (mp == -mm) ys_val = cmplx(sqrthalf,0.0)
   end select
 end if
 
end subroutine ys

!!***
