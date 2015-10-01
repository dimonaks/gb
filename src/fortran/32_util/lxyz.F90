!{\src2tex{textfont=tt}}
!!****f* ABINIT/lxyz.F90
!! NAME
!! lxyz
!!
!! FUNCTION
!! computes the matrix element <Yl'm'|L_idir|Ylm>
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
!!   complex(dpc) :: lidir
!! 
!! NOTES
!! Ylm is the standard complex-valued spherical harmonic, idir is the direction in space of L
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

subroutine lxyz(lp,mp,idir,ll,mm,lidir)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: idir,ll,lp,mm,mp
 complex(dpc),intent(out) :: lidir

!Local variables ---------------------------------------

! *********************************************************************

 lidir = czero
 
 if (lp == ll) then 
   select case (idir)
     case (1) ! Lx
       if (mp == mm + 1) lidir = cmplx(0.5*sqrt(dble((ll-mm)*(ll+mm+1))),0.0)
       if (mp == mm - 1) lidir = cmplx(0.5*sqrt(dble((ll+mm)*(ll-mm+1))),0.0)
     case (2) ! Ly
       if (mp == mm + 1) lidir = cmplx(0.0,-0.5*sqrt(dble((ll-mm)*(ll+mm+1))))
       if (mp == mm - 1) lidir = cmplx(0.0,+0.5*sqrt(dble((ll+mm)*(ll-mm+1))))
     case (3) ! Lz
       if (mp == mm) lidir = cmplx(mm,0.0)
   end select
 end if

end subroutine lxyz

!!***
