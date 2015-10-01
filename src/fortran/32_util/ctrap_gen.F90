!{\src2tex{textfont=tt}}
!!****f* ABINIT/ctrap_gen
!! NAME
!! ctrap_gen
!!
!! FUNCTION
!! Do corrected trapezoidal integral on a given (generalized) grid
!! This routine interfaces ctrap.f (integral on a regular grid)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  func(:)=integrand values
!!
!! OUTPUT
!!  intg=resulting integral by corrected trapezoid
!!
!! PARENTS
!!
!! CHILDREN
!!      ctrap
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ctrap_gen(intg,func,radmesh)

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util, except_this_one => ctrap_gen
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: intg
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: func(radmesh%int_meshsz)

!Local variables-------------------------------
!scalars
 integer :: msz
!arrays
 real(dp),allocatable :: func2(:)

! *************************************************************************

 msz=radmesh%int_meshsz

 if (radmesh%mesh_type==1) then

   call ctrap(msz,func,radmesh%rstep,intg)

 else if (radmesh%mesh_type==2) then

   allocate(func2(msz))
   func2(1:msz)=func(1:msz)*radmesh%radfact(1:msz)
   call ctrap(msz,func2,radmesh%lstep,intg)
   deallocate(func2)

 else if (radmesh%mesh_type==3) then

!  From r1 to end
   allocate(func2(msz))
   func2(2:msz)=func(2:msz)*radmesh%rad(2:msz)
   call ctrap(msz-1,func2(2:msz),radmesh%lstep,intg)
   deallocate(func2)
!  From 0 to r1
   intg = intg &
&   + (func(1)+func(2))*radmesh%rstep*0.5d0

 else if (radmesh%mesh_type==4) then

   allocate(func2(msz))
   func2(1:msz)=func(1:msz)*radmesh%radfact(1:msz)
   call ctrap(msz,func2,radmesh%lstep,intg)
   deallocate(func2)

 end if

end subroutine ctrap_gen
!!***
