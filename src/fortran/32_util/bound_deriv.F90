!{\src2tex{textfont=tt}}
!!****f* ABINIT/bound_deriv
!! NAME
!! bound_deriv
!!
!! FUNCTION
!! Computes derivatives of a function a boundaries of interval (first and last derivative)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  func(n)= array containing function
!!  mesh <type(pawrad_type)>= radial mesh and related data
!!  nn= size of intervall
!!
!! OUTPUT
!!  yp1,ypn= derivatives of func at r(1) and r(n)
!!
!! PARENTS
!!      pawdij0,pawkij,psp7in
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine bound_deriv(func,mesh,nn,yp1,ypn)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments----------------------
 integer, intent(in) :: nn
 real(dp), intent(in) :: func(nn)
 real(dp), intent(out) :: yp1,ypn
 type(pawrad_type),intent(in) :: mesh

!*************************************************************************

 if (mesh%radfact(1)>zero) then
   yp1=1._dp/12._dp/mesh%stepint/mesh%radfact(1) &
&   *(-25._dp*func(1)+48._dp*func(2)-36._dp*func(3)+16._dp*func(4)-3._dp*func(5))
 else
   yp1=(func(2)-func(1))/(mesh%rad(2)-mesh%rad(1))
 end if
 ypn=1._dp/12._dp/mesh%stepint &
& *( 3._dp*func(nn-4)-16._dp*func(nn-3)+36._dp*func(nn-2)-48._dp*func(nn-1) &
& +25._dp*func(nn))/mesh%radfact(nn)

end subroutine bound_deriv

!!***
