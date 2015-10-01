!{\src2tex{textfont=tt}}
!!****f* ABINIT/ratint
!! NAME
!! ratint
!!
!! FUNCTION
!! Bulirsch-Stoer rational interpolation, from Numerical Recipes in Fortran
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! integer :: npts number of points to use to compute the interpolating function
!! real(dp) :: xpt point at which to perform the interpolation or extrapolation
!! read(dp) :: xin(npts), yin(npts) the x and y(x) values for the interpolation
!!
!! OUTPUTS
!! real(dp) :: ypt, yerr: the interpolated point and error estimate
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ratint(npts,xin,xpt,yin,yerr,ypt)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
 real(dp),intent(in) :: xpt
 real(dp),intent(out) :: yerr,ypt
!arrays
 real(dp),intent(in) :: xin(npts),yin(npts)

!Local variables-------------------------------
!scalars
 integer :: ipts,mpts,ns
 real(dp) :: ddr,h1,hh,tt,ww
!arrays
 real(dp) :: cc(npts),dd(npts)

! *************************************************************************
 
 ns = 1
 hh=abs(xpt-xin(1))
 do ipts=1, npts
   h1=abs(xpt-xin(ipts))
   if (h1<tol8) then
     ypt=yin(ipts)
     yerr=tol8
     return
   else if (h1 < hh) then
     ns=ipts
     hh=h1
   end if
   cc(ipts)=yin(ipts)
   dd(ipts)=yin(ipts)+tol16
 end do
 ypt=yin(ns)
 ns=ns-1
 do mpts=1,npts-1
   do ipts=1,npts-mpts
     ww=cc(ipts+1)-dd(ipts)
     h1=xin(ipts+mpts)-xpt
     tt=(xin(ipts)-xpt)*dd(ipts)/h1
     ddr=tt-cc(ipts+1)
     ddr=ww/ddr
     dd(ipts)=cc(ipts+1)*ddr
     cc(ipts)=tt*ddr
   end do
   if (2*ns < npts-mpts) then
     yerr=cc(ns+1)
   else
     yerr=dd(ns)
     ns=ns-1
   end if
   ypt=ypt+yerr
 end do
 return

end subroutine ratint
!!***
