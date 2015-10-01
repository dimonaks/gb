!{\src2tex{textfont=tt}}
!!****f* ABINIT/ylm_cmplx
!! NAME
!! ylm_cmplx
!!
!! FUNCTION
!! Calculate all (complex) spherical harmonics for lx<=4
!! 
!! COPYRIGHT
!!  Copyright (C) 2008-2010 ABINIT group (drh, TRangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! lx= quantum numbers.
!! xx= cartesian coordinate in the x direction
!! yy= cartesian coordinate in the y direction
!! zz= cartesian coordinate in the z direction
!!
!! cartesian coordinates
!! OUTPUT
!!  ylm((lx+1)*(lx+1)) complex spherical harmonics for all l<=lx and all
!!                     possible values of m.
!!
!! SIDE EFFECTS
!!
!! NOTES
!! we are supressing the so-called Condon-Shortley phase
!!
!! PARENTS
!!      mlwfovlp_proj
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ylm_cmplx(lx,ylm,xx,yy,zz)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lx
 real(dp),intent(in) :: xx,yy,zz
!arrays
 complex(dpc),intent(out) :: ylm((lx+1)*(lx+1))

!Local variables-------------------------------
!scalars
 integer :: ii,l1,m1,nc,nn
 real(dp) :: dc,dl,dm,ds,rr,rrs,rs,sq2,w,x,xs,ya,yi,yr
!arrays
 real(dp) :: cosa(lx+1),fact(2*(lx+1)),plm(lx+2,lx+2),qlm(lx+2,lx+2),sgn(lx+1)
 real(dp) :: sina(lx+1)

! *************************************************************************
 
!DEBUG
!write (std_out,*) ' ylm_cmplx : enter'
!ENDDEBUG

!normalization coefficients
 sq2=sqrt(2.0d0)
 fact(1)=1.0d0
 do ii=2,2*lx+1
   fact(ii)=(ii-1)*fact(ii-1)
 end do
 do l1=1,lx+1
   sgn(l1)=(-1.d0)**(l1-1)
   do m1=1,l1
     qlm(l1,m1)=sqrt((2*l1-1)*fact(l1-m1+1)/&
&     (four_pi*fact(l1+m1-1)))
   end do
 end do

!legendre polynomials
 rs=xx**2 + yy**2 + zz**2
 if(rs > tol8) then
   xs=zz**2/rs
   x=zz/sqrt(rs)
   w=sqrt(abs(1.0d0 - xs))
 else
   x=0.0d0

   w=1.0d0
 end if
 plm(1,1)=1.0d0
 plm(2,1)=x
 plm(2,2)=w
 plm(3,2)=3.0d0*x*w
 do m1=1,lx
   dm=m1-1
   if(m1 > 1) then
     plm(m1+1,m1)=x*plm(m1,m1) + 2*dm*w*plm(m1,m1-1)
   end if
   if(m1 < lx) then
     do l1=m1+2,lx+1
       dl=l1-1
       plm(l1,m1)=((2*dl-1)*x*plm(l1-1,m1)&
&       - (dl+dm-1)*plm(l1-2,m1))/(dl-dm)
     end do
   end if
   plm(m1+1,m1+1)=(2*dm+1)*w*plm(m1,m1)
 end do

!azimuthal angle phase factors
 rrs=xx**2 + yy**2
 if(rrs > tol8) then
   rr=sqrt(rrs)
   dc=xx/rr
   ds=yy/rr
 else
   dc=1.0d0
   ds=0.0d0
 end if
 cosa(1)=1.0d0
 sina(1)=0.0d0
 do m1=2,lx+1
   cosa(m1)=dc*cosa(m1-1) - ds*sina(m1-1)
   sina(m1)=ds*cosa(m1-1) + dc*sina(m1-1)
 end do

!combine factors
 do l1=1,lx+1
   do m1=2,l1
     nn=(l1-1)**2 + (l1-1) + (m1-1) + 1
     nc=(l1-1)**2 + (l1-1) - (m1-1) + 1
!    note that we are supressing the so-called Condon-Shortley phase
!    ya=sgn(m1)*qlm(l1,m1)*plm(l1,m1)
     ya=qlm(l1,m1)*plm(l1,m1)
     yr=ya*cosa(m1)
     yi=ya*sina(m1)
     ylm(nc)=sgn(m1)*cmplx(yr,-yi)
     ylm(nn)=cmplx(yr,yi)
   end do
 end do
 do l1=1,lx+1
   nn=(l1-1)**2 + (l1-1) + 1
   ya=qlm(l1,1)*plm(l1,1)
   ylm(nn)=cmplx(ya,0.d0)
 end do

end subroutine ylm_cmplx
!!***
