!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_special_funcs
!! NAME
!! m_special_funcs
!!
!! FUNCTION
!! This module contains routines and functions used to 
!! evaluate special functions frequently needed in Abinit.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2010 ABINIT group (MG,MT,FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_special_funcs

 use defs_basis
 use m_errors

 implicit none

 private

 public :: jbessel_4spline   ! Spherical Bessel functions and derivatives employing a polynomial approximation for q->0
 public :: ylmc              ! Complex Spherical harmonics for l<=3.
 public :: ylmcd             ! First derivative of complex Ylm wrt theta and phi up to l<=3
 public :: factorial         ! Calculates N! returning a real.
 public :: binomcoeff        ! Binominal coefficient n!/(n-k)!
 public :: laguerre          ! Laguerre Polynomial(x,n,a). 
 public :: RadFnH            ! Atomic radial function(r,n,l,Z).
 public :: iradfnh           ! Norm of atomic radial function(a,b,n,l,Z).  

CONTAINS  !===========================================================
!!***

!!****f* m_special_funcs/jbessel_4spline
!! NAME
!!  jbessel_4spline
!!
!! FUNCTION
!!  Compute spherical Bessel functions and derivatives. 
!!  A polynomial approximation is employed for q-->0.
!!  
!! INPUT
!!  ll=l-order of the Bessel function
!!  tol=tolerance below which a Polynomial approximation is employed
!!   both for jl and its derivative (if required)
!!  order=1 if only first derivative is requested
!!        2 if first and second derivatives are requested
!!  xx=where to compute j_l
!!
!! OUTPUT
!!  bes=Spherical Bessel function j_l at xx
!!  besp= first derivative of j_l at xx (only if order>=1)
!!
!! TODO 
!! Remove inline definitions, they are obsolete in F2003
!!
!! PARENTS
!!      m_paw_pwij,psp7nl
!!
!! CHILDREN
!!
!! SOURCE

subroutine jbessel_4spline(bes,besp,ll,order,xx,tol)
!Arguments ---------------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util
!End of the abilint section

 integer,intent(in) :: ll,order
 real(dp),intent(in) :: xx,tol
 real(dp),intent(out) :: bes,besp

!Local variables ---------------------------------------
!scalars
 real(dp) :: bespp
 real(dp) :: arg,bes0a,bes0ap,bes0b,bes0bp,bes1a,bes1ap,bes1b,bes1bp
 real(dp) :: bes2a,bes2ap,bes2b,bes2bp,bes3a,bes3ap,bes3b,bes3bp
 character(len=500) :: msg
! *********************************************************************

! === l=0,1,2 and 3 spherical Bessel functions (and derivatives) ===
 bes0a(arg)=1.0_dp-arg**2/6.0_dp*(1.0_dp-arg**2/20.0_dp)
 bes0b(arg)=sin(arg)/arg
 bes1a(arg)=(10.0_dp-arg*arg)*arg/30.0_dp
 bes1b(arg)=(sin(arg)-arg*cos(arg))/arg**2
 bes2a(arg)=arg*arg/15.0_dp-arg**4/210.0_dp
 bes2b(arg)=((3.0_dp-arg**2)*sin(arg)-3.0_dp*arg*cos(arg))/arg**3
 bes3a(arg)=arg*arg*arg/105.0_dp-arg**5/1890.0_dp+arg**7/83160.0_dp
 bes3b(arg)=(15.0_dp*sin(arg)-15.0_dp*arg*cos(arg)-6.0_dp*arg**2*sin(arg)+arg**3*cos(arg))/arg**4
 bes0ap(arg)=(-10.0_dp+arg*arg)*arg/30.0_dp
 bes0bp(arg)=-(sin(arg)-arg*cos(arg))/arg**2
 bes1ap(arg)=(10.0_dp-3.0_dp*arg*arg)/30.0_dp
 bes1bp(arg)=((arg*arg-2.0_dp)*sin(arg)+2.0_dp*arg*cos(arg))/arg**3
 bes2ap(arg)=(1.0_dp-arg*arg/7.0_dp)*2.0_dp*arg/15.0_dp
 bes2bp(arg)=((4.0_dp*arg*arg-9.0_dp)*sin(arg)+(9.0_dp-arg*arg)*arg*cos(arg))/arg**4
 bes3ap(arg)=(1.0_dp/35-arg*arg/378.0_dp+arg**4/11880.0_dp)*arg*arg
 bes3bp(arg)=((-60.0_dp+27.0_dp*arg*arg-arg**4)*sin(arg)+(60.0_dp*arg-7.0_dp*arg**3)*cos(arg))/arg**5

 ! This is to test jbessel calculation without polynomial approximation for q-->0.
 ! call jbessel(bes,besp,bespp,ll,order,xx)
 ! RETURN

 if (order>2) then 
   MSG_ERROR("Wrong order in jbessel")
 end if

 select case (ll)
 case (0)
   if (xx<TOL) then
     bes=bes0a(xx)
     if (order>=1) besp=bes0ap(xx)
   else
     bes=bes0b(xx)
     if (order>=1) besp=bes0bp(xx)
   end if

 case (1)
  if (xx<TOL) then
    bes=bes1a(xx)
    if (order>=1) besp=bes1ap(xx)
  else
    bes=bes1b(xx)
    if (order>=1) besp=bes1bp(xx)
  end if

 case (2)
   if (xx<TOL) then
     bes=bes2a(xx)
     if (order>=1) besp=bes2ap(xx)
   else
     bes=bes2b(xx)
     if (order>=1) besp=bes2bp(xx)
   end if

 case (3)
   if (xx<TOL) then
     bes=bes3a(xx)
     if (order>=1) besp=bes3ap(xx)
   else
     bes=bes3b(xx)
     if (order>=1) besp=bes3bp(xx)
   end if

 case (4:)
   call jbessel(bes,besp,bespp,ll,order,xx)

 case default
   write(msg,'(a,i4)')' wrong value for ll = ',ll
   MSG_BUG(msg)
 end select

end subroutine jbessel_4spline
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/ylmc
!! NAME
!! ylmc
!!
!! FUNCTION
!!  Return a complex spherical harmonic with l <= 3
!!
!! INPUTS
!!  il=angular quantum number
!!  im=magnetic quantum number
!!  kcart=vector in cartesian coordinates defining the value of \theta and \psi
!!   where calculate the spherical harmonic
!!
!! OUTPUT
!!  ylm= spherical harmonic
!!
!! NOTES
!!  Note the use of double precision complex.
!!  Case l>3 not implemented.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ylmc(il,im,kcart)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,im
 complex(dpc) :: ylmc
!arrays
 real(dp),intent(in) :: kcart(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: LMAX=3
 real(dp),parameter :: PPAD=tol8
 real(dp) :: cosphi,costh,costhreephi,costwophi,r,rxy,sinphi,sinth,sinthreephi,sintwophi
 !complex(dpc) :: new_ylmc
 character(len=500) :: msg

! *************************************************************************

 if (ABS(im)>ABS(il)) then
   write(msg,'(3(a,i0))')' m is,',im,' however it should be between ',-il,' and ',il
   MSG_ERROR(msg)
 end if

 ylmc = czero

 r=SQRT(kcart(1)**2+kcart(2)**2+kcart(3)**2)
 if (r<PPAD) r=r+PPAD
 !$if (r<tol10) RETURN

 rxy=SQRT(kcart(1)**2+kcart(2)**2)
 if (rxy<PPAD)rxy=r+PPAD
!
! Determine theta and phi 
 costh= kcart(3)/r

#if 1
 ! old buggy coding
 sinth= rxy/r
 cosphi= kcart(1)/rxy
 sinphi= kcart(2)/rxy
#else
 sinth=sqrt(abs((one-costh)*(one+costh))) ! abs is needed to prevent very small negative arg
 cosphi=one
 sinphi=zero
 if (sinth>tol10) then
   cosphi=kcart(1)/(r*sinth)
   sinphi=kcart(2)/(r*sinth)
 end if
#endif

 costwophi= two*cosphi**2 - one
 sintwophi= two*sinphi*cosphi
 costhreephi=cosphi*costwophi-sinphi*sintwophi
 sinthreephi=cosphi*sintwophi+sinphi*costwophi

 select case (il)

 case (0)
  ylmc= one/SQRT(four_pi)

 case (1)
  if (ABS(im)==0) then
   ylmc = SQRT(three/(four_pi))*costh
  else if (ABS(im)==1) then
   ylmc = -SQRT(three/(eight*pi))*sinth*CMPLX(cosphi,sinphi)
  else 
   MSG_ERROR("wrong im")
  end if

 case (2)
  if (ABS(im)==0) then
   ylmc = SQRT(5.d0/(16.d0*pi))*(three*costh**2-one)
  else if (ABS(im)==1) then
   ylmc = -SQRT(15.d0/(8.d0*pi))*sinth*costh*cmplx(cosphi,sinphi)
  else if (ABS(im)==2) then
   ylmc = SQRT(15.d0/(32.d0*pi))*(sinth)**2*CMPLX(costwophi,sintwophi)
  else 
   MSG_ERROR("wrong im")
  end if

 case (3)
  if (ABS(im)==0) then
   ylmc= SQRT(7.d0/(16.d0*pi))*(5.d0*costh**3 -3.d0*costh)
  else if (ABS(im)==1) then
   ylmc= -SQRT(21.d0/(64.d0*pi))*sinth*(5.d0*costh**2-one)*CMPLX(cosphi,sinphi)
  else if (ABS(im)==2) then
   ylmc= SQRT(105.d0/(32.d0*pi))*sinth**2*costh*CMPLX(costwophi,sintwophi)
  else if (ABS(im)==3) then
   ylmc=-SQRT(35.d0/(64.d0*pi))*sinth**3*CMPLX(costhreephi,sinthreephi)
  else 
   MSG_ERROR("wrong im")
  end if

 case default
  !write(msg,'(a,i6,a,i6)')' The maximum allowed value for l is,',LMAX,' however l=',il
  !MSG_ERROR(msg)
 end select
!
!=== Treat the case im < 0 ===
 if (im < 0) ylmc=(-one)**(im)*CONJG(ylmc)

 ! FIXME: Use the piece of code below as it works for arbitrary (l,m)
 ! the implementation above is buggy when the vector is along z! 
 ! 
#if 0
! Remember the expression of complex spherical harmonics:
! $Y_{lm}(\theta,\phi)=sqrt{{(2l+1) over (4\pi)} {fact(l-m)/fact(l+m)} } P_l^m(cos(\theta)) e^{i m\phi}$
  new_ylmc = SQRT((2*il+1)*factorial(il-ABS(im))/(factorial(il+ABS(im))*four_pi)) * &
&   ass_leg_pol(il,ABS(im),costh) * CMPLX(cosphi,sinphi)**ABS(im)
  if (im<0) new_ylmc=(-one)**(im)*CONJG(new_ylmc)

  if (ABS(new_ylmc-ylmc)>tol6) then
    !MSG_WARNING("Check new_ylmc")
    !write(*,*)"il,im,new_ylmc, ylmc",il,im,new_ylmc,ylmc
    !write(*,*)"fact",SQRT((2*il+1)*factorial(il-ABS(im))/(factorial(il+ABS(im))*four_pi))
    !write(*,*)"costh,sinth,ass_leg_pol",costh,sinth,ass_leg_pol(il,ABS(im),costh) 
    !write(*,*)"cosphi,sinphi,e^{imphi}",cosphi,sinphi,CMPLX(cosphi,sinphi)**ABS(im)
  end if
  ylmc = new_ylmc
#endif

end function ylmc
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/ylmcd
!! NAME
!! ylmcd
!!
!! FUNCTION
!!  Computes dth and dphi, the first derivatives of complex Ylm as a function of 
!!  th and phi (the angles of the spherical coordinates)
!!  It works for all spherical harmonics with l <= 3
!!
!! INPUTS
!!  il=angular quantum number
!!  im=magnetic quantum number
!!  kcart=cartesian coordinates of the vector where the first derivatives of Ylm are evaluated
!!
!! OUTPUT
!!  dth =derivative of Y_lm with respect to \theta
!!  dphi=derivative of Y_lm with respect to \phi
!!
!! NOTES
!!  Note the use of double precision complex.
!!  Case l>3 not implemented.
!!
!! PARENTS
!!      m_commutator_vkbr
!!
!! CHILDREN
!!
!! SOURCE

subroutine ylmcd(il,im,kcart,dth,dphi)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,im
 complex(dpc),intent(out) :: dphi,dth
!arrays
 real(dp),intent(in) :: kcart(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: LMAX=3
 real(dp),parameter :: PPAD=tol8
 real(dp) :: cosphi,costh,costhreephi,costwophi,r,rxy,sinphi,sinth,sinthreephi,sintwophi
 character(len=500) :: msg

! *************************************************************************

 if (ABS(im)>ABS(il))then
   write(msg,'(3(a,i0))')' m is,',im,' however it should be between ',-il,' and ',il
   MSG_ERROR(msg)
 end if

 dphi=czero; dth=czero

 r=SQRT(kcart(1)**2+kcart(2)**2+kcart(3)**2)
 if (r<PPAD) r=r+PPAD
 !$if (r<tol10) RETURN

 rxy=SQRT(kcart(1)**2+kcart(2)**2)
 if (rxy<PPAD) rxy=r+PPAD

! Determine theta and phi 
 costh= kcart(3)/r
#if 1
 ! old buggy coding
 sinth= rxy/r
 cosphi= kcart(1)/rxy
 sinphi= kcart(2)/rxy
#else
 sinth=sqrt(abs((one-costh)*(one+costh))) ! abs is needed to prevent very small negative arg
 cosphi=one
 sinphi=zero
 if (sinth>tol10) then
   cosphi=kcart(1)/(r*sinth)
   sinphi=kcart(2)/(r*sinth)
 end if
#endif

 costwophi= two*cosphi**2 - one
 sintwophi= two*sinphi*cosphi
 costhreephi=cosphi*costwophi-sinphi*sintwophi
 sinthreephi=cosphi*sintwophi+sinphi*costwophi

 select case (il)

 case (0)
   dth  = czero 
   dphi = czero

 case (1)
   if (ABS(im)==0) then
     dth= -SQRT(three/(four_pi))*sinth
     dphi= czero
   else if (abs(im)==1) then
     dth= -SQRT(3.d0/(8.d0*pi))*costh*CMPLX(cosphi,sinphi)
     dphi=-SQRT(3.d0/(8.d0*pi))*sinth*CMPLX(-sinphi,cosphi)
   end if

 case (2)
   if (ABS(im)==0) then
     dth= -SQRT(5.d0/(16.d0*pi))*6.d0*costh*sinth
     dphi= czero
   else if (ABS(im)==1) then
     dth=  -SQRT(15.d0/(8.d0*pi))*(costh**2-sinth**2)*CMPLX(cosphi,sinphi)
     dphi= -SQRT(15.d0/(8.d0*pi))*costh*sinth*(0.d0,1.d0)*CMPLX(cosphi,sinphi)
   else if (abs(im)==2) then
     dth  = SQRT(15.d0/(32.d0*pi))*2.d0*costh*sinth*CMPLX(costwophi,sintwophi)
     dphi = SQRT(15.d0/(32.d0*pi))*sinth**2*(0.d0,2.d0)*CMPLX(costwophi,sintwophi)
   end if

 case (3)
   if (ABS(im)==0) then
     dth = SQRT(7.d0/(16*pi))*(-15.d0*costh**2*sinth + 3.d0**sinth)
     dphi= czero
   else if (ABS(im)==1) then
     dth= -SQRT(21.d0/(64.d0*pi))*CMPLX(cosphi,sinphi)*(5.d0*costh**3-costh-10.d0*sinth**2*costh)
     dphi=-SQRT(21.d0/(64.d0*pi))*sinth*(5.d0*costh**2-1)*(0.d0,1.d0)*CMPLX(cosphi,sinphi)
   else if (ABS(im)==2) then
     dth =SQRT(105.d0/(32.d0*pi))*(2.d0*sinth*costh**2-sinth**3)*CMPLX(costwophi,sintwophi)
     dphi=SQRT(105.d0/(32*pi))*sinth**2*costh*(0.d0,2.d0)*CMPLX(costwophi,sintwophi)
   else if (abs(im)==3) then
     dth =-SQRT(35.d0/(64.d0*pi))*3.d0*sinth**2*costh*CMPLX(costhreephi,sinthreephi)
     dphi= SQRT(35.d0/(64.d0*pi))*sinth**3*(0.d0,3.d0)*CMPLX(costhreephi,sinthreephi)
   end if

 case default
   write(msg,'(2(a,i0))')' The maximum allowed value for l is,',LMAX,' however, l=',il
   MSG_ERROR(msg)
 end select
!
!=== Treat the case im < 0 ===
 if (im<0) then
   dth = (-one)**(im)*CONJG(dth)
   dphi= (-one)**(im)*CONJG(dphi)
 end if

end subroutine ylmcd
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/factorial
!! NAME
!! factorial
!!
!! FUNCTION
!! Calculates N!. Returns a (dp) real.
!!
!! INPUTS
!!   nn=number to use
!!
!! OUTPUT
!!   factorial= n! (real)
!!
!! PARENTS
!!      gaunt,setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

function factorial(nn)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: factorial

!Local variables ---------------------------------------
!scalars
 integer :: ii
 real(dp) :: ff

! *********************************************************************

 ff=one
 do ii=2,nn
   ff=ff*ii
 end do

 factorial=ff

end function factorial
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/binomcoeff
!! NAME
!! factorial
!!
!! FUNCTION
!! Calculates n!/( k!* (n-k)!). Returns a real (dp)
!!
!! INPUTS
!!   nn=number to use
!!
!! OUTPUT
!!   binomcoeff= n!/( k!* (n-k)!)  (real dp)
!!
!! PARENTS
!!      
!!
!! CHILDREN
!!
!! SOURCE

function binomcoeff(n,k)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: n,k
 real(dp) :: binomcoeff

!Local variables ---------------------------------------
!scalars

! *********************************************************************

 binomcoeff=factorial(n)/(factorial(k)*factorial(n-k))

end function binomcoeff
!!***


!----------------------------------------------------------------------

!!****f* m_special_funcs/laguerre
!! NAME
!! laguerre
!!
!! FUNCTION
!! Laguerre(x,n,a). Returns a (dp) real.
!!
!! INPUTS
!!   x position
!!   n order of laguerre polynomial
!!   a 
!!
!! OUTPUT
!!   Laguerre(x,n,a) (dp)
!!
!! PARENTS
!!   
!!
!! CHILDREN
!!   factorial
!!
!! SOURCE

function laguerre(x,n,a) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in),optional :: n,a
 real(dp)                    :: laguerre
 real(dp),intent(in)         :: x


!Local variables ---------------------------------------
!scalars
 integer :: ii, nn, aa
 

!Local variables ---------------------------------------
!arrays
 real(dp),allocatable :: ff(:)

! *********************************************************************

 if (present(n)) then
   nn=n
 else
   nn=1
 end if

 if (present(a)) then
   aa=a
 else
   aa=0
 end if
 allocate(ff(nn+1))
 ff=0.0_dp
 ff=(/ (binomcoeff(nn+aa,nn-ii)*((-1.0_dp)*x)**ii/factorial(ii) ,ii=0,nn) /)
 laguerre=sum(ff)
 
 deallocate(ff)

end function laguerre
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/RadFnH
!! NAME
!! RadFnH
!!
!! FUNCTION
!!  RadFnH(r,n,l,Z) radial function of atomic wavefunction with nuclear charge Z. 
!!  for quantum number n, and l. 
!!  Default: Fe 3d function. Returns a (dp) real.
!!
!! INPUTS
!!   r radius
!!   n principal quantum number
!!   l quantum number
!!
!! OUTPUT
!!  RadFnH(r,n,l,Z) (dp)
!!
!! PARENTS
!!   
!!
!! CHILDREN
!!  Laguerre
!!  factorial
!!
!! SOURCE

function RadFnH(r,n,l,Z)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in),optional :: n,l
 real(dp) :: RadFnH
 real(dp),intent(in) :: r
 real(dp),intent(in),optional :: Z

!Local variables ---------------------------------------
!scalars
 integer   :: nn,ll
 real(dp)  :: ff,rr,ZZ

! *********************************************************************
 
 if (present(n)) then
   nn=n
 else
   nn=3
 end if

 if (present(l)) then
   ll=l
 else
   ll=2
 end if

 if (present(Z)) then
   ZZ=Z
 else
   ZZ=28.0_dp
 end if
 
 rr=ZZ*r/nn
 ff=exp(log(ZZ*1.0_dp)*(3.0_dp/2.0_dp))*2/nn**2
 ff=ff*sqrt(factorial(nn-ll-1)/factorial(nn+ll))*(2*rr)**ll
 RadFnH=ff*exp(-1*rr)*laguerre(2*rr,nn-ll-1,2*ll+1)

end function RadFnH
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/IRadFnH(a,b,n,l,Z)
!! NAME
!! IRadFnH
!!
!! FUNCTION
!!  IRadFnH(a,b,n,l,Z): Integral of radial function of atomic wavefunction between a and b.
!!  recursive programming using simpson's rule 
!!  iteration depth of m=8 corresponds to relative error of 10^(-12).
!!
!! INPUTS
!!   a lower limit for integration
!!   b upper limit for integration
!!   n principal quantum number
!!   l quantum number
!    Z nuclear charge
!!
!! OUTPUT
!!  IRadFnH(a,b,n,l,Z) (dp)
!!
!! PARENTS
!!   
!!
!! CHILDREN
!!  Laguerre
!!  factorial
!!  RadFnH
!!
!! SOURCE

recursive function IRadFnH(a,b,n,l,Z,m) result(x)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in),optional  :: n,l,m
 real(dp),intent(in):: a
 real(dp),intent(in),optional :: b,Z

!Local variables ---------------------------------------
!scalars
 integer   :: nn,ll,mm
 real(dp)  :: h,bb,ZZ,x

! *********************************************************************
 
 if (present(n)) then
   nn=n
 else
   nn=3
 end if

 if (present(l)) then
   ll=l
 else
   ll=2
 end if

 if (present(Z)) then
   ZZ=Z
 else
   ZZ=28
 end if

 if (present(b)) then
   bb=b
 else
   bb=100.0_dp
 end if

 if (present(m)) then
   mm=m
 else
   mm=0
 end if

 h=(bb-a)/2.0_dp
 if (mm<8) then
  !h=2*h/exp(1.0_dp)
  x=IRadFnH(a,a+h,nn,ll,ZZ,mm+1)+IRadFnH(a+h,bb,nn,ll,ZZ,mm+1)
 else 
  x=RadFnH(a,nn,ll,ZZ)**2*a**2+4.0_dp*RadFnH(a+h,nn,ll,ZZ)**2*(a+h)**2
  x=h/3.0_dp*(x+RadFnH(bb,nn,ll,ZZ)**2*bb**2)
 end if 

end function IRadFnH
!!***

!----------------------------------------------------------------------

END MODULE m_special_funcs
!!***
