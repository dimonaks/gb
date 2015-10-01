!{\src2tex{textfont=tt}}
!!****f* ABINIT/initylmr
!! NAME
!! initylmr
!!
!! FUNCTION
!! Calculate the real spherical harmonics Ylm (and gradients)
!! over a set of (r) vectors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotential
!!  normchoice=0  the input rr vectors are normalized
!!            =1  the norm of the input vector is in nrm() array
!!  nrm(npts) = Depending of normchoice, this array contains
!!              either the weight of the point or the norm of rr.
!!  npts = number of rr vectors
!!  option= 1=compute Ylm(R), 2=compute Ylm(R) and dYlm/dRi (cartesian derivatives),
!!          3=compute Ylm(R), dYlm/dRi and d2Ylm/dRidRj (cartesian derivatives)
!!  rr(3,npts)=  vectors for which ylmr have to be calculated
!!               For each point of the spherical mesh, gives the
!!               coordinates of the corresponding point.
!!
!! OUTPUT
!!  if (option=1, 2 or 3)
!!    ylm(mpsang*mpsang,npts)     = real spherical harmonics for each r point
!!  if (option=2 or 3)
!!    ylmr_gr(1:3,mpsang*mpsang,npts)= gradients of real spherical harmonics
!!  if (option=3)
!!    ylmr_gr(4:9,mpsang*mpsang,npts)= first and second gradients of real spherical harmonics
!!
!! NOTES
!! Remember the expression of complex spherical harmonics:
!! $Y_{lm}(%theta ,%phi)=sqrt{{(2l+1) over (4 %pi)} {fact(l-m) over fact(l+m)} } P_l^m(cos(%theta)) func e^{i m %phi}$
!! Remember the expression of real spherical harmonics as linear combination of imaginary spherical harmonics:
!! $Yr_{lm}(%theta ,%phi)=(Re{Y_{l-m}}+(-1)^m Re{Y_{lm}})/sqrt{2}
!! $Yr_{l-m}(%theta ,%phi)=(Im{Y_{l-m}}-(-1)^m Im{Y_{lm}})/sqrt{2}
!!
!! PARENTS
!!      debug_tools,denfgr,mlwfovlp_ylmfar,pawfrnhat_recipspace,pawgylm,pawinit
!!      pawmkaewf,pspnl_operat_rec,qijb,qijb_bk,qijb_kk,smatrix_pawinit
!!
!! CHILDREN
!!      plm_coeff
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initylmr(mpsang,normchoice,npts,nrm,option,rr,ylmr,ylmr_gr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_28_numeric_noabirule
 use interfaces_32_util, except_this_one => initylmr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang,normchoice,npts,option
!arrays
 real(dp),intent(in) :: nrm(npts),rr(3,npts)
 real(dp),intent(out) :: ylmr(mpsang*mpsang,npts)
 real(dp),optional,intent(out) :: ylmr_gr(3*(option/2)+6*(option/3),mpsang*mpsang,npts)

!Local variables ------------------------------
!scalars
 integer :: dimgr,ilang,inpt,l0,ll,mm
 real(dp) :: cphi,ctheta,fact,onem,rnorm,sphi,stheta,work1,work2,ylmcst,ylmcst2
 logical :: compute_ylm,compute_ylm2gr,compute_ylmgr
!arrays
 real(dp) :: dphi(3),dtheta(3),iphase(mpsang-1),rphase(mpsang-1)
 real(dp),allocatable :: blm(:,:)

!************************************************************************

!What has to be computed ?
 compute_ylm   = (option==1.or.option==2.or.option==3)
 compute_ylmgr =((             option==2.or.option==3).and.present(ylmr_gr))
 compute_ylm2gr=((                          option==3).and.present(ylmr_gr))
 dimgr=3*(option/2)+6*(option/3)

!Initialisation of spherical harmonics
 if (compute_ylm  ) ylmr   (:  ,1:npts)=zero
 if (compute_ylmgr) ylmr_gr(:,:,1:npts)=zero

!Special case for l=0
 if (compute_ylm  ) ylmr(1,1:npts)=1._dp/sqrt(four_pi)
 if (compute_ylmgr) ylmr_gr(1:dimgr,1,1:npts)=zero
 if (mpsang>1) then

!  Loop over all rr
   do inpt=1,npts

!    Load module of rr
     rnorm=one
     if (normchoice==1) rnorm=nrm(inpt)

!    Continue only for r<>0

     if (rnorm>tol10) then

!      Determine theta and phi
       cphi=one
       sphi=zero
       ctheta=rr(3,inpt)/rnorm
!      MM030519 : abs is needed to prevent very small negative arg
       stheta=sqrt(abs((one-ctheta)*(one+ctheta)))
       if (stheta>tol10) then
         cphi=rr(1,inpt)/(rnorm*stheta)
         sphi=rr(2,inpt)/(rnorm*stheta)
       end if
       do mm=1,mpsang-1
         rphase(mm)=dreal(dcmplx(cphi,sphi)**mm)
         iphase(mm)=aimag(dcmplx(cphi,sphi)**mm)
       end do

!      Determine gradients of theta and phi
       if (compute_ylmgr) then
         dtheta(1)=ctheta*cphi
         dtheta(2)=ctheta*sphi
         dtheta(3)=-stheta
         dphi(1)=-sphi
         dphi(2)=cphi
         dphi(3)=zero
       end if

!      COMPUTE Ylm(R)
       if (compute_ylm) then
!        Loop over angular momentum l
         do ilang=2,mpsang
           ll=ilang-1
           l0=ll**2+ll+1
           fact=1._dp/real(ll*(ll+1),dp)
           ylmcst=sqrt(real(2*ll+1,dp)/four_pi)
!          Special case m=0
           ylmr(l0,inpt)=ylmcst*ass_leg_pol(ll,0,ctheta)
!          Compute for m>0
           onem=one
           do mm=1,ll
             onem=-onem
             work1=ylmcst*sqrt(fact)*onem*ass_leg_pol(ll,mm,ctheta)*sqrt(2._dp)
             ylmr(l0+mm,inpt)=work1*rphase(mm)
             ylmr(l0-mm,inpt)=work1*iphase(mm)
             if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
           end do ! End loop over m
         end do  ! End loop over l
       end if

!      COMPUTE dYlm/dRi
       if (compute_ylmgr) then
!        Loop over angular momentum l
         do ilang=2,mpsang
           ll=ilang-1
           l0=ll**2+ll+1
           fact=1._dp/real(ll*(ll+1),dp)
           ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/rnorm
!          Special case m=0
           work1=ylmcst*plm_dtheta(ll,0,ctheta)
           ylmr_gr(1:3,l0,inpt)=work1*dtheta(1:3)
!          Compute for m>0
           onem=one
           do mm=1,ll
             onem=-onem
             work1=ylmcst*sqrt(fact)*onem*plm_dtheta(ll,mm,ctheta)*sqrt(2._dp)
             work2=ylmcst*sqrt(fact)*onem*plm_dphi  (ll,mm,ctheta)*sqrt(2._dp)
             ylmr_gr(1:3,l0+mm,inpt)=rphase(mm)*work1*dtheta(1:3)-iphase(mm)*work2*dphi(1:3)
             ylmr_gr(1:3,l0-mm,inpt)=iphase(mm)*work1*dtheta(1:3)+rphase(mm)*work2*dphi(1:3)
             if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
           end do ! End loop over m
         end do  ! End loop over l
       end if

!      COMPUTE d2Ylm/dRidRj
       if (compute_ylm2gr) then
         allocate(blm(5,mpsang*mpsang))
         call plm_coeff(blm,mpsang,ctheta)

!        Loop over angular momentum l
         do ilang=2,mpsang
           ll=ilang-1
           l0=ll**2+ll+1
           fact=1._dp/real(ll*(ll+1),dp)
           ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/(rnorm**2)
!          Special case m=0
           ylmr_gr(4,l0,inpt)=ylmcst*(-blm(3,l0)*sphi*sphi+blm(4,l0)*cphi*cphi)
           ylmr_gr(5,l0,inpt)=ylmcst*(-blm(3,l0)*cphi*cphi+blm(4,l0)*sphi*sphi)
           ylmr_gr(6,l0,inpt)=ylmcst*blm(1,l0)
           ylmr_gr(7,l0,inpt)=ylmcst*blm(2,l0)*sphi
           ylmr_gr(8,l0,inpt)=ylmcst*blm(2,l0)*cphi
           ylmr_gr(9,l0,inpt)=ylmcst*(blm(3,l0)+blm(4,l0))*sphi*cphi
!          Compute for m>0
           onem=one
           do mm=1,ll
             onem=-onem;ylmcst2=ylmcst*sqrt(fact)*sqrt(two)
             ylmr_gr(4,l0+mm,inpt)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*rphase(mm)-&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
             ylmr_gr(4,l0-mm,inpt)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*iphase(mm)+&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
             ylmr_gr(5,l0+mm,inpt)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*rphase(mm)+&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
             ylmr_gr(5,l0-mm,inpt)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*iphase(mm)-&
&             blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
             ylmr_gr(6,l0+mm,inpt)=ylmcst2*blm(1,l0+mm)*rphase(mm)
             ylmr_gr(6,l0-mm,inpt)=ylmcst2*blm(1,l0+mm)*iphase(mm)
             ylmr_gr(7,l0+mm,inpt)=ylmcst2*(blm(2,l0+mm)*sphi*rphase(mm)+&
&             mm*iphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(7,l0-mm,inpt)=ylmcst2*(blm(2,l0+mm)*sphi*iphase(mm)-&
&             mm*rphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(8,l0+mm,inpt)=ylmcst2*(blm(2,l0+mm)*cphi*rphase(mm)-&
&             mm*iphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(8,l0-mm,inpt)=ylmcst2*(blm(2,l0+mm)*cphi*iphase(mm)+&
&             mm*rphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
             ylmr_gr(9,l0+mm,inpt)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*rphase(mm)-&
&             blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*iphase(mm))
             ylmr_gr(9,l0-mm,inpt)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*iphase(mm)+&
&             blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*rphase(mm))
             if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
           end do ! End loop over m
         end do  ! End loop over l
         deallocate(blm)
       end if

!      End condition r<>0
     end if

!    End loop over rr
   end do

!  End condition l<>0
 end if

end subroutine initylmr
!!***
