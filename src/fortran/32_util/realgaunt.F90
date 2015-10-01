!{\src2tex{textfont=tt}}
!!****f* ABINIT/realgaunt
!! NAME
!! realgaunt
!!
!! FUNCTION
!! This routine compute "real Gaunt coefficients", i.e. gaunt
!! coefficients according to "real spherical harmonics"
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  l_max= max. value of ang. momentum l+1;  Gaunt coeffs up to
!!          [(2*l_max-1,m),(l_max,m),(l_max,m)] are computed
!!
!! OUTPUT
!!  gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)=
!!          selection rules for Gaunt coefficients
!!          if Gaunt coeff. is zero, gntselect=0
!!          if Gaunt coeff. is non-zero, gntselect is the index of
!!                           the coeff. in realgnt(:) array
!!  ngnt= number of non-zero Gaunt coefficients
!!  realgnt((2*l_max-1)**2*l_max**4)= non-zero real Gaunt coefficients
!!
!! PARENTS
!!      paw_mkrhox,pawinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine realgaunt(l_max,ngnt,gntselect,realgnt)

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util, except_this_one => realgaunt
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: l_max
 integer,intent(out) :: ngnt
!arrays
 integer,intent(out) :: gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)
 real(dp),intent(out) :: realgnt((2*l_max-1)**2*(l_max)**4)

!Local variables ------------------------------
!scalars
 integer :: ilm1,ilm2,ilmp1,k0lm1,klm1,l1,l2,ll,lp1,m1,m2,mm,mm1,mm2,mm3,mp1
 real(dp) :: c11,c12,c21,c22,c31,c32,fact,realgnt_tmp
!arrays
 integer,allocatable :: ssgn(:)
 type(coeff3_type), allocatable :: coeff(:)

!************************************************************************

!Compute matrix cc where Sl=cc*Yl (Sl=real sph. harm.)
!------------------------------------------------
 allocate(coeff(4*l_max-3))
 do ll=1,4*l_max-3
   allocate (coeff(ll)%value(2,2*ll-1,2*ll-1))
   coeff(ll)%value(:,:,:)=zero
   coeff(ll)%value(1,ll,ll)=one
   do mm=1,ll-1
     coeff(ll)%value(1,ll+mm,ll+mm)= (-1._dp)**mm/sqrt(2._dp)
     coeff(ll)%value(1,ll-mm,ll+mm)= ( 1._dp)    /sqrt(2._dp)
     coeff(ll)%value(2,ll+mm,ll-mm)=-(-1._dp)**mm/sqrt(2._dp)
     coeff(ll)%value(2,ll-mm,ll-mm)= ( 1._dp)    /sqrt(2._dp)
   end do
 end do

 allocate(ssgn(l_max**2));ssgn(:)=1
 if (l_max>0) then
   do l1=1,l_max-1
     ilm1=1+l1**2+l1
     do m1=-l1,-1
       ssgn(ilm1+m1)=-1
     end do
   end do
 end if

 ngnt=0

!Loop on (lp1,mp1)
!------------------------------------------------
 do lp1=0,l_max-1
   do mp1=-lp1,lp1
     ilmp1=1+lp1**2+lp1+mp1
     k0lm1=ilmp1*(ilmp1-1)/2

!    Loop on (l1,m1)<=(lp1,mp1)
!    ------------------------------------------------
     do l1=0,l_max-1
       do m1=-l1,l1
         ilm1=1+l1**2+l1+m1

         if (ilm1<=ilmp1) then

           klm1=k0lm1+ilm1
           gntselect(:,klm1)=0

!          Loop on (l2,m2)
!          ------------------------------------------------
           do l2=abs(l1-lp1),l1+lp1,2
             do m2=-l2,l2
               ilm2=1+l2**2+l2+m2

!              Real Gaunt coeffs selection rules
!              ------------------------------------------------
               if ((l2<=l1+lp1).and.&
&               (((m1== mp1).and.((m2==0).or.(m2==2*abs(mp1)))).or.&
&               ((m1==-mp1).and.(m2==-abs(m1)-abs(mp1))).or.&
&               ((abs(m1)/=(abs(mp1)).and.&
&               ((m2==ssgn(ilm1)*ssgn(ilmp1)*   (abs(m1)+abs(mp1))).or.&
&               (m2==ssgn(ilm1)*ssgn(ilmp1)*abs(abs(m1)-abs(mp1)))&
               ))))) then

!                Compute selected real Gaunt coefficient
!                ------------------------------------------------
                 realgnt_tmp=zero
                 do mm1=-l1,l1
                   c11=coeff(l1+1)%value(1,l1+mm1+1,l1+m1+1)
                   c12=coeff(l1+1)%value(2,l1+mm1+1,l1+m1+1)
                   do mm2= -lp1,lp1
                     c21=coeff(lp1+1)%value(1,lp1+mm2+1,lp1+mp1+1)
                     c22=coeff(lp1+1)%value(2,lp1+mm2+1,lp1+mp1+1)
                     do mm3= -l2,l2
                       c31=coeff(l2+1)%value(1,l2+mm3+1,l2+m2+1)
                       c32=coeff(l2+1)%value(2,l2+mm3+1,l2+m2+1)
                       fact=c11*c21*c31  -  c12*c22*c31&
&                       -c11*c22*c32  -  c12*c21*c32
                       if((abs(fact)>=tol12).and.(mm3==-mm2-mm1)) &
&                       realgnt_tmp=realgnt_tmp+fact*(-1)**mm2 &
&                       *gaunt(l2,mm3,l1,mm1,lp1,-mm2)
                     end do
                   end do
                 end do

!                Count additional non-zero real Gaunt coeffs
!                ------------------------------------------------
                 if (abs(realgnt_tmp)>=tol12) then
                   ngnt=ngnt+1
                   gntselect(ilm2,klm1)=ngnt
                   realgnt(ngnt)=realgnt_tmp/sqrt(four_pi)
                 end if

!                End loops
!                ------------------------------------------------
               end if
             end do
           end do
         end if
       end do
     end do
   end do
 end do

!Deallocate memory
!------------------------------------------------
 do ll=1,4*l_max-3
   deallocate (coeff(ll)%value)
 end do
 deallocate(coeff,ssgn)

end subroutine realgaunt

!!***
