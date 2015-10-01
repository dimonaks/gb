!{\src2tex{textfont=tt}}
!!****f* ABINIT/mat_slm2ylm
!! NAME
!! mat_slm2ylm
!!
!! FUNCTION
!! For a given angular momentum lcor, change a matrix  of dimension (2*lcor+1)
!! from the Slm to the Ylm basis if option==1 or from Ylm to Slm if !option==2
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!  mat_inp_c= Input matrix
!!  ndij= ndij = 4
!!  option= -1  Change matrix from Slm to Ylm basis
!!           1  Change matrix from Ylm to Slm basis
!!  optspin=  1  Spin up are first
!!            2  Spin dn are first
!!  prtvol=printing volume
!!
!! OUTPUT
!!  mat_inp_c= Output matrix in Ylm or Slm basis according to option
!!
!! NOTES
!!  usefull only in ndij==4
!!
!! PARENTS
!!      setnoccmmp,pawlsylm
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mat_slm2ylm(lcor,mat_inp_c,mat_out_c,ndij,option,optspin,prtvol)

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ndij,lcor,option,optspin,prtvol
!arrays
 complex(dpc) :: mat_inp_c(2*lcor+1,2*lcor+1,ndij),mat_out(2*lcor+1,2*lcor+1,ndij),mat_out_c(2*lcor+1,2*lcor+1,ndij)
 character(len=9),parameter :: dspinc(6)=(/"up       ","down     ","up-up    ","down-down","up-dn    ","dn-up    "/)! optspin 1
 character(len=9),parameter :: dspinc2(6)=(/"up       ","down     ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)! optspin 2

!Local variables ---------------------------------------
!scalars
 integer :: jm,ii,jj,ll,mm,ispden,im,im1,im2
 real(dp),parameter :: invsqrt2=one/sqrt2
 character(len=500) :: message
 real(dp) :: onem
 complex(dpc) :: tmp2
!arrays
 complex(dpc),allocatable :: slm2ylm(:,:)

! *********************************************************************

 if(abs(prtvol)>2) then
   write(message,'(3a)') ch10, "   mat_slm2ylm"
   call wrtout(std_out,message,'COLL')
 end if
 if(ndij/=4) then
   write(message,'(3a)') ch10,"ndij/=4: BUG"
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(abs(prtvol)>2) then
   if(option==1.or.option==3) then
     write(message,'(3a)') ch10,"matrix in Slm basis is changed into Ylm basis"
     call wrtout(std_out,message,'COLL')
   else if(option==2.or.option==4) then
     write(message,'(3a)') ch10,"matrix in Ylm basis is changed into Slm basis"
     call wrtout(std_out,message,'COLL')
   else
     write(message,'(3a)') ch10,"option=/1 or 2 or 3 or 4 : BUG"
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
 ll=lcor
 allocate(slm2ylm(2*ll+1,2*ll+1));slm2ylm=czero
 mat_out=zero
 mat_out_c=czero
 do im=1,2*ll+1
   mm=im-ll-1;jm=-mm+ll+1
   onem=dble((-1)**mm)
   if (mm> 0) then
     slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
     slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
   end if
   if (mm==0) then
     slm2ylm(im,im)=cone
   end if
   if (mm< 0) then
     slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
     slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
   end if
 end do
 if(abs(prtvol)>2) then
   do ispden=1,ndij
     if(optspin==1) then
       if(option==1.or.option==3)&
&       write(message,'(3a)') ch10,&
&       "Input matrix in the Slm basis for component ",trim(dspinc(ispden+2*(ndij/4)))
       if(option==2.or.option==3)&
&       write(message,'(3a)') ch10,&
&       "Input matrix in the Ylm basis for component ",trim(dspinc(ispden+2*(ndij/4)))
     else
       if(option==1.or.option==3)&
&       write(message,'(3a)') ch10,&
&       "Input matrix in the Slm basis for component ",trim(dspinc2(ispden+2*(ndij/4)))
       if(option==2.or.option==3)&
&       write(message,'(3a)') ch10,&
&       "Input matrix in the Ylm basis for component ",trim(dspinc2(ispden+2*(ndij/4)))
     end if
     call wrtout(std_out,message,'COLL')
     do im1=1,lcor*2+1
       write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
&       (mat_inp_c(im1,im2,ispden),im2=1,lcor*2+1)
       call wrtout(std_out,message,'COLL')
     end do
   end do
 end if
 do ispden=1,ndij
   do jm=1,2*ll+1
     do im=1,2*ll+1
       tmp2=czero
       do ii=1,2*ll+1
         do jj=1,2*ll+1
           if(option==1) then
             tmp2=tmp2+mat_inp_c(ii,jj,ispden)*(slm2ylm(im,ii))*CONJG(slm2ylm(jm,jj))
           else if(option==2) then
             tmp2=tmp2+mat_inp_c(ii,jj,ispden)*CONJG(slm2ylm(ii,im))*(slm2ylm(jj,jm))
           end if
         end do
       end do
       mat_out_c(im,jm,ispden)=tmp2
     end do
   end do
 end do ! ispden
 do ii=1,2*ll+1
   do jj=1,2*ll+1
     mat_out(ii,jj,1)=real(mat_out_c(ii,jj,1))
     mat_out(ii,jj,2)=real(mat_out_c(ii,jj,2))
     mat_out(ii,jj,3)=real(mat_out_c(ii,jj,3))
     mat_out(ii,jj,4)=aimag(mat_out_c(ii,jj,3))
!    check that n_{m,m'}^{alpha,beta}=conjg(n_{m',m"}^{beta,alpha}).
     if((abs(aimag(mat_out_c(ii,jj,3))+aimag(mat_out_c(jj,ii,4))).ge.0.0001).or. &
&     (abs(real(mat_out_c(ii,jj,3))-real(mat_out_c(jj,ii,4))).ge.0.0001)) then
       write(message,'(2a,4f10.4)') ch10,&
&       "prb with mat_out_c BUG",mat_out_c(ii,jj,3),mat_out_c(ii,jj,4)
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end do
 end do
 deallocate(slm2ylm)
 end subroutine mat_slm2ylm

!!***
