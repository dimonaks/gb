!{\src2tex{textfont=tt}}
!!****f* ABINIT/mat_mlms2jmj
!! NAME
!! mat_mlms2jmj
!!
!! FUNCTION
!! For a given angular momentum lcor, change a matrix  of dimension 2(2*lcor+1)
!! from the Ylm basis to the J,M_J basis if option==1
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum
!!  ndij= ndij = 4
!!  option=  1 matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis
!!           2 matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis
!!  optspin=  1  Spin up are first
!!            2  Spin dn are first
!!  prtvol=printing volume
!!
!! SIDE EFFECTS
!!  mat_mlms= Input/Ouput matrix in the Ylm basis, size of the matrix is (2*lcor+1,2*lcor+1,ndij)
!!  mat_jmj= Input/Output matrix in the J,M_J basis, size is 2*(2*lcor+1),2*(2*lcor+1)
!!
!! NOTES
!!  usefull only in ndij==4
!!
!! PARENTS
!!      pawlsylm,pawprt,setnoccmmp
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
subroutine mat_mlms2jmj(lcor,mat_mlms,mat_jmj,ndij,option,optspin,prtvol)

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
 complex(dpc),intent(inout) :: mat_mlms(2*lcor+1,2*lcor+1,ndij),mat_jmj(2*(2*lcor+1),2*(2*lcor+1))

!Local variables ---------------------------------------
!scalars
 character(len=9),parameter :: dspinold(6)=(/"up       ","down     ","up-up    ","down-down","up-dn    ","dn-up    "/)
 character(len=9),parameter :: dspin(6)=(/"dn       ","up       ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)
 integer ::ii,im,im1,im2,ispden,jc1,jc2,jj,jm,ll,ml1,ml2,ms1,ms2
 real(dp),parameter :: invsqrt2=one/sqrt2
 complex(dpc),allocatable :: mat_mlms2(:,:)
 integer, allocatable :: ind_msml(:,:)
 complex(dpc),allocatable :: mat_tmp(:,:,:)
 character(len=500) :: message
 real(dp) :: invsqrt2lp1
 complex(dpc) :: tmp2
! *************************************************************************
!arrays
 complex(dpc),allocatable :: mlms2jmj(:,:)
 real(dp) :: xj,xmj

!*********************************************************************
 if(ndij/=4) then
   write(message,'(3a)') ch10,"ndij/=4: BUG"
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(option==1) then
   write(message,'(3a)') ch10,&
&   "matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis"
   call wrtout(std_out,message,'COLL')
 else if(option==2) then
   write(message,'(3a)') ch10,&
&   "matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis"
   call wrtout(std_out,message,'COLL')
 else
   write(message,'(3a)') ch10,"option=/1 and =/2: BUG"
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(option==1) then
   if(optspin==2) then
     if(abs(prtvol)>2)&
&     write(message,'(3a)') ch10,"assume spin dn is the first in the array"
   else if (optspin==1) then
     allocate(mat_tmp(2*lcor+1,2*lcor+1,ndij))
     if(abs(prtvol)>2)&
&     write(message,'(3a)') ch10,"change array in order that spin dn is the first in the array"
     mat_tmp(:,:,1)=mat_mlms(:,:,2)
     mat_tmp(:,:,2)=mat_mlms(:,:,1)
     mat_tmp(:,:,3)=mat_mlms(:,:,4)
     mat_tmp(:,:,4)=mat_mlms(:,:,3)
     mat_mlms(:,:,:)=mat_tmp(:,:,:)
     deallocate(mat_tmp)
   else
     write(message,'(3a)') ch10,"optspin=/1 and =/2: BUG"
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   call wrtout(std_out,message,'COLL')
 end if

 if(option==1.and.abs(prtvol)>2) then
   do ispden=1,ndij
     write(message,'(3a)') ch10,&
&     "Input matrix in the Ylm basis for component ",trim(dspin(ispden+2*(ndij/4)))
     call wrtout(std_out,message,'COLL')
     do im1=1,lcor*2+1
       write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&       (mat_mlms(im1,im2,ispden),im2=1,lcor*2+1)
       call wrtout(std_out,message,'COLL')
     end do
   end do
 end if ! option==1
!--------------- Built indices + allocations
 ll=lcor
 allocate(mlms2jmj(2*(2*ll+1),2*(2*ll+1)));mlms2jmj=czero
 allocate(ind_msml(2,-ll:ll))
 allocate(mat_mlms2(2*(2*lcor+1),2*(2*lcor+1)))
 mlms2jmj=czero
 jc1=0
 do ms1=1,2
   do ml1=-ll,ll
     jc1=jc1+1
     ind_msml(ms1,ml1)=jc1
!    write(6,*) "ind_msml",ms1,ml1,ind_msml(ms1,ml1)
   end do
 end do
!--------------- Change representation of input matrix for ndij==4
 if(option==1) then
   jc1=0
   do ms1=1,2
     do ml1=1,2*ll+1
       jc1=jc1+1
       jc2=0
       do ms2=1,2
         do ml2=1,2*ll+1
           jc2=jc2+1
           if(ms1==ms2) mat_mlms2(jc1,jc2)=mat_mlms(ml1,ml2,ms1)
           if(ms1<ms2) mat_mlms2(jc1,jc2)=mat_mlms(ml1,ml2,3)
           if(ms1>ms2) mat_mlms2(jc1,jc2)=mat_mlms(ml1,ml2,4)
         end do
       end do
     end do
   end do
   if(abs(prtvol)>1) then
     write(message,'(3a)') ch10,"Input matrix in the lms basis for all component"
     call wrtout(std_out,message,'COLL')
     do im1=1,2*(lcor*2+1)
       write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))')&
&       (mat_mlms2(im1,im2),im2=1,2*(lcor*2+1))
       call wrtout(std_out,message,'COLL')
     end do
   end if
 end if  ! option==1
!--------------- built mlms2jmj
!do jj=ll,ll+1    ! the physical value of j are ll-0.5,ll+0.5
!xj(jj)=jj-0.5
 if(ll==0)then
   write(message,'(a,a,a)')' mat_mlms2jmj : BUG - ',ch10,&
&   ' ll should not be equal to zero '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 jc1=0
 invsqrt2lp1=one/sqrt(float(2*lcor+1))
 do jj=ll,ll+1
   xj=float(jj)-half
!  do xj=float(ll)-0.5,float(ll)+0.5,1    !! NON F90 STANDARD
   do jm=-jj,jj-1
     xmj=float(jm)+half
!    do xmj=-xj,xj,1   !! NON F90 STANDARD
!    xmj(jm)=jm+0.5
!    write(6,*) "indices",xj,xmj,jc1
     jc1=jc1+1
     if(nint(xj+0.5)==ll+1) then
       if(nint(xmj+0.5)==ll+1)  then
         mlms2jmj(ind_msml(2,ll),jc1)=1.0   !  J=L+0.5 and m_J=L+0.5
!        write(6,*) "1",jc1,ind_msml(2,ll)
       else if(nint(xmj-0.5)==-ll-1) then
         mlms2jmj(ind_msml(1,-ll),jc1)=1.0   !  J=L+0.5 and m_J=-L-0.5
!        write(6,*) "2",jc1,ind_msml(1,-ll)
       else
         mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
         mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
!        write(6,*) "3",jc1,ind_msml(1,nint(xmj+0.5))
!        write(6,*) "4",jc1,ind_msml(2,nint(xmj-0.5))
       end if
     end if
     if(nint(xj+0.5)==ll) then
       mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
       mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=-invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
!      write(6,*) "5",jc1,ind_msml(1,nint(xmj+0.5))
!      write(6,*) "6",jc1,ind_msml(2,nint(xmj-0.5))
     end if
   end do
 end do
 if(abs(prtvol)>2) then
   write(message,'(3a)') ch10,"Matrix to go from |M_L,M_S> to |J,M_J>"
   call wrtout(std_out,message,'COLL')
   do im1=1,2*(lcor*2+1)
     write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mlms2jmj(im1,im2),im2=1,2*(lcor*2+1))
     call wrtout(std_out,message,'COLL')
   end do
 end if

 do jm=1,2*(2*ll+1)
   do im=1,2*(2*ll+1)
     tmp2=czero
     do ii=1,2*(2*ll+1)
       do jj=1,2*(2*ll+1)
         if(option==1) then
           tmp2=tmp2+mat_mlms2(ii,jj)*CONJG(mlms2jmj(ii,im))*(mlms2jmj(jj,jm))
         else if(option==2) then
           tmp2=tmp2+mat_jmj(ii,jj)*(mlms2jmj(im,ii))*CONJG(mlms2jmj(jm,jj)) ! inv=t*
         end if
       end do
     end do
     if(option==1) then
       mat_jmj(im,jm)=tmp2
     else if(option==2) then
       mat_mlms2(im,jm)=tmp2
     end if
   end do
 end do
 if(option==1.and.abs(prtvol)>=1) then
   write(message,'(3a)') ch10," Matrix in the J,M_J basis"
   call wrtout(std_out,message,'COLL')
   do im1=1,2*(lcor*2+1)
     write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mat_jmj(im1,im2),im2=1,2*(lcor*2+1))
     call wrtout(std_out,message,'COLL')
   end do
 else if(option==2.and.abs(prtvol)>=1) then
   write(message,'(3a)') ch10," Matrix in the m_s m_l basis"
   call wrtout(std_out,message,'COLL')
   do im1=1,2*(lcor*2+1)
     write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mat_mlms2(im1,im2),im2=1,2*(lcor*2+1))
     call wrtout(std_out,message,'COLL')
   end do
   jc1=0
   do ms1=1,2
     do ml1=1,2*ll+1
       jc1=jc1+1
       jc2=0
       do ms2=1,2
         do ml2=1,2*ll+1
           jc2=jc2+1
           if(ms1==ms2) mat_mlms(ml1,ml2,ms1)=mat_mlms2(jc1,jc2)
           if(ms1<ms2) mat_mlms(ml1,ml2,3)=mat_mlms2(jc1,jc2)
           if(ms1>ms2) mat_mlms(ml1,ml2,4)=mat_mlms2(jc1,jc2)
         end do
       end do
     end do
   end do
 end if
 deallocate(mlms2jmj,ind_msml)
 deallocate(mat_mlms2)

 end subroutine mat_mlms2jmj
!!***
