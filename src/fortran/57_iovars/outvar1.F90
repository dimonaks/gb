!{\src2tex{textfont=tt}}
!!****f* ABINIT/outvar1
!! NAME
!! outvar1
!!
!! FUNCTION
!! Echo variables between acell and natom (by alphabetic order)
!! for the ABINIT code.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  jdtset_(0:ndtset_alloc)=actual index of the dataset (equal to dtsets(:)%jdtset)
!!  mxgw_nqlwl=maximal value of input nqptdm for all the datasets
!!  mxmband=maximum number of bands
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnatpawu=maximal value of number of atoms on which +U is applied for all the datasets
!!  mxnatsph=maximal value of input natsph for all the datasets
!!  mxnatvshift=maximal value of input natvshift for all the datasets
!!  mxnimage=maximal value of input nimage for all the datasets
!!  mxnkptgw=maximal value of input nkptgw for all the datasets
!!  mxnkpt=maximal value of input nkpt for all the datasets
!!  mxnqptdm=maximal value of input nqptdm for all the datasets
!!  mxnsppol=maximal value of input nsppol for all the datasets
!!  mxnsym=maximum number of symmetries
!!  mxntypat=maximum number of type of atoms
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set. Use for most dimensioned arrays.
!!  npsp=number of pseudopotentials
!!  nqptdm=the number of q vectors provided by the user to calculate DM in GW
!!  prtvol_glob= if 0, minimal output volume, if 1, no restriction.
!!  response= 1 if response variables must be output, 0 otherwise.
!!  response_(0:ndtset_alloc)= 1 if response variables must be output, 0 otherwise,
!!   for different datasets
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!! The lines of code needed to output the defaults are preserved
!! (see last section of the routine, but are presently disabled)
!!
!!  Note that acell, occ, rprim, xred and vel might have been modified by the
!!  computation, so that their values if choice=1 or choice=2 will differ.
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!      prttagm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outvar1 (choice,dtsets,iout,istatr,istatshft,&
& jdtset_,mxgw_nqlwl,mxmband,mxnatom,mxnatsph,mxnatvshift,mxnimage,mxnkptgw,mxnkpt,mxnqptdm,mxnsppol,mxnsym,mxntypat,&
& mxnatpawu,ndtset,ndtset_alloc,npsp,prtvol_glob,response,response_,results_out,usepaw)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_57_iovars, except_this_one => outvar1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,iout,istatr,istatshft,mxmband,mxgw_nqlwl,mxnatom,mxnatsph,mxnatvshift
 integer,intent(in) :: mxnimage,mxnkpt,mxnkptgw,mxnqptdm,mxnsppol,mxnsym,mxntypat,mxnatpawu,ndtset
 integer,intent(in) :: ndtset_alloc,npsp,prtvol_glob,response,usepaw
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc),response_(ndtset_alloc)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
 character(len=*), parameter :: format01110 ="(1x,a16,1x,(t20,8i8) )"
 character(len=*), parameter :: format01150 ="(1x,a16,1x,(t20,3es16.8))"
 character(len=*), parameter :: format01150a="(1x,a16,a,1x,(t20,3es16.8))"
 character(len=*), parameter :: format01155 ="(1x,a16,1x,(t20,10i5))"
 character(len=*), parameter :: format01155a="(1x,a16,a,1x,(t20,10i5))"
 character(len=*), parameter :: format01160 ="(1x,a16,1x,(t20,3es18.10)) "
 character(len=*), parameter :: format01160a="(1x,a16,a,1x,(t20,3es18.10)) "
 character(len=*), parameter :: format01170 ="(1x,a16,1x,(t20,5f11.6)) "
 character(len=*), parameter :: format01170a="(1x,a16,a,1x,(t20,5f11.6)) "
!scalars
 integer,parameter :: nkpt_max=50
 integer :: allowed,first,iatom,idtset,ii,iimage,ikpt,iscf,istatr_defo
 integer :: istatshft_defo,jdtset,kptopt,marr
 integer :: multi_natom,multi_natfix,multi_natfixx
 integer :: multi_natfixy,multi_natfixz,multi_natsph,multi_natvshift,multi_nberry,multi_nimage,multi_nkpt
 integer :: multi_nkptgw,multi_nqptdm,multi_nshiftk,multi_ntypalch
 integer :: multi_nsppol,multi_gw_nqlwl
 integer :: multi_ntypat,natfix,natfixx,natfixy,natfixz,natom,natsph
 integer :: ndtset_kptopt,nkpt_eff,nimage,nqpt
 integer :: nshiftk,ntypalch,ntypat,tnkpt
 real(dp) :: cpus,kpoint
 character(len=2) :: appen
 character(len=4) :: stringimage
 character(len=16) :: keywd
!arrays
 integer,allocatable :: iatfixio_(:,:),iatfixx_(:,:),iatfixy_(:,:)
 integer,allocatable :: iatfixz_(:,:),intarr(:,:),istwfk_2(:,:)
 integer,allocatable :: jdtset_kptopt(:),natfix_(:),natfixx_(:),natfixy_(:)
 integer,allocatable :: natfixz_(:)
 real(dp),allocatable :: dprarr(:,:)
 character(len=8),allocatable :: strimg(:)

! *************************************************************************

!DEBUG
!write(6,*)' outvar1 : enter '
!ENDDEBUG
!
!Must treat separately the translation of iatfix from the internal
!representation to the input/output representation
 allocate(natfix_(0:ndtset_alloc),iatfixio_(mxnatom,0:ndtset_alloc))
 allocate(natfixx_(0:ndtset_alloc),iatfixx_(mxnatom,0:ndtset_alloc))
 allocate(natfixy_(0:ndtset_alloc),iatfixy_(mxnatom,0:ndtset_alloc))
 allocate(natfixz_(0:ndtset_alloc),iatfixz_(mxnatom,0:ndtset_alloc))
 natfix_(0:ndtset_alloc)=0 ; iatfixio_(:,0:ndtset_alloc)=0
 natfixx_(0:ndtset_alloc)=0 ; iatfixx_(:,0:ndtset_alloc)=0
 natfixy_(0:ndtset_alloc)=0 ; iatfixy_(:,0:ndtset_alloc)=0
 natfixz_(0:ndtset_alloc)=0 ; iatfixz_(:,0:ndtset_alloc)=0
 do idtset=1,ndtset_alloc
!  DEBUG
!  write(6,*)' outvar1 : iatfix_ for idtset= ',idtset
!  ENDDEBUG
   do iatom=1,dtsets(idtset)%natom
!    First look whether the atom is fixed along the three directions
     if( dtsets(idtset)%iatfix(1,iatom)+ &
&     dtsets(idtset)%iatfix(2,iatom)+ &
&     dtsets(idtset)%iatfix(3,iatom)   ==3 )then
       natfix_(idtset)=natfix_(idtset)+1
!      DEBUG
!      write(6,*)' outvar1: iatom,natfix_(idtset)=',iatom,natfix_(idtset)
!      ENDDEBUG
       iatfixio_(natfix_(idtset),idtset)=iatom
     else
!      Now examine each direction, one at a time
       if( dtsets(idtset)%iatfix(1,iatom) ==1)then
         natfixx_(idtset)=natfixx_(idtset)+1
         iatfixx_(natfixx_(idtset),idtset)=iatom
       end if
       if( dtsets(idtset)%iatfix(2,iatom) ==1)then
         natfixy_(idtset)=natfixy_(idtset)+1
         iatfixy_(natfixy_(idtset),idtset)=iatom
       end if
       if( dtsets(idtset)%iatfix(3,iatom) ==1)then
         natfixz_(idtset)=natfixz_(idtset)+1
         iatfixz_(natfixz_(idtset),idtset)=iatom
       end if
     end if
   end do
!  DEBUG
!  write(6,*)' natfix ...'
!  write(6,*)natfix_(idtset),natfixx_(idtset),natfixy_(idtset),natfixz_(idtset)
!  ENDDEBUG
 end do

 allocate(strimg(mxnimage))
 do iimage=1,mxnimage
   if(iimage<10)then
     write(stringimage,'(i1)')iimage
   else if(iimage<100)then
     write(stringimage,'(i2)')iimage
   else if(iimage<1000)then
     write(stringimage,'(i3)')iimage
   else if(iimage<10000)then
     write(stringimage,'(i4)')iimage
   end if
   strimg(iimage)='_'//trim(stringimage)//'img'
 end do
 strimg(1)=''

!Maximal size of dprarr and intarr arrays
 marr=max(3*mxnatom,3*mxnkptgw,mxnkpt*mxnsppol*mxmband,3*mxnkpt,npsp,mxntypat,3*mxnqptdm,3*mxgw_nqlwl,&
& 9*mxnsym,mxnatsph,mxnimage,mxnatvshift*mxnsppol*mxnatom)
 allocate(dprarr(marr,0:ndtset_alloc))
 allocate(intarr(marr,0:ndtset_alloc))

!Set up dimensions : determine whether these are different for different
!datasets.

 multi_gw_nqlwl=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%gw_nqlwl/=dtsets(idtset)%gw_nqlwl)multi_gw_nqlwl=1
   end do
 end if

 multi_nqptdm=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nqptdm/=dtsets(idtset)%nqptdm)multi_nqptdm=1
   end do
 end if

 multi_natfix=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfix_(1)/=natfix_(idtset))multi_natfix=1
   end do
 end if
 if(multi_natfix==0)natfix=natfix_(1)

 multi_natfixx=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixx_(1)/=natfixx_(idtset))multi_natfixx=1
   end do
 end if
 if(multi_natfixx==0)natfixx=natfixx_(1)

 multi_natfixy=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixy_(1)/=natfixy_(idtset))multi_natfixy=1
   end do
 end if
 if(multi_natfixy==0)natfixy=natfixy_(1)

 multi_natfixz=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixz_(1)/=natfixz_(idtset))multi_natfixz=1
   end do
 end if
 if(multi_natfixz==0)natfixz=natfixz_(1)

 multi_natom=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%natom/=dtsets(idtset)%natom)multi_natom=1
   end do
 end if
 if(multi_natom==0)natom=dtsets(1)%natom

 multi_natsph=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%natsph/=dtsets(idtset)%natsph)multi_natsph=1
   end do
 end if
 if(multi_natsph==0)natsph=dtsets(1)%natsph

 multi_natvshift=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%natvshift/=dtsets(idtset)%natvshift)multi_natvshift=1
   end do
 end if


 multi_nberry=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nberry/=dtsets(idtset)%nberry)multi_nberry=1
   end do
 end if

 multi_nimage=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nimage/=dtsets(idtset)%nimage)multi_nimage=1
   end do
 end if
 if(multi_nimage==0)nimage=dtsets(1)%nimage

 multi_nkptgw=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nkptgw/=dtsets(idtset)%nkptgw)multi_nkptgw=1
   end do
 end if

 multi_nkpt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nkpt/=dtsets(idtset)%nkpt)multi_nkpt=1
   end do
 end if

 multi_nshiftk=0
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   first=0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>=1)then
       if(first==0)then
         first=1
         nshiftk=dtsets(idtset)%nshiftk
       else
         if(nshiftk/=dtsets(idtset)%nshiftk)multi_nshiftk=1
       end if
     end if
   end do
 end if

 multi_nsppol=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nsppol/=dtsets(idtset)%nsppol)multi_nsppol=1
   end do
 end if

 multi_ntypalch=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%ntypalch/=dtsets(idtset)%ntypalch)multi_ntypalch=1
   end do
 end if
 if(multi_ntypalch==0)ntypalch=dtsets(1)%ntypalch

 multi_ntypat=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%ntypat/=dtsets(idtset)%ntypat)multi_ntypat=1
   end do
 end if
 if(multi_ntypat==0)ntypat=dtsets(1)%ntypat

!Print each variable, one at a time

 intarr(1,:)=dtsets(:)%accesswff
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'accesswff','INT')

 if(mxnimage==1)then
   do idtset=0,ndtset_alloc
     dprarr(1:3,idtset)=results_out(idtset)%acell(:,1)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,ndtset_alloc,'acell','LEN')
 else
   do idtset=1,ndtset_alloc
     do iimage=1,dtsets(idtset)%nimage
       keywd='acell'//trim(strimg(iimage))
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01160a)trim(keywd),appen,results_out(idtset)%acell(:,iimage)
       else
         write(iout,format01160)trim(keywd),results_out(idtset)%acell(:,iimage)
       end if
     end do
   end do
 end if


!DEBUG
!write(6,*)' outvar1 : before algalch'
!ENDDEBUG
!
!algalch
 if(multi_ntypalch==0)then
   do idtset=0,ndtset_alloc
     intarr(1:ntypalch,idtset)=dtsets(idtset)%algalch(1:ntypalch)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypalch,ndtset_alloc,'algalch','INT')
 else
   do idtset=1,ndtset_alloc
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01155a)'algalch',appen,dtsets(idtset)%algalch(1:dtsets(idtset)%ntypalch)
     else
       write(iout,format01155)'algalch',dtsets(idtset)%algalch(1:dtsets(idtset)%ntypalch)
     end if
   end do
 end if

!atvshift
 if(usepaw>0)then
   if(multi_natvshift==0 .and. multi_nsppol==0 .and. multi_natom==0)then
     if(dtsets(1)%natvshift/=0)then
       do idtset=0,ndtset_alloc
         dprarr(1:dtsets(1)%natvshift*dtsets(1)%nsppol*mxnatpawu,idtset)=&
&         reshape(dtsets(idtset)%atvshift(1:dtsets(1)%natvshift,1:dtsets(1)%nsppol,1:mxnatpawu),&
&         (/dtsets(1)%natvshift*dtsets(1)%nsppol*mxnatpawu/) )
       end do
       call prttagm(dprarr,intarr,iout,jdtset_,-5,marr, &
&       dtsets(1)%natvshift*dtsets(1)%nsppol*mxnatpawu,ndtset_alloc,'atvshift','DPR')
     end if
   else
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%natvshift/=0)then
         if(ndtset>0)then
           jdtset=jdtset_(idtset)
           if(jdtset<10)write(appen,'(i1)')jdtset
           if(jdtset>=10)write(appen,'(i2)')jdtset
           write(iout,format01170a)'atvshift',appen,&
&           dtsets(idtset)%atvshift(1:dtsets(idtset)%natvshift,1:dtsets(idtset)%nsppol,1:mxnatpawu)
         else
           write(iout,format01170)'atvshift',&
&           dtsets(idtset)%atvshift(1:dtsets(idtset)%natvshift,1:dtsets(idtset)%nsppol,1:mxnatpawu)
         end if
       end if
     end do
   end if
 end if

 dprarr(1,:)=dtsets(:)%alpha
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'alpha','DPR')

!amu
 if(multi_ntypat==0)then
   do idtset=0,ndtset_alloc
     dprarr(1:ntypat,idtset)=dtsets(idtset)%amu(1:ntypat)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'amu','DPR')
 else
   do idtset=1,ndtset_alloc
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01160a)'amu',appen,dtsets(idtset)%amu(1:dtsets(idtset)%ntypat)
     else
       write(iout,format01160)'amu',dtsets(idtset)%amu(1:dtsets(idtset)%ntypat)
     end if
   end do
 end if

 intarr(1,:)=dtsets(:)%awtr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'awtr','INT')

 intarr(1,:)=dtsets(:)%bandpp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bandpp','INT')

 intarr(1,:)=dtsets(:)%berryopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'berryopt','INT')

 intarr(1,:)=dtsets(:)%bdberry(1)
 intarr(2,:)=dtsets(:)%bdberry(2)
 intarr(3,:)=dtsets(:)%bdberry(3)
 intarr(4,:)=dtsets(:)%bdberry(4)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,4,ndtset_alloc,'bdberry','INT')

 if(response==1) then
   intarr(1,:)=dtsets(:)%bdeigrf
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bdeigrf','INT')
 end if


!bdgw
 if(multi_nkptgw==0)then
!  Required if for pathscale  to avoid failure with -ff_bounds_check
   if (dtsets(1)%nkptgw > 0) then
     do idtset=0,ndtset_alloc
       intarr(1:2*dtsets(1)%nkptgw*dtsets(1)%nsppol,idtset)=&
&       reshape(dtsets(idtset)%bdgw(1:2,1:dtsets(1)%nkptgw,1:dtsets(1)%nsppol),(/2*dtsets(1)%nkptgw*dtsets(1)%nsppol/) )
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2*dtsets(1)%nkptgw*dtsets(1)%nsppol,ndtset_alloc,'bdgw','INT')
   end if
 else
   do idtset=1,ndtset_alloc
     if(dtsets(idtset)%nkptgw>0)then
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01155a)'bdgw',appen,dtsets(idtset)%bdgw(1:2,1:dtsets(idtset)%nkptgw,1:dtsets(idtset)%nsppol)
       else
         write(iout,format01155)'bdgw',dtsets(idtset)%bdgw(1:2,1:dtsets(idtset)%nkptgw,1:dtsets(idtset)%nsppol)
       end if
     end if
   end do
 end if

 dprarr(1,:)=dtsets(:)%bfield(1)
 dprarr(2,:)=dtsets(:)%bfield(2)
 dprarr(3,:)=dtsets(:)%bfield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'bfield','DPR')

 dprarr(1,:)=dtsets(:)%bmass
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'bmass','DPR')

 dprarr(1,:)=dtsets(:)%boxcenter(1)
 dprarr(2,:)=dtsets(:)%boxcenter(2)
 dprarr(3,:)=dtsets(:)%boxcenter(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'boxcenter','DPR')

 dprarr(1,:)=dtsets(:)%boxcutmin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'boxcutmin','DPR')

 if (usepaw==1) then

   dprarr(1,:)=dtsets(:)%bxctmindg
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'bxctmindg','DPR')

 end if

 dprarr(1,:)=dtsets(:)%charge
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'charge','DPR')

 intarr(1,:)=dtsets(:)%chkexit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'chkexit','INT')

 intarr(1,:)=dtsets(:)%chksymbreak
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'chksymbreak','INT')

 if(dtsets(1)%cpus>one)then
   cpus=dtsets(1)%cpus
   write(iout,'(1x,a16,1x,1p,t20,g10.2,t25,a)') 'cpus',cpus,'(seconds)'
   write(iout,'(1x,a16,1x,1p,t20,g10.2,t25,a)') 'cpum',cpus/60.0_dp,'(minutes)'
   write(iout,'(1x,a16,1x,1p,t20,g10.2,t25,a)') 'cpuh',cpus/3600.0_dp,'(hours)'
 end if

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%corecs(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'corecs','DPR')

 intarr(1,:)=dtsets(:)%delayperm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'delayperm','INT')

!densty
 if(multi_ntypat==0)then
   do idtset=0,ndtset_alloc
!    Only one component of densty is used until now
     dprarr(1:ntypat,idtset)=dtsets(idtset)%densty(1:ntypat,1)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,ndtset_alloc,'densty','DPR')
 else
   do idtset=1,ndtset_alloc
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01160a)'densty',appen,dtsets(idtset)%densty(1:dtsets(idtset)%ntypat,1)
     else
       write(iout,format01160)'densty',dtsets(idtset)%densty(1:dtsets(idtset)%ntypat,1)
     end if
   end do
 end if

 dprarr(1,:)=dtsets(:)%diecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diecut','ENE')

 dprarr(1,:)=dtsets(:)%diegap
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diegap','ENE')

 dprarr(1,:)=dtsets(:)%dielam
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dielam','DPR')

 dprarr(1,:)=dtsets(:)%dielng
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dielng','LEN')

 dprarr(1,:)=dtsets(:)%diemac
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diemac','DPR')

 dprarr(1,:)=dtsets(:)%diemix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diemix','DPR')

 if (any(dtsets(1:ndtset_alloc)%diemixmag/=dtsets(1:ndtset_alloc)%diemix)) then
   dprarr(1,:)=dtsets(:)%diemixmag
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'diemixmag','DPR')
 end if

 dprarr(1,:)=dtsets(:)%dilatmx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dilatmx','DPR')

 dprarr(1,:)=dtsets(:)%dosdeltae
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dosdeltae','ENE')

 dprarr(1,:)=dtsets(:)%dtion
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'dtion','DPR')

!dynimage
 if(multi_nimage==0)then
   if(nimage/=0)then
     do idtset=0,ndtset_alloc
       intarr(1:nimage,idtset)=dtsets(idtset)%dynimage(1:nimage)
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nimage,&
&     ndtset_alloc,'dynimage','INT')
   end if
 else
   do idtset=1,ndtset_alloc
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01155a)'dynimage',appen,dtsets(idtset)%dynimage(1:dtsets(idtset)%nimage)
     else
       write(iout,format01155)'dynimage',appen,dtsets(idtset)%dynimage(1:dtsets(idtset)%nimage)
     end if
   end do
 end if

 dprarr(1,:)=dtsets(:)%ecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecut','ENE')

 dprarr(1,:)=dtsets(:)%ecuteps
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecuteps','ENE')

 dprarr(1,:)=dtsets(:)%ecutsigx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecutsigx','ENE')

 dprarr(1,:)=dtsets(:)%ecutwfn
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecutwfn','ENE')

 dprarr(1,:)=dtsets(:)%ecutsm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'ecutsm','ENE')

 dprarr(1,:)=dtsets(:)%effmass
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'effmass','DPR')

 dprarr(1,:)=dtsets(:)%efield(1)
 dprarr(2,:)=dtsets(:)%efield(2)
 dprarr(3,:)=dtsets(:)%efield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'efield','DPR')

 intarr(1,:)=dtsets(:)%enunit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'enunit','INT')

 dprarr(1,:)=dtsets(:)%eshift
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'eshift','ENE')

 dprarr(1,:)=dtsets(:)%esmear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'esmear','ENE')

 dprarr(1,:)=dtsets(:)%exchmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'exchmix','DPR')

!etotal
 if(choice==2)then
   do idtset=1,ndtset_alloc
     do iimage=1,dtsets(idtset)%nimage
       if(dtsets(idtset)%dynimage(iimage)==1)then
         keywd='etotal'//trim(strimg(iimage))
         iscf=dtsets(idtset)%iscf
         if(iscf>0 .or. iscf==-3)then
           if(ndtset>0)then
             jdtset=jdtset_(idtset)
             if(jdtset<10)write(appen,'(i1)')jdtset
             if(jdtset>=10)write(appen,'(i2)')jdtset
             write(iout,format01160a)trim(keywd),appen,results_out(idtset)%etotal(iimage)
           else
             write(iout,format01160)trim(keywd),results_out(idtset)%etotal(iimage)
           end if
         end if
       end if
     end do
   end do
 end if

 intarr(1,:)=dtsets(:)%exchn2n3d
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'exchn2n3d','INT')

 intarr(1,:)=dtsets(:)%pawfatbnd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'pawfatbnd','INT')

 intarr(1,:)=dtsets(:)%ngfft(7)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fftalg','INT')

 intarr(1,:)=dtsets(:)%ngfft(8)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fftcache','INT')

 intarr(1,:)=dtsets(:)%fftgw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fftgw','INT')

 intarr(1,:)=dtsets(:)%fft_opt_lob
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'fft_opt_lob','INT')

!force
 if(choice==2)then
   do idtset=1,ndtset_alloc
     do iimage=1,dtsets(idtset)%nimage
       if(dtsets(idtset)%dynimage(iimage)==1)then
         keywd='fcart'//trim(strimg(iimage))
         iscf=dtsets(idtset)%iscf
         if(iscf>0)then
           if(ndtset>0)then
             jdtset=jdtset_(idtset)
             if(jdtset<10)write(appen,'(i1)')jdtset
             if(jdtset>=10)write(appen,'(i2)')jdtset
             write(iout,format01160a)trim(keywd),appen,&
&             results_out(idtset)%fcart(:,1:dtsets(idtset)%natom,iimage)
           else
             write(iout,format01160)trim(keywd),&
&             results_out(idtset)%fcart(:,1:dtsets(idtset)%natom,iimage)
           end if
         end if
       end if
     end do
   end do
 end if

 dprarr(1,:)=dtsets(:)%fixmom
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'fixmom','DPR')

 dprarr(1,:)=dtsets(:)%fxcartfactor
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'fxcartfactor','DPR')

 dprarr(1,:)=dtsets(:)%rhoqpmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'rhoqpmix','DPR')

 dprarr(1,:)=dtsets(:)%freqremax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqremax','ENE')

 dprarr(1,:)=dtsets(:)%freqspmax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqspmax','ENE')

 dprarr(1,:)=dtsets(:)%freqsusin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqsusin','DPR')

 dprarr(1,:)=dtsets(:)%freqsuslo
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'freqsuslo','DPR')

 dprarr(1,:)=dtsets(:)%friction
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'friction','DPR')

 intarr(1,:)=dtsets(:)%frzfermi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'frzfermi','INT')

 intarr(1,:)=dtsets(:)%getbseig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getbseig','INT')

 intarr(1,:)=dtsets(:)%getcell
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getcell','INT')

 intarr(1,:)=dtsets(:)%getddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getddk','INT')

 intarr(1,:)=dtsets(:)%getden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getden','INT')

 intarr(1,:)=dtsets(:)%getqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getqps','INT')

 intarr(1,:)=dtsets(:)%getscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getscr','INT')

 intarr(1,:)=dtsets(:)%getsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getsuscep','INT')

 intarr(1,:)=dtsets(:)%getkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getkss','INT')

 intarr(1,:)=dtsets(:)%getocc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getocc','INT')

 intarr(1,:)=dtsets(:)%getvel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getvel','INT')

 intarr(1,:)=dtsets(:)%getwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getwfk','INT')

 intarr(1,:)=dtsets(:)%getwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getwfq','INT')

 intarr(1,:)=dtsets(:)%getxcart
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getxcart','INT')

 intarr(1,:)=dtsets(:)%getxred
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'getxred','INT')

 intarr(1,:)=dtsets(:)%get1den
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'get1den','INT')

 intarr(1,:)=dtsets(:)%get1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'get1wf','INT')

 intarr(1,:)=dtsets(:)%gwcalctyp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwcalctyp','INT')

 intarr(1,:)=dtsets(:)%gwcomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwcomp','INT')

 dprarr(1,:)=dtsets(:)%gwencomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwencomp','ENE')

 intarr(1,:)=dtsets(:)%gwgamma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwgamma','INT')

 intarr(1,:)=dtsets(:)%gw_nqlwl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gw_nqlwl','INT')

 intarr(1,:)=dtsets(:)%gw_sctype
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gw_sctype','INT')

 intarr(1,:)=dtsets(:)%gw_EET
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gw_EET','INT')

 intarr(1,:)=dtsets(:)%gw_EET_nband
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gw_EET_nband','INT')

 intarr(1,:)=dtsets(:)%gw_sigxcore
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gw_sigxcore','INT')

 dprarr(1,:)=dtsets(:)%gw_toldfeig
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'gw_toldfeig','ENE')

 intarr(1,:)=dtsets(:)%gwmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwmem','INT')

 intarr(1,:)=dtsets(:)%gw_nstep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gw_nstep','INT')

 intarr(1,:)=dtsets(:)%gwpara
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwpara','INT')

 intarr(1,:)=dtsets(:)%gwrpacorr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'gwrpacorr','INT')

!gw_qlwl
 if(multi_gw_nqlwl==0)then
!  Required if for pathscale to avoid failure with -ff_bounds_check
   if (dtsets(1)%gw_nqlwl > 0) then
     do idtset=0,ndtset_alloc
       dprarr(1:3*dtsets(1)%gw_nqlwl,idtset) = &
&       reshape(dtsets(idtset)%gw_qlwl(1:3,1:dtsets(1)%gw_nqlwl),(/3*dtsets(1)%gw_nqlwl/) )
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*dtsets(1)%gw_nqlwl,ndtset_alloc,'gw_qlwl','DPR')
   end if
 else
   do idtset=1,ndtset_alloc
     if(dtsets(idtset)%gw_nqlwl>0)then
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01150a)'gw_qlwl',appen,dtsets(idtset)%gw_qlwl(1:3,1:dtsets(idtset)%gw_nqlwl)
       else
         write(iout,format01150)'gw_qlwl',dtsets(idtset)%gw_qlwl(1:3,1:dtsets(idtset)%gw_nqlwl)
       end if
     end if
   end do
 end if

!iatfix
 if(multi_natfix==0)then
   if(natfix/=0)then
     intarr(1:natfix,0:ndtset_alloc)=iatfixio_(1:natfix,0:ndtset_alloc)
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfix,&
&     ndtset_alloc,'iatfix','INT')
   end if
 else
   do idtset=1,ndtset_alloc
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01155a)'iatfix',appen,iatfixio_(1:natfix_(idtset),idtset)
     else
       write(iout,format01155)'iatfix',iatfixio_(1:natfix_(idtset),idtset)
     end if
   end do
 end if

!iatfixx
 if(multi_natfixx==0)then
   if(natfixx/=0)then
     intarr(1:natfixx,0:ndtset_alloc)=iatfixx_(1:natfixx,0:ndtset_alloc)
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfixx,&
&     ndtset_alloc,'iatfixx','INT')
   end if
 else
   do idtset=1,ndtset_alloc
     if(natfixx_(idtset)>0)then
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01155a)'iatfixx',appen,iatfixx_(1:natfixx_(idtset),idtset)
       else
         write(iout,format01155)'iatfixx',iatfixx_(1:natfixx_(idtset),idtset)
       end if
     end if
   end do
 end if

!iatfixy
 if(multi_natfixy==0)then
   if(natfixy/=0)then
     intarr(1:natfixy,0:ndtset_alloc)=iatfixy_(1:natfixy,0:ndtset_alloc)
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfixy,&
&     ndtset_alloc,'iatfixy','INT')
   end if
 else
   do idtset=1,ndtset_alloc
     if(natfixy_(idtset)>0)then
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01155a)'iatfixy',appen,iatfixy_(1:natfixy_(idtset),idtset)
       else
         write(iout,format01155)'iatfixy',iatfixy_(1:natfixy_(idtset),idtset)
       end if
     end if
   end do
 end if

!iatfixz
 if(multi_natfixz==0)then
   if(natfixz/=0)then
     intarr(1:natfixz,0:ndtset_alloc)=iatfixz_(1:natfixz,0:ndtset_alloc)
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,natfixz,&
&     ndtset_alloc,'iatfixz','INT')
   end if
 else
   do idtset=1,ndtset_alloc
     if(natfixz_(idtset)>0)then
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01155a)'iatfixz',appen,iatfixz_(1:natfixz_(idtset),idtset)
       else
         write(iout,format01155)'iatfixz',iatfixz_(1:natfixz_(idtset),idtset)
       end if
     end if
   end do
 end if

!iatsph   need to be printed only if there is some occurence of prtdos==3 or pawfatbnd
 do idtset=1,ndtset_alloc
   if(dtsets(idtset)%prtdos==3.or.dtsets(idtset)%pawfatbnd>0)then
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01155a)'iatsph',appen,dtsets(idtset)%iatsph(1:dtsets(idtset)%natsph)
     else
       write(iout,format01155)'iatsph',dtsets(idtset)%iatsph(1:dtsets(idtset)%natsph)
     end if
   end if
 end do

 intarr(1,:)=dtsets(:)%iboxcut
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iboxcut','INT')

 intarr(1,:)=dtsets(:)%icutcoul
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'icutcoul','INT')

 intarr(1,:)=dtsets(:)%icoulomb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'icoulomb','INT')

 intarr(1,:)=dtsets(:)%idyson
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'idyson','INT')

 intarr(1,:)=dtsets(:)%ieig2rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ieig2rf','INT')

 intarr(1,:)=dtsets(:)%ikhxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ikhxc','INT')

 intarr(1,:)=dtsets(:)%imgmov
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'imgmov','INT')

 intarr(1,:)=dtsets(:)%inclvkb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'inclvkb','INT')

 intarr(1,:)=dtsets(:)%intxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'intxc','INT')

 intarr(1,:)=dtsets(:)%intexact
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'intexact','INT')

 intarr(1,:)=dtsets(:)%ionmov
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ionmov','INT')

 intarr(1,:)=dtsets(:)%iextrapwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iextrapwf','INT')

 intarr(1,:)=dtsets(:)%iprcch
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprcch','INT')

 intarr(1,:)=dtsets(:)%iprcel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprcel','INT')

 intarr(1,:)=dtsets(:)%iprctfvw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprctfvw','INT')

 intarr(1,:)=dtsets(:)%iprcfc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iprcfc','INT')

 intarr(1,:)=dtsets(:)%irdbseig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdbseig','INT')

 intarr(1,:)=dtsets(:)%irdddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdddk','INT')

 intarr(1,:)=dtsets(:)%irdkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdkss','INT')

 intarr(1,:)=dtsets(:)%irdqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdqps','INT')

 intarr(1,:)=dtsets(:)%irdscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdscr','INT')

 intarr(1,:)=dtsets(:)%irdsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdsuscep','INT')

 intarr(1,:)=dtsets(:)%irdden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdden','INT')

 intarr(1,:)=dtsets(:)%irdwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdwfk','INT')

 intarr(1,:)=dtsets(:)%irdwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'irdwfq','INT')

 intarr(1,:)=dtsets(:)%ird1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ird1wf','INT')

 intarr(1,:)=dtsets(:)%iscf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'iscf','INT')

 intarr(1,:)=dtsets(:)%isecur
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'isecur','INT')

 istatr_defo=49
 if(istatr/=istatr_defo)write(iout,format01110) 'istatr',istatr

 istatshft_defo=1
 if(istatshft/=istatshft_defo)write(iout,format01110) 'istatshft',istatshft

!istwfk (must first restore the default istwf=0 for non-allowed k points)
 allocate(istwfk_2(mxnkpt,0:ndtset_alloc))
 do idtset=1,ndtset_alloc
   nqpt=dtsets(idtset)%nqpt
   do ikpt=1,dtsets(idtset)%nkpt
     allowed=1
     do ii=1,3
!      kpoint=dtsets(idtset)%kptns(ii)
       kpoint=dtsets(idtset)%kpt(ii,ikpt)/dtsets(idtset)%kptnrm
       if(nqpt/=0 .and. response_(idtset)==0)&
&       kpoint=kpoint+dtsets(idtset)%qptn(ii)
       if(abs(kpoint)>1.d-10 .and. abs(kpoint-0.5_dp)>1.e-10_dp )&
&       allowed=0
     end do
     if(allowed==0)then
       istwfk_2(ikpt,idtset)=0
     else
       istwfk_2(ikpt,idtset)=dtsets(idtset)%istwfk(ikpt)
     end if
   end do
 end do

!DEBUG
!write(6,*)' outvar1 '
!write(6,*)istwfk_2(:,:)
!ENDDEBUG
!
 if(multi_nkpt==0)then
!  Might restrict the number of k points to be printed
   tnkpt=0
   nkpt_eff=dtsets(1)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
     nkpt_eff=nkpt_max
     tnkpt=1
   end if
   intarr(1:nkpt_eff,0)=0
   intarr(1:nkpt_eff,1:ndtset_alloc)=istwfk_2(1:nkpt_eff,1:ndtset_alloc)
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nkpt_eff,&
&   ndtset_alloc,'istwfk','INT')
   if(tnkpt==1 .and. sum(istwfk_2(1:nkpt_eff,1:ndtset_alloc))/=0 ) &
&   write(iout,'(23x,a)' ) &
&   'outvar1 : prtvol=0, do not print more k-points.'

 else
   do idtset=1,ndtset_alloc
     tnkpt=0
     nkpt_eff=dtsets(idtset)%nkpt
     if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
       nkpt_eff=nkpt_max
       tnkpt=1
     end if
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       if(sum(istwfk_2(1:nkpt_eff,idtset))/=0)then
         write(iout,format01155a)'istwfk',appen,istwfk_2(1:nkpt_eff,idtset)
       end if
     else
       if(sum(istwfk_2(1:nkpt_eff,idtset))/=0)then
         write(iout,format01155)'istwfk',istwfk_2(1:nkpt_eff,idtset)
       end if
     end if
     if(tnkpt==1 .and. sum(istwfk_2(1:nkpt_eff,idtset))/=0 ) &
&     write(iout,'(23x,a)' ) &
&     'outvar1 : prtvol=0, do not print more k-points.'
   end do
 end if
 deallocate(istwfk_2)

 intarr(1,:)=dtsets(:)%ixc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ixc','INT')

 intarr(1,:)=dtsets(:)%ixcpositron
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ixcpositron','INT')

 intarr(1,:)=dtsets(:)%vdw_xc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'vdw_xc','INT')

 if (ndtset > 0) write(iout,format01155) 'jdtset',jdtset_(1:ndtset)

 intarr(1,:)=dtsets(:)%jellslab
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'jellslab','INT')

!kberry
 if(multi_nberry==0)then
   if(dtsets(1)%nberry/=0)then
     do idtset=0,ndtset_alloc
       intarr(1:3*dtsets(1)%nberry,idtset)=&
&       reshape( dtsets(idtset)%kberry(1:3,1:dtsets(1)%nberry), (/3*dtsets(1)%nberry/) )
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*dtsets(1)%nberry,&
&     ndtset_alloc,'kberry','INT')
   end if
 else
   do idtset=1,ndtset_alloc
     if(dtsets(idtset)%nberry>0)then
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01155a)&
&         'kberry',appen,dtsets(idtset)%kberry(1:3,1:dtsets(idtset)%nberry)
       else
         write(iout,format01155)&
&         'kberry',dtsets(idtset)%kberry(1:3,1:dtsets(idtset)%nberry)
       end if
     end if
   end do
 end if


!kpt
 if(multi_nkpt==0)then
!  Might restrict the number of k points to be printed
   tnkpt=0
   nkpt_eff=dtsets(1)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
     nkpt_eff=nkpt_max
     tnkpt=1
   end if
   do idtset=0,ndtset_alloc
     dprarr(1:3*nkpt_eff,idtset)=&
&     reshape(dtsets(idtset)%kpt(1:3,1:nkpt_eff),(/3*nkpt_eff/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*nkpt_eff,&
&   ndtset_alloc,'kpt','DPR')
   if(tnkpt==1) write(iout,'(23x,a)' ) &
&   '       outvar1 : prtvol=0, do not print more k-points.'
 else
   do idtset=1,ndtset_alloc
     tnkpt=0
     nkpt_eff=dtsets(idtset)%nkpt
     if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
       nkpt_eff=nkpt_max
       tnkpt=1
     end if
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01150a)'kpt',appen,dtsets(idtset)%kpt(1:3,1:nkpt_eff)
     else
       write(iout,format01150)'kpt',dtsets(idtset)%kpt(1:3,1:nkpt_eff)
     end if
     if(tnkpt==1) write(iout,'(23x,a)' ) &
&     'outvar1 : prtvol=0, do not print more k-points.'
   end do
 end if

!kptgw
 if(multi_nkptgw==0)then
!  Required if for pathscale to avoid failure with -ff_bounds_check
   if (dtsets(1)%nkptgw > 0) then
     do idtset=0,ndtset_alloc
       dprarr(1:3*dtsets(1)%nkptgw,idtset)=&
&       reshape(dtsets(idtset)%kptgw(1:3,1:dtsets(1)%nkptgw),(/3*dtsets(1)%nkptgw/) )
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3*dtsets(1)%nkptgw,&
&     ndtset_alloc,'kptgw','DPR')
   end if
 else
   do idtset=1,ndtset_alloc
     if(dtsets(idtset)%nkptgw>0)then
       if(ndtset>0)then
         jdtset=jdtset_(idtset)
         if(jdtset<10)write(appen,'(i1)')jdtset
         if(jdtset>=10)write(appen,'(i2)')jdtset
         write(iout,format01150a)'kptgw',appen,dtsets(idtset)%kptgw(1:3,1:dtsets(idtset)%nkptgw)
       else
         write(iout,format01150)'kptgw',dtsets(idtset)%kptgw(1:3,1:dtsets(idtset)%nkptgw)
       end if
     end if
   end do
 end if

 dprarr(1,:)=dtsets(:)%kptnrm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'kptnrm','DPR')

 dprarr(1,:)=dtsets(:)%kptrlen
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'kptrlen','DPR')

 intarr(1,:)=dtsets(:)%kptopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'kptopt','INT')

!kptrlatt
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   ndtset_kptopt=0
   intarr(1:9,0)=reshape( dtsets(0)%kptrlatt(:,:) , (/9/) )
   allocate(jdtset_kptopt(0:ndtset_alloc))
!  Define the set of datasets for which kptopt>0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>0)then
       ndtset_kptopt=ndtset_kptopt+1
       jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
       intarr(1:9,ndtset_kptopt)=reshape( dtsets(idtset)%kptrlatt(:,:) , (/9/) )
     end if
   end do
   if(ndtset_kptopt>0)then
     call prttagm(dprarr,intarr,iout,jdtset_kptopt,3,marr,9,&
&     ndtset_kptopt,'kptrlatt','INT')
   end if
   deallocate(jdtset_kptopt)
 end if

 intarr(1,:)=dtsets(:)%kssform
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'kssform','INT')

 intarr(1,:)=dtsets(:)%ldgapp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'ldgapp','INT')

 intarr(1,:)=dtsets(:)%localrdwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'localrdwf','INT')

 dprarr(1,:)=dtsets(:)%mdftemp
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'mdftemp','DPR')

 dprarr(1,:)=dtsets(:)%mditemp
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'mditemp','DPR')

 dprarr(1,:)=dtsets(:)%mdwall
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'mdwall','LEN')

 intarr(1,:)=dtsets(:)%mffmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'mffmem','INT')

!mixalch
 if(multi_ntypalch==0)then
   do idtset=0,ndtset_alloc
     dprarr(1:dtsets(idtset)%npspalch*ntypalch,idtset)=&
&     reshape(dtsets(idtset)%mixalch(1:dtsets(idtset)%npspalch,1:ntypalch),&
&     (/dtsets(idtset)%npspalch*ntypalch/))
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,dtsets(1)%npspalch*ntypalch,ndtset_alloc,'mixalch','DPR')
 else
   do idtset=1,ndtset_alloc
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01160a)'mixalch',appen,&
&       dtsets(idtset)%mixalch(1:dtsets(idtset)%npspalch,1:dtsets(idtset)%ntypalch)
     else
       write(iout,format01160)'mixalch',&
&       dtsets(idtset)%mixalch(1:dtsets(idtset)%npspalch,1:dtsets(idtset)%ntypalch)
     end if
   end do
 end if

!DEBUG
!write(6,*)' outvar1 : after mixalch '
!ENDDEBUG
!

 intarr(1,:)=dtsets(:)%maxnsym
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,1,ndtset_alloc,'maxnsym','INT')

 intarr(1,:)=dtsets(:)%mkmem
 call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,ndtset_alloc,'mkmem','INT')

 if(response==1)then
   intarr(1,:)=dtsets(:)%mkqmem
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,ndtset_alloc,'mkqmem','INT')
 end if

 if(response==1)then
   intarr(1,:)=dtsets(:)%mk1mem
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,ndtset_alloc,'mk1mem','INT')
 end if

 intarr(1,:)=dtsets(:)%mqgrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'mqgrid','INT')
 if (usepaw==1) then
   intarr(1,:)=dtsets(:)%mqgriddg
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'mqgriddg','INT')
 end if

 intarr(1,:)=natfix_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfix','INT')

 intarr(1,:)=natfixx_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfixx','INT')

 intarr(1,:)=natfixy_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfixy','INT')

 intarr(1,:)=natfixz_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natfixz','INT')

 intarr(1,:)=dtsets(0:ndtset_alloc)%natom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'natom','INT')

!natsph   need to be printed only if there is some occurence of prtdos==3 or
!pawfatbnd>0
 do idtset=1,ndtset_alloc
   if(dtsets(idtset)%prtdos==3.or.dtsets(idtset)%pawfatbnd>0)then
     if(ndtset>0)then
       jdtset=jdtset_(idtset)
       if(jdtset<10)write(appen,'(i1)')jdtset
       if(jdtset>=10)write(appen,'(i2)')jdtset
       write(iout,format01155a)'natsph',appen,dtsets(idtset)%natsph
     else
       write(iout,format01155)'natsph',dtsets(idtset)%natsph
     end if
   end if
 end do

!BEGIN VARIABLES FOR @Bethe-Salpeter
 intarr(1,:)=dtsets(:)%bs_algorithm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bs_algorithm','INT')

 intarr(1,:)=dtsets(:)%bs_haydock_niter
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bs_haydock_niter','INT')

 intarr(1,:)=dtsets(:)%bs_exchange_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bs_exchange_term','INT')

 intarr(1,:)=dtsets(:)%bs_coulomb_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bs_coulomb_term','INT')

 intarr(1,:)=dtsets(:)%bs_calctype
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bs_calctype','INT')

 intarr(1,:)=dtsets(:)%bs_coupling
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,ndtset_alloc,'bs_coupling','INT')

 dprarr(1,:)=dtsets(:)%bs_haydock_tol
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,ndtset_alloc,'bs_haydock_tol','DPR')

 do idtset=0,ndtset_alloc
   intarr(1:2,idtset)=dtsets(idtset)%bs_eh_basis_set(1:2)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,ndtset_alloc,'bs_eh_basis_set','INT')

 do idtset=0,ndtset_alloc
   dprarr(1:3,idtset)=dtsets(idtset)%bs_eh_cutoff(1:3)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'bs_eh_cutoff','ENE')

 do idtset=0,ndtset_alloc
   dprarr(1:3,idtset)=dtsets(idtset)%bs_freq_mesh(1:3)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,ndtset_alloc,'bs_freq_mesh','ENE')
!END VARIABLES FOR @Bethe-Salpeter.

 deallocate(dprarr,intarr)
 deallocate(natfix_,iatfixio_)
 deallocate(natfixx_,iatfixx_)
 deallocate(natfixy_,iatfixy_)
 deallocate(natfixz_,iatfixz_)

!DEBUG
!write(6,*)' outvar1 : end of subroutine '
!if(.true.)stop
!ENDDEBUG
!
end subroutine outvar1
!!***
