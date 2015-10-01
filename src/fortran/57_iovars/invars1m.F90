!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars1m
!! NAME
!! invars1m
!!
!! FUNCTION
!! Initialisation phase : prepare the main input subroutine call by
!! reading all the NO MULTI variables, as well as the dimensions
!! needed for allocating the input arrays in abinit.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG, MKV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  mpi_enreg=informations about MPI parallelization
!!  msym=default maximal number of symmetries
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnimage=maximal value of input nimage for all the datasets
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!               one data set.
!!  string*(*)=string of characters containing all input variables and data
!!  zion_max=maximal valence charge over all psps
!!
!! OUTPUT
!!  bravais_(11,0:ndtset_alloc)=characteristics of Bravais lattice (see symlatt.F90)
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  mband_upper_(0:ndtset_alloc)=list of mband_upper values
!!  mxlpawu=maximal value of input lpawu for all the datasets
!!  mxmband_upper=maximal value of input nband for all the datasets
!!  mxnatpawu=maximal value of number of atoms on which +U is applied for all the datasets
!!  mxnatsph=maximal value of input natsph for all the datasets
!!  mxnatvshift=maximal value of input natvshift for all the datasets
!!  mxnconeq=maximal value of input nconeq for all the datasets
!!  mxnkptgw=maximal value of input nkptgw for all the datasets
!!  mxgw_nqlwl=maximal value of input gw_nqlwl for all the datasets
!!  mxnnos=maximal value of input nnos for all the datasets
!!  mxnqptdm=maximal value of input nqptdm for all the datasets
!!  mxnkpt=maximal value of input nkpt for all the datasets
!!  mxnspinor=maximal value of input nspinor for all the datasets
!!  mxnsppol=maximal value of input nsppol for all the datasets
!!  mxnsym=maximum number of symmetries
!!
!! SIDE EFFECTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here (see invars1.f for more details on the
!!   initialized records)
!!
!! TODO
!!  MG: What about using modules to store the maximum dimensions as global variables?
!!
!! PARENTS
!!      m_ab6_invars_f90
!!
!! CHILDREN
!!      indefo1,invars1,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine invars1m(bravais_,dmatpuflag,dtsets,iout,lenstr,mband_upper_,mpi_enreg,&
& msym,mxgw_nqlwl,mxlpawu,mxmband_upper,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxnconeq,&
& mxnkpt,mxnkptgw,mxnimage,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,&
& ndtset,ndtset_alloc,string,zion_max)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_57_iovars, except_this_one => invars1m
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,msym,mxnatom,mxnimage,ndtset,ndtset_alloc
 integer,intent(out) :: dmatpuflag,mxgw_nqlwl,mxlpawu,mxmband_upper,mxnatpawu,mxnatsph
 integer,intent(out) :: mxnatvshift,mxnconeq,mxnkpt,mxnkptgw,mxnnos
 integer,intent(out) :: mxnqptdm,mxnspinor,mxnsppol,mxnsym
 real(dp),intent(in) :: zion_max
 character(len=*),intent(inout) :: string
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(out) :: bravais_(11,0:ndtset_alloc)
 integer,intent(out) :: mband_upper_(0:ndtset_alloc)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer :: i1,idtset,ierr,ii,jdtset,lpawu,mband_upper
!integer :: iatom,itypat
 real(dp):: weight
 character(len=500) :: message
!arrays
 integer :: bravais(11)
 integer,allocatable :: symafm_(:,:)
 integer,allocatable :: symrel_(:,:,:,:)
 real(dp),allocatable :: tnons_(:,:,:)
 integer,allocatable :: symafm(:),symrel(:,:,:)
 real(dp),allocatable :: tnons(:,:)

!******************************************************************

!Here, allocation of the three arrays that depend on msym.
!DEBUG
!write(*,*)'invars1m start, mxnimage=',mxnimage
!ENDDEBUG


 allocate(symrel_(3,3,msym,0:ndtset_alloc))
 allocate(symafm_(msym,0:ndtset_alloc))
 allocate(tnons_(3,msym,0:ndtset_alloc))
 allocate(symafm(msym),symrel(3,3,msym),tnons(3,msym))

!Set up default values (note that the default acell, amu
!mkmem, mkmem1,mkqmem, and nkpt must be overcome

!DEBUG
!write(*,*)'invars1m: call indefo1'
!ENDDEBUG

 do idtset=0,ndtset_alloc
   call indefo1(dtsets(idtset))
 end do

!natom and nimage are already initialized in invars0
 dtsets(0)%natom=-1
 dtsets(0)%nimage=1

 bravais_(:,0)=0
 symafm_(:,0)=1
 symrel_(:,:,:,0)=0
 symrel_(1,1,:,0)=1 ; symrel_(2,2,:,0)=1 ; symrel_(3,3,:,0)=1
 tnons_(:,:,0)=0.0_dp
 ierr=0

!Loop on datasets
 do idtset=1,ndtset_alloc
   if(mpi_enreg%paral_compil==0.and.dtsets(idtset)%paral_kgb==1)then
     dtsets(idtset)%paral_kgb=0
     write(message, '(8a)' )ch10,&
&     ' abinit : WARNING -',ch10,&
&     '  When ABINIT is compiled without MPI flag,',ch10,&
&     '  setting paral_kgb/=0 is useless. paral_kgb has been reset to 0.',ch10,&
&     '  Action : modify compilation option or paral_kgb in the input file.'
     call wrtout(std_out,  message,'COLL')
   end if
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0

   write(message, '(a,a,i6)') ch10,&
&   ' invars1m : enter jdtset=',jdtset
   call wrtout(std_out,message,'COLL')

!  Input default values

   bravais(:)=bravais_(:,0)
   symafm(:)=symafm_(:,0)
   symrel(:,:,:)=symrel_(:,:,:,0)
   tnons(:,:)=tnons_(:,:,0)
   if(dtsets(idtset)%paral_kgb/=0) then
     allocate (mpi_enreg%keywp(5,250),mpi_enreg%trialproc(2))
     mpi_enreg%keywp=0;mpi_enreg%trialproc=0
   end if

   call invars1(bravais,dtsets(idtset),iout,jdtset,lenstr,&
&   mband_upper,mpi_enreg,msym,string,symafm,symrel,tnons,zion_max)

   bravais_(:,idtset)=bravais(:)
   mband_upper_ (idtset)=mband_upper
   symafm_(:,idtset)=symafm(:)
   symrel_(:,:,:,idtset)=symrel(:,:,:)
   tnons_(:,:,idtset)=tnons(:,:)

!  DEBUG
!  write(6,*)' invars1m : idtset, dtsets(idtset)%tolsym=',idtset, dtsets(idtset)%tolsym
!  write(6,*)' dtsets(idtset)%paral_kgb=',dtsets(idtset)%paral_kgb
!  ENDDEBUG

   if(dtsets(idtset)%paral_kgb/=0) then

     if(mpi_enreg%me==0)then

       write(message, '(a,a,a,i4,a,i4,a,a,a)' ) ch10,&
&       ' WARNING in invars1m',&
&       ' For dataset=', idtset,&
&       '  a possible choice for less than ', mpi_enreg%trialproc(1) ,' processors is:',ch10,&
&       '  nproc    npkpt    npband     npfft    bandpp    weight'
       call wrtout(std_out,  message,'COLL')
       do ii=1,250

         if(mpi_enreg%keywp(1,ii)/=0) then
           weight=mpi_enreg%keywp(3,ii)/mpi_enreg%keywp(4,ii)/4.d0
           if ((mpi_enreg%keywp(3,ii)<=50).and.(mpi_enreg%keywp(4,ii)==1)) weight=1
           write(message, '(a,i4,a,i4,a,i4,a,i4,a,i4,a,f8.2)' ) &
&           '  ',mpi_enreg%keywp(1,ii), '    ',mpi_enreg%keywp(2,ii), '     ',mpi_enreg%keywp(3,ii), '      ',&
&           mpi_enreg%keywp(4,ii), '       ',mpi_enreg%keywp(5,ii),'       ',weight
           call wrtout(std_out,  message,'COLL')
         end if

       end do

       if(mpi_enreg%trialproc(2)==1) then
         write(message,'(a,a,a,a)' )ch10,&
&         ' invars1m :  launch a parallel version of ABINIT with a number of processor among the above list,',&
&         ' and the associated input variables npkpt, npband, npfft and bandpp. ',&
&         ' The optimal weight is close to 1.'
         call wrtout(std_out,  message,'COLL')
       end if
     end if

     ierr=ierr+mpi_enreg%trialproc(2)
     deallocate(mpi_enreg%keywp,mpi_enreg%trialproc)
   end if

 end do

 if(ierr/=0)then
   call leave_new('COLL')
 end if

 mxmband_upper =maxval(mband_upper_ (1:ndtset_alloc))

 dmatpuflag=0;mxnatpawu=0;mxlpawu=0
 mxnatsph=dtsets(1)%natsph
 mxnatvshift=dtsets(1)%natvshift
 mxnconeq=dtsets(1)%nconeq
 mxgw_nqlwl = dtsets(1)%gw_nqlwl
 mxnkptgw=dtsets(1)%nkptgw
 mxnkpt  =dtsets(1)%nkpt
 mxnnos  =dtsets(1)%nnos
 mxnqptdm=dtsets(1)%nqptdm
 mxnspinor=dtsets(1)%nspinor
 mxnsppol=dtsets(1)%nsppol

!Get MAX dimension over datasets
 do ii=1,ndtset_alloc
   mxnatsph=max(dtsets(ii)%natsph,mxnatsph)
   mxnconeq=max(dtsets(ii)%nconeq,mxnconeq)
   mxgw_nqlwl = max(dtsets(ii)%gw_nqlwl,mxgw_nqlwl)
   mxnkptgw=max(dtsets(ii)%nkptgw,mxnkptgw)
   mxnkpt  =max(dtsets(ii)%nkpt,mxnkpt)
   mxnnos  =max(dtsets(ii)%nnos,mxnnos)
   mxnqptdm=max(dtsets(ii)%nqptdm,mxnqptdm)
   mxnspinor=max(dtsets(ii)%nspinor,mxnspinor)
   mxnsppol=max(dtsets(ii)%nsppol,mxnsppol)
   if (dtsets(ii)%usepawu>0) then
     if (dtsets(ii)%usedmatpu/=0) dmatpuflag=1
     lpawu=maxval(dtsets(ii)%lpawu(:))
     mxlpawu=max(lpawu,mxlpawu)
     dtsets(ii)%natpawu=count(dtsets(ii)%lpawu(dtsets(ii)%typat((/(i1,i1=1,dtsets(ii)%natom)/)))/=-1)
     mxnatpawu=max(dtsets(ii)%natpawu,mxnatpawu)
     if (dtsets(ii)%macro_uj/=0) dtsets(ii)%natvshift=lpawu*2+1
   end if
   mxnatvshift=max(dtsets(ii)%natvshift,mxnatvshift)
 end do

 deallocate(symafm,symrel,tnons)

!mxnsym=maxval(dtsets(1:ndtset_alloc)%nsym) ! This might not work properly with HP compiler
 mxnsym=dtsets(1)%nsym
 do idtset=1,ndtset_alloc
   mxnsym=max(dtsets(idtset)%nsym,mxnsym)
 end do

 do idtset=0,ndtset_alloc
   allocate(dtsets(idtset)%atvshift(mxnatvshift,mxnsppol,mxnatpawu))
   allocate(dtsets(idtset)%bdgw(2,mxnkptgw,mxnsppol))
   allocate(dtsets(idtset)%dmatpawu(2*mxlpawu+1,2*mxlpawu+1,max(mxnsppol,mxnspinor),mxnatpawu*dmatpuflag))
   allocate(dtsets(idtset)%gw_qlwl(3,mxgw_nqlwl))
   allocate(dtsets(idtset)%kpt(3,mxnkpt))
   allocate(dtsets(idtset)%kptgw(3,mxnkptgw))
   allocate(dtsets(idtset)%kptns(3,mxnkpt))
   allocate(dtsets(idtset)%iatsph(mxnatsph))
   allocate(dtsets(idtset)%dynimage(mxnimage))
   allocate(dtsets(idtset)%istwfk(mxnkpt))
   allocate(dtsets(idtset)%nband(mxnkpt*mxnsppol))
   allocate(dtsets(idtset)%occ_orig(mxmband_upper*mxnkpt*mxnsppol))
   allocate(dtsets(idtset)%qmass(mxnnos))
   allocate(dtsets(idtset)%qptdm(3,mxnqptdm))
   allocate(dtsets(idtset)%symafm(mxnsym))
   allocate(dtsets(idtset)%symrel(3,3,mxnsym))
   allocate(dtsets(idtset)%tnons(3,mxnsym))
   allocate(dtsets(idtset)%wtatcon(3,mxnatom,mxnconeq))
   allocate(dtsets(idtset)%wtk(mxnkpt))
   dtsets(idtset)%symrel(:,:,:)=symrel_(:,:,1:mxnsym,idtset)
   dtsets(idtset)%symafm(:)    =symafm_(1:mxnsym,idtset)
   dtsets(idtset)%tnons (:,:)  =tnons_ (:,1:mxnsym,idtset)
 end do
!DEBUG
!write(std_out,*)'invars1m dtsets(:)%natvshift',dtsets(:)%natvshift
!ENDDEBUG
 deallocate(symafm_,symrel_,tnons_)

!DEBUG
!write(6,*)' invars1m : exit'
!write(6,*)' invars1m : dtsets(0)%vel_orig(:,1,1)=',dtsets(0)%vel_orig(:,1,1)
!write(6,*)' mxnatvshift,mxnsppol,mxnatpawu=',mxnatvshift,mxnsppol,mxnatpawu
!dtsets(0)%atvshift(:,:,:)=zero
!write(6,*)' succeeded to zero atvshift0 '
!ENDDEBUG

!stop
!ENDDEBUG

end subroutine invars1m
!!***
