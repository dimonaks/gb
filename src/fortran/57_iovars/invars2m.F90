!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars2m
!! NAME
!! invars2m
!!
!! FUNCTION
!! Initialisation phase - main input routine.
!! Big loop on the datasets :
!! - for each of the datasets, write one line about the crystallographic data
!! - call invars2, that read the eventual single dataset input values ;
!! - compute mgfft,mpw,nfft,... for this data set ;
!! - compute quantities for the susceptibility computation
!! - compute the memory needs for this data set.
!!  *** At the output of this routine, all the dtsets input variables are known ***
!! The content of dtsets should not be modified anymore afterwards.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bravais_(11,0:ndtset_alloc)=characteristics of Bravais lattice
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  mband_upper_(0:ndtset_alloc)=list of mband_upper values
!!  mpi_enreg=informations about MPI parallelization
!!  msym=default maximal number of symmetries
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!  string*(*)=character string containing all the input data.
!!   Initialized previously in instrng.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!!   initialized previously.
!!
!! SIDE EFFECTS
!!  dtpawuj=<type macro_uj_type>The potential shift and the occupation for DFT+U
!           calculations.
!!
!! NOTES
!! The outputs of this routine are the values of input variables,
!! their default value is stored at the index 0 of the last dimension
!! of their multi-dataset representation.
!!
!! PARENTS
!!      m_ab6_invars_f90
!!
!! CHILDREN
!!      clnmpi_fft,distrb2,getdim_nloc,getmpw,getng,initmpi_fft,invars2
!!      mati3inv,memorf,memory,metric,mkrdim,prtspgroup,setshells,symq3,wrtout
!!      wvl_atoms_data_set,wvl_memory
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine invars2m(bravais_,dtsets,iout,lenstr,&
&  mband_upper_,mpi_enreg,msym,ndtset,ndtset_alloc,npsp,pspheads,string)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_56_recipspace
 use interfaces_57_iovars, except_this_one => invars2m
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,msym,ndtset,ndtset_alloc,npsp
 character(len=*),intent(in) :: string
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: bravais_(11,0:ndtset_alloc),mband_upper_(0:ndtset_alloc)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables -------------------------------
!scalars
 integer :: cplex,exchn2n3d,getcell,idtset,ii,intxc,iprcch,iprcel
 integer :: iscf,isym,jdtset,lmnmax
 integer :: lmnmax_eff,lmnmaxso,lnmax,lnmax_eff,lnmaxso,mband,mband_upper
 integer :: me_fft,mffmem,mgfft,mgfftdg,mgfftdiel,mgfftf,mkmem,mpsang,mpspso
 integer :: mpssoang,mpw,mpw_k,mqgrid_ff,mqgrid_vl,n1xccc,natom,nfft,nfftdg
 integer :: nfftdiel,nfftf,nkpt,nproc_fft,npulayit,npwdiel,nqpt,nspden,nspinor
 integer :: nsppol,nsym,ntypat,occopt,optddk,optforces,optphon,optstress
 integer :: optstrs,paral_fft,pawcpxocc,pawmixdg,pawnhatxc,pawspnorb,pawstgylm,prtvol,ptgroupma,response
 integer :: spgroup,timrev,usepaw,useylm,xclevel
 real(dp) :: diecut,dilatmx,ecut,ecut_eff,ecutdg_eff,ecutsus,ucvol
 character(len=500) :: message
!arrays
 integer :: bravais(11),mkmems(3),ngfft(18),ngfftdg(18),ngfftdiel(18)
 integer :: ngfftf(18),nloalg(5),ngfftc(3)
 integer,allocatable :: istwfk(:),nband(:),symq(:,:,:),symrec(:,:,:)
 integer,allocatable :: symrel(:,:,:)
 real(dp) :: genafm(3),gmet(3,3),gprimd(3,3),kpt_diel(3),qphon(3),rmet(3,3)
 real(dp) :: rprimd(3,3)
 real(dp),allocatable :: kpt_with_shift(:,:),zionpsp(:)

!*************************************************************************

!PATCH invars2m PATCH init mpw_k
 mpw_k=0

 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
!  Space group output
   bravais(:)=bravais_(:,idtset)
   spgroup   =dtsets(idtset)%spgroup
   ptgroupma =dtsets(idtset)%ptgroupma
   genafm(:) =dtsets(idtset)%genafm(:)
!  DEBUG
!  write(6,*)' invars2m : ptgroupma=',ptgroupma
!  ENDDEBUG
   call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)

!  Get actual dimensions
   mband_upper  =mband_upper_(idtset)
   nkpt  =dtsets(idtset)%nkpt
   nsppol=dtsets(idtset)%nsppol
   ntypat=dtsets(idtset)%ntypat

!  Initialize ngfftc to the initial guess for the coarse mesh
   ngfftc(:) = 2

!  Allocate arrays
   allocate(istwfk(nkpt))
   allocate(kpt_with_shift(3,nkpt),zionpsp(npsp))
   zionpsp(:)=pspheads(1:npsp)%zionpsp

!  jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0 ! repetition of line 122
   usepaw=dtsets(idtset)%usepaw

!  Here, nearly all the remaining input variables are initialized
   call invars2(bravais,dtsets(idtset),iout,jdtset,lenstr,&
&   mband_upper,msym,npsp,pspheads,string,usepaw,zionpsp)

!  We set the internal wvl variables.
   if (dtsets(idtset)%usewvl == 1) then
     call wvl_atoms_data_set(dtsets(idtset)%acell_orig(1:3,1), dtsets(idtset))
   end if

   deallocate(zionpsp)

!  Compute mgfft,mpw,nfft for this data set
   call mkrdim(dtsets(idtset)%acell_orig(1:3,1),dtsets(idtset)%rprim_orig(1:3,1:3,1),rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!  For GW or BSE calculations, we only use (npwwfn|nshwfn|ecutwfn) G-vectors read from the KSS file,
!  therefore the FFT box for the density should be defined according to ecut=ecutwfn.
   if ( ANY(dtsets(idtset)%optdriver == (/RUNL_SCREENING,RUNL_SIGMA,RUNL_SCGW,RUNL_BSE/)) ) then

     call setshells(dtsets(idtset)%ecutwfn,dtsets(idtset)%npwwfn,dtsets(idtset)%nshwfn,&
&     dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'wfn',ucvol)

!    MG: Hack to avoid portability problems under gfortran and g95:
!    getng and getmpw are indeed quite sensitive if ecut is small
!    and, in the GW tests, mpw and ngfft might depend on the compiler used.
!    the problem shows up if we use npwwfn instead of ecutwfn, a good
!    reason for removing nshwfn and npwwfn!
     dtsets(idtset)%ecutwfn=dtsets(idtset)%ecutwfn-tol14
!    MG: This is a kind of a hack, but the problem is ecutwfn that is too much redundant!
     dtsets(idtset)%ecut=dtsets(idtset)%ecutwfn

!    Close the shell for (W|chi0)
     call setshells(dtsets(idtset)%ecuteps,dtsets(idtset)%npweps,dtsets(idtset)%nsheps,&
&     dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'eps',ucvol)

!    Close the shell for the exchange term.
     call setshells(dtsets(idtset)%ecutsigx,dtsets(idtset)%npwsigx,dtsets(idtset)%nshsigx,&
&     dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'sigx',ucvol)

   end if ! (SIGMA|SCREENING|SCGW|BSE)

   ecut     =dtsets(idtset)%ecut
   dilatmx  =dtsets(idtset)%dilatmx
   ngfft(:) =dtsets(idtset)%ngfft(:)
   istwfk(:)=dtsets(idtset)%istwfk(1:nkpt)
   nsym     =dtsets(idtset)%nsym

   allocate(symrel(3,3,nsym))
   symrel(:,:,1:nsym)=dtsets(idtset)%symrel(:,:,1:nsym)
   ecut_eff=ecut*dilatmx**2

!  MPIWF : here, set up the complete ngfft, containing the information
!  for the parallelisation of the FFT
!  Default values for sequentiel case
   paral_fft=0 ; mpi_enreg%paral_fft=0
   nproc_fft=1;mpi_enreg%nproc_fft=1
   me_fft=0;mpi_enreg%me_fft=0
   mpi_enreg%fft_option_lob=1

   if(mpi_enreg%paral_compil_fft==1 .and. dtsets(idtset)%usewvl == 0)then
!    Compute the actual values
     paral_fft=1           ! parallelisation over FFT
!    Should fill the values of nproc_fft and me_fft by MPI calls
     mband=maxval(dtsets(idtset)%nband(1:nkpt*nsppol))
     allocate(nband(nkpt*nsppol))
     nband(1:nkpt*nsppol)=dtsets(idtset)%nband(1:nkpt*nsppol)
     allocate(mpi_enreg%proc_distrb(nkpt,mband,nsppol))
     call distrb2(mband,nband,nkpt,nsppol,mpi_enreg)
     dtsets(idtset)%mgfft=0
     call initmpi_fft(dtsets(idtset),mpi_enreg)
     nproc_fft=mpi_enreg%nproc_fft
     me_fft=mpi_enreg%me_fft
     mpi_enreg%paral_fft=paral_fft
     deallocate(mpi_enreg%proc_distrb,nband)
   end if

   if (usepaw==1) then
     write(message,'(2a)') ch10,' getng is called for the coarse grid:'
     call wrtout(std_out,message,'COLL')
   end if
   call getng(dtsets(idtset)%boxcutmin,ecut_eff,gmet,me_fft,mgfft,nfft,&
&   ngfft,nproc_fft,nsym,mpi_enreg%fft_option_lob,paral_fft,symrel)
   dtsets(idtset)%ngfft(:)=ngfft(:)
   dtsets(idtset)%mgfft=mgfft
   dtsets(idtset)%nfft=nfft
   if (mpi_enreg%paral_compil_fft==1) then
!    creation of arrays for FFT parallelization here because in initmpi_fft mgfft is not known
     allocate(mpi_enreg%nplanes_fft(dtsets(idtset)%nkpt))
     allocate(mpi_enreg%ind_fft_planes(dtsets(idtset)%nkpt,ngfft(2)))
   end if

   kpt_with_shift(:,:)=dtsets(idtset)%kpt(:,1:nkpt)/dtsets(idtset)%kptnrm
   nqpt=dtsets(idtset)%nqpt
   response=0
   exchn2n3d=dtsets(idtset)%exchn2n3d
   nproc_fft=ngfft(10) ; me_fft=ngfft(11)
   if(dtsets(idtset)%rfddk/=0 .or. dtsets(idtset)%rfelfd/=0 .or. &
&   dtsets(idtset)%rfmgfd/=0 .or.  dtsets(idtset)%rfphon/=0 .or. &
&   dtsets(idtset)%rfstrs/=0 .or.  dtsets(idtset)%rfuser/=0    ) response=1
   if(response/=0)then
!    This value of mpw is used in the first part of respfn.f
     call getmpw(ecut_eff,exchn2n3d,gmet,istwfk,kpt_with_shift,mpi_enreg,mpw_k,nkpt)
   end if
   qphon(:)=zero
   if(nqpt/=0)then
     qphon(:)=dtsets(idtset)%qptn(:)
     kpt_with_shift(1,:)=kpt_with_shift(1,:)+qphon(1)
     kpt_with_shift(2,:)=kpt_with_shift(2,:)+qphon(2)
     kpt_with_shift(3,:)=kpt_with_shift(3,:)+qphon(3)
   end if
   if (dtsets(idtset)%usewvl == 0) then
     call getmpw(ecut_eff,exchn2n3d,gmet,istwfk,kpt_with_shift,mpi_enreg,mpw,nkpt)
   else
     mpw = 0
   end if

!  The dimensioning, in the RF case, should be done only with mpw,
!  but mpw is used in the first part of respfn.f, and should at least
!  be equal to mpw_k . The chosen way to code is not optimal, only convenient :
!  it leads to a small waste of memory.
   if(response/=0 .and. mpw_k>mpw)mpw=mpw_k

!  Compute mgfftdiel,npwdiel,nfftdiel for this data set
   iprcel=dtsets(idtset)%iprcel
   iscf=dtsets(idtset)%iscf
   if((modulo(iprcel,100)>=20 .and.modulo(iprcel,100) < 71).or. iscf==-1)then

!    Get diecut, and the fft grid to be used for the susceptibility computation
     diecut=abs(dtsets(idtset)%diecut)
     if( dtsets(idtset)%diecut < zero )then
       ecutsus=ecut
     else
       ecutsus= ( sqrt(ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
     end if
!    Beware, for the dielectric matrix fftalg=ngfftdiel(7) is default here
     ngfftdiel(1:3)=0 ; ngfftdiel(7)=101 ; ngfftdiel(8:18)=dtsets(idtset)%ngfft(8:18)
     if(iscf==-1)ngfftdiel(7)=102
     ecut_eff=ecutsus*dilatmx**2
     call getng(dtsets(idtset)%boxcutmin,ecut_eff,gmet,me_fft,mgfftdiel,nfftdiel,ngfftdiel,&
&     nproc_fft,nsym,mpi_enreg%fft_option_lob,paral_fft,symrel)
!    Compute the size of the dielectric matrix : npwdiel
     kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
     ecut_eff=diecut*dilatmx**2
     call getmpw(ecut_eff,exchn2n3d,gmet,(/1/),kpt_diel,mpi_enreg,npwdiel,1)

   else

     npwdiel=1 ; mgfftdiel=1 ; nfftdiel=1 ; ngfftdiel(1:8)=1

   end if

   dtsets(idtset)%ngfft(:)=ngfft(:)

!  In case of PAW, compute fine FFT parameters
   if (usepaw==1) then
     ecutdg_eff=dtsets(idtset)%pawecutdg*dtsets(idtset)%dilatmx**2
     ngfftdg(:)=dtsets(idtset)%ngfftdg(:)
     write(message,'(2a)') ch10,' getng is called for the fine grid:'
     call wrtout(std_out,message,'COLL')
!    Start with the coarse mesh as an initial guess for the fine mesh
!    This ensures that the fine mesh will not be any coarser than the coarse mesh in each dimension
     ngfftc(:) = ngfft(1:3)
     call getng(dtsets(idtset)%bxctmindg,ecutdg_eff,gmet,me_fft,mgfftdg,&
&     nfftdg,ngfftdg,nproc_fft,nsym,mpi_enreg%fft_option_lob,paral_fft,symrel,ngfftc)
     dtsets(idtset)%ngfftdg(:)=ngfftdg(:)
     dtsets(idtset)%mgfftdg=mgfftdg
     dtsets(idtset)%nfftdg=nfftdg
   end if
!  Transfer other data needed to compute the memory
   getcell=dtsets(idtset)%getcell
   prtvol=dtsets(idtset)%prtvol
   intxc=dtsets(idtset)%intxc
   xclevel=dtsets(idtset)%xclevel
   useylm=dtsets(idtset)%useylm
   optforces=dtsets(idtset)%optforces
   if (dtsets(idtset)%toldff>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%tolrff>tol16.and.optforces==0) optforces=1
   if (dtsets(idtset)%ionmov>tol16.and.optforces==0) optforces=1
   optstress=dtsets(idtset)%optstress
   optddk=0;optphon=0;optstrs=0
   if (dtsets(idtset)%rfddk>0) optddk=1
   if (dtsets(idtset)%rfelfd>0.or.dtsets(idtset)%rf1elfd>0.or.&
&   dtsets(idtset)%rf2elfd>0.or.dtsets(idtset)%rf3elfd>0) optddk=1
   if (dtsets(idtset)%rfmgfd>0) optddk=1
   if (dtsets(idtset)%rfphon>0.or.dtsets(idtset)%rf1phon>0.or.&
&   dtsets(idtset)%rf2phon>0.or.dtsets(idtset)%rf3phon>0) optphon=1
   if (dtsets(idtset)%rfstrs>0) optstrs=1

   allocate(nband(nkpt*nsppol))
   nband(1:nkpt*nsppol)=dtsets(idtset)%nband(1:nkpt*nsppol)
   mband=maxval(nband(1:nkpt*nsppol))
   dtsets(idtset)%mband=mband

!  mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! Likely problems with the HP compiler
!  n1xccc=maxval(pspheads(1:npsp)%xccc)
   mpsang=1
   n1xccc=pspheads(1)%xccc
   do ii=1,npsp
     mpsang=max(pspheads(ii)%lmax+1,mpsang)
     n1xccc=max(pspheads(ii)%xccc,n1xccc)
   end do

!  Determine the maximum number of projectors, for the set of pseudo atom
   call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtsets(idtset)%mixalch,&
&   npsp,dtsets(idtset)%npspalch,ntypat,dtsets(idtset)%ntypalch,pspheads)

   nspinor=dtsets(idtset)%nspinor

!  Treatment of the effect of using a spin-orbit part
!  Warning : mpspso is different for each dataset.
   mpspso=1
   do ii=1,npsp
     if(nspinor/=1)then
       if(pspheads(ii)%pspso/=0)then
         if(dtsets(idtset)%so_psp(ii)/=0)then
           mpspso=2
         end if
       end if
     end if
   end do
!  In case of no spin-orbit
   if(mpspso==1)then
     mpssoang=mpsang ; lmnmax_eff =lmnmax; lnmax_eff =lnmax
   else ! spin-orbit will be used
     mpssoang=2*mpsang-1 ; lmnmax_eff =lmnmaxso ; lnmax_eff =lnmaxso
   end if
!  lmnmax is not used if the Ylm are not used
   if (useylm==0) lmnmax_eff =lnmax_eff

   iprcch=dtsets(idtset)%iprcch
   nloalg(:)=dtsets(idtset)%nloalg(:)
   mffmem=dtsets(idtset)%mffmem
   mqgrid_ff=dtsets(idtset)%mqgrid
   if (usepaw==0) mqgrid_vl=dtsets(idtset)%mqgrid
   if (usepaw==1) mqgrid_vl=dtsets(idtset)%mqgriddg
   natom=dtsets(idtset)%natom
   npulayit=dtsets(idtset)%npulayit
   nspden=dtsets(idtset)%nspden
   occopt=dtsets(idtset)%occopt
   pawcpxocc=dtsets(idtset)%pawcpxocc
   pawmixdg=dtsets(idtset)%pawmixdg
   pawnhatxc=dtsets(idtset)%pawnhatxc
   pawspnorb=dtsets(idtset)%pawspnorb
   pawstgylm=dtsets(idtset)%pawstgylm
   if (usepaw==0) then
     mgfftf=mgfft;nfftf=nfft;ngfftf(:)=ngfft(:)
   else
     mgfftf=mgfftdg;nfftf=nfftdg;ngfftf(:)=ngfftdg(:)
   end if

!  Compute the memory needs for this data set.
   if(response==0)then
     if (dtsets(idtset)%usewvl == 0) then
       mkmem=dtsets(idtset)%mkmem

       call memory(n1xccc,getcell,idtset,dtsets(idtset)%icoulomb,intxc,dtsets(idtset)%ionmov,iout,iprcch,&
&       iprcel,iscf,jdtset,lmnmax_eff,lnmax_eff,mband,mffmem,mgfft,mgfftdiel,mgfftf,mkmem,&
&       mpi_enreg,mpsang,mpssoang,mpw,mqgrid_ff,mqgrid_vl,natom,nband,nfft,nfftdiel,nfftf,&
&       ngfft,ngfftdiel,ngfftf,nkpt,nloalg,npsp,npulayit,npwdiel,nspden,nspinor,&
&       nsppol,nsym,ntypat,occopt,optforces,1,optstress,pawcpxocc,pawmixdg,&
&       pawnhatxc,pawspnorb,pawstgylm,prtvol,pspheads,&
&       dtsets(idtset)%typat,ucvol,usepaw,useylm,xclevel)
     else
       if (mpi_enreg%me == 0) then
         call wvl_memory(dtsets(idtset), idtset, mpi_enreg, npsp, 1, pspheads)
       end if
     end if
   else
!    Compute the value of cplex, for which one needs symrec
     allocate(symq(4,2,nsym),symrec(3,3,nsym))
     do isym=1,nsym
       call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
     end do
     call symq3(nsym,qphon,symq,symrec,timrev)
     cplex=2-timrev
     deallocate(symq,symrec)
     mkmems(1)=dtsets(idtset)%mkmem
     mkmems(2)=dtsets(idtset)%mkqmem
     mkmems(3)=dtsets(idtset)%mk1mem
     call memorf(cplex,n1xccc,getcell,idtset,intxc,iout,iprcel,&
&     iscf,jdtset,lmnmax_eff,lnmax_eff,mband,mffmem,mgfft,&
&     mkmems,mpi_enreg,mpsang,mpssoang,mpw,mqgrid_ff,natom,nband,nfft,&
&     ngfft,nkpt,nloalg,nspden,nspinor,nsppol,nsym,&
&     ntypat,occopt,optddk,optphon,1,optstrs,prtvol,useylm,xclevel)
   end if
!  Copy input values from file

   dtsets(idtset)%mpw=mpw

!  Deallocate temporary arrays (when they will really be temporary !)
   deallocate(istwfk,kpt_with_shift,nband)
   deallocate(symrel)

   if(mpi_enreg%paral_compil_fft==1)then
     call clnmpi_fft(mpi_enreg)
   end if
 end do

!DEBUG
!write(6,*)' invars2m : exit '
!ENDDEBUG

end subroutine invars2m
!!***
