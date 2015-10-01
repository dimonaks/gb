!{\src2tex{textfont=tt}}
!!****f* ABINIT/indefo
!! NAME
!! indefo
!!
!! FUNCTION
!! Initialisation phase : defaults values for most input variables
!! (some are initialized earlier)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG,MM,FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are given a default value here.
!!   The dataset with number 0 should be the reference default value
!!   in the remaining of the code.
!!
!! NOTES
!! The outputs of this routine are the defaults values of input
!! variables, stored at the index 0 of the last dimension of their
!! multi-dataset representation.
!!
!! PARENTS
!!      m_ab6_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine indefo(dtsets,ndtset_alloc)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_gwdefs

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc
!arrays
 type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------------
!scalars
 integer :: idtset,ii,jdtset

!******************************************************************

 DBG_ENTER("COLL")

!
!Set up default values. All variables to be output in outvars.f
!should have a default, even if a nonsensible one can be
!chosen to garantee print in that routine.


!These variables have already been initialized, for idtset/=0
 dtsets(0)%kptrlatt(1:3,1:3)=0
 dtsets(0)%ptgroupma=0
 dtsets(0)%spgroup=0
 dtsets(0)%shiftk(:,:)=half
 dtsets(0)%tolsym=tol8
 dtsets(0)%znucl(:)=zero

!WARNING : set default in all datasets, including idtset=0 !!!
!Use alphabetic order

 do idtset=0,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset
!  A
   dtsets(idtset)%algalch(:)=0
   dtsets(idtset)%alpha=one
!  Note that this default value might be overriden for specific datasets later, in case of MPI_IO
   dtsets(idtset)%accesswff=IO_MODE_FORTRAN
   dtsets(idtset)%atvshift(:,:,:)=zero
   dtsets(idtset)%awtr=1
!  B
   dtsets(idtset)%bdberry(1:4)=0
   dtsets(idtset)%bdeigrf=-1
   dtsets(idtset)%bdgw=0
   dtsets(idtset)%bmass=ten
   dtsets(idtset)%boxcenter(1:3)=half
   dtsets(idtset)%boxcutmin=two
   dtsets(idtset)%brvltt=0
   dtsets(idtset)%bxctmindg=two
!  C
   dtsets(idtset)%charge=zero
   dtsets(idtset)%chkexit=0
   dtsets(idtset)%chksymbreak=1
   dtsets(idtset)%corecs(:) = zero
!  D
   dtsets(idtset)%delayperm=0
   dtsets(idtset)%diecut=2.2_dp
   dtsets(idtset)%dielng=1.0774841_dp
   dtsets(idtset)%diemac=1.0d6
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%diemix=one
   else
     dtsets(idtset)%diemix=0.7_dp
   end if
   dtsets(idtset)%diemixmag=dtsets(idtset)%diemix
   dtsets(idtset)%diegap=0.1_dp
   dtsets(idtset)%dielam=half
   dtsets(idtset)%dilatmx=one
   dtsets(idtset)%dmatpuopt=2
   if (size(dtsets(idtset)%dmatpawu,4)>0) dtsets(idtset)%dmatpawu=-10._dp
   dtsets(idtset)%dmatudiag=0
   dtsets(idtset)%dmft_iter=0
   dtsets(idtset)%dmft_nwli=0
   dtsets(idtset)%dmft_nwlo=0
   dtsets(idtset)%dmft_mxsf=0.3_dp
   dtsets(idtset)%dmft_rslf=0
   dtsets(idtset)%dmft_solv=4
   dtsets(idtset)%dmftbandi=0
   dtsets(idtset)%dmftbandf=0
   dtsets(idtset)%dmftcheck=0
   dtsets(idtset)%dosdeltae=zero
   dtsets(idtset)%dtion=100.0_dp
   dtsets(idtset)%dynimage=1
!  E
   dtsets(idtset)%ecut=-one
   dtsets(idtset)%ecuteps=zero
   dtsets(idtset)%ecutsigx=zero ! The true default value is ecut . This is defined in invars2.F90
   dtsets(idtset)%ecutsm=zero
   dtsets(idtset)%ecutwfn=zero ! The true default value is ecut . This is defined in invars2.F90
   dtsets(idtset)%effmass=one
   dtsets(idtset)%enunit=0
   dtsets(idtset)%eshift=zero
   dtsets(idtset)%esmear=0.04_dp
   dtsets(idtset)%exchn2n3d=0
   dtsets(idtset)%exchmix=quarter
!  F
   dtsets(idtset)%fftgw=21
   dtsets(idtset)%fft_opt_lob=0
   dtsets(idtset)%fixmom=-99.99_dp
   dtsets(idtset)%freqsusin=one
   dtsets(idtset)%freqsuslo=one
   dtsets(idtset)%optfreqsus=2
   dtsets(idtset)%freqremax=zero
   dtsets(idtset)%freqspmax=zero
   dtsets(idtset)%friction=0.001_dp
   dtsets(idtset)%frzfermi=0
   dtsets(idtset)%fxcartfactor=one ! Should be adjusted to the H2 conversion factor
!  G
   dtsets(idtset)%gw_nstep =30
   dtsets(idtset)%gwgamma =0
   if ( dtsets(idtset)%gw_nqlwl > 0 ) then
     dtsets(idtset)%gw_qlwl(:,:)=zero
     dtsets(idtset)%gw_qlwl(1,1)=0.00001_dp
     dtsets(idtset)%gw_qlwl(2,1)=0.00002_dp
     dtsets(idtset)%gw_qlwl(3,1)=0.00003_dp
   end if
   dtsets(idtset)%gw_EET=-1
   dtsets(idtset)%gw_EET_nband=-1
   dtsets(idtset)%gw_sigxcore=0
   dtsets(idtset)%gw_sctype =  GWSC_one_shot
   dtsets(idtset)%gw_toldfeig=0.1/Ha_eV
   dtsets(idtset)%getbseig=0
   dtsets(idtset)%getcell =0
   dtsets(idtset)%getddk  =0
   dtsets(idtset)%getden  =0
   dtsets(idtset)%getkss  =0
   dtsets(idtset)%getocc  =0
   dtsets(idtset)%getqps  =0
   dtsets(idtset)%getscr  =0
   dtsets(idtset)%getsuscep=0
   dtsets(idtset)%getvel  =0
   dtsets(idtset)%getwfk  =0
   dtsets(idtset)%getwfq  =0
   dtsets(idtset)%getxcart=0
   dtsets(idtset)%getxred =0
   dtsets(idtset)%get1den =0
   dtsets(idtset)%get1wf  =0
   dtsets(idtset)%gwcalctyp=0
   dtsets(idtset)%gwcomp=0
   dtsets(idtset)%gwencomp=2.0_dp
   dtsets(idtset)%gwmem=11
   dtsets(idtset)%gwpara=1
   dtsets(idtset)%gwrpacorr=0

!  I
   if(dtsets(idtset)%natsph/=0) then
!    do not use iatsph(:) but explicit boundaries
!    to avoid to read to far away in the built array (/ ... /)
     dtsets(idtset)%iatsph(1:dtsets(idtset)%natsph)=(/ (ii,ii=1,dtsets(idtset)%natsph) /)
   else
     dtsets(idtset)%iatsph(:)=0
   end if
   dtsets(idtset)%iboxcut=0
   dtsets(idtset)%idyson=1
   dtsets(idtset)%icutcoul=3
   dtsets(idtset)%ieig2rf=0
   dtsets(idtset)%iextrapwf=0
   dtsets(idtset)%ikhxc=0
   dtsets(idtset)%inclvkb=1
   dtsets(idtset)%intexact=0
   dtsets(idtset)%intxc=0
   dtsets(idtset)%imgmov=0
   dtsets(idtset)%ionmov=0
   dtsets(idtset)%iprcch=2
   dtsets(idtset)%iprcel=0
   dtsets(idtset)%iprctfvw=0
   dtsets(idtset)%iprcfc=0
   dtsets(idtset)%irdbseig=0
   dtsets(idtset)%irdddk=0
   dtsets(idtset)%irdden=0
   dtsets(idtset)%irdkss=0
   dtsets(idtset)%irdqps=0
   dtsets(idtset)%irdscr=0
   dtsets(idtset)%irdsuscep=0
   dtsets(idtset)%irdwfk=0
   dtsets(idtset)%irdwfq=0
   dtsets(idtset)%ird1wf=0
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%iscf=7
   else
     dtsets(idtset)%iscf=17
   end if
   dtsets(idtset)%isecur=0
   dtsets(idtset)%istatr = -1
   dtsets(idtset)%istatshft = -1
   dtsets(idtset)%istwfk(:)=0
   dtsets(idtset)%ixc=1
   dtsets(idtset)%ixcpositron=1
!  J
   dtsets(idtset)%jpawu(:)=zero
!  K
!  We should not do that, kberry can be smaller than 1:20
!  dtsets(idtset)%kberry(1:3,1:20)=0
   dtsets(idtset)%kberry(1:3,:)=0
   dtsets(idtset)%kpt(:,:)=zero
   dtsets(idtset)%kptgw(:,:)=zero
   dtsets(idtset)%kptnrm=one
   dtsets(idtset)%kptopt=1
   dtsets(idtset)%kptrlen=30.0_dp
   dtsets(idtset)%kssform=1
!  L
   dtsets(idtset)%ldgapp=0
!  dtsets(idtset)%lexexch(:)=-1 initalized in invars1
   dtsets(idtset)%localrdwf=1
!  dtsets(idtset)%lpawu(:)=-1 initalized in invars1

!  M
   dtsets(idtset)%mband = -1
   dtsets(idtset)%mdftemp=300.0_dp
   dtsets(idtset)%mditemp=300.0_dp
   dtsets(idtset)%mdwall=10000_dp
   dtsets(idtset)%mffmem=1
   dtsets(idtset)%mgfft = -1
   dtsets(idtset)%mgfftdg = -1
   dtsets(idtset)%mixalch(:,:)=zero
   dtsets(idtset)%mpw = -1
   dtsets(idtset)%mqgrid=3001
   dtsets(idtset)%mqgriddg=3001
!  N
   dtsets(idtset)%natrd = -1
   dtsets(idtset)%nband(:)=0
   dtsets(idtset)%nbandsus=-1
   dtsets(idtset)%nbdblock=1
   dtsets(idtset)%nbdbuf=0
   dtsets(idtset)%nberry=1
   dtsets(idtset)%nbandkss=0
   dtsets(idtset)%nctime=0
   dtsets(idtset)%ndtset = -1
   dtsets(idtset)%ndynimage = 1
   dtsets(idtset)%npwkss=0
   dtsets(idtset)%ndyson=-1
   dtsets(idtset)%nfft = -1
   dtsets(idtset)%nfftdg = -1

   dtsets(idtset)%nfreqim=0
   dtsets(idtset)%nfreqre=0
   dtsets(idtset)%nfreqsp=0
   dtsets(idtset)%nfreqsus=0
   dtsets(idtset)%diismemory=8
!  dtsets(idtset)%npband=1
!  dtsets(idtset)%npfft=1

   dtsets(idtset)%npulayit=7
!  dtsets(idtset)%bandpp=1   !initialized in invars1m

!  ngfft is a special case
   dtsets(idtset)%ngfft(1:8)=0

!  fftalg=ngfft(7) is machine-dependent.
   dtsets(idtset)%ngfft(7)=112
#if defined FC_FUJITSU
!  For the vpp fujitsu, it is better to have the last digit at 1.
   dtsets(idtset)%ngfft(7)=111
#elif defined HAVE_FFT_ASL
!  For the NEC computer, the library FFT routine is (presently) faster than the routines from Stefan
   dtsets(idtset)%ngfft(7)=200
#elif defined HAVE_FFT_FFTW3
   dtsets(idtset)%ngfft(7)=312
#endif

   dtsets(idtset)%ngfft(8)=16
!  fftcache=ngfft(8) is machine-dependent.
#if defined FC_FUJITSU || defined HAVE_FFT_ASL || defined FC_HITACHI
   dtsets(idtset)%ngfft(8)=8096    ! Large value for the vector machines
#elif defined i386
   dtsets(idtset)%ngfft(8)=256     ! Was optimized for my PII 450MHz
#endif

!  DEBUG
!  dtsets(idtset)%ngfft(7)=401
!  ENDDEBUG
   if (dtsets(idtset)%paral_kgb>0.and.idtset>0) then
     dtsets(idtset)%ngfft(7)    = 401  ! fftalg      = 401
   end if

   dtsets(idtset)%ngfftdg(:)=dtsets(idtset)%ngfft(:)
   dtsets(idtset)%nline=4

!  nloalg is also a special case
   dtsets(idtset)%nloalg(1)=4
   dtsets(idtset)%nloalg(2)=4
   dtsets(idtset)%nloalg(3)=199
   dtsets(idtset)%nloalg(4)=10
   dtsets(idtset)%nloalg(5)=dtsets(idtset)%usepaw
#if defined i386
!  nloalg(2) is machine-dependent, the default is 4, but 8 is better for the P6 and hp .
   dtsets(idtset)%nloalg(2)=8
#elif defined HAVE_FFT_ASL
   dtsets(idtset)%nloalg(1)=2    ! Important to use opernl2 for this vector machine
#endif

   if (dtsets(idtset)%paral_kgb>0.and.idtset>0) then
     dtsets(idtset)%intxc       = 0
     dtsets(idtset)%fft_opt_lob = 2
   end if

   dtsets(idtset)%nnsclo=0
   dtsets(idtset)%nomegasf=100
   dtsets(idtset)%nomegasrd=9
   dtsets(idtset)%nomegasi=12
   dtsets(idtset)%noseinert=1.0d5
   dtsets(idtset)%npspalch=0
   dtsets(idtset)%npweps=0
   dtsets(idtset)%npwsigx=0
   dtsets(idtset)%npwwfn=0
   dtsets(idtset)%nqpt=0
   dtsets(idtset)%nscforder=16
   dtsets(idtset)%nsheps=0
   dtsets(idtset)%nshiftk=1
   dtsets(idtset)%nshsigx=0
   dtsets(idtset)%nshwfn=0

!  dtsets(idtset)%nsppol  defined in invars1m
!  dtsets(idtset)%nspinor defined in invars1m
!  dtsets(idtset)%nspden  defined in invars1m

   dtsets(idtset)%nstep=30
   dtsets(idtset)%ntime=0
   dtsets(idtset)%ntimimage=1
   dtsets(idtset)%ntypalch=0
   dtsets(idtset)%ntyppure = -1
   dtsets(idtset)%nwfshist=0
!  O
   dtsets(idtset)%occopt=1
   dtsets(idtset)%occ_orig(:)=zero
   dtsets(idtset)%omegasrdmax=1.0_dp/Ha_eV  ! = 1eV
   dtsets(idtset)%omegasimax=50/Ha_eV
   dtsets(idtset)%optcell=0
   dtsets(idtset)%optforces=1
   dtsets(idtset)%optstress=1
   dtsets(idtset)%optnlxccc=1
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%ortalg=2
   else
     dtsets(idtset)%ortalg=-2
   end if
!  P
   dtsets(idtset)%pawcpxocc=1
   dtsets(idtset)%pawecutdg=-one
   dtsets(idtset)%pawfatbnd=0
   dtsets(idtset)%pawlcutd=10
   dtsets(idtset)%pawlmix=10
   dtsets(idtset)%pawmixdg=0
   dtsets(idtset)%pawnhatxc=1
   dtsets(idtset)%pawntheta=12
   dtsets(idtset)%pawnphi=13
   dtsets(idtset)%pawnzlm=1
   dtsets(idtset)%pawoptmix=0
   dtsets(idtset)%pawovlp=5._dp
   dtsets(idtset)%pawprtden=0
   dtsets(idtset)%pawprtdos=0
   dtsets(idtset)%pawprtvol=0
   dtsets(idtset)%pawprtwf=0
   dtsets(idtset)%pawstgylm=1
   dtsets(idtset)%pawujat=1
   dtsets(idtset)%pawujrad=20.0_dp
   dtsets(idtset)%pawujv=0.1_dp/Ha_eV
   dtsets(idtset)%pawusecp=1
   dtsets(idtset)%pawxcdev=1
   dtsets(idtset)%pitransform=0
   dtsets(idtset)%ptcharge(:) = zero
   dtsets(idtset)%positron=0
   dtsets(idtset)%posnstep=50
   dtsets(idtset)%posocc=one
   dtsets(idtset)%postoldfe=0.000001_dp
   dtsets(idtset)%postoldff=zero
   dtsets(idtset)%ppmodel=1
   dtsets(idtset)%ppmfrq=zero
   dtsets(idtset)%prepanl=0
   dtsets(idtset)%prepgkk=0
   dtsets(idtset)%prtbbb=0
   dtsets(idtset)%prtbltztrp=0
   dtsets(idtset)%prtcif=0
   dtsets(idtset)%prtcml=0
   dtsets(idtset)%prtcs=0
   dtsets(idtset)%prtden=1
   dtsets(idtset)%prtdensph=0
   dtsets(idtset)%prtdipole=0
   dtsets(idtset)%prtdos=0
   dtsets(idtset)%prtdosm=0
   dtsets(idtset)%prtefg=0
   dtsets(idtset)%prteig=1
   dtsets(idtset)%prtelf=0
   dtsets(idtset)%prtfc=0
   dtsets(idtset)%prtfsurf=0
   dtsets(idtset)%prtgden=0
   dtsets(idtset)%prtgeo=0
   dtsets(idtset)%prtgkk=0
   dtsets(idtset)%prtkden=0
   dtsets(idtset)%prtkpt = -1
   dtsets(idtset)%prtlden=0
   dtsets(idtset)%prtnabla=0
   dtsets(idtset)%prtnest=0
   dtsets(idtset)%prtposcar=0
   dtsets(idtset)%prtpot=0
   dtsets(idtset)%prtspcur=0
   dtsets(idtset)%prtstm=0
   dtsets(idtset)%prtvha=0
   dtsets(idtset)%prtvhxc=0
   dtsets(idtset)%prtvxc=0
   dtsets(idtset)%prtvol=0
   dtsets(idtset)%prtwant=0
   dtsets(idtset)%prtwf=1
   dtsets(idtset)%prtxangst = 1
   dtsets(idtset)%prtxcart = 1
   dtsets(idtset)%prtxred = 1
   dtsets(idtset)%prtxml = 0
   dtsets(idtset)%prt1dm=0
!  Q
   dtsets(idtset)%qmass(:)=ten
   dtsets(idtset)%qprtrb(1:3)=0
   dtsets(idtset)%qpt(1:3)=zero
   dtsets(idtset)%qptdm(:,:)=zero
   dtsets(idtset)%qptn(1:3)=zero
   dtsets(idtset)%qptnrm=one
   dtsets(idtset)%quadmom(:) = zero
!  R
   dtsets(idtset)%ratsph(:)=two
   dtsets(idtset)%recefermi=zero
   dtsets(idtset)%recgratio=1
   dtsets(idtset)%recnpath=500
   dtsets(idtset)%recnrec=10
   dtsets(idtset)%recrcut=zero
   dtsets(idtset)%recptrott=0
   dtsets(idtset)%rectesteg=0
   dtsets(idtset)%rectolden=zero
   dtsets(idtset)%rcut=zero
   dtsets(idtset)%rdmnb=0
   dtsets(idtset)%restartxf=0
   dtsets(idtset)%rfasr=0
   dtsets(idtset)%rfatpol(1:2)=1
   dtsets(idtset)%rfddk=0
   dtsets(idtset)%rfdir(1:3)=0
   dtsets(idtset)%rfelfd=0
   dtsets(idtset)%rfmgfd=0
   dtsets(idtset)%rfmeth=1
   dtsets(idtset)%rfphon=0
   dtsets(idtset)%rfstrs=0
   dtsets(idtset)%rfuser=0

   dtsets(idtset)%rf1atpol(1:2)=1
   dtsets(idtset)%rf1dir(1:3)=0
   dtsets(idtset)%rf1elfd=0
   dtsets(idtset)%rf1phon=0

   dtsets(idtset)%rf2atpol(1:2)=1
   dtsets(idtset)%rf2dir(1:3)=0
   dtsets(idtset)%rf2elfd=0
   dtsets(idtset)%rf2phon=0

   dtsets(idtset)%rf3atpol(1:2)=1
   dtsets(idtset)%rf3dir(1:3)=0
   dtsets(idtset)%rf3elfd=0
   dtsets(idtset)%rf3phon=0

   dtsets(idtset)%rhoqpmix=one

!  S
   dtsets(idtset)%scphon_supercell(:)=1
   dtsets(idtset)%scphon_temp=zero
   dtsets(idtset)%sciss=zero
   dtsets(idtset)%signperm=1
   dtsets(idtset)%slabwsrad=zero
   dtsets(idtset)%smdelta=0
   dtsets(idtset)%soenergy=zero
   dtsets(idtset)%spbroad=0.1
   dtsets(idtset)%spgaxor = -1
   dtsets(idtset)%spgorig = -1
   dtsets(idtset)%spmeth=0
   dtsets(idtset)%spnorbscl=one
   dtsets(idtset)%stmbias=zero
   dtsets(idtset)%strfact=100.0_dp
   dtsets(idtset)%strprecon=one
   dtsets(idtset)%strtarget(1:6)=zero
   dtsets(idtset)%supercell(:)=1
   dtsets(idtset)%suskxcrs=0
   dtsets(idtset)%symchi=0
   dtsets(idtset)%symsigma=0
!  T
   dtsets(idtset)%td_maxene=zero
   dtsets(idtset)%td_mexcit=0
   dtsets(idtset)%tfkinfunc=0
   dtsets(idtset)%timopt = -1
   dtsets(idtset)%tl_nprccg = 30
   dtsets(idtset)%tl_radius = real(0, dp)
   dtsets(idtset)%tphysel=zero
   dtsets(idtset)%toldfe=zero
   dtsets(idtset)%toldff=zero
   dtsets(idtset)%tolimg=5.0d-5
   dtsets(idtset)%tolrff=zero
   dtsets(idtset)%tolmxf=5.0d-5
   dtsets(idtset)%tolvrs=zero
   dtsets(idtset)%tolwfr=zero
   dtsets(idtset)%tsmear=0.04_dp
!  U
   dtsets(idtset)%upawu(:)=zero
   dtsets(idtset)%usekden=0
   dtsets(idtset)%userec=0
!  dtsets(idtset)%usedmft=0 initialized in invars1
!  dtsets(idtset)%usedmatpu=0 initialized in invars1
!  dtsets(idtset)%usepawu=0 initialized in invars1
   dtsets(idtset)%usexcnhat=-1
   dtsets(idtset)%useylm=0
!  V
   dtsets(idtset)%vacnum = -1
   dtsets(idtset)%vcutgeo(:)=zero
   dtsets(idtset)%vdw_xc = 0
   dtsets(idtset)%vdw_nwan(:)=0
   dtsets(idtset)%vdw_supercell(:)=0
   dtsets(idtset)%vis=100.0_dp
   dtsets(idtset)%vprtrb(1:2)=zero
!  W
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%wfoptalg=0
   else
     dtsets(idtset)%wfoptalg=10
   end if
   if (dtsets(idtset)%paral_kgb>0.and.idtset>0) then
     dtsets(idtset)%wfoptalg=4          ! compulsory for paral_kgb paralization:overwrite
   end if
   dtsets(idtset)%wtatcon(:,:,:)=zero
   dtsets(idtset)%wtk=one
   dtsets(idtset)%wvl_cpmult  = 10._dp
   dtsets(idtset)%wvl_crmult  = 6._dp
   dtsets(idtset)%wvl_fpmult  = 10._dp
   dtsets(idtset)%wvl_frmult  = 10._dp
   dtsets(idtset)%wvl_hgrid   = 0.5_dp
   dtsets(idtset)%wvl%n(:)    = -1
   dtsets(idtset)%wvl%ni(:)   = -1
   dtsets(idtset)%wvl%ntot    = -1
   dtsets(idtset)%wvl%h(:)    = -1._dp
   dtsets(idtset)%wvl_nprccg  = 10
   dtsets(idtset)%w90iniprj   = 1
   dtsets(idtset)%w90prtunk   = 0

!  X
   dtsets(idtset)%xclevel  = 0
!  Y
!  Z
   dtsets(idtset)%zcut=3.67493260d-03  ! = 0.1eV
   dtsets(idtset)%ziontypat(:)=zero

!  BEGIN VARIABLES FOR @Bethe-Salpeter
   dtsets(idtset)%bs_algorithm    =1
   dtsets(idtset)%bs_haydock_niter=100
   dtsets(idtset)%bs_exchange_term=1
   dtsets(idtset)%bs_coulomb_term=1
   dtsets(idtset)%bs_calctype=1
   dtsets(idtset)%bs_coupling=1

   dtsets(idtset)%bs_haydock_tol=0.01_dp

   dtsets(idtset)%bs_eh_basis_set=(/0,0/)
   dtsets(idtset)%bs_eh_cutoff = (/smallest_real,greatest_real,greatest_real/)
   dtsets(idtset)%bs_freq_mesh=(/zero,zero,0.01_dp/Ha_eV/)
!  END VARIABLES FOR @Bethe-Salpeter.

 end do

 DBG_EXIT("COLL")

end subroutine indefo
!!***
