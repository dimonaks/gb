!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars2
!!
!! NAME
!! invars2
!!
!! FUNCTION
!! Initialize variables for the ABINIT code, for one particular
!! dataset, characterized by jdtset.
!! Note : some parameters have already been read in invars0 and invars1,
!! and were used to dimension the arrays needed here.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprim in the axes
!!               of the conventional bravais lattice (*2 if center/=0)
!! iout=unit number for echoed output
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! mband=maximum number of bands for any k-point
!! msym=default maximal number of symmetries
!! npsp=number of pseudopotentials
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!! string*(*)=character string containing all the input data.
!!  Initialized previously in instrng.
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! zionpsp(npsp)=valence charge of each type of atom (coming from the psp files,
!!   except when invars2 is used in newsp, in which case zion is set to zero)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  dtset=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!! All rest of arguments given alphabetically from acell (length scales)
!! to wtk (k-point weights), used to control running of the main routine.
!! See abinit_help for definitions of these variables.
!! These values become defined by being read from string,
!! that contains all information from the input file,
!! in a compressed, standardized, format
!! At the input, they already contain a default value.
!!
!! NOTES
!!
!! PARENTS
!!      invars2m,m_ab6_invars_f90,newsp
!!
!! CHILDREN
!!      atmdata,chkneu,inkpts,intagm,invacuum,leave_new,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine invars2(bravais,dtset,iout,jdtset,lenstr,&
& mband,msym,npsp,pspheads,string,usepaw,zionpsp)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

 !$use m_gwdefs,      only: gw_sctype_from_name

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_57_iovars, except_this_one => invars2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,jdtset,lenstr,mband,msym,npsp,usepaw
 character(len=*),intent(in) :: string
 type(dataset_type),intent(inout) :: dtset
!arrays
 integer,intent(in) :: bravais(11)
 real(dp),intent(in) :: zionpsp(npsp)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer :: bantot,berryopt,dmatsize,getocc,iatom,ii,ikpt,index,ionmov,intimage
 integer :: iprcch,ipsp,iscf,itypat,jj,kptopt,lpawu,marr,natom,nband1,nberry
 integer :: niatcon,nimage,nkpt,npspalch,nqpt,nspinor,nsppol,nsym,ntypalch,ntypat,ntyppure
 integer :: occopt,occopt_tmp,response,tfband,tnband,tread,tread_alt,tread_key
 integer :: itol, itol_gen, ds_input
 real(dp) :: amu_default,areaxy,charge,fband,kptrlen,nelectjell,rcov
 real(dp) :: rhoavg,sumalch,zelect,zval
 real(dp) :: toldfe_, tolrff_, toldff_, tolwfr_, tolvrs_
 character(len=2) :: symbol
 character(len=30) :: token
 character(len=500) :: message
 character(len=fnlen) :: keyw
!arrays
 integer :: vacuum(3)
 integer,allocatable :: iatcon(:),natcon(:)
 integer,allocatable :: intarr(:) !Dummy arguments for subroutine 'intagm' to parse input file
 real(dp) :: tsec(2)
 real(dp),allocatable :: mass_psp(:)
 real(dp),allocatable :: dprarr(:) !Dummy arguments for subroutine 'intagm' to parse input file

! *************************************************************************

!DEBUG
!write(6,*)' invars2 : enter '
!write(6,*)string(1:lenstr)
!ENDDEBUG

 call timab(191,1,tsec)

!Compute the maximum size of arrays intarr and dprarr
 natom=dtset%natom
 nkpt=dtset%nkpt
 nspinor=dtset%nspinor
 nsppol=dtset%nsppol
 ntypat=dtset%ntypat
 nimage=dtset%nimage
 dmatsize=0
 if (dtset%usepawu>0.and.dtset%usedmatpu/=0) then
   do iatom=1,natom
     lpawu=dtset%lpawu(dtset%typat(iatom))
     if (lpawu/=-1) dmatsize=dmatsize+nsppol*nspinor*(2*lpawu+1)**2
   end do
 end if
 marr=max(3*natom,nkpt*nsppol*mband,2*dtset%nkptgw*nsppol,dmatsize,&
& 3*nkpt,npsp,ntypat,9*msym,60,3*dtset%nconeq*natom,nimage,3*dtset%nqptdm,&
& dtset%natvshift*nsppol*natom)
 allocate(intarr(marr),dprarr(marr))

!----------------------------------------------------------------------------

!****   Read parameters which set remaining array dimensions ****

!Note : some parameters have already been read in invars0 and invars1
!Also, some checking is needed here.

!Read ngfft(1), ngfft(2), and ngfft(3),
!then ngfft(7)=fftalg and ngfft(8)=fftcache.
!Read ngfftdg(1:3)

 token = 'ngfft'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ngfft(1:3)=intarr(1:3)
 token = 'fftalg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%ngfft(7)=intarr(1)
   if (usepaw==1) dtset%ngfftdg(7)=intarr(1)
 end if
 token = 'fftcache'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%ngfft(8)=intarr(1)
   if (usepaw==1) dtset%ngfftdg(8)=intarr(1)
 end if

 token = 'fft_opt_lob'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%fft_opt_lob=intarr(1)

 token = 'ngfftdg'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ngfftdg(1:3)=intarr(1:3)

 token = 'ng'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  The use of the "ng" input variable is forbidden since version 1.8.',ch10,&
&   '  Action : take "ng" out of your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!token = 'npfft'
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
!if(tread==1) dtset%npfft=intarr(1)

!token = 'npband'
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
!if(tread==1) dtset%npband=intarr(1)

!token = 'npkpt'
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
!if(tread==1) dtset%npkpt=intarr(1)

!token = 'bandpp'
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
!if(tread==1) dtset%bandpp=intarr(1)

 token = 'mqgrid'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%mqgrid=intarr(1)
 token = 'mqgriddg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%mqgriddg=intarr(1)

!Make special arrangements to check nband: may be a scalar
!(for occopt=0, 1 or 3, 4, 5, 6, 7) or may be an array (for occopt=2)

 token = 'occopt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%occopt=intarr(1)
 occopt=dtset%occopt

 token = 'gwcalctyp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gwcalctyp=intarr(1)

 token = 'gwcomp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gwcomp=intarr(1)

 token = 'gwencomp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%gwencomp=dprarr(1)

 if (usepaw==1) then
   token = 'gw_sigxcore'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%gw_sigxcore=intarr(1)
 end if

 token = 'gwmem'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gwmem=intarr(1)

 keyw='gw_sctype_name'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),keyw,tread_key,'KEY')
 if (tread_key==1.and..FALSE.) then ! the call below leads to an undefined reference!
!  $dtset%gw_sctype=gw_sctype_from_name(keyw)
 else
   token = 'gw_sctype'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%gw_sctype=intarr(1)
 end if

 token = 'gw_nstep'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gw_nstep=intarr(1)

 token = 'gw_toldfeig'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%gw_toldfeig=dprarr(1)

 token = 'rhoqpmix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%rhoqpmix=dprarr(1)

 token = 'nfreqim'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nfreqim=intarr(1)

 token = 'freqremax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%freqremax=dprarr(1)

 token = 'nfreqre'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nfreqre=intarr(1)
 token = 'nfreqsp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nfreqsp=intarr(1)
 token = 'freqspmax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%freqspmax=dprarr(1)

 token = 'gwrpacorr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gwrpacorr=intarr(1)

 token = 'supercell'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%supercell(1:3)=intarr(1:3)
 token = 'scphon_supercell'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%scphon_supercell(1:3)=intarr(1:3)

!RESPFN integer input variables (needed here to get the value of response
!Presently, rfmeth is not used.
!Warning : rfddk,rfelfd,rfmgfd,rfphon,rfstrs and rfuser are also read in invars1
 token = 'rfasr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfasr=intarr(1)
 token = 'rfatpol'
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfatpol(1:2)=intarr(1:2)
 token = 'rfdir'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfdir(1:3)=intarr(1:3)
 token = 'rfddk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfddk=intarr(1)
 token = 'rfelfd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfelfd=intarr(1)
 token = 'rfmgfd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfmgfd=intarr(1)
 token = 'rfmeth'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfmeth=intarr(1)
 token = 'rfphon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfphon=intarr(1)
 token = 'rfstrs'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfstrs=intarr(1)
 token = 'rfuser'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rfuser=intarr(1)

 response=0
 if(dtset%rfddk/=0 .or. dtset%rfelfd/=0 .or. dtset%rfmgfd/=0 .or. &
& dtset%rfphon/=0 .or. dtset%rfstrs/=0 .or. dtset%rfuser/=0 ) response=1

!NONLINEAR integer input variables (same definition as for rfarr)
!Presently, rf?asr, rf?meth,rf?strs and rf?thrd are not used
 token = 'rf1atpol'
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf1atpol(1:2)=intarr(1:2)
 token = 'rf1dir'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf1dir(1:3)=intarr(1:3)
 token = 'rf1elfd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf1elfd=intarr(1)
 token = 'rf1phon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf1phon=intarr(1)

 token = 'rf2atpol'
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf2atpol(1:2)=intarr(1:2)
 token = 'rf2dir'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf2dir(1:3)=intarr(1:3)
 token = 'rf2elfd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf2elfd=intarr(1)
 token = 'rf2phon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf2phon=intarr(1)

 token = 'rf3atpol'
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf3atpol(1:2)=intarr(1:2)
 token = 'rf3dir'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf3dir(1:3)=intarr(1:3)
 token = 'rf3elfd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf3elfd=intarr(1)
 token = 'rf3phon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rf3phon=intarr(1)

 response=0
 if(dtset%rfddk/=0 .or. dtset%rfelfd/=0 .or. &
& dtset%rfmgfd/=0 .or. dtset%rfphon/=0 .or. &
& dtset%rfstrs/=0 .or. dtset%rfuser/=0 .or. &
& dtset%rf1elfd/=0 .or. dtset%rf1phon/=0 .or. &
& dtset%rf2elfd/=0 .or. dtset%rf2phon/=0 .or. &
& dtset%rf3elfd/=0 .or. dtset%rf3phon/=0 ) response=1

 token = 'prepanl'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prepanl=intarr(1)

 token = 'prepgkk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prepgkk=intarr(1)

!real(dp) input variables
 token = 'boxcutmin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%boxcutmin=dprarr(1)
 token = 'charge'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%charge=dprarr(1)
 token = 'dosdeltae'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%dosdeltae=dprarr(1)
 token = 'dtion'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%dtion=dprarr(1)
 token = 'ecut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%ecut=dprarr(1)
 token = 'esmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%esmear=dprarr(1)
 token = 'scphon_temp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%scphon_temp=dprarr(1)
 token = 'sciss'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%sciss=dprarr(1)
 token = 'tsmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%tsmear=dprarr(1)
 token = 'vis'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%vis=dprarr(1)
 token = 'ecutsm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%ecutsm=dprarr(1)
 token = 'exchmix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%exchmix=dprarr(1)
 token = 'dilatmx'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%dilatmx=dprarr(1)
 token = 'fxcartfactor'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%fxcartfactor=dprarr(1)
 token = 'strfact'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%strfact=dprarr(1)
 token = 'freqsusin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%freqsusin=dprarr(1)
 token = 'freqsuslo'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%freqsuslo=dprarr(1)
 token = 'optfreqsus'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%optfreqsus=intarr(1)
 token = 'effmass'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%effmass=dprarr(1)

 token = 'mditemp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%mditemp=dprarr(1)
 token = 'mdftemp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1)then
   dtset%mdftemp=dprarr(1)
 else
!  Default mdftemp is mditemp
   dtset%mdftemp=dtset%mditemp
 end if
!Recursion input variables

 token = 'tfkinfunc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%tfkinfunc=intarr(1)
 token = 'recgratio'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%recgratio=intarr(1)
 token = 'recefermi'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%recefermi=dprarr(1)
 token = 'recnpath'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%recnpath=intarr(1)
 token = 'recnrec'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%recnrec=intarr(1)
 token = 'recrcut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
 if(tread==1) dtset%recrcut=dprarr(1)
 token = 'recptrott'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%recptrott=intarr(1)
 token = 'rectesteg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%rectesteg=intarr(1)
 token = 'rectolden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%rectolden=dprarr(1)
!Constant NPT Molecular Dynamics Input variables
 token = 'noseinert'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%noseinert=dprarr(1)
 token = 'bmass'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%bmass=dprarr(1)
 token = 'qmass'
 call intagm(dprarr,intarr,jdtset,marr,dtset%nnos,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%qmass(:)=dprarr(1:dtset%nnos)
 token = 'tphysel'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%tphysel=dprarr(1)
 token = 'strtarget'
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%strtarget(1:6)=dprarr(1:6)
 token = 'vcutgeo'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%vcutgeo(1:3)=dprarr(1:3)

 token = 'friction'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%friction=dprarr(1)
 token = 'mdwall'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
 if(tread==1) dtset%mdwall=dprarr(1)
 token = 'fixmom'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%fixmom=dprarr(1)
 token = 'eshift'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%eshift=dprarr(1)
 token = 'boxcenter'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%boxcenter(1:3)=dprarr(1:3)
 token = 'ecuteps'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%ecuteps=dprarr(1)

 token = 'ecutsigx'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) then
   dtset%ecutsigx=dprarr(1)
!else
!  if(dtset%optdriver==RUNL_SIGMA .or. dtset%optdriver==RUNL_SCGW) dtset%ecutsigx=dtset%ecut
 endif
 token = 'ecutwfn'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) then
   dtset%ecutwfn=dprarr(1)
 else
   if(dtset%optdriver==RUNL_SCREENING .or. dtset%optdriver==RUNL_SIGMA &
&       .or. dtset%optdriver==RUNL_SCGW) dtset%ecutwfn=dtset%ecut
 endif
 token = 'omegasimax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%omegasimax=dprarr(1)
 token = 'omegasrdmax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%omegasrdmax=dprarr(1)
 token = 'soenergy'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%soenergy=dprarr(1)
 token = 'spbroad'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%spbroad=dprarr(1)
 token = 'stmbias'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%stmbias=dprarr(1)

 token = 'awtr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%awtr=intarr(1)
 token = 'inclvkb'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%inclvkb=intarr(1)
 token = 'nomegasf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nomegasf=intarr(1)
 token = 'spmeth'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%spmeth=intarr(1)
 token= 'symchi'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%symchi=intarr(1)

 token = 'getscr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getscr=intarr(1)
 token = 'gwgamma'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gwgamma=intarr(1)
 token = 'irdsuscep'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdsuscep=intarr(1)
 token = 'irdscr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdscr=intarr(1)
 token = 'nomegasi'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nomegasi=intarr(1)
 token = 'ppmodel'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ppmodel=intarr(1)
 token = 'symsigma'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%symsigma=intarr(1)

 token = 'fftgw'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%fftgw=intarr(1)
 token = 'getsuscep'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getsuscep=intarr(1)
 token = 'getqps'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getqps=intarr(1)
 token = 'getkss'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getkss=intarr(1)
 token='gwpara'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gwpara=intarr(1)
 token = 'irdqps'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdqps=intarr(1)
 token = 'irdkss'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdkss=intarr(1)
 token = 'ppmfrq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%ppmfrq=dprarr(1)
 token = 'rcut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
 if(tread==1) dtset%rcut=dprarr(1)
 token = 'zcut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%zcut=dprarr(1)

!q-points for long wave-length limit.
 if (dtset%gw_nqlwl>0) then
   token = 'gw_qlwl'
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%gw_nqlwl,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtset%gw_qlwl(1:3,1:dtset%gw_nqlwl) = reshape(dprarr(1:3*dtset%gw_nqlwl),(/3,dtset%gw_nqlwl/))
 end if

!@bethe_salpeter

!integer scalars
 token = 'bs_algorithm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bs_algorithm=intarr(1)

 token = 'bs_haydock_niter'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bs_haydock_niter=intarr(1)

 token = 'bs_exchange_term'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bs_exchange_term=intarr(1)
 
 token = 'bs_calctype'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bs_calctype=intarr(1)

 token = 'bs_coulomb_term'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bs_coulomb_term=intarr(1)

 token = 'bs_coupling'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bs_coupling=intarr(1)

 token = 'irdbseig'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdbseig=intarr(1)

 token = 'getbseig'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getbseig=intarr(1)

!integer arrays
 token = 'bs_eh_basis_set'
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bs_eh_basis_set(1:2)=intarr(1:2)

!real arrays
 token = 'bs_eh_cutoff'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%bs_eh_cutoff(1:3)=dprarr(1:3)

 token = 'bs_freq_mesh'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%bs_freq_mesh(1:3)=dprarr(1:3)

 token = 'bs_haydock_tol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%bs_haydock_tol=dprarr(1)

 token = 'gw_EET'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gw_EET=intarr(1)

 token = 'gw_EET_nband'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gw_EET_nband=intarr(1)

 if(dtset%optdriver==RUNL_RDM) then
   token = 'rdmnb'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%rdmnb=intarr(1)
 end if !dtset%optdriver==7

 token = 'ntypalch'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ntypalch=intarr(1)

 ntypalch=dtset%ntypalch
 if(ntypalch>ntypat)then
   write(message, '(6a,i4,a,i4,a,a)' ) ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  The input variable ntypalch must be smaller than ntypat, while they are',ch10,&
&   '  ntypalch=',dtset%ntypalch,', and ntypat=',ntypat,ch10,&
&   '  Action : check ntypalch vs ntypat in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 ntyppure=ntypat-ntypalch
 dtset%ntyppure=ntyppure
 npspalch=npsp-ntyppure
!DEBUG
!write(6,*)' invars2m : npspalch=',npspalch
!ENDDEBUG
 dtset%npspalch=npspalch
 if(npspalch<0)then
   write(message, '(4a,i4,2a,i4,a,a)' ) ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  The number of available pseudopotentials, npsp=',npsp,ch10,&
&   '  is smaller than the requested number of types of pure atoms, ntyppure=',ntyppure,ch10,&
&   '  Action : check ntypalch versus ntypat and npsp in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(ntypalch>0)then
   token = 'algalch'
   call intagm(dprarr,intarr,jdtset,marr,ntypalch,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%algalch(1:ntypalch)=intarr(1:ntypalch)
   token = 'mixalch'
   call intagm(dprarr,intarr,jdtset,marr,npspalch*ntypalch,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtset%mixalch(1:npspalch,1:ntypalch)=&
&   reshape(dprarr(1:npspalch*ntypalch),(/npspalch,ntypalch/))
   do itypat=1,ntypalch
     sumalch=sum(dtset%mixalch(:,itypat))
     if(abs(sumalch-one)>tol10)then
       write(message, '(4a,i4,2a,f8.2,4a)' ) ch10,&
&       ' invars2: ERROR -',ch10,&
&       '  For the alchemical atom number',itypat,ch10,&
&       '  the sum of the pseudopotential coefficients is',sumalch,ch10,&
&       '  while it should be one.',ch10,&
&       '  Action : check the content of the input variable mixalch.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if
   end do
 end if

!Compute ziontypat
!When the pseudo-atom is pure, simple copy
 if(ntyppure>0)then
   do itypat=1,ntyppure
     dtset%ziontypat(itypat)=zionpsp(itypat)
   end do
 end if
!When the pseudo-atom is alchemical, must make mixing
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     dtset%ziontypat(itypat)=zero
     do ipsp=ntyppure+1,npsp
       dtset%ziontypat(itypat)=dtset%ziontypat(itypat) &
&       +dtset%mixalch(ipsp-ntyppure,itypat-ntyppure)*zionpsp(ipsp)
     end do
   end do
 end if

!DEBUG
!write(6,*)' invars2 : dtset%ziontypat=',dtset%ziontypat(1:ntypat)
!write(6,*)' invars2 : zionpsp=',zionpsp(1:npsp)
!write(6,*)' invars2 : ntyppure,npspalch,ntypalch=',ntyppure,npspalch,ntypalch
!write(6,*)' invars2 : mixalch=',dtset%mixalch(1:npspalch,1:ntypalch)
!ENDDEBUG

 charge=dtset%charge

 if (occopt==0 .or. occopt==1 .or. (occopt>=3 .and. occopt<=7) ) then

   token = 'nband'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tnband,'INT')
   if(tnband==1) then
     nband1=intarr(1)
   else
!    Default value in the metallic case, or in the insulating case
     fband=0.5_dp
     if(occopt==1)fband=0.125_dp
     if (dtset%usewvl == 1) fband = zero
     token = 'fband'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tfband,'DPR')
     if(tfband==1)then
       fband=dprarr(1)
       write(message, '(a,es16.8,a)' )&
&       ' invars2: read the value of fband=',fband,' from input file.'
     else
       write(message, '(a,es16.8)' )&
&       ' invars2: take the default value of fband=',fband
     end if
     call wrtout(std_out,message,'COLL')
!    First compute the total valence charge
     zval=0.0_dp
     do iatom=1,natom
       zval=zval+dtset%ziontypat(dtset%typat(iatom))
     end do
     zelect=zval-charge
!    Then select the minimum number of bands, and add the required number
!    Note that this number might be smaller than the one computed
!    by a slightly different formula in invars1
     nband1=dtset%nspinor*&
&     ((ceiling(zelect-1.0d-10)+1)/2 + ceiling( fband*natom - 1.0d-10 ))
   end if

!  Set nband to same input number for each k point and spin
!  where nband1 is the eventual input, computed value, or default
   do ikpt=1,nkpt*nsppol
     dtset%nband(ikpt)=nband1
   end do

 else if (occopt==2) then
!  Give nband explicitly for each k point and spin

   token = 'nband'
   call intagm(dprarr,intarr,jdtset,nkpt*nsppol,nkpt*nsppol,string(1:lenstr),&
&   token,tnband,'INT')
   if(tnband==1) dtset%nband(1:nkpt*nsppol)=intarr(1:nkpt*nsppol)

 else
   write(message, '(a,a,a,a,i8,a,a,a)' ) ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  occopt=',occopt,' not allowed.',ch10,&
&   '  Action : correct your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!DEBUG
!write(6,*)' invars2: after read nband '
!write(6,*)dtset%nband(1:nkpt*nsppol)
!stop
!ENDDEBUG

!----------------------------------------------------------------------------

!****   Read other parameters  ****

!All checking should be done in chkinp.f


 token = 'chksymbreak'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%chksymbreak=intarr(1)

!Get array
 token = 'getocc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getocc=intarr(1)
 getocc=dtset%getocc
 token = 'getwfk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getwfk=intarr(1)
 token = 'getxcart'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getxcart=intarr(1)
 token = 'getxred'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getxred=intarr(1)
 token = 'getden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getden=intarr(1)
 token = 'getcell'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getcell=intarr(1)
 token = 'getwfq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getwfq=intarr(1)
 token = 'get1wf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%get1wf=intarr(1)
 token = 'getddk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getddk=intarr(1)
 token = 'getvel'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%getvel=intarr(1)
 token = 'get1den'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%get1den=intarr(1)

 token = 'getwf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  The use of the "getwf" variable is forbidden since version 2.0.',ch10,&
&   '  Action : replace "getwf" by "getwfk" in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 token = 'accesswff'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%accesswff=intarr(1)

 token = 'dynimage'
 call intagm(dprarr,intarr,jdtset,marr,nimage,string(1:lenstr),token,tread,'INT')
 if(tread==1)then
   dtset%dynimage(1:nimage)=intarr(1:nimage)
 end if
 dtset%ndynimage=0
 do ii=1,dtset%nimage
   if(dtset%dynimage(ii)/=0)dtset%ndynimage=dtset%ndynimage+1
 end do

 token = 'enunit'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%enunit=intarr(1)
 token = 'exchn2n3d'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%exchn2n3d=intarr(1)
 token = 'iboxcut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%iboxcut=intarr(1)
 token = 'icutcoul'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%icutcoul=intarr(1)
 token = 'imgmov'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%imgmov=intarr(1)
 token = 'ionmov'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ionmov=intarr(1)
 token = 'intxc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%intxc=intarr(1)
 token = 'iextrapwf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%iextrapwf=intarr(1)
 token = 'iprcch'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%iprcch=intarr(1)
 else
   if (dtset%iextrapwf==1) dtset%iprcch=6
 end if
 token = 'iprcel'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%iprcel=intarr(1)
 token = 'iprctfvw'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%iprctfvw=intarr(1)
 token = 'iprcfc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%iprcfc=intarr(1)
 token = 'irdden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdden=intarr(1)
 token = 'irdwfk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdwfk=intarr(1)
 token = 'iscf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%iscf=intarr(1)
 else if (dtset%optdriver==1.and.dtset%iscf>=10) then
   dtset%iscf=dtset%iscf-10
 end if
 token = 'isecur'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%isecur=intarr(1)

!Reading ixc must be immediately followed by reading xcname
 token = 'ixc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ixc=intarr(1)
 keyw='xcname'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),keyw,tread_key,'KEY')
 if(tread_key==1)then
   if(tread==1)then
     write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&     ' invars1: ERROR -',ch10,&
&     '  ixc and xcname cannot be specified simultaneously',ch10,&
&     '  for the same dataset.',ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   else
!    Note that xcname is a 'key' variable : its value is stored in keyw at output of intagm
     if(trim(keyw)=='PW92')dtset%ixc=7
     tread=1
   end if
 end if
 dtset%xclevel=0
 if( ( 1<=dtset%ixc .and. dtset%ixc<=10).or.(30<=dtset%ixc .and. dtset%ixc<=39) )dtset%xclevel=1 ! LDA
 if( (11<=dtset%ixc .and. dtset%ixc<=19).or.(23<=dtset%ixc .and. dtset%ixc<=29) )dtset%xclevel=2 ! GGA
 if( 20<=dtset%ixc .and. dtset%ixc<=22 )dtset%xclevel=3 ! ixc for TDDFT kernel tests

 token = 'vdw_xc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%vdw_xc=intarr(1)

 token = 'positron'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%positron=intarr(1)
 token = 'ixcpositron'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ixcpositron=intarr(1)
 token = 'posnstep'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%posnstep=intarr(1)
 token = 'posocc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%posocc=dprarr(1)

 token = 'frzfermi'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%frzfermi=intarr(1)
 token = 'nqpt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nqpt=intarr(1)
 token = 'ieig2rf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ieig2rf=intarr(1)
 token = 'smdelta'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%smdelta=intarr(1)
 token = 'bdeigrf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bdeigrf=intarr(1)
 token = 'restartxf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%restartxf=intarr(1)
 token = 'optcell'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%optcell=intarr(1)
 token = 'irdwfq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdwfq=intarr(1)
 token = 'ird1wf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ird1wf=intarr(1)
 token = 'irdddk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%irdddk=intarr(1)
 token = 'kptopt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%kptopt=intarr(1)
 token = 'chkexit'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%chkexit=intarr(1)
 token = 'ikhxc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ikhxc=intarr(1)
 token = 'nbdbuf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1)then
   dtset%nbdbuf=intarr(1)
 else
   if(response/=1 .and. dtset%iscf<0)dtset%nbdbuf=2*dtset%nspinor
   if(response==1 .and. 3<=occopt .and. occopt<=7 )dtset%nbdbuf=2*dtset%nspinor
 end if
 token = 'localrdwf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%localrdwf=intarr(1)
 token = 'optforces'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%optforces=intarr(1)
 token = 'optstress'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%optstress=intarr(1)
 token = 'optnlxccc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%optnlxccc=intarr(1)
 token = 'nberry'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nberry=intarr(1)
 token = 'bdberry'
 call intagm(dprarr,intarr,jdtset,marr,2*nsppol,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%bdberry(1)=intarr(1); dtset%bdberry(2)=intarr(2)
   if(nsppol==2)then
     dtset%bdberry(3)=intarr(3); dtset%bdberry(4)=intarr(4)
   end if
 end if
 token = 'delayperm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%delayperm=intarr(1)
 token = 'signperm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%signperm=intarr(1)

 token = 'nbandkss'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nbandkss=intarr(1)

 token = 'npwkss'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npwkss=intarr(1)

 token = 'wfoptalg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%wfoptalg=intarr(1)
 token = 'nbdblock'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nbdblock=intarr(1)
 token = 'kssform'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%kssform=intarr(1)
 token = 'td_mexcit'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%td_mexcit=intarr(1)
 token = 'npweps'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npweps=intarr(1)
 token = 'npulayit'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npulayit=intarr(1)
 token = 'diismemory'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%diismemory=intarr(1)
 token = 'nwfshist'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nwfshist=intarr(1)

 token = 'npwsigx'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npwsigx=intarr(1)

 token = 'npwwfn'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npwwfn=intarr(1)

 token = 'nscforder'
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), token, tread, 'INT')
 if (tread == 1) dtset%nscforder = intarr(1)

 token = 'nsheps'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nsheps=intarr(1)

 token = 'nshsigx'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nshsigx=intarr(1)

 token = 'nshwfn'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nshwfn=intarr(1)

 token = 'nomegasrd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nomegasrd=intarr(1)

 token = 'corecs'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'DPR')
 if(tread==1)then
   dtset%corecs(1:ntypat)=dprarr(1:ntypat)
 end if

 token = 'corecs'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'DPR')
 if(tread==1)then
   dtset%corecs(1:ntypat)=dprarr(1:ntypat)
 end if

 token = 'ptcharge'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'DPR')
 if(tread==1)then
   dtset%ptcharge(1:ntypat)=dprarr(1:ntypat)
 end if

 token = 'quadmom'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'DPR')
 if(tread==1)then
   dtset%quadmom(1:ntypat)=dprarr(1:ntypat)
 end if

 token = 'pawcpxocc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%pawcpxocc=intarr(1)
 else if (dtset%nspinor==2.and.(dtset%usepawu>=1.or.dtset%usedmft>0)) then
   dtset%pawcpxocc=2
 else if (dtset%pawspnorb>0.and.(dtset%kptopt<=0.or.dtset%kptopt>=3)) then
   if (dtset%optdriver/=0.or.dtset%ionmov<6.or.dtset%iscf<10) dtset%pawcpxocc=2
 end if
 token = 'pawecutdg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) then
   dtset%pawecutdg=dprarr(1)
!  else                           ! MT 2006, nov 28th: this defaut was "dangerous"
!  dtset%pawecutdg=two*dtset%ecut
 end if
 token = 'pawlcutd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawlcutd=intarr(1)
 token = 'pawlmix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawlmix=intarr(1)
 token = 'pawmixdg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%pawmixdg=intarr(1)
 else if (dtset%npfft>1.and.usepaw==1) then
   dtset%pawmixdg=1
 end if
 token = 'pawnhatxc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawnhatxc=intarr(1)
 token = 'pawntheta'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawntheta=intarr(1)
 token = 'pawnphi'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawnphi=intarr(1)
 token = 'pawnzlm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawnzlm=intarr(1)
 token = 'pawoptmix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%pawoptmix=intarr(1)
 else
   if (usepaw==1.and.dtset%iscf<10.and.(dtset%positron==1.or.dtset%positron<0)) dtset%pawoptmix=1
 end if
 token = 'pawovlp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%pawovlp=dprarr(1)
 token = 'pawprtden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%pawprtden=intarr(1)
   write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&   ' invars2 : COMMENT -',ch10,&
&   '  The pawprtden variable is deprecated and will soon be removed,',ch10,&
&   '  use prtden instead next time.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if
 token = 'pawprtdos'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawprtdos=intarr(1)
 token = 'pawprtvol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawprtvol=intarr(1)
 token = 'pawprtwf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawprtwf=intarr(1)
 token = 'pawstgylm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawstgylm=intarr(1)
 token = 'pawusecp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawusecp=intarr(1)
 token = 'pawxcdev'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawxcdev=intarr(1)
 token = 'spnorbscl'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%spnorbscl=dprarr(1)

 if (dtset%usedmft>0) then
   token = 'dmft_iter'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmft_iter=intarr(1)
   token = 'dmft_mxsf'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtset%dmft_mxsf=dprarr(1)
   token = 'dmft_nwli'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmft_nwli=intarr(1)
   token = 'dmft_nwlo'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmft_nwlo=intarr(1)
   token = 'dmft_rslf'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmft_rslf=intarr(1)
   token = 'dmft_solv'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmft_solv=intarr(1)
   token = 'dmftbandi'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmftbandi=intarr(1)
   token = 'dmftbandf'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmftbandf=intarr(1)
   token = 'dmftcheck'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmftcheck=intarr(1)
 end if

 if (dtset%usepawu>0.or.dtset%usedmft>0) then
   token = 'upawu'
   call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'ENE')
   if(tread==1)then
     dtset%upawu(1:ntypat)=dprarr(1:ntypat)
   end if
   token = 'jpawu'
   call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'ENE')
   if(tread==1)then
     dtset%jpawu(1:ntypat)=dprarr(1:ntypat)
   end if
   token = 'dmatpuopt'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmatpuopt=intarr(1)
   token = 'dmatudiag'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%dmatudiag=intarr(1)

   if (dtset%usedmatpu/=0) then
     token = 'dmatpawu'
     if (dmatsize>0) then
       call intagm(dprarr,intarr,jdtset,marr,dmatsize,string(1:lenstr),token,tread,'DPR')
       if(tread==1) then
         ii=1;jj=0
         do iatom=1,natom
           lpawu=dtset%lpawu(dtset%typat(iatom))
           if (lpawu/=-1) then
             dtset%dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:nsppol*nspinor,ii)= &
&             reshape(dprarr(jj+1:jj+nsppol*nspinor*(2*lpawu+1)**2),(/2*lpawu+1,2*lpawu+1,nsppol*nspinor/))
             ii=ii+1;jj=jj+nsppol*nspinor*(2*lpawu+1)**2
           end if
         end do
       else
         write(message, '(6a)' )ch10,&
&         ' invars2: ERROR -',ch10,&
&         '  When LDA/GGA+U is activated and usedmatpu/=0, dmatpawu MUST be defined.',ch10,&
&         '  Action : add dmatpawu keyword in input file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end if
   end if
 end if

 if (dtset%macro_uj>0) then

   token = 'pawujat'
   dtset%pawujat=minval(pack((/ (ii,ii=1,dtset%natom) /) ,dtset%lpawu(dtset%typat(:))>0))
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%pawujat=intarr(1)

   token = 'pawujrad'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
   if(tread==1) dtset%pawujrad=dprarr(1)

   token = 'pawujv'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
   if(tread==1) dtset%pawujv=dprarr(1)

 end if

 token = 'pawfatbnd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawfatbnd=intarr(1)

 token = 'bxctmindg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) then
   dtset%bxctmindg=dprarr(1)
 else
   dtset%bxctmindg=dtset%boxcutmin
 end if

 token = 'usexcnhat'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%usexcnhat=intarr(1)

 token = 'useylm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%useylm=intarr(1)
   if ((usepaw==1).and.(dtset%useylm==0)) then
     write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  Pseudopotential file is PAW format (pspcod=7) while',ch10,&
&     '  input variable "useylm" has the incompatible value 0 !',ch10,&
&     '  Action : change psp format or "useylm" value in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if ((dtset%tfkinfunc==2).and.(dtset%useylm==0)) then
     write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  You are using recursion method (tfkinfunc=2)  while',ch10,&
&     '  input variable "useylm" has the incompatible value 0 !',ch10,&
&     '  Action : change  "useylm" value in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
 if (usepaw==1) dtset%useylm=1
 if(dtset%tfkinfunc==2) then
   dtset%useylm=1
   dtset%userec=1
 end if
 token = 'ireadc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  The use of the "ireadc" variable is forbidden since version 2.0.',ch10,&
&   '  Action : replace "ireadc" by "irdwfk" in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 ionmov=dtset%ionmov ; iprcch=dtset%iprcch ; iscf=dtset%iscf ; nqpt=dtset%nqpt
 kptopt=dtset%kptopt; nberry=dtset%nberry ; berryopt=dtset%berryopt

!Dielectric real(dp) input variables
!Reading of diemix/diemixmag must be inserted after iprcel
 token = 'diecut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%diecut=dprarr(1)
!Special treatment if iscf==-1
 if(iscf==-1) dtset%diecut=four*dtset%ecut
 token = 'dielng'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
 if(tread==1) dtset%dielng=dprarr(1)
 token = 'diemac'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%diemac=dprarr(1)

 token = 'diemix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) then
   dtset%diemix=dprarr(1)
 else
   if (mod(dtset%iprcel,100)>19) dtset%diemix=one
 end if

 token = 'diemixmag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) then
   dtset%diemixmag=dprarr(1)
 else
   if (dtset%iscf<10.or.dtset%nspden==1.or.dtset%iprcel==0.or.(dtset%iprcel>70.and.dtset%iprcel<80)) then
     dtset%diemixmag=dtset%diemix
   else
     dtset%diemixmag=-dtset%diemix
   end if
 end if

 token = 'diegap'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%diegap=dprarr(1)
 token = 'dielam'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%dielam=dprarr(1)

 token = 'td_maxene'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%td_maxene=dprarr(1)

!ACFD input variables
 token = 'idyson'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%idyson=intarr(1)
 token = 'ndyson'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ndyson=intarr(1)
 token = 'intexact'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%intexact=intarr(1)
 token = 'nbandsus'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nbandsus=intarr(1)
 token = 'ldgapp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ldgapp=intarr(1)
 token = 'suskxcrs'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%suskxcrs=intarr(1)

 if((iscf==5.or.iscf==6) .and. ionmov==4 .and. iprcch/=3 )then
   iprcch=3
   dtset%iprcch=iprcch
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' invars2: COMMENT -',ch10,&
&   '  When ionmov==4 and iscf==5 or 6, iprcch must be 3.',ch10,&
&   '  Set iprcch to 3.'
   call wrtout(std_out,message,'COLL')
 end if


 token = 'mffmem'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%mffmem=intarr(1)

!Set default values of 0 for occupation numbers and k pt wts
 bantot=0

!nkpt and nband must be defined to execute following loop
 if ( tnband == 1 ) then
   do ikpt=1,nkpt*nsppol
     do ii=1,dtset%nband(ikpt)
       bantot=bantot+1
       dtset%occ_orig(bantot)=0.0_dp
     end do
   end do
 end if

!NON DEVLOP CASE
 token = 'nloalg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%nloalg(1)=mod(intarr(1),10)
   dtset%nloalg(5)=intarr(1)/10
 end if
!DEVLOP CASE. Note that the fifth component is not used now.
!token = 'nloalg'
!call intagm(dprarr,intarr,jdtset,marr,5,string(1:lenstr),token,tread,'INT')
!if(tread==1) then
!dtset%nloalg(1)=intarr(1) ; dtset%nloalg(2)=intarr(2)
!dtset%nloalg(3)=intarr(3) ; dtset%nloalg(4)=intarr(4) ; dtset%nloalg(4)=intarr(5)
!end if
!ENDDEVLOP

!LOOP variables
 token = 'nline'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nline=intarr(1)
 token = 'nnsclo'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nnsclo=intarr(1)
 token = 'nstep'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nstep=intarr(1)
 token = 'ntime'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ntime=intarr(1)
 token = 'ntimimage'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%ntimimage=intarr(1)
 token = 'nfreqsus'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nfreqsus=intarr(1)

 token = 'nctime'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nctime=intarr(1)

 token = 'ortalg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%ortalg=intarr(1)
 else if (dtset%wfoptalg>=10 .and. dtset%ortalg>0) then
   dtset%ortalg=-dtset%ortalg
 end if

 token = 'pitransform'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pitransform=intarr(1)

!Print variables
 token = 'prtcs'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtcs=intarr(1)
 token = 'prtefg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtefg=intarr(1)
 token = 'prteig'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prteig=intarr(1)
 token = 'prtfc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtfc=intarr(1)
 token = 'prtvol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtvol=intarr(1)
 token = 'prtden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   if (intarr(1)<=0) then
     dtset%prtden=intarr(1)
   else
     dtset%prtden=1
     dtset%pawprtden=intarr(1)-1
   end if
 end if
 token = 'prtdipole'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtdipole=intarr(1)



 token = 'prtdos'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtdos=intarr(1)
 token = 'prtdosm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtdosm=intarr(1)
 token = 'prtelf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtelf=intarr(1)
 token = 'prtfsurf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtfsurf=intarr(1)
 token = 'prtgden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtgden=intarr(1)
 token = 'prtgeo'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtgeo=intarr(1)
 token = 'prtkden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtkden=intarr(1)
 token = 'prtlden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtlden=intarr(1)
 token = 'prtcif'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtcif=intarr(1)
 token = 'prtcml'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtcml=intarr(1)
 token = 'prtnabla'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtnabla=intarr(1)
 token = 'prtnest'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtnest=intarr(1)
 token = 'prtpot'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtpot=intarr(1)
 token = 'prtspcur'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtspcur=intarr(1)
!if(dtset%prtspcur /= 0) stop "Error: spincurrent is not implemented yet"
 token = 'prtstm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtstm=intarr(1)
 token = 'prt1dm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prt1dm=intarr(1)
 if(dtset%nspden==4 .and. dtset%prt1dm==1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  prt1dm is not implemented for nspden 4 yet', ch10,&
&   '  Action : change  "prt1dm=0" value in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 token = 'prtvha'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtvha=intarr(1)
 token = 'prtvhxc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtvhxc=intarr(1)
 token = 'prtvxc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtvxc=intarr(1)
 token = 'prtwant'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtwant=intarr(1)
 token = 'prtwf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%prtwf=intarr(1)
   if ((dtset%tfkinfunc==2).and.(dtset%prtwf==1)) then
     write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  You are using recursion method (tfkinfunc=2)  while',ch10,&
&     '  input variable "prtwf" has the incompatible value 1 !',ch10,&
&     '  Action : change  "prtwf=0" value in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
 if (dtset%tfkinfunc==2) dtset%prtwf=0
 token = 'prtbbb'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtbbb=intarr(1)
 token = 'prtgkk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtgkk=intarr(1)
 token = 'prtxml'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtxml=intarr(1)
 token = 'prtxangst'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtxangst=intarr(1)
 token = 'prtxcart'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtxcart=intarr(1)
 token = 'prtxred'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%prtxred=intarr(1)

 token = 'qprtrb'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%qprtrb(:)=intarr(1:3)

 token = 'strprecon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%strprecon=dprarr(1)

!WVL - Wavelets related values
 token = 'wvl_hgrid'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%wvl_hgrid=dprarr(1)
 dtset%wvl%h(:) = dtset%wvl_hgrid
 token = 'wvl_crmult'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%wvl_crmult=dprarr(1)
 token = 'wvl_frmult'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%wvl_frmult=dprarr(1)
 token = 'wvl_cpmult'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%wvl_cpmult=dprarr(1)
 token = 'wvl_fpmult'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%wvl_fpmult=dprarr(1)
 token = 'wvl_nprccg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%wvl_nprccg=intarr(1)
 token = 'tl_nprccg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%tl_nprccg=intarr(1)
 token = 'tl_radius'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%tl_radius=dprarr(1)

!Wannier90 interface related variables
!w90iniprj
 token = 'w90iniprj'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%w90iniprj=intarr(1)
   if ( usepaw == 0 .and. ( dtset%w90iniprj /= 2  .and. dtset%w90iniprj/=0 .and. dtset%w90iniprj /= 1 )) then
     write(message, '(4a,i4,2a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  w90iniprj should be set to 0, 1 or 2, however, it was ',&
&     dtset%w90iniprj,ch10,&
&     '  Action : check the values of w90iniprj.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if ( usepaw == 1 .and. ( dtset%w90iniprj < 2 .or. dtset%w90iniprj>6 ) &
&   .and. ( dtset%w90iniprj /= 1 .and. dtset%w90iniprj/=0 )) then
     write(message, '(4a,i4,2a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  In the PAW case, the only valid values for w90iniprj are 0, 1, 2, 5 and 6 however, it was ',&
&     dtset%w90iniprj,ch10,&
&     '  Action : check the values of w90iniprj.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
!w90prtunk
 token = 'w90prtunk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%w90prtunk=intarr(1)
   if ( dtset%w90prtunk < 0 ) then
     write(message, '(7a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  w90prtunk should be greater or equal to zero, however, it was ',&
&     dtset%w90prtunk,ch10, '  Action : check the values of w90prtunk.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

!Wannier90 - GW quasiparticle interface
 if(dtset%prtwant==3) then
   token = 'getqps'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%getqps=intarr(1)
   token = 'irdqps'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%irdqps=intarr(1)
 end if

!van der Waals with Wannier functions (Silvestrelli's approach) 
 if(dtset%vdw_xc==10) then
   token = 'vdw_nwan'
   call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),token,tread,'INT')
   if(tread==1)then
     dtset%vdw_nwan(1)=intarr(1)
     dtset%vdw_nwan(2)=intarr(2)
     if ( dtset%vdw_nwan(2) < 0 ) then
       write(message, '(6a)' )ch10,&                                         
&       ' invars2: ERROR -',ch10,&                                           
&       ' second component of vdw_nwan should be greater or equal to zero',& 
&       ch10, '  Action : check the components of vdw_nwan.'    
       call wrtout(std_out,message,'COLL')                                   
       call leave_new('COLL')                                                
     end if                                                                  
   end if
   token = 'vdw_supercell'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
   if(tread==1)then
     dtset%vdw_supercell(1)=intarr(1)
     dtset%vdw_supercell(2)=intarr(2)
     dtset%vdw_supercell(3)=intarr(3)
     if (any (dtset%vdw_supercell(:) < 0) ) then
       write(message, '(6a)' )ch10,&                                           
&       ' invars2: ERROR -',ch10,&                                             
&       ' all components of vdw_supercell should be greater or equal to zero',&     
&       ch10, '  Action : check the components of vdw_supercell.'                   
       call wrtout(std_out,message,'COLL')                                     
       call leave_new('COLL')                                                  
     end if                                                                   
   end if
 end if

!Tolerance variables
 token = 'tolimg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%tolimg=dprarr(1)
 token = 'tolmxf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%tolmxf=dprarr(1)

!find which tolXXX are defined generically and for this jdtset
 tolwfr_=zero
 toldfe_=zero
 toldff_=zero
 tolrff_=zero
 tolvrs_=zero
 itol=0
 itol_gen=0
 token = 'tolwfr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolwfr_=dprarr(1)
     itol_gen=itol_gen+1
   else
     dtset%tolwfr=dprarr(1)
     itol=itol+1
   end if
 end if
 token = 'toldff'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldff_=dprarr(1)
     itol_gen=itol_gen+1
   else
     dtset%toldff=dprarr(1)
     itol=itol+1
   end if
 end if
 token = 'tolrff'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolrff_=dprarr(1)
     itol_gen=itol_gen+1
   else
     dtset%tolrff=dprarr(1)
     itol=itol+1
   end if
 end if
 token = 'toldfe'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldfe_=dprarr(1)
     itol_gen=itol_gen+1
   else
     dtset%toldfe=dprarr(1)
     itol=itol+1
   end if
 end if
 token = 'tolvrs'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     tolvrs_=dprarr(1)
     itol_gen=itol_gen+1
   else
     dtset%tolvrs=dprarr(1)
     itol=itol+1
   end if
 end if

!check for multiple definitions of tolXXX for the present dataset
 if (itol > 1 .or. itol_gen > 1) then
   write(message, '(6a)' )ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  Only one of the tolXXX variables may be defined at once.',ch10,&
&   '  Action : check values of tolvrs, toldfe, tolrff, tolwfr, and toldff.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!if no value is given for jdtset, use defaults
 if (itol == 0 .and. itol_gen == 1) then
   dtset%tolwfr=tolwfr_
   dtset%toldfe=toldfe_
   dtset%toldff=toldff_
   dtset%tolrff=tolrff_
   dtset%tolvrs=tolvrs_
 end if

!Tolerance variables for electrons-positron tc-dft
 toldfe_=zero;toldff_=zero
 itol=0;itol_gen=0
 token = 'postoldff'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldff_=dprarr(1)
     itol_gen=itol_gen+1
   else
     dtset%postoldff=dprarr(1)
     itol=itol+1
   end if
   dtset%postoldfe=zero
 end if
 token = 'postoldfe'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'ENE',ds_input)
 if(tread==1) then
   if (ds_input == 0) then
     toldfe_=dprarr(1)
     itol_gen=itol_gen+1
   else
     dtset%postoldfe=dprarr(1)
     itol=itol+1
   end if
 end if
 if (itol > 1.or.itol_gen >1) then
   write(message, '(6a)' )ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  Only one of the postolXXX variables may be defined at once.',ch10,&
&   '  Action : check values of posoldfe and postoldff.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (itol==0.and.itol_gen==1) then
   dtset%postoldfe=toldfe_
   dtset%postoldff=toldff_
 end if

 token = 'usekden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%usekden=intarr(1)


 token = 'vprtrb'
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),token,tread,'ENE')
 if(tread==1) dtset%vprtrb(:)=dprarr(1:2)

 if((dtset%pawfatbnd>0.or.dtset%prtdos==3) .and. dtset%natsph>0)then
   token = 'iatsph'
   call intagm(dprarr,intarr,jdtset,marr,dtset%natsph,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%iatsph(1:dtset%natsph)=intarr(1:dtset%natsph)
 end if

 token = 'prtdensph'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%prtdensph=intarr(1)
 else if (dtset%nsppol==1.and.dtset%nspden==2) then
   dtset%prtdensph=1
 end if

 token = 'ratsph'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'LEN')
 if(tread==1.and.(usepaw==0.or.dtset%prtdos/=3.or.dtset%pawprtdos/=2)) then
   dtset%ratsph(1:ntypat)=dprarr(1:ntypat)
 else if ((dtset%prtdensph==1.or.dtset%prtdos==3).and.usepaw==1) then
   do itypat=1,ntypat
     dtset%ratsph(itypat)=pspheads(itypat)%pawheader%rpaw
   end do
 end if

!Initialize atvshift
 if(dtset%natvshift>0)then
   token = 'atvshift'
   call intagm(dprarr,intarr,jdtset,marr,dtset%natvshift*nsppol*dtset%natpawu,string(1:lenstr),&
&   token,tread,'ENE')
   if(tread==1) dtset%atvshift(1:dtset%natvshift,1:nsppol,1:dtset%natpawu)=&
&   reshape(dprarr(1:dtset%natvshift*nsppol*dtset%natpawu),(/dtset%natvshift,nsppol,dtset%natpawu/))
 end if

!Initialize wtatcon
 if(dtset%nconeq>0)then

!  Read and check natcon
   allocate(natcon(dtset%nconeq))
   token = 'natcon'
   call intagm(dprarr,intarr,jdtset,marr,dtset%nconeq,string(1:lenstr),token,tread,'INT')
   if(tread==1)then
     natcon(:)=intarr(1:dtset%nconeq)
   else
     write(message, '(6a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  When nconeq is positive, natcon MUST be defined.',ch10,&
&     '  Action : check the values of nconeq and natcon.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   do ii=1,dtset%nconeq
     if(natcon(ii)<0)then
       write(message, '(a,a,a,a,a,a,i4,a,i4,a,a)' )ch10,&
&       ' invars2: ERROR -',ch10,&
&       '  All the components of natcon must be greater than 0.',ch10,&
&       '  The component',ii,' is equal to ',natcon(ii),ch10,&
&       '  Action : check the values of natcon.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end do
   niatcon=sum(natcon(:))

!  Read and check iatcon
   allocate(iatcon(niatcon))
   token = 'iatcon'
   call intagm(dprarr,intarr,jdtset,marr,niatcon,string(1:lenstr),token,tread,'INT')
   if(tread==1)then
     iatcon(:)=intarr(1:niatcon)
   else
     write(message, '(6a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  When nconeq is positive, natcon MUST be defined.',ch10,&
&     '  Action : check the values of nconeq and natcon.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   do ii=1,niatcon
     if(iatcon(ii)<0)then
       write(message, '(a,a,a,a,a,a,i4,a,i4,a,a)' )ch10,&
&       ' invars2: ERROR -',ch10,&
&       '  All the components of iatcon must be greater than 0.',ch10,&
&       '  The component',ii,' is equal to ',iatcon(ii),ch10,&
&       '  Action : check the values of iatcon.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end do

!  Read wtatcon, and unfold it.
   token = 'wtatcon'
   call intagm(dprarr,intarr,jdtset,marr,3*niatcon,string(1:lenstr),token,tread,'DPR')
   if(tread/=1)then
     write(message, '(6a)' )ch10,&
&     ' invars2: ERROR -',ch10,&
&     '  When nconeq is positive, wtatcon MUST be defined.',ch10,&
&     '  Action : check the values of nconeq and wtatcon.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   index=0
   do ii=1,dtset%nconeq
     do jj=1,natcon(ii)
       dtset%wtatcon(1:3,iatcon(jj+index),ii)=dprarr(1+3*(jj+index-1):3+3*(jj+index-1))
     end do
     index=index+natcon(ii)
   end do

   deallocate(iatcon,natcon)

 end if

!Initialize the list of k and q points, as well as wtk and istwfk
 intimage=1 ; if(nimage>2)intimage=(1+nimage)/2
 call invacuum(jdtset,lenstr,natom,dtset%rprimd_orig(1:3,1:3,intimage),string,vacuum,&
& dtset%xred_orig(1:3,1:natom,intimage))

!DEBUG
!write(6,*)' invars2: before inkpts, vacuum=',vacuum(:)
!stop
!ENDDEBUG

!In case of a Berryphase calculation, put response = 1
!in order to set istwfk = 1 at all k-points

 if (abs(dtset%berryopt) > 0) response = 1

 nsym=dtset%nsym
 call inkpts(bravais,dtset%chksymbreak,iout,iscf,dtset%istwfk(1:nkpt),jdtset,&
& dtset%kpt(:,1:nkpt),kptopt,dtset%kptnrm,dtset%kptrlatt,kptrlen,&
& lenstr,nsym,nkpt,nqpt,dtset%nshiftk,&
& nsym,occopt,dtset%qpt,dtset%qptnrm,response,&
& dtset%rprimd_orig(1:3,1:3,intimage),&
& dtset%shiftk,string,dtset%symafm(1:nsym),&
& dtset%symrel(:,:,1:nsym),vacuum,dtset%wtk(1:nkpt))

 dtset%kptrlen=kptrlen
!The fact that qptnrm is positive, non-zero, has been checked in inkpts.
 dtset%qptn(:)=dtset%qpt(:)/dtset%qptnrm

 dtset%kptns(:,1:nkpt)=dtset%kpt(:,1:nkpt)/dtset%kptnrm
 if(nqpt>=1 .and. dtset%optdriver/=1)then
   dtset%kptns(1,1:nkpt)=dtset%kptns(1,1:nkpt)+dtset%qptn(1)
   dtset%kptns(2,1:nkpt)=dtset%kptns(2,1:nkpt)+dtset%qptn(2)
   dtset%kptns(3,1:nkpt)=dtset%kptns(3,1:nkpt)+dtset%qptn(3)
 end if

 if(dtset%nkptgw>0) then

   token = 'bdgw'
   call intagm(dprarr,intarr,jdtset,marr,2*dtset%nkptgw*dtset%nsppol,string(1:lenstr),token,tread,'INT')
   if(tread==1) &
&   dtset%bdgw(1:2,1:dtset%nkptgw,1:dtset%nsppol)=reshape(intarr(1:2*dtset%nkptgw*dtset%nsppol),(/2,dtset%nkptgw,dtset%nsppol/))

   token = 'kptgw'
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%nkptgw,string(1:lenstr),&
&   token,tread,'DPR')
   if(tread==1) dtset%kptgw(1:3,1:dtset%nkptgw)=&
&   reshape(dprarr(1:3*dtset%nkptgw),(/3,dtset%nkptgw/))

 end if

!RS
 if(dtset%nqptdm>0) then
   token = 'qptdm'
   call intagm(dprarr,intarr,jdtset,marr,3*dtset%nqptdm,string(1:lenstr),&
&   token,tread,'DPR')
   if(tread==1) dtset%qptdm(1:3,1:dtset%nqptdm)=&
&   reshape(dprarr(1:3*dtset%nqptdm),(/3,dtset%nqptdm/))
 end if
!end RS

!Only read occ if (iscf >0 or iscf=-1 or iscf=-3) and (occopt==0 or occopt==2)
 if  (iscf>0.or.iscf==-1.or.iscf==-3)  then
   if (occopt==2 .and. getocc==0) then
!    Read occ(nband(kpt)*nkpt*nsppol) explicitly
     write(message, '(a)' )&
&     ' invars2: reading occ(nband*nkpt*nsppol) explicitly'
     call wrtout(std_out,message,'COLL')
     token = 'occ'
     call intagm(dprarr,intarr,jdtset,marr,bantot,&
&     string(1:lenstr),token,tread,'DPR')
     if(tread==1) dtset%occ_orig(1:bantot)=dprarr(1:bantot)
   else if(occopt==0) then
!    Read usual occupancy--same for all k points and spins
     token = 'occ'
     call intagm(dprarr,intarr,jdtset,marr,dtset%nband(1),&
&     string(1:lenstr),token,tread,'DPR')
     if(tread==1) dtset%occ_orig(1:dtset%nband(1))=dprarr(1:dtset%nband(1))
!    Fill in full occ array using input values for each k and spin
!    (make a separate copy for each k point and spin)
     do ikpt=1,nkpt*nsppol
       dtset%occ_orig(1+(ikpt-1)*dtset%nband(1):ikpt*dtset%nband(1))=&
&       dtset%occ_orig(1:dtset%nband(1))
     end do
   end if
 end if

 if(dtset%jellslab/=0)then

   token = 'slabwsrad'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
   if(tread==1) dtset%slabwsrad=dprarr(1)

!  Update number of electrons taking into account jellium
!  Suppose that the cell has z axis perpendicular to x and y axes.
!  This will be checked later
   areaxy=abs(dtset%rprimd_orig(1,1,1)*dtset%rprimd_orig(2,2,1)-dtset%rprimd_orig(1,2,1)*dtset%rprimd_orig(2,1,1))
   rhoavg=three/(four_pi*dtset%slabwsrad**3)
   nelectjell=areaxy*(dtset%slabzend-dtset%slabzbeg)*rhoavg
   charge=charge-nelectjell

 end if

!Initialize occ if occopt==1 or 3 ... 7,
!while if getocc/=0, make a fake initialization
!If iscf>0, check the charge of the system, and compute nelect.
 occopt_tmp=occopt
 if(getocc/=0)occopt_tmp=1
 call chkneu(charge,dtset,occopt_tmp)

!Initialize Berry phase vectors
!Should check that nberry is smaller than 20
 if(berryopt>0 .and. nberry>0)then
   token = 'kberry'
   call intagm(dprarr,intarr,jdtset,marr,3*nberry,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%kberry(1:3,1:nberry)=reshape(intarr(1:3*nberry), (/3,nberry/))
 end if

!Find the default mass
 allocate(mass_psp(npsp))
 do ipsp=1,npsp
   call atmdata(amu_default,rcov,symbol,dtset%znucl(ipsp))
   mass_psp(ipsp)=amu_default
 end do
!When the pseudo-atom is pure, simple copy
 if(ntyppure>0)then
   do itypat=1,ntyppure
     dtset%amu(itypat)=mass_psp(itypat)
   end do
 end if
!When the pseudo-atom is alchemical, must make mixing
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     dtset%amu(itypat)=zero
     do ipsp=ntyppure+1,npsp
       dtset%amu(itypat)=dtset%amu(itypat)+dtset%mixalch(ipsp-ntyppure,itypat-ntyppure)*mass_psp(ipsp)
     end do
   end do
 end if
 deallocate(mass_psp)

 token = 'amu'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%amu(1:ntypat)=dprarr(1:ntypat)
 token = 'densty'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread,'DPR')
 if(tread==1) dtset%densty(1:ntypat,1)=dprarr(1:ntypat)


 token = 'so_psp'
 call intagm(dprarr,intarr,jdtset,marr,npsp,string(1:lenstr),token,tread,'INT')
 if(tread==1)then
   dtset%so_psp(1:npsp)=intarr(1:npsp)
 end if

!This is to be removed after a few version (from 5.4)
 token = 'so_typat'
 call intagm(dprarr,intarr,jdtset,marr,ntypat,string(1:lenstr),token,tread_alt,'INT')
 if(tread_alt==1)then
   write(message, '(6a)' ) ch10,&
&   ' invars2: ERROR -',ch10,&
&   '  The use of the input variable so_typat has been discontinued in v5.4.',ch10,&
&   '  Action : please switch to so_psp (dimensioned to npsp - with different meanings of 0 and 1.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 deallocate(intarr,dprarr)

 call timab(191,2,tsec)

end subroutine invars2
!!***
