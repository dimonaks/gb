!{\Src2tex{textfont=tt}}
!!****f* ABINIT/invars1
!! NAME
!! invars1
!!
!! FUNCTION
!! Initialize the dimensions needed to allocate the input arrays
!! for one dataset characterized by jdtset, by
!! taking from string the necessary data.
!! Perform some preliminary checks and echo these dimensions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, AR, MKV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iout=unit number of output file
!!  jdtset=number of the dataset looked for
!!  lenstr=actual length of string
!!  mpi_enreg=informations about MPI parallelization
!!  msym=default maximal number of symmetries
!!  zion_max=maximal valence charge over all psps
!!
!! OUTPUT
!!  mband_upper=estimation of the maximum number of bands for any k-point
!!
!! SIDE EFFECTS
!! Input/Output (the default value is given in the calling routine)
!!  dtset=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!!   initialized, while some others will still be initialized later.
!!   The list of records of dtset initialized in the present routine is:
!!   acell_orig,densty,iatfix,kptopt,kptrlatt,
!!   mkmem,mkqmem,mk1mem,natsph,natvshift,nconeq,nkptgw,nkpt,
!!   nshiftk,nqptdm,optdriver,
!!   rprim_orig,rprimd_orig,shiftk,
!!   spgroup,spinat,typat,vel_orig,xred_orig
!!  bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!!  symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,1:msym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!!  string*(*)=string of characters containing all input variables and data
!!
!! NOTES
!! Must set up the geometry of the system, needed to compute
!! k point grids in an automatic fashion.
!! Treat separately mband_upper, since
!! fband, charge and zion_max must be known for being able to initialize it.
!!
!! Defaults are provided in the calling routine.
!! Defaults are also provided here for the following variables :
!! mband_upper, occopt, fband, charge
!! They should be kept consistent with defaults of the same variables
!! provided to the invars routines.
!!
!! PARENTS
!!      invars1m,newsp
!!
!! CHILDREN
!!      atmdata,distrb2,finddistrproc,ingeo,initmpi_grid,inkpts,intagm,inupper
!!      invacuum,leave_new,mkrdim,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine invars1(bravais,dtset,iout,jdtset,lenstr,mband_upper,mpi_enreg,msym,&
& string,symafm,symrel,tnons,zion_max)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_42_parser
 use interfaces_51_manage_mpi
 use interfaces_57_iovars, except_this_one => invars1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,jdtset,lenstr,msym
 integer,intent(out) :: mband_upper
 real(dp),intent(in) :: zion_max
 character(len=*),intent(inout) :: string
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
!arrays
 integer,intent(inout) :: bravais(11),symafm(msym),symrel(3,3,msym)
 real(dp),intent(inout) :: tnons(3,msym)

!Local variables-------------------------------
 character :: blank=' ',string1
!scalars
 integer :: chksymbreak,found,iband,ii,ikpt,iimage,index_blank,index_lower
 integer :: index_typsymb,index_upper,ipsp,iscf,isppol,intimage,leave,marr
 integer :: me,mkmem,natom,nband_k,nbdblock,nkpt,nkpt_me
 integer :: npsp,nqpt,nspink,nspinor,nsppol,occopt,response
 integer :: rfddk,rfelfd,rfmgfd,rfphon,rfstrs,rfuser
 integer :: tfband,tnband,tnbdblock,tread,tread_alt,twfoptalg
 integer :: wfoptalg
 real(dp) :: charge,dummy,fband,kptnrm,kptrlen,qptnrm,rcov
 character(len=2) :: string2,symbol
 character(len=30) :: token
 character(len=500) :: message
!arrays
 integer :: useri(5),vacuum(3)
 integer,allocatable :: intarr(:),istwfk(:),nband(:)
 real(dp) :: acell(3),qpt(3),userr(5),rprim(3,3)
 real(dp),allocatable :: dprarr(:),kpt(:,:),reaalloc(:),vel(:,:),wtk(:),xred(:,:)
 character(len=6) :: nm_mkmem(3)

!************************************************************************

!DEBUG
!write(6,*)' invars1 : enter '
!write(6,*)' jdtset = ',jdtset
!call flush(6)
!stop
!ENDDEBUG

!Read parameters
 npsp=dtset%npsp
 wfoptalg=dtset%wfoptalg
 nbdblock=dtset%nbdblock
 marr=npsp
 if (npsp < 3) marr = 3
 allocate(intarr(marr),dprarr(marr))

!---------------------------------------------------------------------------

 rfddk=0; rfelfd=0 ; rfmgfd=0; rfphon=0 ; rfstrs=0 ; rfuser=0
 token = 'rfddk'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) rfddk=intarr(1)
 token = 'rfelfd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) rfelfd=intarr(1)
 token = 'rfmgfd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) rfmgfd=intarr(1)
 token = 'rfphon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) rfphon=intarr(1)
 token = 'rfstrs'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) rfstrs=intarr(1)
 token = 'rfuser'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) rfuser=intarr(1)

 response=0
 if(rfddk/=0 .or. rfelfd/=0 .or. rfmgfd/=0 .or. &
& rfphon/=0 .or. rfstrs/=0 .or. rfuser/=0)response=1

 token = 'optdriver'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1)then
   dtset%optdriver=intarr(1)
   if(response==1 .and. dtset%optdriver/=1)then
     write(message,'(a,a,a,a,i3,a,a,a,a,i3,a,i3,a,i3,a,i3,a,i3,a,a,i3,a,a,a,a)' )ch10,&
&     ' invars1 : ERROR -',ch10,&
&     '  The input variable optdriver=',dtset%optdriver,ch10,&
&     '  This is in conflict with the values of the other input variables,',ch10,&
&     '  rfphon=',rfphon,'  rfddk=',rfddk,'  rfelfd=',rfelfd,'  rfmgfd=',rfmgfd,'  rfstrs=',rfstrs,ch10,&
&     '  rfuser=',rfuser,ch10,&
&     '  Action : check the values of optdriver, rfphon, rfddk, rfelfd, rfmgfd, rfstrs',ch10,&
&     '   and rfuser in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 else
!  If optdriver was not read, while response=1, set optdriver to 1
   if(response==1)dtset%optdriver=1
 end if

!---------------------------------------------------------------------------
!For now, waiting express parallelisation for recursion 
 token = 'tfkinfunc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%tfkinfunc=intarr(1)

!---------------------------------------------------------------------------

 natom=dtset%natom

!No default value for znucl
 token = 'znucl'
 call intagm(dprarr,intarr,jdtset,marr,npsp,string(1:lenstr),token,tread,'DPR')
 if(tread==1)then
   dtset%znucl(1:npsp)=dprarr(1:npsp)
 end if
 if(tread/=1)then
   write(message, '(a,a,a,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  The array znucl MUST be initialized in the input file,',&
&   ' while this is not done.',&
&   ch10,'  Action : initialize znucl in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Special treatment of _TYPAX (from a CML file), taking into account
!the fact that znucl does NOT depend on the dataset
!Examine all occurences of '_TYPAX'

 do
   index_typsymb=index(string(1:lenstr),'_TYPAX')
   if(index_typsymb==0)exit
!  Replace '_TYPAX' by '_TYPAT'
   string(index_typsymb:index_typsymb+5)='_TYPAT'
   index_upper=index_typsymb+5
!  Must start from the first blank after the tag (including possible dtset_char)
   index_upper=index(string(index_upper:lenstr),blank)+index_upper-1
   index_lower=index_upper

!  Examine all atoms (the end of the symbol string is delimited by a XX )
   do
     index_blank=index(string(index_upper:lenstr),blank)+index_upper-1
     string2=string(index_blank+1:index_blank+2)
     if(string2=="XX")exit
     found=0
!    Find the matching symbol
     do ipsp=1,npsp
       call atmdata(dummy,rcov,symbol,dtset%znucl(ipsp))
       call inupper(symbol)
       call inupper(string2)

!      DEBUG
!      write(6,'(a)')' invars1 : before test, trim(adjustl(symbol)),trim(adjustl(string2))'
!      write(6,'(5a)' )'"',trim(adjustl(symbol)),'","',trim(adjustl(string2)),'"'
!      ENDDEBUG

       if(trim(adjustl(symbol))==trim(adjustl(string2)))then
         found=1
         index_upper=index_blank+1
!        Cannot deal properly with more that 9 psps
         if(ipsp>=10)then
           write(message,'(4a)' )ch10,&
&           ' invars1 : BUG -',ch10,&
&           '  Need to use a pseudopotential with number larger than 9. Not allowed yet.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if

!        DEBUG
!        write(6,*)' invars1 : found ipsp=',ipsp
!        ENDDEBUG

         write(string1,'(i1)')ipsp
         string(index_lower:index_lower+1)=blank//string1
         index_lower=index_lower+2
       end if
     end do ! ipsp
!    if not found ...
     if(found==0)then
       write(message,'(9a)' )ch10,&
&       ' invars1 : ERROR -',ch10,&
&       '  Did not find matching pseudopotential for CML atomic symbol,',ch10,&
&       '  with value ',string2,ch10,&
&       '  Action : check that the atoms required by the CML file correspond to one psp file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end do ! Loop on atoms
!  One should find blanks after the last significant type value
   string(index_lower:index_blank+2)=blank
 end do ! loop to identify _TYPAX

!DEBUG
!write(6,*)' invars1 : trim(string)=',trim(string)
!ENDDEBUG

!---------------------------------------------------------------------------

!Here, set up quantities that are related to geometrical description
!of the system (acell,rprim,xred), as well as
!initial velocity(vel), and spin of atoms (spinat),
!the symmetries (symrel,symafm, and tnons)
!and the list of fixed atoms (iatfix,iatfixx,iatfixy,iatfixz).
!Arrays have already been
!dimensioned thanks to the knowledge of msym and mxnatom

 useri(1)=dtset%useria
 useri(2)=dtset%userib
 useri(3)=dtset%useric
 useri(4)=dtset%userid
 useri(5)=dtset%userie
 userr(1)=dtset%userra
 userr(2)=dtset%userrb
 userr(3)=dtset%userrc
 userr(4)=dtset%userrd
 userr(5)=dtset%userre

!ji: We need to read the electric field before calling ingeo
!****** Temporary ******

 token = 'berryopt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')

 if(tread==1) dtset%berryopt=intarr(1)
 if (dtset%berryopt == 4 .or. dtset%usewvl == 1) then
   token = 'efield'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
   if (tread==1) dtset%efield(1:3) = dprarr(1:3)
 end if
 if (dtset%berryopt == 5) then
   token = 'bfield'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
   if (tread==1) dtset%bfield(1:3) = dprarr(1:3)
 end if

!We need to know nsppol/nspinor/nspden before calling ingeo
 nsppol=dtset%nsppol
 token = 'nsppol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) nsppol=intarr(1)

 token = 'SpinPolarized'  !Alternate SIESTA definition of nsppol
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread_alt,'LOG')
 if(tread_alt==1)then
   if(tread==1)then
     write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&     ' invars1: ERROR -',ch10,&
&     '  nsppol and SpinPolarized cannot be specified simultaneously',ch10,&
&     '  for the same dataset.',ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,  message,'COLL')
     leave=1
   else
!    Note that SpinPolarized is a logical input variable
     nsppol=1
     if(intarr(1)==1)nsppol=2
     tread=1
   end if
 end if
 dtset%nsppol=nsppol

 token = 'nspinor'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nspinor=intarr(1)

!Has to read pawspnorb now, in order to adjust nspinor
 token = 'pawspnorb'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%pawspnorb=intarr(1)
 if (dtset%usepaw>0.and.dtset%pawspnorb>0) dtset%nspinor=max(2,dtset%nspinor)
 nspinor=dtset%nspinor

 token = 'nspden'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   dtset%nspden=intarr(1)
 else
   dtset%nspden=dtset%nsppol
 end if

 allocate(vel(3,natom),xred(3,natom))
 intimage=2 ; if(dtset%nimage==1)intimage=1
 do ii=1,dtset%nimage+1
   iimage=ii
   if(dtset%nimage==1 .and. ii==2)exit
   if(dtset%nimage==2 .and. ii==3)exit
   if(dtset%nimage>=2 .and. ii==intimage)cycle ! Will do the intermediate reference image at the last reading
   if(dtset%nimage>=2 .and. ii==dtset%nimage+1)iimage=intimage 

   write(message,'(a,i6)')' invars1 : treat image number',iimage
   call wrtout(std_out,message,'COLL')

!  Need to reset nsym to default value for each image
   dtset%nsym=0

!  Call ingeo for each image in turn, with the possible default values
   acell=dtset%acell_orig(1:3,iimage)
   rprim=dtset%rprim_orig(1:3,1:3,iimage)
   vel=dtset%vel_orig(1:3,1:natom,iimage)
   xred=dtset%xred_orig(1:3,1:natom,iimage)
   call ingeo(acell,dtset%berryopt,bravais,dtset%efield,dtset%genafm(1:3),dtset%iatfix(1:3,1:natom),&
&   dtset%icoulomb,iimage,iout,jdtset,dtset%jellslab,lenstr,msym,natom,dtset%nimage,dtset%nspden,dtset%nsppol,&
&   dtset%nsym,dtset%ntypat,dtset%pawspnorb,dtset%ptgroupma,&
&   rprim,dtset%slabzbeg,dtset%slabzend,dtset%spgroup,dtset%spinat(1:3,1:natom),&
&   string,symafm,dtset%symmorphi,symrel,tnons,dtset%tolsym,dtset%typat(1:natom),vel,xred)
   dtset%acell_orig(1:3,iimage)=acell
   dtset%rprim_orig(1:3,1:3,iimage)=rprim
   dtset%vel_orig(1:3,1:natom,iimage)=vel
   dtset%xred_orig(1:3,1:natom,iimage)=xred
   call mkrdim(dtset%acell_orig(1:3,iimage),dtset%rprim_orig(1:3,1:3,iimage),dtset%rprimd_orig(1:3,1:3,iimage))
 end do
 deallocate(vel,xred)

!Examine whether there is some vacuum space in the unit cell
 call invacuum(jdtset,lenstr,natom,dtset%rprimd_orig(1:3,1:3,intimage),string,vacuum,&
& dtset%xred_orig(1:3,1:natom,intimage))

!DEBUG
!write(6,*)' invars1: before inkpts, vacuum=',vacuum(:)
!stop
!ENDDEBUG

!---------------------------------------------------------------------------

!Set up k point grid number
!First, get additional information
 dtset%kptopt=1
 token = 'kptopt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%kptopt=intarr(1)
 iscf=5
 token = 'iscf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) iscf=intarr(1)


 dtset%natsph=dtset%natom
 token = 'natsph'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%natsph=intarr(1)

 token = 'natvshift'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%natvshift=intarr(1)

 token = 'nconeq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nconeq=intarr(1)

!For the time being (v4.3) keep the opportunity to use the old input variable name
 token = 'ngwpt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nkptgw=intarr(1)
 token = 'nkptgw'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nkptgw=intarr(1)
 if (dtset%nkptgw<0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input nkptgw must be >= 0, but was ',dtset%nkptgw,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Number of points for long wavelength limit.
 token = 'gw_nqlwl'
!Default is dtset%gw_nqlwl=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%gw_nqlwl=intarr(1)
 if (dtset%gw_nqlwl<0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input gw_nqlwl must be > 0, but was ',dtset%gw_nqlwl,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Read number of k-points (if specified)
 nkpt=0
 if(dtset%kptopt==0)nkpt=1
 token = 'nkpt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) nkpt=intarr(1)
 dtset%nkpt=nkpt

 if (nkpt<0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input nkpt must be >= 0, but was ',nkpt,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

 if (nkpt==0 .and. dtset%kptopt==0 ) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  For kptopt=0, input nkpt must be > 0, but was ',nkpt,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Will compute the actual value of nkpt, if needed. Otherwise,
!test that the value of nkpt is OK, if kptopt/=0
!Set up dummy arrays istwfk, kpt, wtk

 if(nkpt/=0 .or. dtset%kptopt/=0)then
   allocate(istwfk(nkpt),kpt(3,nkpt),wtk(nkpt))
!  Here, occopt is also a dummy argument
   occopt=1 ; dtset%nshiftk=1 ; dtset%kptrlatt(:,:)=0
   kptrlen=20.0_dp ; wtk(:)=1.0_dp 
   dtset%shiftk(:,:)=half

   nqpt=0
   token = 'nqpt'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) nqpt=intarr(1)

   chksymbreak=1
   token = 'chksymbreak'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) chksymbreak=intarr(1)

!  Use the first image to predict k points,
!  except if an intermediate image is available
   intimage=1 ; if(dtset%nimage>2)intimage=(1+dtset%nimage)/2
   call inkpts(bravais,chksymbreak,iout,iscf,istwfk,jdtset,&
&   kpt,dtset%kptopt,kptnrm,dtset%kptrlatt,kptrlen,lenstr,&
&   msym,nkpt,nqpt,dtset%nshiftk,dtset%nsym,occopt,&
&   qpt,qptnrm,response,dtset%rprimd_orig(1:3,1:3,intimage),&
&   dtset%shiftk,string,symafm,symrel,vacuum,wtk)

   deallocate(istwfk,kpt,wtk)
!  nkpt has been computed, as well as the k point grid, if needed
   dtset%nkpt=nkpt
 end if

!dtset%nqptdm=0
 token = 'nqptdm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nqptdm=intarr(1)

 if (dtset%nqptdm<0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input nqptdm must be >= 0, but was ',dtset%nqptdm,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!---------------------------------------------------------------------------

 token = 'nnos'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%nnos=intarr(1)

!Perform the first checks

 leave=0

!Check that nkpt is greater than 0
 if (nkpt<=0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  After inkpts, nkpt must be > 0, but was ',nkpt,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,  message,'COLL')
   leave=1
 end if

!Check that nsppol is 1 or 2
 if (nsppol/=1 .and. nsppol/=2) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input nsppol must be 1 or 2, but was ',nsppol,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   leave=1
 end if

!Check that nspinor is 1 or 2
 if (nspinor/=1 .and. nspinor/=2) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input nspinor must be 1 or 2, but was ',nspinor,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   leave=1
 end if

!Check that nspinor and nsppol are not 2 together
 if (nsppol==2 .and. nspinor==2) then
   write(message, '(8a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  nspinor and nsappol cannot be 2 together !',ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   leave=1
 end if

!Here, leave if an error has been detected earlier
 if(leave==1) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' invars1 : WARNING -',ch10,&
&   '  Other errors might be present in the input file. '
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if


!Now, take care of mband_upper

 mband_upper=1
 occopt=1
 fband=0.5_dp

 token = 'occopt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) occopt=intarr(1)

!Also read fband, that is an alternative to nband. The default
!is different for occopt==1 and for metallic occupations.
 if(occopt==1)fband=0.125_dp
 token = 'fband'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tfband,'DPR')
 if(tfband==1)fband=dprarr(1)

!fband cannot be used when occopt==0 or occopt==2
 if(tfband==1 .and. (occopt==0 .or. occopt==2) )then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' invars1 : ERROR - ',ch10,&
&   '  fband cannot be used if occopt==0 or occopt==2 ',ch10,&
&   '  Action : correct your input file, suppress fband, or change occopt.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 allocate(nband(nkpt*nsppol))
 tnband=0

 if (occopt==0 .or. occopt==1 .or. (occopt>=3 .and. occopt<=7) ) then
   token = 'nband'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tnband,'INT')
!  Note : mband_upper is initialized, not nband
   if(tnband==1) mband_upper=intarr(1)

   if(tfband==1 .and. tnband==1)then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' invars1 : ERROR -',ch10,&
&     '  fband and nband cannot be used together. ',ch10,&
&     '  Action : correct your input file, suppress either fband or nband.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  In case nband was not read, use fband, either read, or the default,
!  to provide an upper limit for mband_upper
   if(tnband==0)then

!    Checks whether the call was not issued by newsp, in which case zion_max=0.0_dp
     if(abs(zion_max) < 1.0d-10 )then
       write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&       ' invars1 : ERROR -',ch10,&
&       '  Found that nband in not defined in the input file,',ch10,&
&       '  while the running code is newsp. This is not allowed.',ch10,&
&       '  Action : correct your input file : define nband. '
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     charge=0.0_dp
     token = 'charge'
     call intagm(dprarr,intarr,jdtset,marr,1,&
&     string(1:lenstr),token,tread,'DPR')
     if(tread==1) charge=dprarr(1)

!    Only take into account negative charge, to compute maximum number of bands
     if(charge > 0.0_dp)charge=0.0_dp

     mband_upper=nspinor*((nint(zion_max)*natom+1)/2 - floor(charge/2.0_dp)&
&     + ceiling(fband*natom-1.0d-10))

     nband(:)=mband_upper

!    DEBUG
!    write(6,*)' invars1 : zion_max,natom,fband,mband_upper '
!    write(6,*)zion_max,natom,fband,mband_upper
!    ENDDEBUG

   end if

   nband(:)=mband_upper

 else if (occopt==2) then
   allocate(reaalloc(nkpt*nsppol))
   token = 'nband'
   call intagm(reaalloc,nband,jdtset,nkpt*nsppol,nkpt*nsppol,&
&   string(1:lenstr),token,tnband,'INT')
   if(tnband==1)then
     do ikpt=1,nkpt*nsppol
       if (nband(ikpt)>mband_upper) mband_upper=nband(ikpt)
     end do
   end if
   deallocate(reaalloc)
 else
   write(message, '(a,a,a,a,i8,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  occopt=',occopt,' is not an allowed value.',ch10,&
&   '  Action : correct your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Check that mband_upper is greater than 0
 if (mband_upper<=0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Maximal nband must be > 0, but was ',mband_upper,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

 token = 'wfoptalg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr), &
& token,twfoptalg,'INT')
 if (twfoptalg == 1)wfoptalg=intarr(1)

 token = 'nbdblock'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr), &
& token,tnbdblock,'INT')
 if(tnbdblock ==1)nbdblock=intarr(1)

!paral_kgb is now read in invars0

 token = 'npimage'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npimage=intarr(1)

 token = 'npkpt'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npkpt=intarr(1)

 token = 'npfft'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npfft=intarr(1)

 token = 'npband'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%npband=intarr(1)

 token = 'bandpp'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%bandpp=intarr(1)

 if(dtset%paral_kgb/=0.and.(dtset%bandpp/=1.or.dtset%npband/=1.or.&
& dtset%npfft/=1.or.dtset%npkpt/=1))then
   if(dtset%npkpt*dtset%npfft*dtset%npband > mpi_enreg%nproc )then
     write(message,'(10a)') ch10,&
&     ' abinit : WARNING -',ch10,&
&     '  The product of npkpt, npfft and npband is bigger than the number of processors.',ch10,&
&     '  The user-defined values of npkpt, npfft or npband will be modified,',ch10,&
&     '  in order to bring this product below nproc .',ch10,&
&     '  At present, only a very simple algorithm is used ...'
     call wrtout(std_out,message,'COLL')
     if(dtset%npkpt*dtset%npband <= mpi_enreg%nproc)then
       dtset%npfft=1
       write(message,'(4a)') ch10,&
&       ' abinit : WARNING -',ch10,&
&       '  Set npfft to 1'
       call wrtout(std_out,message,'COLL')
     else if(dtset%npkpt <= mpi_enreg%nproc)then
       dtset%npfft=1
       dtset%npband=1
       write(message,'(4a)') ch10,&
&       ' abinit : WARNING -',ch10,&
&       '  Set npfft and npband to 1'
       call wrtout(std_out,message,'COLL')
     else
       dtset%npfft=1
       dtset%npband=1
       dtset%npkpt=1
       write(message,'(4a)') ch10,&
&       ' abinit : WARNING -',ch10,&
&       '  Set npfft, npband, and npkpt to 1'
       call wrtout(std_out,message,'COLL')
     end if
   else if(dtset%npkpt*dtset%npfft*dtset%npband < mpi_enreg%nproc )then
     write(message,'(5a)') ch10,&
&     ' abinit : STOP -',ch10,&
&     ' The number of processor must not be greater than npfft*npband*npkpt',&
&     ' when npfft or npkpt or npband are chosen manually in the input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
 call finddistrproc(dtset,mband_upper,mpi_enreg)
!Initialize mpi_enreg for this dataset
!Prior to this line, the content of mpi_enreg was dataset-independent.
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_respfn=0
 mpi_enreg%paral_level=0
 mpi_enreg%paralbd=0
 mpi_enreg%mode_para="n"
 mpi_enreg%ngroup_respfn=1
 mpi_enreg%me_respfn=-1

 
!--IF CUDA AND RECURSION: FOR NOW CUDA ONLY BAND PARALLELISATION

 if(dtset%tfkinfunc==2 .and. mpi_enreg%nproc/=1)then
   dtset%npband = dtset%npband*dtset%npkpt*dtset%npfft
   dtset%npkpt = 1
   dtset%npfft = 1
   write(message, '(8a,i6,a)' ) ch10,&
&   ' invars1: WARNING -',ch10,&
&   '  If HAVE_GPU_CUDA and recursion are used ',ch10,&
&   '  only the band parallelisation is active, we set:',ch10,&
&   '  npfft= 1, npkpt= 1, npband=',dtset%npband,' .'
   call wrtout(std_out,  message,'COLL')
 end if


 mpi_enreg%nproc_kpt=dtset%npkpt
 mpi_enreg%nproc_fft=dtset%npfft
 mpi_enreg%nproc_band=dtset%npband

 me=mpi_enreg%me

 if(dtset%paral_kgb/=0)then
   mpi_enreg%paral_compil_fft=1
   mpi_enreg%nproc_kpt=dtset%npkpt
   mpi_enreg%nproc_fft=dtset%npfft
   mpi_enreg%nproc_band=dtset%npband
   call initmpi_grid(mpi_enreg)
 else
   mpi_enreg%paral_compil_fft=0
   mpi_enreg%nproc_kpt=mpi_enreg%nproc
   mpi_enreg%nproc_fft=1
   mpi_enreg%nproc_band=1
   mpi_enreg%me_kpt=mpi_enreg%me
   call initmpi_grid(mpi_enreg)
 end if

!Compute the maximal number of k points treated per processor
 nkpt_me=nkpt
 if(mpi_enreg%paral_compil==1 .and. dtset%usewvl == 0) then
!  Determine who I am
!  Define k-points distribution
!  Note that nkpt_me may differ from processor to processor
!  This fact will NOT be taken into account when
!  the memory needs will be evaluated in the subroutine memory.
!  Also, the reduction of k points due to symmetry in RF calculations
!  is NOT taken into account. This should be changed later ...
   allocate(mpi_enreg%proc_distrb(nkpt,mband_upper,nsppol))
   if(response==0)then

     if (mod(wfoptalg,10) == 1) then

       if (nkpt >= mpi_enreg%nproc_kpt) then
         mpi_enreg%paralbd=0
       else
         if (nbdblock == 1) then
           mpi_enreg%paralbd=0
         else
           mpi_enreg%paralbd=nbdblock
         end if
       end if

     else  ! wfoptalg==0,2,3,4

       mpi_enreg%paralbd=0

     end if ! wftoptalg

!    Set paral_compil_respfn and ngroup_respfn so distrb2 calculates the right proc_distrb for the actual dtset
     mpi_enreg%paral_compil_respfn=dtset%paral_rf
     mpi_enreg%ngroup_respfn=dtset%ngroup_rf
     call distrb2(mband_upper,nband,nkpt,nsppol,mpi_enreg)
     mpi_enreg%paral_compil_respfn=0
     mpi_enreg%ngroup_respfn=-1
     nkpt_me=0
     do ikpt=1,nkpt
       if (mpi_enreg%paralbd == 0) then
         if(minval(abs(mpi_enreg%proc_distrb(ikpt,1,1:nsppol)-me))==0)&
&         nkpt_me=nkpt_me+1
       else !paralbd > 1
         do isppol=1,nsppol
           nband_k=nband(ikpt+(isppol-1)*nkpt)
           if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))==0)&
&           nkpt_me=nkpt_me+1
         end do
       end if ! paralbd
     end do ! ikpt=1,nkpt

   else ! response==1

!    Wrongly assumes that the number of elements of the
!    k-point sets of the two spin polarizations is the maximal
!    value of one of these k-point sets ...
!    This is to be corrected when RF is implemented
!    for spin-polarized case.
     nkpt_me=0
     mpi_enreg%paralbd=1

!    Set paral_compil_respfn and ngroup_respfn so distrb2 calculates the right proc_distrb for the actual dtset
     mpi_enreg%paral_compil_respfn=dtset%paral_rf
     mpi_enreg%ngroup_respfn=dtset%ngroup_rf
     call distrb2(mband_upper,nband,nkpt,nsppol,mpi_enreg)
     mpi_enreg%paral_compil_respfn=0
     mpi_enreg%ngroup_respfn=-1

     do isppol=1,nsppol
       nspink=0
       do ikpt=1,nkpt
         do iband=1,nband(ikpt+(isppol-1)*nkpt)
           if(mpi_enreg%proc_distrb(ikpt,iband,isppol)==me)then
             nspink=nspink+1
             exit
           end if
         end do ! iband
       end do ! ikpt
       if(nspink>nkpt_me)nkpt_me=nspink
     end do ! isppol
!    If the number of bands was estimated, there might be a side effect
!    when the definitive number of bands is known. k points
!    might be attributed to different processors than the present
!    proc_distrb describes. At most, the number of k points could
!    increase by 1 ...
     if(tnband==0)nkpt_me=nkpt_me+1
!    In any case, the maximal number of k points is nkpt
     if(nkpt_me>nkpt)nkpt_me=nkpt
   end if
   deallocate(mpi_enreg%proc_distrb)
 end if

!Take care of mkmems. Use the generic name -mkmem- for mkmem as well as mkqmem
!and mk1mem.
 nm_mkmem(1)='mkmem '
 nm_mkmem(2)='mkqmem'
 nm_mkmem(3)='mk1mem'

 do ii=1,3

!  Read in mkmem here if it is in the input file
   if(ii==1)then
     token = 'mkmem'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   else if(ii==2)then
     token = 'mkqmem'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   else if(ii==3)then
     token = 'mk1mem'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   end if

!  Note that mkmem is used as a dummy variable, representing mkmem as well
!  as mkqmem, and mk1mem.
   if(tread==1) then
     mkmem=intarr(1)
     if (mkmem<0) then
!      mkmem is unreasonable; must be zero or positive
       write(message, '(a,a,a,a,a,a,a,a,i6,a,a,a,a)' ) ch10,&
&       ' invars1: WARNING -',ch10,&
&       '  ',nm_mkmem(ii),' must be positive or nul but ',&
&       nm_mkmem(ii),' =',mkmem,ch10,&
&       '  Use default ',nm_mkmem(ii),' = nkpt .'
       call wrtout(std_out,message,'COLL')
       mkmem=nkpt
     end if

   else

!    mkmem was not set in the input file so default to incore solution
     write(message, '(a,a,a,a,a,a)' ) &
&     ' invars1: ',nm_mkmem(ii),' undefined in the input file.',&
&     ' Use default ',nm_mkmem(ii),' = nkpt'
     call wrtout(std_out,message,'COLL')
     mkmem=nkpt
   end if

!  Check whether nkpt distributed on the processors <= mkmem;
!  if so then may run entirely in core,
!  avoiding i/o to disk for wavefunctions and kg data.
!  mkmem/=0 to avoid i/o; mkmem==0 to use disk i/o for nkpt>=1.
   if (nkpt_me<=mkmem .and. mkmem/=0 ) then
     write(message, '(a,i5,a,a,a,i5,a)' ) &
&     ' invars1: With nkpt_me=',nkpt_me,&
&     ' and ',nm_mkmem(ii),' = ',mkmem,', ground state wf handled in core.'
     call wrtout(std_out,message,'COLL')
     if(nkpt_me<mkmem .and. nkpt_me/=0)then
       write(message, '(a,a,a)' ) &
&       ' Resetting ',nm_mkmem(ii),' to nkpt_me to save memory space.'
       mkmem=nkpt_me
       call wrtout(std_out,message,'COLL')
     end if
   else if(mkmem/=0)then
     write(message, '(a,i5,a,a,a,i5,a,a,a,a,a)' ) &
&     ' invars1: With nkpt_me=',nkpt_me,&
&     ' and ',nm_mkmem(ii),' = ',mkmem,&
&     ' ground state wf require disk i/o.',ch10,&
&     ' Resetting ',nm_mkmem(ii),' to zero to save memory space.'
     mkmem=0
     call wrtout(std_out,message,'COLL')
   end if

   if(ii==1)dtset%mkmem=mkmem
   if(ii==2)dtset%mkqmem=mkmem
   if(ii==3)dtset%mk1mem=mkmem

!  check condition for no empty processors
!  mkmem_eff = ceiling(dtset%nkpt/mpi_enreg%nproc_kpt)
   if (mpi_enreg%nproc_kpt - floor(nsppol*dtset%nkpt*one/mkmem) &
&   >= mkmem) then
     write (message,'(4a,I6,3a,I6,6a,I6,2a)') ch10, &
&     ' invars1 : WARNING -',ch10, &
&     '   your number of k-points (',dtset%nkpt,') will not distribute correctly:',ch10, &
&     '   with the current number of processors (',mpi_enreg%nproc_kpt,').',ch10,&
&     '   you will leave some empty and this is forbidden.',ch10,ch10, &
&     ' ACTION: you can reduce number of processors to ', floor(nsppol*dtset%nkpt*one/mkmem),ch10, &
&     '   without losing speed'
     call wrtout(std_out,message,'COLL')
!    call leave_new('COLL')
   end if

!  End the loop on the three possiblities mkmem, mkqmem, mk1mem.
 end do


!---------------------------------------------------------------------------
!Some PAW+DMFT keywords
 token = 'usedmft'
 dtset%usedmft=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%usedmft=intarr(1)


!Some PAW+U keywords
 token = 'usepawu'
 dtset%usepawu=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%usepawu=intarr(1)
 
 if (dtset%usepawu > 0 .and. dtset%usedmft > 0) then
   write(message, '(10a)' ) ch10,&
&   ' invars1: WARNING -',ch10,&
&   '  usedmft and usepawu are both activated ',ch10,&
&   '  This is not an usual calculation:',ch10,&
&   '  usepawu will be put to 10:',ch10,&
&   '   LDA+U potential and energy will be put to zero'
   call wrtout(06,  message,'COLL')
 end if
 if ( dtset%usedmft > 0 .and. dtset%usepawu >= 0 ) dtset%usepawu = 10

 dtset%usedmatpu=0
 dtset%lpawu(1:dtset%ntypat)=-1

 if (dtset%usepawu>0) then
   token = 'lpawu'
   call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%lpawu(1:dtset%ntypat)=intarr(1:dtset%ntypat)

   token = 'usedmatpu'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%usedmatpu=intarr(1)

 end if

!Some PAW+Exact exchange keywords
 token = 'useexexch'
 dtset%useexexch=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%useexexch=intarr(1)


 dtset%lexexch(1:dtset%ntypat)=-1

 if (dtset%useexexch>0) then
   token = 'lexexch'
   call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtset%lexexch(1:dtset%ntypat)=intarr(1:dtset%ntypat)
 end if

!Macro_uj (determination of U in PAW+U), governs also allocation of atvshift
 dtset%macro_uj = 0
 token = 'macro_uj'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%macro_uj=intarr(1)

 dtset%normpawu(1:dtset%ntypat)=0
 token = 'normpawu'
 call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),token,tread,'INT')
 if(tread==1) dtset%normpawu(1:dtset%ntypat)=intarr(1:dtset%ntypat)

 deallocate(nband)
 deallocate(intarr,dprarr)

!DEBUG
!write(6,*)' invars1 : exit '
!call flush(6)
!stop
!ENDDEBUG

end subroutine invars1
!!***
