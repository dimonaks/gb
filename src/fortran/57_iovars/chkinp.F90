!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkinp
!! NAME
!! chkinp
!!
!! FUNCTION
!! Check consistency of input data against itself.
!! Please : use the alphabetic order
!! Please : use the routines chkint_eq, chkint_ne, chkint_ge, chkint_le, and chkdpr
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, MKV, DRH, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for output file
!!  mpi_enreg=informations about MPI parallelization
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      chkdpr,chkgrp,chkint,chkint_eq,chkint_ge,chkint_le,chkint_ne,chkorthsy
!!      dtsetcopy,dtsetfree,leave_new,metric,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chkinp(dtsets,iout,mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads)

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_errors
#if defined HAVE_ETSF_IO
 use etsf_io
#endif

 use m_numeric_tools,  only : iseven

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_53_abiutil
 use interfaces_57_iovars, except_this_one => chkinp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,ndtset,ndtset_alloc,npsp
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 logical :: twvl
 integer :: bantot,fftalg,ia,iatom,ib,iband,idtset,ierr,ii,iimage,ikpt,ilang,intimage,ierrgrp
 integer :: ipsp,isppol,isym,itypat,jdtset,maxiatsph,maxidyn,mband,miniatsph,minidyn,mod10,mu,natom
 integer :: nfft,nfftdg,nkpt,nspden,nspinor,nsppol,optdriver,response,usepaw,usewvl
 real(dp) :: delta,sumalch,sumocc,ucvol,wvl_hgrid,zatom
 character(len=500) :: message,msg
 type(dataset_type) :: dt
!arrays
 integer :: cond_values(3),nprojmax(0:3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: frac(:,:)
 character(len=32) :: cond_string(3)
 character(len=32) :: input_name

! *************************************************************************

 DBG_ENTER("COLL")

!Print machine precision (other machine parameters are computed
!in the dlamch function, see Lapack library)
 write(message,'(a,a,1p,e24.16)' ) ch10,&
& ' chkinp: machine precision is ',epsilon(0.0_dp)
 call wrtout(std_out,  message,'COLL')

!Some initialisations
 ierr=0
 cond_string(1:3)=' '
 cond_values(1:3)=(/0,0,0/)

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc

   jdtset=dtsets(idtset)%jdtset
   if(ndtset==0)jdtset=0

   if(jdtset/=0)then
     write(message, '(a,a,a,i2,a)' ) ch10,&
&     ' chkinp: Checking input parameters for consistency,',&
&     ' jdtset=',jdtset,'.'
   else
     write(message, '(a,a)' ) ch10,&
&     ' chkinp: Checking input parameters for consistency.'
   end if
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  Will test directly on the dataset "dt"
   call dtsetCopy(dt, dtsets(idtset))

!  Copy or initialize locally a few input dataset values
   fftalg   =dt%ngfft(7)
   natom    =dt%natom
   nkpt     =dt%nkpt
   nspden   =dt%nspden
   nspinor  =dt%nspinor
   nsppol   =dt%nsppol
   optdriver=dt%optdriver
   usepaw   =dt%usepaw
   usewvl   =dt%usewvl
   intimage=1 ; if(dtsets(idtset)%nimage>2)intimage=2
   rprimd(:,:)=dtsets(idtset)%rprimd_orig(:,:,intimage)    ! For the purpose of checking symmetries
   response=0
   if(dt%rfelfd/=0 .or. dt%rfmgfd/=0 .or. dt%rfphon/=0 .or. dt%rfstrs/=0 .or. dt%rfddk/=0  )response=1
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!  =====================================================================================================
!  Check the values of variables, using alphabetical order
!  PLEASE : use the routines chkint_eq, chkint_ne, chkint_ge, chkint_le, chkdpr

!  accesswff
!  Must be one of 0, 1, 3
   call chkint_eq(0,0,cond_string,cond_values,ierr,'accesswff',dt%accesswff,3,&
&   (/IO_MODE_FORTRAN,IO_MODE_MPI,IO_MODE_ETSF/),iout)
!  However, if mpi_io is not enabled, must be one of 0, 3.
   if(mpi_enreg%paral_compil_mpio==0)then
     cond_string(1)='enable_mpi_io' ;  cond_values(1)=0
!    Make sure that accesswff is 0 or 3
     call chkint_eq(1,1,cond_string,cond_values,ierr,'accesswff',dt%accesswff,2,(/IO_MODE_FORTRAN,IO_MODE_ETSF/),iout)
   end if

!  amu
!  Check that atomic masses are > 0 if ionmov = 1
   if (dt%ionmov==1) then
     do itypat=1,dt%ntypat
       cond_string(1)='ionmov' ; cond_values(1)=1
       write(input_name,'(a4,i1,a1)')'amu(',itypat,')'
       call chkdpr(1,1,cond_string,cond_values,ierr,input_name,dt%amu(itypat),1,tol8,iout)
     end do
   end if

!  bdberry
   if(dt%berryopt>0 .and. dt%berryopt/=4 .and. dt%berryopt/=5 .and. dt%nberry>0)then
     cond_string(1)='berryopt' ; cond_values(1)=dt%berryopt
     cond_string(2)='nberry'   ; cond_values(2)=dt%nberry
     do ii=1,2*nsppol
       write(input_name,'(a4,i1,a1)')'bdberry(',ii,')'
       call chkint_ge(2,2,cond_string,cond_values,ierr,input_name,dt%bdberry(ii),1,iout)
     end do
!    bdberry(2) must be greater than bdberry(1)
     cond_string(3)='bdberry(1)'   ; cond_values(3)=dt%bdberry(1)
     call chkint_ge(3,3,cond_string,cond_values,ierr,'bdberry(2)',dt%bdberry(2),dt%bdberry(1),iout)
     if(nsppol==2)then
!      bdberry(4) must be greater than bdberry(3)
       cond_string(3)='bdberry(3)'   ; cond_values(3)=dt%bdberry(3)
       call chkint_ge(3,3,cond_string,cond_values,ierr,'bdberry(4)',dt%bdberry(4),dt%bdberry(3),iout)
     end if
!    Make sure all nband(nkpt) are >= bdberry
     do isppol=1,nsppol
       do ikpt=1,nkpt
         if (dt%nband(ikpt+(isppol-1)*nkpt)<=dt%bdberry(2*isppol)) then
           cond_string(1)='ikpt'   ; cond_values(1)=ikpt
           cond_string(2)='isppol' ; cond_values(2)=isppol
           cond_string(3)='nband'  ; cond_values(3)=dt%nband(ikpt+(isppol-1)*nkpt)
           call chkint_le(0,3,cond_string,cond_values,ierr,&
&           'bdberry',dt%bdberry(2*isppol),dt%nband(ikpt+(isppol-1)*nkpt),iout)
           if(ierr==1)exit
         end if
       end do
     end do
   end if

!  berryopt
!  berryopt must be between -5 or -3 to +5
   call chkint_eq(0,0,cond_string,cond_values,ierr,'berryopt',dt%berryopt,10,(/-5,-3,-2,-1,0,1,2,3,4,5/),iout)
!  berryopt not allowed when nspinor/=1
   if(nspinor/=1)then
     cond_string(1)='nspinor' ; cond_values(1)=nspinor
     call chkint_eq(1,1,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/0/),iout)
   end if
!  berryopt must be positive when mkmem==0
   if(dt%mkmem==0)then
     cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
     call chkint_ge(1,1,cond_string,cond_values,ierr,'berryopt',dt%berryopt,0,iout)
   end if
!  berryopt must be positive when occopt==1
   if(dt%occopt/=1)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkint_ge(1,1,cond_string,cond_values,ierr,'berryopt',dt%berryopt,0,iout)
   end if
!  berryopt cannot be 4 or 5 when toldfe, toldff and tolrff are zero or negative
   if ((dt%toldfe < tiny(one)).and.(dt%toldff < tiny(one)).and.(dt%tolrff < tiny(one))) then
     cond_string(1)='toldfe' ; cond_values(1)=dt%toldfe
     cond_string(2)='toldff' ; cond_values(2)=dt%toldff
     cond_string(3)='tolrff' ; cond_values(3)=dt%tolrff
     call chkint_ne(3,3,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/4/),iout)
     call chkint_ne(3,3,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/5/),iout)
   end if

!  Non-zero berryopt and usepaw==1 cannot be done unless response==0
!  Non-zero berryopt and usepaw==1 cannot be done unless nsppol=1
!  Non-zero berryopt and usepaw==1 cannot be done unless nspinor=1
   if (usepaw==1.and.dt%berryopt/=0) then
     cond_string(1)='usepaw' ; cond_values(1)=1
     cond_string(2)='berryopt' ; cond_values(2)=dt%berryopt
     call chkint_eq(1,2,cond_string,cond_values,ierr,'response',response,1,(/0/),iout)
     call chkint_eq(1,2,cond_string,cond_values,ierr,'nsppol',dt%nsppol,1,(/1/),iout)
     call chkint_eq(1,2,cond_string,cond_values,ierr,'nspinor',dt%nspinor,1,(/1/),iout)
   end if

!  Non-zero berryopt and usepaw==1 and kptopt/=3 cannot be done unless symmorphi=0 (that is,
!  nonsymmorphic symmetries dont work yet
   if (usepaw==1.and.dt%berryopt/=0.and.dt%kptopt/=3) then
     cond_string(1)='usepaw'; cond_values(1)=1
     cond_string(2)='berryopt'; cond_values(2)=dt%berryopt
     cond_string(3)='kptopt'; cond_values(3)=dt%kptopt
     call chkint_eq(1,3,cond_string,cond_values,ierr,'symmorphi',dt%symmorphi,1,(/0/),iout)
     if (dt%symmorphi /= 0) then
       write(message, '(6a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '   under the current conditions, symmorphi must be zero! ',ch10,&
&       '  Action : set symmorphi to zero in the input file.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  Non-zero berryopt and usepaw==1 and kpt // requires nproc to be a divisor of nkpt
   if (usepaw==1.and.dt%berryopt/=0.and.mpi_enreg%nproc>1) then
     if (mod(dt%nkpt,mpi_enreg%nproc) /= 0) then
       write(message, '(6a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '   For berryopt /= 0 with PAW in parallel, nproc must be a divisor of nkpt ',ch10,&
&       '  Action : change number of processes or kpts such that nproc divides nkpt evenly '
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  boxcutmin
   if(response==1)then
     cond_string(1)='response' ; cond_values(1)=1
     call chkdpr(1,1,cond_string,cond_values,ierr,'boxcutmin',dt%boxcutmin,0,two,iout)
   end if

!  chksymbreak
   call chkint_eq(0,0,cond_string,cond_values,ierr,'chksymbreak',dt%chksymbreak,2,(/0,1/),iout)
   if(dt%chksymbreak==1)then
!    Check the values of tnons
     do isym=1,dt%nsym
       do ii=1,3
         delta=dt%tnons(ii,isym)*eight
         if(delta-nint(delta)>tol6)then
           delta=dt%tnons(ii,isym)*three*four
           if(delta-nint(delta)>tol6)then
             write(message, '(6a,i4,2a,9i3,2a,3es16.6,4a)' ) ch10,&
&             ' chkinp: ERROR -',ch10,&
&             '   Chksymbreak=1 . A potentially symmetry-breaking value of tnons has been observed :',ch10,&
&             '   for the symmetry number ',isym,ch10,&
&             '   symrel is ',dt%symrel(1:3,1:3,isym),ch10,&
&             '   tnons is ',dt%tnons(1:3,isym),ch10,&
&             '   Please, read the description of the input variable chksymbreak,',ch10,&
&             '   then, if you feel confident, you might switch it to zero, or consult with the forum.'
             call wrtout(iout,message,'COLL')
             call wrtout(std_out,  message,'COLL')
             ierr=ierr+1
           end if
         end if
       end do
     end do
   end if

!  diecut
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
     cond_string(2)='4*ecut' ; cond_values(1)=4*dt%ecut
!    Checks that presently diecut is 4*ecut
     call chkdpr(1,1,cond_string,cond_values,ierr,'diecut',dt%diecut,0,4*dt%ecut,iout)
   end if

!  diemac
   call chkdpr(0,0,cond_string,cond_values,ierr,'diemac',dt%diemac,1,0.01_dp,iout)

!  dmatpuopt
   if (dt%usepawu==1.or.dt%usepawu==2.or.dt%usepawu==3.or.dt%usepawu==10) then
     cond_string(1)='usepawu' ; cond_values(1)=1
     call chkint_eq(0,1,cond_string,cond_values,ierr,'dmatpuopt',dt%dmatpuopt,10,(/1,2,3,4,5,6,7,8,9,10/),iout)
   end if

!  dmatudiag
   if (dt%usepawu==1.or.dt%usepawu==2.or.dt%usepawu==3.or.dt%usepawu==10) then
     cond_string(1)='usepawu' ; cond_values(1)=1
     call chkint_eq(0,1,cond_string,cond_values,ierr,'dmatudiag',dt%dmatudiag,3,(/0,1,2/),iout)
   end if


!  dmftbandi, dmftbandf
   if (dt%usedmft>0) then
     call chkint_eq(0,1,cond_string,cond_values,ierr,'dmftcheck',dt%dmftcheck,4,(/-1,0,1,2/),iout)
     if(dt%dmftcheck/=-1) then
       cond_string(1)='usedmft' ; cond_values(1)=1
       cond_string(2)='mband' ; cond_values(1)=dt%mband
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftbandi',dt%dmftbandi,1,iout)
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftbandf',dt%dmftbandf,1,iout)
       call chkint_le(0,2,cond_string,cond_values,ierr,'dmftbandi',dt%dmftbandi,dt%mband,iout)
       call chkint_le(0,2,cond_string,cond_values,ierr,'dmftbandf',dt%dmftbandf,dt%mband,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_iter',dt%dmft_iter,0,iout)
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_nwlo',dt%dmft_nwlo,1,iout)
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_nwli',dt%dmft_nwli,1,iout)
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_rslf',dt%dmft_rslf,2,(/0,1/),iout)
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_mxsf',dt%dmft_mxsf,1,zero,iout)
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_mxsf',dt%dmft_mxsf,-1,one,iout)
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_solv',dt%dmft_solv,4,(/-1,0,1,2/),iout)
     end if
   end if

!  dosdeltae
   call chkdpr(0,0,cond_string,cond_values,ierr,'dosdeltae',dt%dosdeltae,1,0.0_dp,iout)

!  dynimage between 0 and 1
   maxidyn=maxval(dt%dynimage(:))
   minidyn=minval(dt%dynimage(:))
   call chkint_ge(0,0,cond_string,cond_values,ierr,'dynimage',minidyn,0,iout)
   call chkint_le(0,0,cond_string,cond_values,ierr,'dynimage',maxidyn,1,iout)

!  ecut
!  With planewaves, one must use positive ecut
   if(usewvl==0)then
     if (abs(dt%ecut+1._dp)<tol8) then
       write(message, '(6a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '   The input keyword "ecut" is compulsory !',ch10,&
&       '  Action : add a value for "ecut" in the input file.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     else
       cond_string(1)='usewvl' ; cond_values(1)=usewvl
       call chkdpr(1,1,cond_string,cond_values,ierr,'ecut',dt%ecut,1,tol8,iout)
     end if
   end if

!  pawecutdg (placed here to stop before ngfftdg)
   if (usepaw==1) then
     call chkdpr(1,0,cond_string,cond_values,ierr,'pawecutdg',dt%pawecutdg,1,tol8,iout)
     cond_string(1)='ecut' ; cond_values(1)=dt%ecut
     call chkdpr(1,1,cond_string,cond_values,ierr,'pawecutdg',dt%pawecutdg,1,dt%ecut,iout)
   end if

!  ecuteps
   if( ANY(optdriver==(/RUNL_SCREENING,RUNL_SCGW/)) )then
     call chkdpr(0,0,cond_string,cond_values,ierr,'ecuteps',dt%ecuteps,1,0.0_dp,iout)
     if(dt%fftgw<20 .and. dt%fftgw/=0)then
       if(dt%ecutwfn<dt%ecuteps-tol8)then
         write(message,'(4a,es16.6,a,es16.6,a,6a)')ch10,&
&         ' chkinp : ERROR -',ch10,&
&         '  The values of ecutwfn and ecuteps are ', dt%ecutwfn,' and ',dt%ecuteps,ch10,&
&         '  With fftgw lower than 20, one expect ecuteps to be smaller or equal to ecutwfn.',ch10,&
&         '  Indeed, one is wasting memory without gaining CPU time or accuracy.',ch10,&
&         '  Action : use another value of fftgw (e.g. 21), or adjust ecutwfn with ecuteps.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end if
   end if

!  ecutsigx
!  @MG FIXME reinstate this check, after having rewritten FFT treatment in GW
   if( ANY( optdriver==(/RUNL_SIGMA,RUNL_SCGW/) ) .and..FALSE.)then
     call chkdpr(0,0,cond_string,cond_values,ierr,'ecutsigx',dt%ecutsigx,1,0.0_dp,iout)
     if(dt%fftgw<20)then
       if(dt%ecutwfn<dt%ecutsigx-tol8)then
         write(message,'(4a,es16.6,a,es16.6,a,6a)')ch10,&
&         ' chkinp : ERROR -',ch10,&
&         '  The values of ecutwfn and ecutsigx are ', dt%ecutwfn,' and ',dt%ecutsigx,ch10,&
&         '  With fftgw lower than 20, one expect ecutsigx to be smaller or equal to ecutwfn.',ch10,&
&         '  Indeed, one is wasting memory without gaining CPU time or accuracy.',ch10,&
&         '  Action : use another value of fftgw (e.g. 21), or adjust ecutwfn with ecutsigx.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end if
   end if

   if ( optdriver==RUNL_SCGW) then ! All *eps, *wfn, *sigx variables must be defined.
     if ( ALL( (/dt%npwwfn,dt%nshwfn/) <tol6 ) .and. dt%ecutwfn<tol6 ) then
       msg = ' Only one of the three variables ecutwfn, npwwfn, or nshwfn can be non-zero.'
       MSG_WARNING(msg)
       ierr=ierr+1
     end if
     if ( ALL( (/dt%npweps,dt%nsheps/) <tol6 ).and. dt%ecuteps<tol6 ) then
       msg = ' Only one of the three variables ecuteps, npweps, or nsheps can be non-zero.'
       MSG_WARNING(msg)
       ierr=ierr+1
     end if
     if ( ALL( (/dt%npwsigx,dt%nshsigx/) <tol6 ) .and. dt%ecutsigx<tol6 ) then
       msg = ' Only one of the three variables ecutsigx, npwsigx, or nshsigx can be non-zero.'
       MSG_WARNING(msg)
       ierr=ierr+1
     end if
   end if

   if ( optdriver==RUNL_BSE) then ! Check for BSE calculations that are not implemented.
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nsppol' ,dt%nsppol, 1,(/1/),iout)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',dt%nspinor,1,(/1/),iout)
   end if

   if ( ANY(optdriver==(/RUNL_SCREENING, RUNL_SIGMA, RUNL_SCGW/)) ) then 
!    
!    Check for GW calculations that are not implemented.
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',dt%nspinor,1,(/1/),iout)
!    
!    Avoid wasting CPUs if nsppol==2.
     if (dt%nsppol==2.and..not.iseven(mpi_enreg%nproc).and.mpi_enreg%nproc>1) then
       write(msg,'(3a)') "Spin-polarized GW calculations should be run with an even number of processors ",ch10,&
&       " for achieving an optimal distribution of memory and CPU load. Change the number of processors."
       MSG_WARNING(msg)
       ierr=ierr+1
     end if
   end if

!  ecutsm
   call chkdpr(0,0,cond_string,cond_values,ierr,'ecutsm',dt%ecutsm,1,0.0_dp,iout)
!  With non-zero optcell, one must use non-zero ecutsm
   if(dt%optcell/=0 )then
     cond_string(1)='optcell' ; cond_values(1)=dt%optcell
     call chkdpr(1,1,cond_string,cond_values,ierr,'ecutsm',dt%ecutsm,1,tol8,iout)
   end if

!  ecutwfn <= ecut. This is also needed for the correct evaluation
!  of the Kleynman-Bylander form factors as the spline in Psps% is done with ecut
!  while we need |q+G| up to ecut. enlargement due to the q is already
!  taken into account by enlarging the spline mesh by around 20%.
   if ( ANY(optdriver==(/RUNL_SCREENING,RUNL_SIGMA,RUNL_BSE/)) ) then
     call chkdpr(0,0,cond_string,cond_values,ierr,'ecutwfn',dt%ecuteps,1,0.0_dp,iout)
     if(dt%ecut<dt%ecutwfn-tol8)then
       write(message,'(4a,es16.6,a,es16.6,a,6a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The values of ecut and ecutwfn are ', dt%ecut,' and ',dt%ecutwfn,ch10,&
&       '  One expects ecutwfn to be smaller or equal to ecut.',ch10,&
&       '  Action : adjust ecutwfn with ecut.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  enable_mpi_io
   if(dt%accesswff==IO_MODE_MPI) then
     cond_string(1)='accesswff' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'enable_mpi_io',mpi_enreg%paral_compil_mpio,1,(/1/),iout)
   end if


!  exchmix
   call chkdpr(0,0,cond_string,cond_values,ierr,'exchmix',dt%exchmix,1,0.0_dp,iout)

!  fftgw
   call chkint_eq(0,0,cond_string,cond_values,ierr,'fftgw',dt%fftgw,8,(/00,01,10,11,20,21,30,31/),iout)

!  fixmom
   if(nsppol==1)then
     cond_string(1)='nsppol' ; cond_values(1)=1
     call chkdpr(1,1,cond_string,cond_values,ierr,'fixmom',dt%fixmom,0,-99.99_dp,iout)
   end if
   if(optdriver==RUNL_RESPFN)then
     cond_string(1)='optdriver' ; cond_values(1)=1
     call chkdpr(1,1,cond_string,cond_values,ierr,'fixmom',dt%fixmom,0,-99.99_dp,iout)
   end if
   if(optdriver==RUNL_SUSCEP)then
     cond_string(1)='optdriver' ; cond_values(1)=2
     call chkdpr(1,1,cond_string,cond_values,ierr,'fixmom',dt%fixmom,0,-99.99_dp,iout)
   end if
   if(dt%prtdos==1)then
     cond_string(1)='prtdos' ; cond_values(1)=1
     call chkdpr(1,1,cond_string,cond_values,ierr,'fixmom',dt%fixmom,0,-99.99_dp,iout)
   end if

!  frzfermi
   call chkint_eq(0,0,cond_string,cond_values,ierr,'frzfermi',dt%frzfermi,2,(/0,1/),iout)

!  fxcartfactor
   call chkdpr(0,0,cond_string,cond_values,ierr,'fxcartfactor',dt%fxcartfactor,1,zero,iout)

!  getxred
   if(dt%getxcart/=0)then
     cond_string(1)='getxcart' ; cond_values(1)=dt%getxcart
!    Make sure that dt%getxred is 0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'getxred',dt%getxred,1,(/0/),iout)
   end if

!  gw_sctype
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gw_sctype',dt%gw_sctype,&
&   4,(/GWSC_one_shot,GWSC_only_W,GWSC_only_G,GWSC_both_G_and_W/),iout)

!  gw_sigxcore
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gw_sigxcore',dt%gw_sigxcore,2,(/0,1/),iout)

!  gwcomp
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gwcomp',dt%gwcomp,2,(/0,1/),iout)

   if (dt%gwcomp/=0) then
     if (optdriver==RUNL_SCREENING .and. ( dt%awtr /=1 .or. dt%spmeth /=0 )) then
       write(message,'(6a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  When gwcomp/=0, the Adler-Wiser formula with time-reversal should be used',ch10,&
&       '  Action : set awtr to 1 or/and spmeth to 0'
       call wrtout(iout,message,'COLL')
       call leave_new('COLL')
     end if

!    Extrapolar trick with HF, SEX and COHSEX is meaningless for Sigma
     if(optdriver==RUNL_SIGMA) then
       mod10=MOD(dt%gwcalctyp,10)
       if ( ANY(mod10 == (/SIG_HF, SIG_SEX, SIG_COHSEX/)) ) then
         write(message,'(6a)' )ch10,&
         ' chkinp: ERROR -',ch10,&
         '  gwcomp/=0, is meaningless in the case of HF, SEX or COHSEX calculations. ',ch10,&
         '  Action : set gwcomp to 0 or change gwcalctyp'
         call wrtout(iout,message,'COLL')
         call leave_new('COLL')
       end if
     end if
   end if

!  gwmem
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gwmem',dt%gwmem,4,(/0,1,10,11/),iout)

!  gwpara
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gwpara',dt%gwpara,3,(/0,1,2/),iout)

!  gwrpacorr
   if(dt%gwrpacorr>0) then
     mod10=MOD(dt%gwcalctyp,10)
     if( optdriver /= RUNL_SCREENING ) then
       write(message,'(6a)' )ch10,&
       ' chkinp: ERROR -',ch10,&
       '  gwrpacorr>0 can only be used when calculating the screening',ch10,&
       '  Action : set gwrpacorr to 0 or optdriver to 3'
       call wrtout(iout,message,'COLL')
       call leave_new('COLL')
     end if
     if( mod10 /= SIG_GW_AC ) then
       write(message,'(6a)' )ch10,&
       ' chkinp: ERROR -',ch10,&
       '  gwrpacorr>0 can only be used with purely imaginary frequencies',ch10,&
       '  Action : set gwrpacorr to 0 or change gwcalctyp'
       call wrtout(iout,message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  iatsph between 1 and natom
   maxiatsph=maxval(dt%iatsph(1:dt%natsph))
   miniatsph=minval(dt%iatsph(1:dt%natsph))
   call chkint_ge(0,0,cond_string,cond_values,ierr,'iatsph',miniatsph,1,iout)
   call chkint_le(0,0,cond_string,cond_values,ierr,'iatsph',maxiatsph,natom,iout)

!  icoulomb
   call chkint_eq(0,0,cond_string,cond_values,ierr,'icoulomb',dt%icoulomb,2,(/0,1/),iout)
   if (dt%nspden > 2) then
     cond_string(1)='nspden' ; cond_values(1)=nspden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'icoulomb',dt%icoulomb,1,(/0/),iout)
   end if

!  ieig2rf
   if(optdriver==RUNL_RESPFN.and.usepaw==1)then
     cond_string(1)='optdriver' ; cond_values(1)=1
     cond_string(2)='usepaw'    ; cond_values(2)=usepaw
     call chkint_eq(1,2,cond_string,cond_values,ierr,'ieig2rf',dt%ieig2rf,1,(/0/),iout)
   end if

!  iextrapwf
   call chkint_eq(0,0,cond_string,cond_values,ierr,'iextrapwf',dt%iextrapwf,2,(/0,1/),iout)
   if (dt%mkmem==0) then
     cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iextrapwf',dt%iextrapwf,2,(/0,1/),iout)
   end if

!  imgmov
   call chkint_eq(0,0,cond_string,cond_values,ierr,'imgmov',dt%imgmov,5,(/0,1,2,9,13/),iout)

!  intxc
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
!    Make sure that dt%intxc is 0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'intxc',dt%intxc,1,(/0/),iout)
   end if
!  TEMPORARY
   if(optdriver==RUNL_RESPFN)then ! Make sure that dt%intxc is 0
     cond_string(1)='optdriver' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'intxc',dt%intxc,1,(/0/),iout)
   end if

!  ionmov
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ionmov',&
&   dt%ionmov,26,&
&   (/0,1,2,3,4,5,6,7,8,9,10,12,13,14,20,30,70,80,90,91,92,95,96,97,98,99/),&
&   iout)
!  When optcell/=0, ionmov must be 2, 3 or 13
   if(dt%optcell/=0)then
     cond_string(1)='optcell' ; cond_values(1)=dt%optcell
!    Make sure that ionmov==2, 3 or 13
     call chkint_eq(1,1,cond_string,cond_values,ierr,'ionmov',dt%ionmov,5,(/2,3,13,97,98/),iout)
   end if

!  iprcch
!  if (abs(dt%iprcch)>=1.and.dt%iscf>=10.and.nspden==4) then
!  write(message,'(8a)')ch10,&
!  &   ' chkinp : ERROR -',ch10,&
!  &   '  When non-collinear magnetism is activated (nspden=4),',ch10,&
!  &   '  iprcch/=0 is not compatible with SCF mixing on density (iscf>=10) !',ch10,&
!  &   '  Action : choose SCF mixing on potential (iscf<10) or change iprcch value.'
!  call wrtout(std_out,message,'COLL')
!  ierr=ierr+1
!  end if
   if(dt%iextrapwf==1)then
     cond_string(1)='iextrapwf' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iprcch',dt%iprcch,2,(/5,6/),iout)
   end if
   if (dt%iprcch<0.and.mod(dt%iprcel,100)>=61.and.(dt%iprcel<71.or.dt%iprcel>79)) then
     cond_string(1)='iprcel';cond_values(1)=dt%iprcel
     call chkint_ge(0,2,cond_string,cond_values,ierr,'iprcch',dt%iprcch,0,iout)
   end if

!  iprcel
   call chkint(0,0,cond_string,cond_values,ierr,'iprcel',dt%iprcel,1,(/0/),1,21,iout)   !  0 or superior to 21
   if(nsppol==2 .and. (dt%occopt>=3 .and. dt%occopt<=7).and.mod(dt%iprcel,10)>49 )then
     write(message,'(8a)')ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  For spin-polarized metallic systems (occopt>3),',ch10,&
&     '  only RPA dielectric matrix can be evaluated) !',ch10,&
&     '  Action : change iprcel value in input file (mod(iprcel,100)<50) !'
     call wrtout(std_out,message,'COLL')
     ierr=ierr+1
   end if

!  iprctfvw
   call chkint_eq(0,0,cond_string,cond_values,ierr,'iprctfvw',dt%iprctfvw,4,(/0,1,2,3/),iout)
   if((nsppol/=1).or.(dt%iscf > 9))then
     cond_string(1)='nsppol' ; cond_values(1)=nsppol
     cond_string(2)='iscf' ; cond_values(2)=dt%iscf
     call chkint_eq(2,2,cond_string,cond_values,ierr,'iprctfvw',dt%iprctfvw,1,(/0/),iout)
   end if


!  iscf
   call chkint_eq(0,0,cond_string,cond_values,ierr,&
&   'iscf',dt%iscf,18,(/-3,-2,-1,1,2,3,4,5,6,7,11,12,13,14,15,16,17,22/),iout)
!  If ionmov==4, iscf must be 2, 12, 5 or 6.
   if(dt%ionmov==4)then
     cond_string(1)='ionmov' ; cond_values(1)=4
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,4,(/2,12,5,6/),iout)
   end if
!  If PAW, iscf cannot be -1, 11
   if (usepaw==1) then
     cond_string(1)='PAW' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,11,(/-3,-2,2,3,4,7,12,13,14,17,22/),iout)
   end if
!  Mixing on density is only allowed for GS calculations
!  or for GW calculations where it is not used.
   if(optdriver/=RUNL_GSTATE .and. ALL(optdriver/=(/RUNL_SCREENING,RUNL_SIGMA,RUNL_BSE,RUNL_SCGW/)) ) then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_le(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,9,iout)
   end if
!  mixing on density is not allowed with some preconditioners
   if (dt%iprctfvw /= 0) then
     cond_string(1)='iprctfvw' ; cond_values(1)=dt%iprctfvw
     call chkint_le(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,9,iout)
   end if
!  When pawoptmix=1 and nspden=4, iscf must be >=10
   if(dt%pawoptmix/=0.and.nspden==4)then
     cond_string(1)='nspden'    ; cond_values(1)=nspden
     cond_string(2)='pawoptmix' ; cond_values(2)=dt%pawoptmix
     call chkint_ge(2,2,cond_string,cond_values,ierr,'iscf',dt%iscf,10,iout)
   end if

!  istwfk
   if(response==1 .and. maxval( abs(dt%istwfk(1:nkpt)-1) ) >0)then
!    Force istwfk to be 1 for RF calculations
!    Other choices cannot be realized yet, because of the ddk perturbation.
     write(message,'(8a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  When response==1, all the components of istwfk must be 1.',ch10,&
&     '  Not yet programmed for time-reversal symmetry.',ch10,&
&     '  Action : set istwfk to 1 for all k-points'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if(dt%nbandkss/=0 .and. dt%kssform/=3 .and. maxval( abs(dt%istwfk(1:nkpt)-1) ) >0)then
     write(message,'(8a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  When nbandkss/=0 and kssform/=3 all the components of istwfk must be 1.',ch10,&
&     '  Not yet programmed for time-reversal symmetry.',ch10,&
&     '  Action : set istwfk to 1 for all k-points'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if(dt%berryopt/=0 .and. maxval(dt%istwfk(:))/=1)then
     write(message,'(8a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  When berryopt/=0, all the components of istwfk must be 1.',ch10,&
&     '  Not yet programmed for time-reversal symmetry.',ch10,&
&     '  Action : set istwfk to 1 for all k-points'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if ((dt%wfoptalg==4.or.dt%wfoptalg==14).and.maxval(dt%istwfk(:)-2)>0) then
     write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  Only the gamma point can use time-reversal and wfoptalg=4 or 14',ch10,&
&     '  Action : put istwfk to 1 or remove k points with half integer coordinates ',ch10,&
&     '  Also contact ABINIT group to say that you need that option.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if ((dt%wfoptalg==4.or.dt%wfoptalg==14).and.minval(dt%istwfk(:)-2)==0.and.fftalg/=401) then
     write(message, '(4a,i3,a,a,a)' ) ch10,&
&     ' chkinp :  ERROR -',ch10,&
&     ' For istwfk=2, the value fftalg= ',fftalg, &
&     ' is not allowed in case of wfoptalg=4 or 14 !', ch10,&
     ' Change if to fftalg=401.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if (dt%accesswff == IO_MODE_ETSF) then    ! ETSF_IO, current limitation to istwfk 1
     do ikpt=1,nkpt
       if(dt%istwfk(ikpt)/=1)then
         write(message,'(8a,I0,a)' ) ch10,&
&         ' chkinp: ERROR -',ch10,&
&         '  When accesswff==3, all the components of istwfk must be 1.',ch10,&
&         '  Not yet programmed for time-reversal symmetry.',ch10,&
&         ' Action : modify value of istwfk to be 1 for all k points in input file.'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,message,'COLL')
!        The following line was wrong : istwfk is used previously, to dimension the
!        number of plane waves, so it cannot be modified here ...
!        (And moreover, no input variable should be modified in chkinp, that is
!        supposed only to do the checking).
!        dt%istwfk(:) = 1
         ierr=ierr+1
         exit
       end if
     end do
   end if

!  ixc
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'ixc',dt%ixc,25,(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,20,21,22,23,24,26,27/),-1,0,iout) ! One of the values, or negative
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
!    Make sure that ixc is 1, 7, 8, 20, 21 or 22 or negative
     call chkint(1,1,cond_string,cond_values,ierr,&
&     'ixc',dt%ixc,6,(/1,7,8,20,21,22/),-1,0,iout)
   end if
   if(response==1)then
     cond_string(1)='response' ; cond_values(1)=1
!    Make sure that ixc is between 0 and 9, or 11, 12, 14, 15, 23 or 24 or negative
     call chkint(1,1,cond_string,cond_values,ierr,&
&     'ixc',dt%ixc,16,(/0,1,2,3,4,5,6,7,8,9,11,12,14,15,23,24/),-1,0,iout)
   end if
   if(nspden/=1)then
     cond_string(1)='nspden' ; cond_values(1)=nspden
!    Make sure that ixc is 0, 1 , the gga, or Fermi-Amaldi, or negative
     call chkint(1,1,cond_string,cond_values,ierr,&
&     'ixc',dt%ixc,17,(/0,1,7,8,9,11,12,13,14,15,16,17,20,23,24,26,27/),-1,0,iout)
   end if

!  ixcpositron
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ixcpositron',dt%ixcpositron,7,(/0,1,11,2,3,31,4/),iout)

!  kptnrm and kpt
!  Coordinates components must be between -1 and 1.
   if(dt%kptnrm<1.0-1.0d-10)then
     write(message, '(a,a,a,a,es22.14,a,a,a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  The input variable kptnrm is',dt%kptnrm,' while it must be >=1.0_dp.',&
&     ch10,'  Action : change the input variable kptnrm.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
   do ikpt=1,nkpt
     do mu=1,3
       if ( abs(dt%kpt(mu,ikpt))> dt%kptnrm*1.0000001_dp ) then
         write(message, '(a,a,a,a,i5,a,a,a,a,3es22.14,a,a,a,a)' ) ch10,&
&         ' chkinp: ERROR -',ch10,&
&         '  For k point number',ikpt,'  the reduced coordinates',ch10,&
&         '  generated by the input variables kpt and kptnrm are',ch10,&
&         dt%kpt(1,ikpt)/dt%kptnrm,dt%kpt(2,ikpt)/dt%kptnrm,dt%kpt(3,ikpt)/dt%kptnrm,ch10,&
&         '  while they must be between -1.0_dp and 1.0_dp (included).',ch10,&
&         '  Action : check kpt and kptnrm in the input file.'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end do

!  jellslab
   call chkint_eq(0,0,cond_string,cond_values,ierr,'jellslab',dt%jellslab,2,(/0,1/),iout)

   if (dt%jellslab==1) then
     if(dt%nimage>1)then
       cond_string(1)='nimage' ; cond_values(1)=dt%nimage
       call chkint_eq(1,1,cond_string,cond_values,ierr,'jellslab',dt%jellslab,1,(/0/),iout)
     end if
!    slabwsrad must be positive
     cond_string(1)='jellslab' ; cond_values(1)=dt%jellslab
     call chkdpr(1,0,cond_string,cond_values,ierr,'slabwsrad',dt%slabwsrad,1,zero,iout)
!    slabzbeg must be positive
     call chkdpr(1,0,cond_string,cond_values,ierr,'slabzbeg',dt%slabzbeg,1,zero,iout)
!    slabzend must be bigger than slabzbeg
     call chkdpr(1,0,cond_string,cond_values,ierr,'slabzend',dt%slabzend,1,dt%slabzbeg,iout)
!    rprimd(3,3) must be bigger than slabzend
     call chkdpr(1,0,cond_string,cond_values,ierr,'rprimd33',rprimd(3,3),1,dt%slabzend,iout)
!    Third real space primitive translation has to be orthogonal to the other ones,
!    actually, for convenience it is useful that rprimd is something like:
!    a  b  0
!    c  d  0
!    0  0  e
     if(abs(rprimd(1,3))+abs(rprimd(2,3))+abs(rprimd(3,1))+abs(rprimd(3,2))>tol12) then
       write(message, '(a,a,a,a,a,a)' ) ch10,&
&       ' chkinp : ERROR - ',ch10,&
&       '  Third real space vector is not orthogonal to the other ones,',ch10,&
&       '  this is needed to use jellium'
       call wrtout(std_out,message,'COLL')
       ierr=ierr+1
     end if

!    Atoms have to be placed in the vacuum space
     do iatom=1,natom
       zatom=(dt%xred_orig(3,iatom,intimage)-anint(dt%xred_orig(3,iatom,intimage)-half+tol6))*rprimd(3,3)
       if(abs(zatom-dt%slabzbeg)<tol8 .or. abs(zatom-dt%slabzend)<tol8) then
         if(dt%znucl(dt%typat(iatom))>tol6) then
           write(message,'(4a,i5,a)') ch10,&
&           ' chkinp : WARNING -',ch10,&
&           ' atom number=',iatom,' lies precisely on the jellium edge !'
           call wrtout(std_out,message,'COLL')
         end if
         cycle
       end if
       if(zatom>dt%slabzbeg .and. zatom<dt%slabzend) then
         write(message, '(a,a,a,a,i5,a)' ) ch10,&
&         ' chkinp : ERROR - ',ch10,&
&         ' atom number=',iatom,' is inside the jellium slab.'
         call wrtout(std_out,message,'COLL')
         ierr=ierr+1
       end if
     end do
   end if

!  kssform
   call chkint_eq(0,0,cond_string,cond_values,ierr,'kssform',dt%kssform,3,(/0,1,3/),iout)

   if (dt%kssform/=0 .and. dt%nbandkss/=0) then ! Check for outkss limitations.
     call wrtout(std_out," Checking if input is consistent with KSS generation",'COLL')
     call chkint_eq(0,0,cond_string,cond_values,ierr,'paral_kgb',dt%paral_kgb,1,(/0/),iout)
     call chkint_eq(0,0,cond_string,cond_values,ierr,'accesswff',dt%accesswff,2,(/IO_MODE_FORTRAN,IO_MODE_ETSF/),iout)
   end if

!  localrdwf
   call chkint_eq(0,0,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,2,(/0,1/),iout)
   if(dt%mkmem==0)then
     cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
     call chkint_eq(1,1,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,1,(/1/),iout)
   end if
   if(dt%mkqmem==0)then
     cond_string(1)='mkqmem' ; cond_values(1)=dt%mkqmem
     call chkint_eq(1,1,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,1,(/1/),iout)
   end if
   if(dt%mk1mem==0)then
     cond_string(1)='mk1mem' ; cond_values(1)=dt%mk1mem
     call chkint_eq(1,1,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,1,(/1/),iout)
   end if

!  macro_uj
   if(dt%macro_uj/=0) then
     if (dt%ionmov/=0) then
       write(message, '(6a,i2,2a,i2,3a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  Determination of U can not be combined with ionic movements.',ch10,&
&       '  Here  ionmov= ',dt%ionmov,ch10,&
&       '  and macro_uj=',dt%macro_uj,'.',ch10,&
&       '  Action: change ionmov in input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     else if (dt%nstep<3) then
       write(message, '(6a,i1,2a,i2,3a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  Determination of U needs at least 3 scf steps:',ch10,&
&       '        nstep = ',dt%nstep,ch10,&
&       '  and macro_uj=',dt%macro_uj,'.',ch10,&
&       '  Action: increase nstep in input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  mffmem
   call chkint_eq(0,0,cond_string,cond_values,ierr,'mffmem',dt%mffmem,2,(/0,1/),iout)

!  mixalch
!  For each type of atom, the sum of the psp components
!  must be one.
   if(dt%ntypalch>0)then
     do itypat=1,dt%ntypalch
       sumalch=sum(dt%mixalch(:,itypat))
       if(abs(sumalch-one)>tol10)then
         if(dt%npspalch<=6)then
           write(message, '(2a,6es12.4)' ) ch10,&
&           ' chkinp : mixalch(:,itypat)=',dt%mixalch(:,itypat)
         end if
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         write(message, '(4a,i4,2a,f8.2,4a)' ) ch10,&
&         ' chkinp: ERROR -',ch10,&
&         '  For the alchemical atom number',itypat,ch10,&
&         '  the sum of the pseudopotential coefficients is',sumalch,ch10,&
&         '  while it should be one.',ch10,&
&         '  Action : check the content of the input variable mixalch.'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if

!  natom
   if(dt%prtgeo>0)then
     cond_string(1)='prtgeo' ; cond_values(1)=dt%prtgeo
     call chkint_le(1,1,cond_string,cond_values,ierr,'natom',natom,9999,iout)
   end if

!  nband
!  Make sure all nband(nkpt) are > 0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       if (dt%nband(ikpt+(isppol-1)*nkpt)<=0) then
         cond_string(1)='ikpt' ; cond_values(1)=ikpt
         cond_string(2)='isppol' ; cond_values(2)=isppol
         call chkint_ge(0,2,cond_string,cond_values,ierr,'nband',dt%nband(ikpt+(isppol-1)*nkpt),1,iout)
       end if
     end do
   end do

   if(mpi_enreg%nproc/=1 .and. nsppol==2 .and. dt%usewvl == 0)then
     do ikpt=1,nkpt
       if (dt%nband(ikpt)/=dt%nband(ikpt+nkpt)) then
         write(message, '(8a,i4,a,2i5,a)' ) ch10,&
&         ' chkinp : in the parallel k point case, for each k point,',ch10,&
&         '  the number of bands in the spin up case must be equal to',ch10,&
&         '  the number of bands in the spin down case.',ch10,&
&         '  This is not the case for the k point number :',ikpt,&
&         '  The number of bands spin up and down are :',dt%nband(ikpt),dt%nband(ikpt+nkpt),&
&         '  Action : change nband, or use the sequential version of ABINIT.'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if

!  nbandkss
!  Must be greater or equal to -1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nbandkss',dt%nbandkss,-1,iout)
!  When ionmov/=0
   if(dt%ionmov/=0 .and. dt%nbandkss/=0)then
     write(message,'(14a)') ch10,&
&     ' chkinp: WARNING -',ch10,&
&     '  Ions (or cell) are allowed to move (ionmov/=0),',ch10,&
&     '  and a _KSS file is requested (nbandkss/=0).',ch10,&
&     '  A _KSS file will be created at each geometry-optimisation step.',ch10,&
&     '  Note that this is time consuming !',ch10,&
&     '  Action : use datasets (one for geometry optimisation,',ch10,&
&     '           one for states output).'
     call wrtout(std_out,message,'COLL')
   end if

!  nbdblock
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,iout)
!  When wfoptalg==0, nbdblock must be 1
   if(mod(dt%wfoptalg,10)==0)then
     cond_string(1)='wfoptalg' ; cond_values(1)=0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,(/1/),iout)
   end if
!  When wfoptalg==2, nbdblock must be 1
   if(dt%wfoptalg==2)then
     cond_string(1)='wfoptalg' ; cond_values(1)=2
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,(/1/),iout)
   end if
!  When wfoptalg==3, nbdblock must be 1, and iscf must be -2
   if(dt%wfoptalg==3)then
     cond_string(1)='wfoptalg' ; cond_values(1)=3
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,(/1/),iout)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,1,(/-2/),iout)
   end if
!  When wfoptalg==4, nbdblock must be a divisor of nband
   if(mod(dt%wfoptalg,10)==4)then
     do isppol=1,nsppol
       do ikpt=1,nkpt
         if(mod(dt%nband(ikpt+(isppol-1)*nkpt),dt%nbdblock)/=0) then
           write(message, '(8a)' ) ch10,&
&           ' chkinp: ERROR -',ch10,&
&           '  For the moment, when wfoptalg=4,',ch10,&
&           '  nband must be a multiple of nbdblock.',ch10,&
&           '  Action : check the value of the input variable nbdblock.'
           call wrtout(iout,message,'COLL')
           call wrtout(std_out,  message,'COLL')
           call leave_new('COLL')
         end if
       end do
     end do
   end if

!  nberry
!  must be between 0 and 20
   if(dt%berryopt/=0)then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'nberry',dt%nberry,0,iout)
     call chkint_le(0,0,cond_string,cond_values,ierr,'nberry',dt%nberry,20,iout)
     if(mpi_enreg%paral_compil==1)then
!      MPI Parallel case
       if ((dt%nberry/=0).and.(dt%berryopt > 0).and.(dt%berryopt /= 4).and.(dt%berryopt /= 5)) then
         write(message,'(a,a,a,a,a,a,a,a,i4,a,a,a)')ch10,&
&         ' chkinp : ERROR -',ch10,&
&         '  Berry phase calculation of polarisation with positive berryopt is not',ch10,&
&         '  allowed in the parallel version of ABINIT.',ch10,&
&         '  So, the value of nberry=',dt%nberry,' is not allowed,',ch10,&
&         '  Action : change berryopt to negative values or change nberry, or use the sequential version.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end if
   end if

!  ndynimage
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'ndynimage',dt%ndynimage,1,iout)

!  nfft and nfftdg
!  Must have nfft<=nfftdg
   if (usepaw==1) then
     nfft  =dt%ngfft(1)  *dt%ngfft(2)  *dt%ngfft(3)
     nfftdg=dt%ngfftdg(1)*dt%ngfftdg(2)*dt%ngfftdg(3)
     cond_string(1)='nfft' ; cond_values(1)=nfft
     call chkint(1,1,cond_string,cond_values,ierr,'nfftdg',nfftdg,1,(/0/),1,nfft,iout) ! Must be 0 or nfft
   end if

!  diismemory
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'diismemory',dt%diismemory,1,iout)

!  nimage
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nimage',dt%nimage,1,iout)
   if (usewvl==1) then
     cond_string(1)='usewvl' ; cond_values(1)=usewvl
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nimage',dt%nimage,1,(/1/),iout)
   end if
   if (optdriver/=RUNL_GSTATE) then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nimage',dt%nimage,1,(/1/),iout)
   end if

!  nkpt
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nkpt',nkpt,1,iout)
!  If prtdos==2 or 3, must be greater or equal to 2
   if(dt%prtdos==2 .or. dt%prtdos==3)then
     cond_string(1)='prtdos' ; cond_values(1)=dt%prtdos
     call chkint_ge(1,1,cond_string,cond_values,ierr,'nkpt',nkpt,2,iout)
   end if
!  Must be smaller than 50 if iscf=-2 (band structure)
!  while prteig=0 and prtvol<2, except if kptopt>0
   if(dt%iscf==-2 .and. dt%prteig==0 .and. dt%prtvol<2 .and. dt%kptopt<=0)then
     cond_string(1)='iscf'   ; cond_values(1)=dt%iscf
     cond_string(2)='prteig' ; cond_values(2)=dt%prteig
     cond_string(3)='prtvol' ; cond_values(3)=dt%prtvol
     call chkint_le(1,3,cond_string,cond_values,ierr,'nkpt',nkpt,50,iout)
   end if

!  npband
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npband',dt%npband,1,iout)

!  npfft
!  Must be greateror equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npfft',dt%npfft,1,iout)
!  If usepaw==1 and pawmixdg==0, npfft must be equal to 1
   if(usepaw==1 .and. dt%pawmixdg==0)then
     cond_string(1)='usepaw  ' ; cond_values(1)=usepaw
     cond_string(2)='pawmixdg' ; cond_values(2)=dt%pawmixdg
     call chkint_eq(1,2,cond_string,cond_values,ierr,'npfft',dt%npfft,1,(/1/),iout)
   end if

!  npimage
!  Must be greater or equal to 1
!  call chkint_ge(0,0,cond_string,cond_values,ierr,'npfft',dt%npfft,1,iout)
!  At present, parallelism over images is not coded ...
   call chkint_eq(0,0,cond_string,cond_values,ierr,'npimage',dt%npimage,1,(/1/),iout)

!  npkpt
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npkpt',dt%npkpt,1,iout)

!  nproj
!  If there is more than one projector for some angular momentum
!  channel of some pseudopotential
   do ilang=0,3
!    nprojmax(ilang)=maxval(pspheads(1:npsp)%nproj(ilang)) ! Likely problems with HP compiler
     nprojmax(ilang)=pspheads(1)%nproj(ilang)
     if(npsp>=2)then
       do ii=2,npsp
         nprojmax(ilang)=max(pspheads(ii)%nproj(ilang),nprojmax(ilang))
       end do
     end if
   end do

   if (maxval(nprojmax(0:3))>1) then
     if (usepaw==0.and.optdriver==RUNL_SCREENING.and.dt%inclvkb/=0) then
       ierr=ierr+1
       write(message,'(4a)')ch10,&
&       ' inclvkb /= 0 not implemented for pseudos with more than one projector per l-channel ',ch10,&
&       ' Use inclvkb == 0 in the input file '
       call wrtout(std_out,message,'COLL')
     end if
   end if

!  npwkss
!  Must be greater or equal to -1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npwkss',dt%npwkss,-1,iout)

!  nqpt
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nqpt',dt%nqpt,2,(/0,1/),iout)

!  nscforder
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nscforder',dt%nscforder,10,(/8,14,16,20,24,30,40,50,60,100/),iout)

!  nspden
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nspden',nspden,3,(/1,2,4/),iout)

   if(nsppol==2)then  !  When nsppol=2, nspden must be 2
     cond_string(1)='nsppol' ; cond_values(1)=2
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspden',nspden,1,(/2/),iout)
   end if
   if(nspden==2 .and. nsppol==1 .and. response==1)then
     write(message,'(16a)')ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  nspden==2 together with nsppol==1 is not allowed',ch10,&
&     '  for response function calculations.',ch10,&
&     '  For antiferromagnetic materials, use nspden==2 and nsppol=2.',ch10,&
&     '  In this case, Shubnikov symmetries will be used to decrease',ch10,&
&     '  the number of perturbations. In a future version, it will also be',ch10,&
&     '  used to decrease the number of spin components (to be coded).',ch10,&
&     '  Action : change nsppol to 1, or check nspden.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if(nspden==4.and.response==1)then
     write(message,'(6a)')ch10,&
&     ' chkinp : ERROR -',ch10,&
&     '  nspden==4 not allowed in response formalism.',ch10,&
&     '  Non collinear magnetism not yet implemented in perturbative treatment.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
!  TR symmetry not allowed for NC magnetism, in the present version
!  (to be investigated further)
   if (nspden==4.and.(dt%kptopt==1.or.dt%kptopt==2)) then
     write(message, '(10a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  When non-collinear magnetism is activated (nspden=4),',ch10,&
&     '  time-reversal symmetry cannot be used in the present',ch10,&
&     '  state of the code (to be checked and validated).',ch10,&
&     '  Action: choose kptopt different from 1 or 2.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     ierr=ierr+1
   end if
!  When iprcch<0 or 3, nspden must be 1 or 2
   if(dt%iprcch<0.or.dt%iprcch==3)then
     cond_string(1)='iprcch' ; cond_values(1)=dt%iprcch
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspden',nspden,2,(/1,2/),iout)
   end if
!  When ionmov=4 and iscf>10, nspden must be 1 or 2
   if(dt%ionmov==4.and.dt%iscf>10)then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     cond_string(1)='iscf' ; cond_values(1)=dt%iprcch
     call chkint_eq(1,2,cond_string,cond_values,ierr,'nspden',nspden,2,(/1,2/),iout)
   end if
!  When iprcel>49, nspden must be 1 or 2
   if(mod(dt%iprcel,100)>49)then
     cond_string(1)='iprcel' ; cond_values(1)=dt%iprcel
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspden',nspden,2,(/1,2/),iout)
   end if

!  nspinor
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nspinor',nspinor,2,(/1,2/),iout)
   if(nspden==2)then !  When nspden=2, nspinor must be 1
     cond_string(1)='nspden' ; cond_values(1)=2
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if

   if(nspden==4)then  !  When nspden=4, nspinor must be 2
     cond_string(1)='nspden' ; cond_values(1)=4
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/2/),iout)
   end if
!  When iscf=-1, nspinor must be 1
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
!    Make sure that nsppol is 1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if
!  spin-orbit is not implemented for the strain perturbation
   if(dt%rfstrs/=0)then
     cond_string(1)='rfstrs' ; cond_values(1)=dt%rfstrs
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if
!  When usepawu=2, nspinor must be 1
   if(dt%usepawu==2)then
     cond_string(1)='usepawu' ; cond_values(1)=2
!    Make sure that nspinor is 1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if

!  nsppol
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nsppol',nsppol,2,(/1,2/),iout)

!  nsym
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nsym',dt%nsym,1,iout)
!  check if nsym=1 in phonon calculation in finite electric field
   if( (response==1) .and. (dt%berryopt==4) ) then
     cond_string(1)='response' ; cond_values(1)=1
     cond_string(2)='berryopt' ; cond_values(2)=4
     call chkint_eq(1,2,cond_string,cond_values,ierr,'nsym',dt%nsym,1,(/1/),iout)
   end if

!  ntime
   call chkint_ge(0,0,cond_string,cond_values,ierr,'ntime',dt%ntime,0,iout)

!  ntimimage
   call chkint_ge(0,0,cond_string,cond_values,ierr,'ntimimage',dt%ntimimage,1,iout)

!  ntypalch
   if (usepaw==1) then
     cond_string(1)='pspcod' ; cond_values(1)=7
     call chkint_eq(1,1,cond_string,cond_values,ierr,'ntypalch',dt%ntypalch,1,(/0/),iout)
   end if

!  occ
!  Do following tests only for occopt==0 or 2, when occupation numbers are needed
   if ((dt%iscf>0.or.dt%iscf==-1.or.dt%iscf==-3) .and. (dt%occopt==0 .or. dt%occopt==2) ) then
!    make sure occupation numbers (occ(n)) were defined:
     sumocc=zero
     bantot=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iband=1,dt%nband(ikpt+(isppol-1)*nkpt)
           bantot=bantot+1
           sumocc=sumocc+dt%occ_orig(bantot)
           if (dt%occ_orig(bantot)<zero) then
             write(message, '(a,a,a,a,2i6,a,e20.10,a,a,a)' )  ch10,&
&             ' chkinp: ERROR -',ch10,&
&             'iband,ikpt=',iband,ikpt,' has negative occ=',dt%occ_orig(bantot),' =>stop',&
&             ch10,'  Action : correct this occupation number in input file.'
             call wrtout(iout,message,'COLL')
             call wrtout(std_out,  message,'COLL')
             call leave_new('COLL')
           end if
         end do
       end do
     end do
     if (sumocc<=1.0d-8) then
       write(message, '(a,a,a,a,1p,e20.10,a,a,a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  Sum of occ=',sumocc, ' =>occ not defined => stop',ch10,&
&       '  Action : correct the array occ in input file.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ierr=ierr+1
     end if
   end if

!  occopt
   call chkint_eq(0,0,cond_string,cond_values,ierr,'occopt',dt%occopt,8,(/0,1,2,3,4,5,6,7/),iout)
!  When prtdos==1, occopt must be between 3 and 7
   if(dt%prtdos==1)then
     write(cond_string(1), "(A)") 'prtdos'
     cond_values(1)=1
!    Make sure that occopt is 3,4,5,6, or 7
     call chkint_eq(1,1,cond_string,cond_values,ierr,'occopt',dt%occopt,5,(/3,4,5,6,7/),iout)
   end if


!  optcell
   call chkint_eq(0,0,cond_string,cond_values,ierr,'optcell',dt%optcell,10,(/0,1,2,3,4,5,6,7,8,9/),iout)
!  With dt%berryopt=4, one must have optcell==0
   if(dt%berryopt==4)then
     cond_string(1)='berryopt' ; cond_values(1)=dt%berryopt
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optcell',dt%optcell,1,(/0/),iout)
   end if

!  optdriver and PAW
   if(usepaw==1)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(1,1,cond_string,cond_values,ierr,&
&     'optdriver',optdriver,6,(/RUNL_GSTATE,RUNL_RESPFN,RUNL_SCREENING,RUNL_SIGMA,RUNL_SCGW,RUNL_BSE/),iout)
   end if
!  Non-linear response calculations
   if(nspinor/=1)then
     cond_string(1)='nspinor' ; cond_values(1)=nspinor
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if
   if(dt%occopt/=1 .and. dt%occopt/=2)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if
   if(dt%mkmem==0)then
     cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if
   if(dt%kptopt/=2 .and. dt%kptopt/=3)then
     cond_string(1)='kptopt' ; cond_values(1)=dt%kptopt
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if
   if ((dt%ixc /= 0).and.(dt%ixc /= 3).and.(dt%ixc /= 7).and.(dt%ixc /= 8)) then
     cond_string(1)='ixc' ; cond_values(1)=dt%ixc
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if

!  optforces
!  When ionmov>0, optforces must be >0
   if(dt%ionmov>0)then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optforces',dt%optforces,2,(/1,2/),iout)
   end if
!  When iscf=22, optforces must be 0 or 2
   if(dt%iscf==22)then
     cond_string(1)='iscf' ; cond_values(1)=dt%iscf
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optforces',dt%optforces,2,(/0,2/),iout)
   end if

!  optstress
!  When optcell>0, optstress must be >0
   if(dt%optcell>0)then
     cond_string(1)='optcell' ; cond_values(1)=dt%optcell
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optstress',dt%optstress,1,(/1/),iout)
   end if

!  paral_kgb
   call chkint_eq(0,0,cond_string,cond_values,ierr,'paral_kgb',dt%paral_kgb,2,(/0,1/),iout)
!  Warning
   if(dt%paral_kgb==1.and.dt%accesswff/=IO_MODE_MPI) then
     write(message,'(14a)' ) ch10,&
&     ' chkinp: WARNING -',ch10,&
&     '  When k-points/bands/FFT parallelism is activated',ch10,&
&     '  (paral_kgb=1), only MPI-IO input/output is allowed !',ch10,&
&     '  accesswff/=1 in your input file',ch10,&
&     '  You will not be able to perform input/output !'
     call wrtout(std_out,message,'COLL')
   end if

!  pawcpxocc
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawcpxocc',dt%pawcpxocc,2,(/1,2/),iout)
     if (dt%usepawu>=1.and.nspinor==2.and.dt%pawcpxocc==1) then
       write(message, '(8a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  When non-collinear magnetism is activated ,',ch10,&
&       '  and LDA+U activated ',ch10,&
&       '  PAW occupancies must be complex !'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ierr=ierr+1
     else if (dt%pawspnorb==1.and.(dt%kptopt==0.or.dt%kptopt>=3).and.dt%pawcpxocc==1) then
       if (optdriver==RUNL_GSTATE.and.dt%iscf<10) then
         write(message, '(14a)' ) ch10,&
&         ' chkinp: ERROR -',ch10,&
&         '  When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&         '  and time-reversal symmetry is broken (kptopt/=1 and kptopt/=2)',ch10,&
&         '  PAW occupancies are complex !',ch10,&
&         '  Their imaginary part is used to evaluate total energy by direct',ch10,&
&         '  scheme, needed here because SCF potential mixing has been chosen (iscf<10).',ch10,&
&         '  Action: put pawcpxocc=2 in input file, or choose SCF density mixing (iscf>=10).'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         ierr=ierr+1
       else if (optdriver==RUNL_GSTATE.and.dt%iscf>=10) then
         write(message, '(14a)' ) ch10,&
&         ' chkinp: WARNING -',ch10,&
&         '  When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&         '  and time-reversal symmetry is broken (kptopt/=1 and kptopt/=2)',ch10,&
&         '  PAW occupancies are complex !',ch10,&
&         '  By setting pawcpxocc=1 in input file, their imaginary part',ch10,&
&         '  is not computed. As a consequence, total energy computed',ch10,&
&         '  is not available. Put pawcpxocc=2 in input file if you want it.'
         call wrtout(std_out,  message,'COLL')
       else
         write(message, '(14a)' ) ch10,&
&         ' chkinp: ERROR -',ch10,&
&         '  When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&         '  and time-reversal symmetry is broken (kptopt/=1 and kptopt/=2)',ch10,&
&         '  PAW occupancies are complex !',ch10,&
&         '  Action: put pawcpxocc=2 in input file to compute their imaginary part.'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         ierr=ierr+1
       end if
     end if
     if (dt%pawspnorb==1.and.dt%kptopt==0) then
       write(message, '(10a)' ) ch10,&
&       ' chkinp: WARNING -',ch10,&
&       '  When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&       '  time-reversal symmetry might be broken.',ch10,&
&       '  Using kptopt=0 might be risky: if (kx,ky,kz) is present in k-points list,',ch10,&
&       '  (-kx,-ky,-kz) (or equivalent) should also be present.'
       call wrtout(std_out,  message,'COLL')
     end if
   end if

!  pawfatbnd
   call chkint_eq(0,0,cond_string,cond_values,ierr,'pawfatbnd',dt%pawfatbnd,3,(/0,1,2/),iout)
   if(usepaw/=1.and.dt%pawfatbnd>1) then
     write(message,'(4a)') ch10,&
&     ' chkinp: WARNING -',ch10,&
&     '  pawfatbnd without PAW is not possible'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if(dt%prtdosm>0.and.dt%pawfatbnd>0)then
     write(message,'(4a)') ch10,&
&     ' chkinp: WARNING -',ch10,&
&     ' pawfatbnd>0  and prtdosm>0 are not compatible '
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  pawlcutd
   if (usepaw==1) then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'pawlcutd',dt%pawlcutd,0,iout)
   end if

!  pawlmix
   if (usepaw==1) then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'pawlmix',dt%pawlmix,0,iout)
   end if

!  pawmixdg
   if (usepaw==1) then
     if(dt%ionmov==4)then
       cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawmixdg',dt%pawmixdg,1,(/1/),iout)
     end if
     if(dt%iscf==5.or.dt%iscf==6.or.dt%iscf==15.or.dt%iscf==16)then
       cond_string(1)='iscf' ; cond_values(1)=dt%iscf
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawmixdg',dt%pawmixdg,1,(/1/),iout)
     end if
   end if

!  pawnhatxc
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawnhatxc',dt%pawnhatxc,2,(/0,1/),iout)
   end if

!  pawnzlm
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawnzlm',dt%pawnzlm,2,(/0,1/),iout)
   end if

!  pawoptmix
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawoptmix',dt%pawoptmix,2,(/0,1/),iout)
   end if

!  pawprtdos
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawprtdos',dt%pawprtdos,3,(/0,1,2/),iout)
   end if

!  pawprtvol
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawprtvol',dt%pawprtvol,7,(/-3,-2,-1,0,1,2,3/),iout)
   end if

!  pawspnorb
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawspnorb',dt%pawspnorb,2,(/0,1/),iout)
     if (dt%pawspnorb==1.and.(dt%kptopt==1.or.dt%kptopt==2)) then
       write(message, '(10a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&       '  time-reversal symmetry is broken; k-points cannot',ch10,&
&       '  be generated using TR-symmetry.',ch10,&
&       '  Action: choose kptopt different from 1 or 2.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ierr=ierr+1
     end if
   end if

!  pawstgylm
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawstgylm',dt%pawstgylm,2,(/0,1/),iout)
   end if

!  pawusecp
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,2,(/0,1/),iout)
     if (dt%mkmem/=0)then
       cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,1,(/1/),iout)
     end if
     if (dt%mk1mem/=0)then
       cond_string(1)='mk1mem' ; cond_values(1)=dt%mk1mem
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,1,(/1/),iout)
     end if
     if (dt%mkqmem/=0)then
       cond_string(1)='mkqmem' ; cond_values(1)=dt%mkqmem
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,1,(/1/),iout)
     end if
   end if

!  pawxcdev
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawxcdev',dt%pawxcdev,3,(/0,1,2/),iout)
   end if

!  pitransform
   call chkint_eq(0,0,cond_string,cond_values,ierr,'pitransform',&
&   dt%pitransform,3,(/0,1,2/),iout)
!  When imgmov is not one of 9 or 13, pitransform must be 0
   if(dt%imgmov/=9 .and. dt%imgmov/=13 )then
     cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
!    Make sure that pitransform=0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'pitransform',dt%pitransform,1,(/0/),iout)
   end if

!  positron
   call chkint_eq(0,0,cond_string,cond_values,ierr,'positron',dt%positron,7,(/-20,-10,-2,-1,0,1,2/),iout)
   if ((dt%positron==2.or.dt%positron<0).and.(dt%ixcpositron==3.or.dt%ixcpositron==31)) then
     if ((dt%ixc<11.or.dt%ixc>17).and.dt%ixc/=23.and.dt%ixc/=26.and.dt%ixc/=27) then
       write(message, '(10a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  For the electronic ground-state calculation in presence of a positron,',ch10,&
&       '  when GGA is selected for electron-positron correlation (ixcpositron=3 or 31),',ch10,&
&       '  electron-electron XC must also be GGA !',ch10,&
&       '  Action: choose another psp file.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ierr=ierr+1
     end if
   end if
   if (dt%positron/=0.and.dt%ionmov==5) then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'positron',dt%positron,1,(/0/),iout)
   end if
   if (dt%positron<0.and.usepaw==0) then
     write(message, '(8a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  You cannot use positron<0 (automatic two-component DFT)',ch10,&
&     '  with norm-conserving pseudopotentials !',ch10,&
&     '  Action: choose PAW.'
     call wrtout(iout   ,message,'COLL')
     call wrtout(std_out,message,'COLL')
     ierr=ierr+1
   end if
   if ((dt%positron==1.or.dt%positron<0).and.dt%iscf<10.and.dt%tolvrs>tiny(one)) then
     write(message, '(10a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  You cannot perform a positronic ground-state calculation (positron=1 or <0)',ch10,&
&     '  using SCF potential mixing (iscf<10) and tolvrs !',ch10,&
&     '  (in that case, the potential is constant)',ch10,&
&     '  Action: change iscf or select another convergence criterion.'
     call wrtout(iout   ,message,'COLL')
     call wrtout(std_out,message,'COLL')
     ierr=ierr+1
   end if

!  posocc
   call chkdpr(0,0,cond_string,cond_values,ierr,'posocc',dt%posocc,-1,one,iout)

!  postoldfe, postoldff
   call chkdpr(0,0,cond_string,cond_values,ierr,'postoldff',dt%postoldff,1,zero,iout)
   if (dt%positron<0) then
     if ( (abs(dt%postoldfe)> tiny(0.0_dp).and.abs(dt%postoldff)> tiny(0.0_dp)).or.&
&     (abs(dt%postoldfe)<=tiny(0.0_dp).and.abs(dt%postoldff)<=tiny(0.0_dp))) then
       write(message,'(8a)' ) ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  One and only one of the input tolerance criteria postoldfe or postoldff',ch10,&
&       '  must differ from zero !',ch10,&
&       '  Action : change postoldfe or postldff in input file.'
       call wrtout(iout   ,message,'COLL')
       call wrtout(std_out,message,'COLL')
       ierr=ierr+1
     end if
     if (abs(dt%postoldff)>tiny(0.0_dp).and.dt%optforces/=1)then
       write(message,'(6a)' ) ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  When postoldff is set to a non-zero value, optforces must be set to 1 !',ch10,&
&       '  Action : change your input file.'
       call wrtout(iout   ,message,'COLL')
       call wrtout(std_out,message,'COLL')
       ierr=ierr+1
     end if
   end if

!  prepanl
!  Cannot prepare a non-linear calculation with ixc is not 3 or 7
   if ((dt%ixc /= 0).and.(dt%ixc /= 3).and.(dt%ixc /= 7).and.(dt%ixc /= 8)) then
     cond_string(1)='ixc' ; cond_values(1)=dt%ixc
     call chkint_ne(1,1,cond_string,cond_values,ierr,'prepanl',dt%prepanl,1,(/1/),iout)
   end if
!  Must have prtden=1 to prepare a nonlinear calculation
   if (dt%prtden /= 1) then
     cond_string(1)='prtden' ; cond_values(1)=dt%prtden
     call chkint_ne(1,1,cond_string,cond_values,ierr,'prepanl',dt%prepanl,1,(/1/),iout)
   end if

!  prtbbb
!  Not allowed for PAW
   if(usepaw==1.and.dt%prtbbb==1)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtbbb',dt%prtbbb,1,(/0/),iout)
   end if



!  prtden
   if (usepaw==1) then
     call chkint_le(0,0,cond_string,cond_values,ierr,'prtden',dt%prtden,6,iout)
   else
     call chkint_le(0,0,cond_string,cond_values,ierr,'prtden',dt%prtden,1,iout)
   end if

!  prtdensph
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'prtdensph',dt%prtdensph,2,(/0,1/),iout)
   end if

!  prtdos
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtdos',dt%prtdos,4,(/0,1,2,3/),iout)

!  prtdosm
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtdosm',dt%prtdosm,2,(/0,1/),iout)
   if(usepaw==1.and.dt%pawprtdos==1)then
     cond_string(1)='pawprtdos' ; cond_values(1)=dt%pawprtdos
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtdosm',dt%prtdosm,1,(/0/),iout)
   end if

!  prtelf
   call chkint_ge(0,0,cond_string,cond_values,ierr,'prtelf',dt%prtkden,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtelf',dt%prtelf,1,(/0/),iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtelf',dt%prtelf,1,(/0/),iout)
   end if

!  prtfsurf
!  only one shift allowed (gamma)
   if (dt%prtfsurf == 1) then

     if (abs(dt%kptrlatt(1,2))+abs(dt%kptrlatt(1,3))+abs(dt%kptrlatt(2,3))+&
&     abs(dt%kptrlatt(2,1))+abs(dt%kptrlatt(3,1))+abs(dt%kptrlatt(3,2)) /= 0 ) then
       ierr=ierr+1
       write(message,'(4a)')ch10,&
&       ' prtfsurf does not work with non-diagonal kptrlatt ', ch10,&
&       ' Action: set nshift 1 and shiftk 0 0 0'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
     if (dt%nshiftk > 1) then
       ierr=ierr+1
       write(message,'(4a)') ch10,&
&       ' prtfsurf does not work with multiple kpt shifts ', ch10, &
&       ' Action: set nshift 1 and shiftk 0 0 0'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
     if (sum(abs(dt%shiftk(:,1:dt%nshiftk))) > tol8) then
       ierr=ierr+1
       write(message,'(4a)')ch10,&
&       ' prtfsurf does not work with non-zero kpt shift ',ch10,&
&       ' Action: set nshift 1 and shiftk 0 0 0'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

!    Occcupations, Fermi level and k weights have to be calculated correctly.
     if (.not.(dt%iscf>1.or.dt%iscf==-3)) then
       ierr=ierr+1
       write(message,'(4a)')ch10,&
&       ' prtfsurf==1 requires either iscf>1 or iscf==-3 ',ch10,&
&       ' Action: change iscf in the input file. '
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

!    Make sure all nband are equal (well it is always enforced for metals)
     if (any(dt%nband(:) /= maxval(dt%nband(:)) )) then
       write(message,'(6a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The number of bands has to be constant for the output of the Fermi surface.',ch10,&
&       '  Action : set all the nbands to the same value in your input file'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

   end if ! prtfsurf==1

!  prtgden
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'prtgden',dt%prtgden,1,(/0/),1,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint(0,1,cond_string,cond_values,ierr,&
&     'prtgden',dt%prtgden,1,(/0/),0,0,iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint(0,1,cond_string,cond_values,ierr,&
&     'prtgden',dt%prtgden,1,(/0/),0,0,iout)
   end if

!  prtkden
   call chkint_ge(0,0,cond_string,cond_values,ierr,'prtkden',dt%prtkden,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtkden',dt%prtkden,1,(/0/),iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtkden',dt%prtkden,1,(/0/),iout)
   end if

!  prtlden
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'prtlden',dt%prtlden,1,(/0/),1,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint(0,1,cond_string,cond_values,ierr,&
&     'prtlden',dt%prtlden,1,(/0/),0,0,iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint(0,1,cond_string,cond_values,ierr,&
&     'prtlden',dt%prtlden,1,(/0/),0,0,iout)
   end if

!  prtstm
   call chkint_ge(0,0,cond_string,cond_values,ierr,'prtstm',dt%prtstm,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%occopt/=7)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%nstep/=1)then
     cond_string(1)='nstep' ; cond_values(1)=dt%nstep
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%ionmov/=0)then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%tolwfr<tol6)then
     cond_string(1)='tolwfr' ; cond_values(1)=0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtden/=0)then
     cond_string(1)='prtden' ; cond_values(1)=dt%prtden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtnabla>0)then
     cond_string(1)='prtnabla' ; cond_values(1)=dt%prtnabla
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtvxc>0)then
     cond_string(1)='prtvxc' ; cond_values(1)=dt%prtvxc
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtvha>0)then
     cond_string(1)='prtvha' ; cond_values(1)=dt%prtvha
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtvhxc>0)then
     cond_string(1)='prtvhxc' ; cond_values(1)=dt%prtvhxc
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if

!  prtwant
#if !defined HAVE_WANNIER90
   if(dt%prtwant==2) then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '         prtwant==2 is only relevant if wannier90 library is linked',ch10,&
&     '         Action: check compilation options'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
#endif

!  prtwf
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtwf',dt%prtwf,3,(/0,1,2/),iout)

!  prtxangst
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtxangst',dt%prtxangst,2,(/0,1/),iout)

!  prtxcart
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtxcart',dt%prtxcart,2,(/0,1/),iout)

!  prtxred
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtxred',dt%prtxred,2,(/0,1/),iout)
   if(dt%prtxangst==0 .and. dt%prtxcart==0)then
     cond_string(1)='prtxangst' ; cond_values(1)=dt%prtxangst
     cond_string(2)='prtxcart' ; cond_values(2)=dt%prtxcart
     call chkint_eq(2,2,cond_string,cond_values,ierr,'prtxred',dt%prtxred,1,(/1/),iout)
   end if

!  ratsph
!  If PAW and (prtdos==3 or dt%prtdensph==1), must be greater than PAW radius
   if(usepaw==1.and.(dt%prtdos==3.or.dt%prtdensph==1))then
     do itypat=1,dt%ntypat
       if (pspheads(itypat)%pawheader%rpaw>dt%ratsph(itypat)) then
         write(message, '(10a,i2,a,f15.12,3a)' ) ch10,&
&         ' chkinp: ERROR -',ch10,&
&         '  Projected DOS/density is required in the framework of PAW !',ch10,&
&         '  The radius of spheres in which DOS/density has to be projected',ch10,&
&         '  must be greater or equal than the (max.) PAW radius !',ch10,&
&         '  Rpaw(atom_type ',itypat,')= ',pspheads(itypat)%pawheader%rpaw,' au',ch10,&
&         '  Action : modify value of ratsph in input file.'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         ierr=ierr+1
       end if
     end do
   end if


!  recgratio
   if (dt%tfkinfunc==2) then
     write(message, '(a,a)' ) ch10,&
&     '=== RECURSION METHOD ==========================================================='
     call wrtout(ab_out,message,'COLL')
     cond_string(1)='tfkinfunc' ; cond_values(1)=2
     call chkint_ge(0,1,cond_string,cond_values,ierr,'recgratio',dt%recgratio,1,iout)
     if(dt%recgratio>1) then
       write(message, '(a,a)' )&
&       '=== Coarse Grid is used in recursion ==========================================='
       call wrtout(ab_out,message,'COLL')
       write(message, '(a,i2,a,a,i2,a,i2,a,i2)' ) 'grid ratio =',dt%recgratio,&
&       ch10,'fine grid =   ',dt%ngfft(1),' ',dt%ngfft(2),' ',dt%ngfft(3)
       call wrtout(ab_out,message,'COLL')
       write(message, '(a,i2,a,i2,a,i2)' ) 'coarse grid = ',&
&       dt%ngfft(1)/dt%recgratio,' ',dt%ngfft(2)/dt%recgratio,' ',dt%ngfft(3)/dt%recgratio
       call wrtout(ab_out,message,'COLL')
     else
       write(message, '(a,i2,a,i2,a,i2)' ) 'fine grid =   ',dt%ngfft(1),' ',dt%ngfft(2),' ',dt%ngfft(3)
       call wrtout(ab_out,message,'COLL')

     end if
   end if


!  rfatpol
   call chkint_ge(0,0,cond_string,cond_values,ierr,'rfatpol(1)',dt%rfatpol(1),1,iout)
   cond_string(1)='natom' ; cond_values(1)=natom
   call chkint_le(1,1,cond_string,cond_values,ierr,'rfatpol(2)',dt%rfatpol(2),natom,iout)

!  rprimd
!  With optcell beyond 4, one has constraints on rprimd.
   if(dt%optcell==4 .or. dt%optcell==7 )then
     cond_string(1)='optcell' ; cond_values(1)=4
     if(dt%optcell==7)cond_values(1)=7
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,2)',rprimd(1,2),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,3)',rprimd(1,3),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,1)',rprimd(2,1),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,1)',rprimd(3,1),0,0.0_dp,iout)
   else if(dt%optcell==5 .or. dt%optcell==8 )then
     cond_string(1)='optcell' ; cond_values(1)=5
     if(dt%optcell==8)cond_values(1)=8
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,1)',rprimd(2,1),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,3)',rprimd(2,3),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,2)',rprimd(1,2),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,2)',rprimd(3,2),0,0.0_dp,iout)
   else if(dt%optcell==6 .or. dt%optcell==9 )then
     cond_string(1)='optcell' ; cond_values(1)=6
     if(dt%optcell==9)cond_values(1)=9
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,1)',rprimd(3,1),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,2)',rprimd(3,2),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,3)',rprimd(1,3),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,3)',rprimd(2,3),0,0.0_dp,iout)
   end if

!  so_psp
   do ipsp=1,npsp
!    Check that so_psp is between 0 and 3
     if ( dt%so_psp(ipsp)<0 .or. dt%so_psp(ipsp)>3 ) then
       write(message, '(a,a,a,a,i3,a,i3,a,a,a,a,a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  so_psp(',ipsp,' ) was input as',&
&       dt%so_psp(ipsp),' .',ch10,&
&       '  Input value must be 0, 1, 2, or 3.',ch10,&
&       ' Action : modify value of so_psp (old name : so_typat) in input file.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ierr=ierr+1
     end if
!    If nspinor=1, the spin-orbit contribution cannot be taken into account
     if ( nspinor==1 .and. (dt%so_psp(ipsp)==2 .or. dt%so_psp(ipsp)==3) ) then
       write(message, '(a,a,a,a,i2,a,i3,a,a,a,a,a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  so_psp(',ipsp,') was input as',&
&       dt%so_psp(ipsp),', with nspinor=1.',ch10,&
&       '  When nspinor=1, so_psp cannot be required to be 2 or 3.',ch10,&
&       '  Action : modify value of so_psp (old name : so_typat) or nspinor in input file.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ierr=ierr+1
     end if
!    
!    TO BE ACTIVATED IN V5.5
!    Send a warning if the spin-orbit contribution cannot be included due to nspinor=1
!    if ( nspinor==1 .and. dt%so_psp(ipsp)==1 .and. pspheads(ipsp)%pspso/=1 ) then
!    write(message, '(a,a,a,a,i2,a,i3,a,a,a,a,a)' ) ch10,&
!    &    ' chkinp: WARNING -',ch10,&
!    &    '  so_psp(',ipsp,') was input as',&
!    &     dt%so_psp(ipsp),', with nspinor=1 and a pseudopotential containing spin-orbit contribution.',ch10,&
!    &    '  However, when nspinor=1, existing psp spin-orbit contribution is not taken into account.'
!    call wrtout(iout,message,'COLL')
!    call wrtout(std_out,  message,'COLL')
!    end if
!    END OF TO BE ACTIVATED
!    
!    If nspden=4, so_psp must be 1
!    if ( dt%so_psp(ipsp)/=1 .and. nspden==4 ) then
!    write(message, '(a,a,a,a,i2,a,i3,7a)' ) ch10,&
!    &    ' chkinp: ERROR -',ch10,&
!    &    '  so_psp(',ipsp,') was input as',&
!    &     dt%so_psp(ipsp),', with nspden=4.',ch10,&
!    &    '  However, non-collinear magnetism is not yet implemented with spin-orbit.',ch10,&
!    &    '  so_psp must be 1 for each pseudopotential type.',ch10,&
!    &    '  Action : modify value of so_psp or nspden in input file.'
!    call wrtout(iout,message,'COLL')
!    call wrtout(std_out,  message,'COLL')
!    ierr=ierr+1
!    end if
   end do

!  stmbias
   cond_string(1)='prtstm' ; cond_values(1)=dt%prtstm
   if(dt%prtstm/=0)then
!    If non-zero prtstm, stmbias cannot be zero : test is positive or zero
     if(dt%stmbias>-tol10)then
!      Then, enforce positive
       call chkdpr(1,1,cond_string,cond_values,ierr,'stmbias',dt%stmbias,1,2*tol10,iout)
     end if
   else
     call chkdpr(1,1,cond_string,cond_values,ierr,'stmbias',dt%stmbias,0,zero,iout)
   end if

!  symafm
   if(nsppol==1 .and. nspden==2)then
!    At least one of the symmetry operations must be antiferromagnetic
     if(minval(dt%symafm(1:dt%nsym))/=-1)then
       write(message, '(8a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  When nsppol==1 and nspden==2, at least one of the symmetry operations',ch10,&
&       '  must be anti-ferromagnetic (symafm=-1), in order to deduce the spin-down density',ch10,&
&       '  from the spin-up density.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       write(message, '(7a)' ) &
&       '  However, it is observed that none of the symmetry operations is anti-ferromagnetic.',ch10,&
&       '  Action : Check the atomic positions, the input variables spinat, symrel, tnons, symafm.',ch10,&
&       '           In case your system is not antiferromagnetic (it might be ferrimagnetic ...),',ch10,&
&       '           you must use nsppol=2 with nspden=2 (the latter being the default when nsppol=2).'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  symrel and tnons
!  Check the point group closure (TODO should check the spatial group closure !!)
   call chkgrp(dt%nsym,dt%symafm,dt%symrel,ierrgrp)
   if (ierrgrp==1) ierr=ierr+1

!  Check the orthogonality of the symmetry operations
!  (lengths and absolute values of scalar products should be preserved)
   call chkorthsy(gprimd,iout,dt%nsym,rmet,rprimd,dt%symrel)

!  symchi
   if (dt%symchi/=0.and.dt%symchi/=1.and.dt%symchi/=2) then
     write(message, '(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  symchi  was input as ',dt%symchi,ch10,&
&     '  Input value must be 0, 1, or 2.',ch10,&
&     ' Action : modify value of symchi in input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     ierr=ierr+1
   end if

!  symsigma
   if (dt%symsigma/=0.and.dt%symsigma/=1.and.dt%symsigma/=2) then
     write(message, '(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&     ' chkinp: ERROR -',ch10,&
&     '  symsigma  was input as',dt%symsigma,ch10,&
&     '  Input value must be 0, 1, or 2.',ch10,&
&     ' Action : modify value of symsigma in input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     ierr=ierr+1
   end if

!  MG now it is possible to perform a GW calculation with non-symmorphic operations if required by the user
!  tnons
   if (dt%symmorphi==0) then
     if(dt%nbandkss/=0)then
       do isym=1,dt%nsym
         if(sum(dt%tnons(:,isym)**2)>tol6)then
           write(message, '(6a,i3,a,3f8.4,3a)' ) ch10,&
&           ' chkinp: ERROR -',ch10,&
&           '  When nbandkss/=0, all the components of tnons must be zero.',ch10,&
&           '  However, for the symmetry operation number ',isym,', tnons =',dt%tnons(:,isym),'.',ch10,&
&           '  Action : use the symmetry finder (nsym=0) with symmorphi==0.'
           call wrtout(iout,message,'COLL')
           call wrtout(std_out,  message,'COLL')
           call leave_new('COLL')
         end if
       end do
     end if
     if (ANY(optdriver==(/RUNL_SCREENING,RUNL_SIGMA/) ))then
       do isym=1,dt%nsym
         if (sum(dt%tnons(:,isym)**2)>tol6) then
           ierr=ierr+1
           write(message, '(6a,i3,a,3f8.4,3a)' ) ch10,&
&           ' chkinp: ERROR -',ch10,&
&           '  When optdriver==RUNL_SCREENING or RUNL_SIGMA, all the components of tnons must be zero.',ch10,&
&           '  However, for the symmetry operation number ',isym,', tnons =',dt%tnons(:,isym),'.',ch10,&
&           '  Action : use the symmetry finder (nsym=0) with symmorphi==0.'
           call wrtout(iout,message,'COLL')
         end if
       end do
     end if
   end if !of symmorphi

!  toldff
   call chkdpr(0,0,cond_string,cond_values,ierr,'toldff',dt%toldff,1,zero,iout)

!  tolimg
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolimg',dt%tolimg,1,zero,iout)

!  tolrff
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolrff',dt%tolrff,1,zero,iout)

!  tolwfr
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolwfr',dt%tolwfr,1,zero,iout)

!  tsmear
   call chkdpr(0,0,cond_string,cond_values,ierr,'tsmear',dt%tsmear,1,zero,iout)
!  Check that tsmear is non-zero positive for metallic occupation functions
   if(3<=dt%occopt .and. dt%occopt<=7)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkdpr(1,1,cond_string,cond_values,ierr,'tsmear',dt%tsmear,1,tol8,iout)
   end if

!  usedmatpu
   if (usepaw==1.and.dt%usepawu==1) then
!    abs(dt%usedmatpu)<=nstep
     cond_string(1)='nstep' ; cond_values(1)=dt%nstep
     call chkint_le(1,1,cond_string,cond_values,ierr,'abs(usedmatpu)',abs(dt%usedmatpu),dt%nstep-1,iout)
!    lpawu must be constant or -1
     if (dt%usedmatpu/=0) then
       do itypat=1,dt%ntypat
         if (dt%lpawu(itypat)/=-1.and.dt%lpawu(itypat)/=maxval(dt%lpawu(:))) then
           write(message, '(6a)' ) ch10,&
&           ' chkinp: ERROR -',ch10,&
&           '  When usedmatpu/=0 (use of an initial density matrix for LDA+U),',ch10,&
&           '  lpawu must be equal for all types of atoms on which +U is applied !'
           call wrtout(iout,message,'COLL')
           call wrtout(std_out,  message,'COLL')
           call leave_new('COLL')
         end if
       end do
     end if
   end if

!  usedmft
   if (dt%usedmft>0) then
     cond_string(1)='usedmft' ; cond_values(1)=1
     call chkint_eq(0,1,cond_string,cond_values,ierr,'usedmft',dt%usedmft,2,(/1,2/),iout)
   end if

!  useexexch and lexexch
!  Local exact-exchange and restrictions
   if(dt%useexexch/=0)then
     cond_string(1)='useexexch' ; cond_values(1)=dt%useexexch
     call chkint_eq(1,1,cond_string,cond_values,ierr,'useexexch',dt%useexexch,1,(/1/),iout)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usepaw',usepaw,1,(/1/),iout)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'pawxcdev',dt%pawxcdev,2,(/1,2/),iout)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'ixc',dt%ixc,2,(/11,23/),iout)
     do itypat=1,dt%ntypat
       cond_string(1)='lexexch' ; cond_values(1)=dt%lexexch(itypat)
       call chkint_eq(1,1,cond_string,cond_values,ierr,'lexexch',dt%lexexch(itypat),5,(/-1,0,1,2,3/),iout)
     end do
   end if

!  usekden
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usekden',dt%usekden,2,(/0,1/),iout)
   if(dt%usekden==0)then
     cond_string(1)='usekden' ; cond_values(1)=dt%usekden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtkden',dt%prtkden,1,(/0/),iout)
   end if
   if(dt%usekden/=0)then
     cond_string(1)='usekden' ; cond_values(1)=dt%usekden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usewvl',usewvl,1,(/0/),iout)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usepaw',usepaw,1,(/0/),iout)
   end if

!  usepawu and lpawu
!  PAW+U and restrictions
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usepawu',dt%usepawu,5,(/0,1,2,3,10/),iout)
   if(dt%usepawu/=0)then
     cond_string(1)='usepawu' ; cond_values(1)=dt%usepawu
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usepaw',usepaw,1,(/1/),iout)
     do itypat=1,dt%ntypat
       cond_string(1)='lpawu' ; cond_values(1)=dt%lpawu(itypat)
       call chkint_eq(1,1,cond_string,cond_values,ierr,'lpawu',dt%lpawu(itypat),5,(/-1,0,1,2,3/),iout)
     end do
     if(dt%pawspnorb>0) then
       write(message,'(6a)' ) ch10,&
&       ' chkinp: WARNING -',ch10,&
&       '  LDA+U+SpinOrbit is still on test ',ch10,&
&       '  (not yet in production)'
       call wrtout(std_out,message,'COLL')
!      call leave_new('COLL')
     end if
   end if

!  useexexch AND usepawu
!  Restriction when use together
   if(dt%useexexch>0.and.dt%usepawu>0)then
     do itypat=1,dt%ntypat
       if (dt%lpawu(itypat)/=dt%lexexch(itypat).and.&
&       dt%lpawu(itypat)/=-1.and.dt%lexexch(itypat)/=-1) then
         write(message, '(8a,i2,3a)' ) ch10,&
&         ' chkinp: ERROR - ',ch10,&
&         '  When PAW+U (usepawu>0) and local exact-exchange (useexexch>0)',ch10,&
&         '  are selected together, they must apply on the same',ch10,&
&         '  angular momentum (lpawu/=lexexch forbidden, here for typat=',&
&         itypat,') !',ch10,'  Action: correct your input file.'
         call wrtout(iout,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if

!  usexcnhat
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usexcnhat',dt%usexcnhat,3,(/-1,0,1/),iout)

!  useylm
   call chkint_eq(0,0,cond_string,cond_values,ierr,'useylm',dt%useylm,2,(/0,1/),iout)
   if (usepaw==1) then
     write(cond_string(1), "(A)") 'pspcod'
     cond_values(1)=7
     call chkint_eq(1,1,cond_string,cond_values,ierr,'useylm',dt%useylm,1,(/1/),iout)
   end if

!  wfoptalg
!  Must be greater or equal to 0
   call chkint_ge(0,0,cond_string,cond_values,ierr,'wfoptalg',dt%wfoptalg,0,iout)
!  wfoptalg==0,1,10,11,4,14,5 if PAW
   if (usepaw==1) then
     write(cond_string(1), "(A)") 'usepaw'
     cond_values(1)=1
     call chkint_eq(0,1,cond_string,cond_values,ierr,'wfoptalg',dt%wfoptalg,7,(/0,1,10,11,4,14,5/),iout)
   end if
   if (fftalg/=400 .and. fftalg/=401) then   ! If fftalg/=400, cannot use wfoptalg=4 or 14 (lopbcg algo)
     write(cond_string(1), "(A)") 'fftalg'
     cond_values(1)=fftalg
     call chkint_eq(0,1,cond_string,cond_values,ierr,'wfoptalg',dt%wfoptalg,6,(/0,1,2,3,10,11/),iout)
   end if

!  wtk
!  Check that no k point weight is < 0:
   do ikpt=1,nkpt
     if (dt%wtk(ikpt)< -tiny(0.0_dp) ) then
       write(message, '(a,a,a,a,i5,a,1p,e12.4,a,a,a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  At k point number',ikpt,'  wtk=',dt%wtk(ikpt),' <0.',ch10,&
&       '  Action : check wtk in input file. Each wtk must be >=0.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if
   end do

!  xred
!  Check that two atoms are not on top of each other
   do iimage=1,dt%nimage
     if(natom>1)then
       allocate(frac(3,natom))
       do ia=1,natom
!        Map reduced coordinate xred(mu,ia) into [0,1)
         frac(1,ia)=dt%xred_orig(1,ia,iimage)-aint(dt%xred_orig(1,ia,iimage))+0.5_dp-sign(0.5_dp,dt%xred_orig(1,ia,iimage))
         frac(2,ia)=dt%xred_orig(2,ia,iimage)-aint(dt%xred_orig(2,ia,iimage))+0.5_dp-sign(0.5_dp,dt%xred_orig(2,ia,iimage))
         frac(3,ia)=dt%xred_orig(3,ia,iimage)-aint(dt%xred_orig(3,ia,iimage))+0.5_dp-sign(0.5_dp,dt%xred_orig(3,ia,iimage))
       end do
       do ia=1,natom-1
         do ib=ia+1,natom
           if( abs(frac(1,ia)-frac(1,ib))<1.0d-6 .and. &
&           abs(frac(2,ia)-frac(2,ib))<1.0d-6 .and. &
&           abs(frac(3,ia)-frac(3,ib))<1.0d-6         ) then
             if(iimage>1)then
               write(message,'(2a,i5)') ch10,&
&               ' The following was observed for image=',iimage
               call wrtout(iout,message,'COLL')
               call wrtout(std_out,message,'COLL')
             end if
             write(message, '(a,a,a,a,i4,a,i4,a,a,a,a,a,a)' ) ch10,&
&             ' chkinp: ERROR - ',ch10,&
&             '  Atoms number',ia,' and',ib,' are located at the same point',&
&             ' of the unit cell',ch10,&
&             '  (periodic images are taken into account).',ch10,&
&             '  Action: change the coordinate of one of these atoms in the input file.'
             call wrtout(iout,message,'COLL')
             call wrtout(std_out,  message,'COLL')
             call leave_new('COLL')
           end if
         end do
       end do
       deallocate(frac)
     end if
   end do

!  znucl
!  Check that znucl and znuclpsp agree
   do ipsp=1,npsp
     if (abs(dt%znucl(ipsp)-pspheads(ipsp)%znuclpsp)> tol12 ) then
       write(message, '(4a,i5,a,es12.4,a,a,es12.4,2a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  For pseudopotential ',ipsp,'  znucl from user input file= ',dt%znucl(ipsp),ch10,&
&       '  while znucl from pseudopotential file=',pspheads(ipsp)%znuclpsp,ch10,&
&       '  Action : check znucl in input file, or check psp file. They must agree.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if
   end do

!  bandFFT
   if(dt%paral_kgb==1) then
     if (mod(dt%wfoptalg,10) /= 4) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of wfoptalg is found to be ',dt%wfoptalg,ch10,&
&       '  This is not allowed in the case of band-FFT parallelization.',ch10,&
&       '  Action : put wfoptalg = 4 or 14 in your input file'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
!    if (nspinor /= 1) then
!    write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
!    &       ' chkinp : ERROR -',ch10,&
!    &       '  The value of nspinor is found to be ',nspinor,ch10,&
!    &       '  This is not allowed in the case of band-FFT parallelization.',ch10,&
!    &       '  Action : put nspinor = 1 in your input file'
!    call wrtout(std_out,message,'COLL')
!    call leave_new('COLL')
!    end if
!    Make sure all nband are equal
     if (any(dt%nband(1:nkpt*nsppol) /= maxval(dt%nband(1:nkpt*nsppol)) )) then
       write(message,'(a,a,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The number of bands have to remain constant in the case of band-FFT parallelization.',ch10,&
&       '  Action : set all the nbands to the same value in your input file'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(maxval(abs(dt%istwfk(1:nkpt)-1)) > 1)then
       write(message,'(8a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  One of the components of istwfk is not equal to 1 or 2.',ch10,&
&       '  Time-reversal symmetry is not yet programmed in the case of band-FFT parallelization.',ch10,&
&       '  Action : set istwfk to 1 or 2 for all k-points'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (dt%mkmem == 0) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of mkmem is found to be ',dt%mkmem,ch10,&
&       '  An out-of-core solution can''t be used in the case of band-FFT parallelization.',ch10,&
&       '  Action : put mkmem = nkpt in your input file'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  WVL - wavelets checks and limitations
   if(dt%usewvl == 1) then
     if (dt%wvl_hgrid <= 0) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of wvl_hgrid is found to be ',dt%wvl_hgrid,ch10,&
&       '  This value is mandatory and must be positive.',ch10,&
&       '  Action : put wvl_hgrid to a positive value in your input file'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (dt%iscf /= 2) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of iscf is found to be ',dt%iscf,ch10,&
&       '  Only simple potential mixing is allowed with wavelets.',ch10,&
&       '  Action : put iscf = 2 in your input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(dt%icoulomb /= 1)then
       write(message,'(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&       ' chkinp: ERROR -',ch10,&
&       '  The value of icoulomb is found to be ',dt%icoulomb,ch10,&
&       '  The real space computation of hartree potential is mandatory with wavelets.',ch10,&
&       '  Action : put icoulomb = 1 in your input file.'
       call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (dt%nsym /= 1) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of nsym is found to be ',dt%nsym,ch10,&
&       '  No symetry operations are allowed for isolated systems.',ch10,&
&       '  Action : put nsym = 1 in your input file'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (dt%optstress > 0) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of optstress is found to be ', dt%optstress, ch10,&
&       '  There is no stress computation available with the wavelet code.',ch10,&
&       '  Action : put optstress = 0 in your input file'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (usepaw == 1) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of usepaw is found to be ',usepaw,ch10,&
&       '  The wavelet computation is not allowed in the framework of PAW.',ch10,&
&       '  Action : use HGH or GTH pseudo-potentials'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (dt%nspden > 2) then
       write(message,'(a,a,a,a,i3,a,a,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of nspden is found to be ', dt%nspden, ch10, &
&       '  The wavelet computation is not allowed with non-colinear spin.',ch10,&
&       '  Action : put nspden = 1 or 2 in your input file'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (dt%nspden /= dt%nsppol) then
       write(message,'(a,a,a,a,i3,a,a,i3,a,a)')ch10,&
&       ' chkinp : ERROR -',ch10,&
&       '  The value of nspden is found to be ', dt%nspden, ch10, &
&       '  and the one of nsppol is found to be ', dt%nsppol, ch10, &
&       '  In wavelet computation, nspden and nsppol must be equal.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
!    We check the consistency of occupation, empty bands are not allowed.
     if (dt%nsppol == 2) then
       mband = dt%nelect
     else
       mband = dt%mband
     end if
     do ii = 1, mband, 1
       if (dt%occ_orig(ii) < tol8) then
         write(message,'(a,a,a,a,f6.4,a,a,a,a,a,a)') ch10,&
&         ' chkinp : ERROR -',ch10,&
&         '  One value of occ is found to be ', dt%occ_orig(ii), ch10, &
&         '  The wavelet computation is not allowed empty bands.',ch10,&
&         '  Action : use occopt = 1 for automatic band filling or', ch10, &
&         '  change occ value in your input file'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if


!  If molecular dynamics or structural optimization is being done
!  (dt%ionmov>0), make sure not all atoms are fixed
!  if (dt%ionmov > 0) then
!  if (natfix == natom) then
!  write(message, '(a,a,a,a,i4,a,i5,a,a,i5,a,a,a,a,a,a)' ) ch10,&
!  &   ' setup1: ERROR -',ch10,&
!  &   '  ionmov is ',dt%ionmov,' and number of fixed atoms is ',natfix,ch10,&
!  &   '  while number of atoms natom is ',natom,'.',ch10,&
!  &   '  Thus all atoms are fixed and option ionmov to move atoms',&
!  &           ' is inconsistent.',ch10,&
!  &   '  Action : change ionmov or natfix and iatfix in input file and resubmit.'
!  call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
!  end if
!  end if

!  Should check that the symmetry operations are consistent with iatfixx,
!  iatfixy and iatfixz (diagonal symmetry operations)

!  Should check values of fftalg

!  rfasr=2 possible only when electric field response is computed.

!  Must have nqpt=1 for rfphon=1

!  ** Here ends the checking section **************************************

   call dtsetfree(dt)

!  End do loop on idtset
 end do

 if (maxval(dtsets(:)%usewvl) > 0) then
   write(message,'(5A)') ch10,&
&   ' chkinp: COMMENT - comparison between wvl_hgrid and ecut',ch10,&
&   '  real-space mesh | eq. Ec around atoms | eq. Ec further from atoms'
   call wrtout(std_out,message,'COLL')
   wvl_hgrid = zero
   twvl = .false.
   do idtset=1,ndtset_alloc
!    Give an indication to the equivalent ecut corresponding to
!    given hgrid.
     if (dtsets(idtset)%usewvl == 1 .and. &
&     wvl_hgrid /= dtsets(idtset)%wvl_hgrid) then
       write(message,'(F11.3,A,F16.1,A,F16.1,A)') &
&       dtsets(idtset)%wvl_hgrid, " bohr  |", &
&       two * pi * pi / (dtsets(idtset)%wvl_hgrid ** 2), " Ht  | ", &
&       half * pi * pi / (dtsets(idtset)%wvl_hgrid ** 2), " Ht"
       call wrtout(std_out,message,'COLL')
       wvl_hgrid = dtsets(idtset)%wvl_hgrid
     end if
     twvl = twvl .or. (dtsets(idtset)%usewvl == 1 .and. dtsets(idtset)%accesswff /= IO_MODE_ETSF)
   end do
   if (twvl) then
     write(message,'(8a)') ch10,&
&     ' chkinp: WARNING -',ch10,&
&     '  Restart files from wavelets in non ETSF format does not follow', ch10, &
&     '  the ABINIT standards.', ch10, &
&     '  Put accesswff to 3 to use ETSF retart files.'
     call wrtout(std_out,message,'COLL')
   end if
 end if

!If there was a problem, then stop.
 if (ierr/=0) then
   write(message,'(4a,i4,a)')ch10,&
&   ' chkinp: ERROR -',ch10,&
&   '  Checking consistency of input data against itself gave ',ierr,' inconsistencies.'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 DBG_EXIT("COLL")

end subroutine chkinp
!!***
