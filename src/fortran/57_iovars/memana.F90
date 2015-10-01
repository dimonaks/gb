!{\src2tex{textfont=tt}}
!!****f* ABINIT/memana
!! NAME
!! memana
!!
!! FUNCTION
!! Analysis of the memory and disk space needed for the job,
!! thanks to the data computed in the calling routine : for each
!! array, the number of blocks of size mpw or nfft bytes, and the
!! additional memory occupation;
!! the list of arrays that are used for each chain.
!!
!! According to the value of the option variable,
!! the routine will eventually try to allocate this amount of memory,
!! and if it fails, estimate the maximum value nfft compatible with
!! the available memory.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cadd(marrays)= count of bytes needed in addition of cmpw, cfftc and cfft.
!!  cfft(marrays) =for each array, count of blocks of size nfft bytes (coarse grid, if PAW)
!!  cfftf(marrays)=for each array, count of blocks of size nfft bytes (fine grid, if PAW)
!!  chain(marrays,nchain)=logical variable, that informs whether an array
!!    belongs to a given chain.
!!  cmpw(marrays)=for each array, count of blocks of size mpw bytes.
!!  dttyp(marrays)=datatype of the array : 4 for integers, 8 for real(dp)
!!  iout=unit number for output of formatted data.
!!  iprcel=govern the choice of preconditioner for the SCF cycle
!!  iscf=governs the choice of SCF algorithm, or non-SCF calculation.
!!  marrays=maximal number of arrays (or group of arrays) to be monitored.
!!  mbcg=number of MB needed for the cg array.
!!  mbdiskpd=number of MB needed to store a density or potential file on disk
!!  mbdiskwf=number of MB needed to store a wavefunction file on disk
!!  mbf_fftgr=number of MB needed for the f_fftgr array.
!!  mbgylm=number of MB needed for the pawfgrtab%gylm array (paw only)
!!  mffmem =governs the number of FFT arrays which are fit in core memory
!!  mkmem =maximum number of k points which can fit in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  natom =number of atoms in unit cell
!!  nchain=number of chains to be used in the estimation of memory.
!!  nfft =(effective) number of FFT grid points (for one processor) (coarse grid, if PAW)
!!  nfftf=(effective) number of FFT grid points (for one processor) (fine grid, if PAW)
!!  occopt=option for occupation numbers. If 3<=occopt<=7, varying occupation
!!  option : if 0 , no test of available memory
!!           if 1 , the routine tries to allocate the estimated memory, for testing
!!                    purposes, and if a failure occurs, the routine stops.
!!           if 2 , like 1, but before stopping, the routine will provide
!!                    an estimation of the available memory.
!!  prtvol=control print volume
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      memorf,memory
!!
!! CHILDREN
!!      leave_new,leave_test,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine memana(cadd,cfft,cfftf,chain,cmpw,dttyp,iout,iprcel,iscf,&
& marrays,mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm,mffmem,&
& mkmem,mpi_enreg,mpw,natom,nchain,nfft,nfftf,occopt,option,prtvol)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,iprcel,iscf,marrays,mffmem,mkmem,mpw,natom,nchain
 integer,intent(in) :: nfft,nfftf,occopt,option,prtvol
 real(dp),intent(in) :: mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: cadd(marrays),cfft(marrays),cfftf(marrays),cmpw(marrays)
 integer,intent(in) :: dttyp(marrays)
 logical,intent(in) :: chain(marrays,nchain)

!Local variables-------------------------------
!scalars
 integer :: biggest,ichain,ier,ier1,ier2,ier3,ier4,ier5,ier6,ier7,ier8,ii
!integer :: jj,kk
 integer :: mu,nmbytes,nquarter_mbytes,quit
 real(dp) :: mbbigarr,mbbiggest
 character(len=500) :: message
!arrays
 integer,allocatable :: cdpfft(:),cdpfftf(:),cdpmpw(:),cintfft(:),cintfftf(:)
 integer,allocatable :: cintmpw(:)
 real(dp),allocatable :: bigarray(:,:),bigarray1(:,:),bigarray2(:,:)
 real(dp),allocatable :: bigarray3(:,:),bigarray4(:,:),bigarray5(:,:)
 real(dp),allocatable :: bigarray6(:,:),bigarray7(:,:),bigarray8(:,:),cdpadd(:)
 real(dp),allocatable :: cintadd(:),mbdpadd(:),mbdpfft(:),mbdpfftf(:)
 real(dp),allocatable :: mbdpmpw(:),mbintadd(:),mbintfft(:),mbintfftf(:)
 real(dp),allocatable :: mbintmpw(:),mbother(:),mbtot(:)

! **************************************************************************

!DEBUG
!write(6,*)' memana : enter'
!write(6,*)' memana : nchain=',nchain
!ENDDEBUG
 allocate(cdpfftf(nchain),cdpfft(nchain),cdpmpw(nchain))
 allocate(cintfftf(nchain),cintfft(nchain),cintmpw(nchain))
 allocate(cdpadd(nchain),cintadd(nchain))
 allocate(mbdpadd(nchain),mbdpfftf(nchain),mbdpfft(nchain),mbdpmpw(nchain))
 allocate(mbintadd(nchain),mbintfftf(nchain),mbintfft(nchain),mbintmpw(nchain))
 allocate(mbother(nchain),mbtot(nchain))

 biggest=0
 mbbiggest=0.0_dp

!For each chain, compute the number of bytes
 do ichain=1,nchain

!  First, the number of integer or real(dp), fft, mpw or add blocks
   cdpmpw(ichain) =sum(cmpw(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintmpw(ichain)=sum(cmpw(:),MASK=(dttyp(:)==4).and.chain(:,ichain))
   cdpfftf(ichain) =sum(cfftf(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintfftf(ichain)=sum(cfftf(:),MASK=(dttyp(:)==4).and.chain(:,ichain))
   cdpfft(ichain) =sum(cfft(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintfft(ichain)=sum(cfft(:),MASK=(dttyp(:)==4).and.chain(:,ichain))
   cdpadd(ichain) =sum(cadd(:),MASK=(dttyp(:)==8).and.chain(:,ichain))
   cintadd(ichain)=sum(cadd(:),MASK=(dttyp(:)==4).and.chain(:,ichain))

!  Compute the corresponding number of Mbytes
   mbdpmpw(ichain) =8*cdpmpw(ichain) *dble(mpw) /1024._dp**2
   mbintmpw(ichain)=4*cintmpw(ichain)*dble(mpw) /1024._dp**2
   mbdpfftf(ichain) =8*cdpfftf(ichain) *dble(nfftf)/1024._dp**2
   mbintfftf(ichain)=4*cintfftf(ichain)*dble(nfftf)/1024._dp**2
   mbdpfft(ichain) =8*cdpfft(ichain) *dble(nfft)/1024._dp**2
   mbintfft(ichain)=4*cintfft(ichain)*dble(nfft)/1024._dp**2
   mbdpadd(ichain) =8*cdpadd(ichain)              /1024._dp**2
   mbintadd(ichain)=4*cintadd(ichain)             /1024._dp**2
   mbother(ichain) =dble(231+6*natom)/1024._dp
   if(3<=occopt .and. occopt<=7)mbother(ichain)=dble(991+natom)/1024._dp

!  Compute the total number of Mbytes
   mbtot(ichain)=mbdpmpw(ichain)+mbintmpw(ichain)&
&   +mbdpfftf(ichain)+mbintfftf(ichain)&
&   +mbdpfft(ichain)+mbintfft(ichain)&
&   +mbdpadd(ichain)+mbintadd(ichain)+mbother(ichain)

!  Select the biggest chain
   if(mbtot(ichain)>mbbiggest)then
     mbbiggest=mbtot(ichain)
     biggest=ichain
   end if
 end do
!When iprcel<20, the biggest chains cannot be number 8 or 9 ...
 if(modulo(iprcel,100)<20 .and. (biggest==8 .or. biggest==9))then
   write(message,'(a,a,a,a,i3,a,a,a)') ch10,&
&   ' memana : BUG -',ch10,&
&   '  The biggest chain is number',biggest,' while iprcel==20.',ch10,&
&   '  This is not allowed.'
   call wrtout(std_out,message,'COLL')
 end if

 write(message, '(a,f11.3,a)' ) &
& 'P This job should need less than                 ',&
& mbbiggest+tol10,' Mbytes of memory. '
 call wrtout(std_out,message,'COLL')
 call wrtout(iout,message,'COLL')

 if(prtvol>=10)then
   if(biggest==1)write(message,'(a)')'P Max. in main chain + fourwf.f '
   if(biggest==2)write(message,'(a)')'P Max. in main chain + nonlop.f + opernl.f '
   if(biggest==3)write(message,'(a)')'P Max. in XC chain '
   if(biggest==4)write(message,'(a)')'P Max. in mkrho chain '
   if(biggest==5)write(message,'(a)')'P Max. in fourdp chain '
   if(biggest==6)write(message,'(a)')'P Max. in parallel k-point chain '
   if(biggest==7)write(message,'(a)')'P Max. in newvtr chain '
   if(biggest==8)write(message,'(a)')'P Max. in suscep chain '
   if(biggest==9)write(message,'(a)')'P Max. in dielmt chain '
   if(biggest==10)write(message,'(a)')'P Max. in tddft chain '
   call wrtout(iout,message,'COLL')

   write(message, '(a,i11,a,f11.3,a)' )&
&   'P',cintmpw(biggest),' blocks of mpw  integer numbers, for',&
&   mbintmpw(biggest)+tol10,' Mbytes. '
   call wrtout(iout,message,'COLL')
   write(message, '(a,i11,a,f11.3,a)' )&
&   'P',cdpmpw(biggest),' blocks of mpw  real(dp)  numbers, for',&
&   mbdpmpw(biggest)+tol10,' Mbytes. '
   call wrtout(iout,message,'COLL')
   if (nfft==nfftf) then
     if(mbintfft(biggest)+mbintfftf(biggest)>0.001)then
       write(message, '(a,i11,a,f11.3,a)' )&
&       'P',cintfft(biggest)+cintfftf(biggest),' blocks of nfft integer numbers, for',&
&       mbintfft(biggest)+mbintfftf(biggest)+tol10,' Mbytes. '
       call wrtout(iout,message,'COLL')
     end if
     write(message, '(a,i11,a,f11.3,a)' )&
&     'P',cdpfft(biggest)+cdpfftf(biggest),' blocks of nfft real(dp)  numbers, for',&
&     mbdpfft(biggest)+mbdpfftf(biggest)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
   else
     if(mbintfftf(biggest)>0.001)then
       write(message, '(a,i11,a,f11.3,a)' )&
&       'P',cintfftf(biggest),' blocks of nfft (fine grid) integer numbers, for',&
&       mbintfftf(biggest)+tol10,' Mbytes. '
       call wrtout(iout,message,'COLL')
     end if
     write(message, '(a,i11,a,f11.3,a)' )&
&     'P',cdpfftf(biggest),' blocks of nfft (fine grid) real(dp)  numbers, for',&
&     mbdpfftf(biggest)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
     if(mbintfft(biggest)>0.001)then
       write(message, '(a,i11,a,f11.3,a)' )&
&       'P',cintfft(biggest),' blocks of nfft (coarse grid) integer numbers, for',&
&       mbintfft(biggest)+tol10,' Mbytes. '
       call wrtout(iout,message,'COLL')
     end if
     write(message, '(a,i11,a,f11.3,a)' )&
&     'P',cdpfft(biggest),' blocks of nfft (coarse grid) real(dp)  numbers, for',&
&     mbdpfft(biggest)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
   end if
   if(mbintadd(biggest)>0.001)then
     write(message, '(a,11x,a,f11.3,a)' )&
&     'P',' Additional     integer numbers, for',mbintadd(biggest)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
   end if
   write(message, '(a,11x,a,f11.3,a)' )&
&   'P',' Additional     real(dp)  numbers, for',mbdpadd(biggest)+tol10,' Mbytes. '
   call wrtout(iout,message,'COLL')
   write(message, '(a,11x,a,f11.3,a)' )&
&   'P',' With residue estimated to be       ',mbother(biggest)+tol10,' Mbytes. '
   call wrtout(iout,message,'COLL')
   write(message, '(a)' )'P'
   call wrtout(iout,message,'COLL')
   write(message, '(a)' )&
&   'P Comparison of the memory needs of different chains'
   call wrtout(iout,message,'COLL')

   write(message, '(a,f11.3,a)' )&
&   'P Main chain + fourwf.f           ',mbtot(1)+tol10,' Mbytes. '
   call wrtout(iout,message,'COLL')
   write(message, '(a,f11.3,a)' )&
&   'P Main chain + nonlop.f + opernl.f',mbtot(2)+tol10,' Mbytes. '
   call wrtout(iout,message,'COLL')

!  The next chains are not defined in the RF case.
   if(nchain>2)then
     write(message, '(a,f11.3,a)' )&
&     'P XC chain                        ',mbtot(3)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
     write(message, '(a,f11.3,a)' )&
&     'P mkrho chain                     ',mbtot(4)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
     write(message, '(a,f11.3,a)' )&
&     'P fourdp chain                    ',mbtot(5)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
     if(mpi_enreg%paral_compil_kpt==1)then
       write(message, '(a,f11.3,a)' )&
&       '- parallel k-point chain          ',mbtot(6)+tol10,' Mbytes. '
       call wrtout(iout,message,'COLL')
     end if
     write(message, '(a,f11.3,a)' )&
&     'P newvtr chain                    ',mbtot(7)+tol10,' Mbytes. '
     call wrtout(iout,message,'COLL')
     if(modulo(iprcel,100)>=20.and.modulo(iprcel,100)<70)then
       write(message, '(a,f11.3,a)' )&
&       'P suscep chain                    ',mbtot(8)+tol10,' Mbytes. '
       call wrtout(iout,message,'COLL')
       write(message, '(a,f11.3,a)' )&
&       'P dielmt chain                    ',mbtot(9)+tol10,' Mbytes. '
       call wrtout(iout,message,'COLL')
     end if
     if(iscf==-1)then
       write(message, '(a,f11.3,a)' )&
&       'P tddft  chain                    ',mbtot(10)+tol10,' Mbytes. '
     end if
   end if ! nchain>2

 end if

!--------------------------------------------------------------------

 write(message, '(a)' ) &
& '  Rough estimation (10% accuracy) of disk space for files :'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message, '(a,f11.3,a,a,f11.3,a)' ) &
& '  WF disk file :',mbdiskwf+tol10,' Mbytes ;',&
& ' DEN or POT disk file :',mbdiskpd+tol10,' Mbytes.'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if(mkmem==0)then
   write(message, '(a)' )&
&   '  mkmem==0 => use of 2 WF temporary disk files'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 if(mffmem==0 .and. iscf>0)then
   if(iscf==1)then
     write(message, '(a,a,a)' )&
&     '  mffmem==0, iscf==1 => use of 1 FFT temporary disk file,',ch10,&
&     '                       5 times bigger than a DEN file.'
   else if(iscf==2.or.iscf==12)then
     write(message, '(a,a,a)' )&
&     '  mffmem==0, iscf==2 => use of 1 FFT temporary disk file,',ch10,&
&     '                       3 times bigger than a DEN file.'
   else if(iscf==3.or.iscf==13)then
     write(message, '(a,a,a)' )&
&     '  mffmem==0, iscf==3 => use of 1 FFT temporary disk file,',ch10,&
&     '                       4 times bigger than a DEN file.'
   else if(iscf==4.or.iscf==14)then
     write(message, '(a,a,a)' )&
&     '  mffmem==0, iscf==4 => use of 1 FFT temporary disk file,',ch10,&
&     '                       6 times bigger than a DEN file.'
   else if(iscf==5)then
     write(message, '(a,a,a)' )&
&     '  mffmem==0, iscf==5 => use of 1 FFT temporary disk file,',ch10,&
&     '                       10 times bigger than a DEN file.'
   else if(iscf==6)then
     write(message, '(a,a,a)' )&
&     '  mffmem==0, iscf==6 => use of 1 FFT temporary disk file,',ch10,&
&     '                       10 times bigger than a DEN file.'
   else if(iscf==7.or.iscf==17)then
     write(message, '(a,a,a)' )&
&     '  mffmem==0, iscf==7 => use of 1 FFT temporary disk file,',ch10,&
&     '                       (2+2*npulayit) times bigger than a DEN file.'
   end if
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Temporary message - estimation of PAW specific data has to be done...
!Have to add the usepaw argument to use this.
!if (usepaw==1) then
!write(message,'(5a)') '  WARNING: You are using PAW formalism;',ch10,&
!&       '           Above estimations do not take PAW',ch10,&
!&       '           specific data into account !'
!call wrtout(iout,message,'COLL')
!call wrtout(std_out,message,'COLL')
!end if

 write(message,'(80a,a)') ('=',mu=1,80),ch10
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

!--------------------------------------------------------------------
!Here, each processor must test its memory, so use
!the PERS mode for error messages, followed by synchronisation

 mbbigarr=max(mbf_fftgr,mbcg,mbgylm)
 if(mbbigarr==mbcg) then
   write(message, '(a,f12.4,a)' ) &
&   ' Biggest array : cg(disk), with',mbcg+tol10,' MBytes.'
 else if (mbbigarr==mbf_fftgr) then
   write(message, '(a,f12.4,a)' ) &
&   ' Biggest array : f_fftgr(disk), with',mbf_fftgr+tol10,' MBytes.'
 else if (mbbigarr==mbgylm)then
   write(message, '(a,f12.4,a)' ) &
&   ' Biggest array : pawfgrtab%gylm(gr), with',mbgylm+tol10,' MBytes.'
 end if
 call wrtout(std_out,message,'COLL')

 quit=0

 if(option>=1)then

!  Test the ability to allocate the biggest array
   nquarter_mbytes=4.0_dp*mbbigarr+1.0_dp
   allocate(bigarray(32*1024,nquarter_mbytes),STAT=ier)
   if(ier/=0)then
     write(message,'(a,a,a,a,f11.3,a,a,a,a,a,a,a)') ch10,&
&     ' memana : ERROR -',ch10,&
&     '  Test failed to allocate an array of',mbbigarr,' Mbytes',ch10,&
&     '  It is not worth to continue ',ch10,&
&     '  Action : modify input variable to fit the available memory,',ch10,&
&     '  or increase limit on maximal array size.'
     call wrtout(std_out,message,'PERS')
     if(option==1)then
       call leave_new('PERS')
     else
       quit=1
     end if
   end if
   if(mpi_enreg%paral_compil_kpt==1)then
!    BEGIN TF_CHANGES
     call leave_test(mpi_enreg)
!    END TF_CHANGES
   end if
   if(quit==0)then
     write(message,'(a,f11.3,a)')&
&     ' memana : allocated an array of',mbbigarr+tol10,&
&     ' Mbytes, for testing purposes. '
     call wrtout(std_out,message,'COLL')
   end if
   if(allocated(bigarray))deallocate(bigarray)

!  Test the ability to allocate the needed total memory : use 8 segments,
!  hoping that the maximal segment size is not so much smaller than the
!  total memory
   nquarter_mbytes=0.5_dp*mbbiggest+1.0_dp
   allocate(bigarray1(32*1024,nquarter_mbytes),STAT=ier1)
   allocate(bigarray2(32*1024,nquarter_mbytes),STAT=ier2)
   allocate(bigarray3(32*1024,nquarter_mbytes),STAT=ier3)
   allocate(bigarray4(32*1024,nquarter_mbytes),STAT=ier4)
   allocate(bigarray5(32*1024,nquarter_mbytes),STAT=ier5)
   allocate(bigarray6(32*1024,nquarter_mbytes),STAT=ier6)
   allocate(bigarray7(32*1024,nquarter_mbytes),STAT=ier7)
   allocate(bigarray8(32*1024,nquarter_mbytes),STAT=ier8)
   if(ier1/=0 .or. ier2/=0 .or. ier3/=0 .or. ier4/=0 .or.&
&   ier5/=0 .or. ier6/=0 .or. ier7/=0 .or. ier8/=0      )then
     write(message,'(a,a,a,a,f11.3,a,a,a,a,a,a,a)') ch10,&
&     ' memana : ERROR -',ch10,&
&     '  Test failed to allocate ',mbbiggest,' Mbytes',ch10,&
&     '  It is not worth to continue ',ch10,&
&     '  Action : modify input variable to fit the available memory.',ch10,&
&     '  or increase limit on available memory.'
     call wrtout(std_out,message,'PERS')
     if(option==1)then
       call leave_new('PERS')
     else
       quit=1
     end if
   end if

   if(quit==0)then
     write(message,'(a,f11.3,a,a,a)')&
&     ' memana : allocated ',mbbiggest,&
&     ' Mbytes, for testing purposes. ',ch10,&
&     ' The job will continue.'
     call wrtout(std_out,message,'COLL')
   end if
   if(allocated(bigarray1))deallocate(bigarray1)
   if(allocated(bigarray2))deallocate(bigarray2)
   if(allocated(bigarray3))deallocate(bigarray3)
   if(allocated(bigarray4))deallocate(bigarray4)
   if(allocated(bigarray5))deallocate(bigarray5)
   if(allocated(bigarray6))deallocate(bigarray6)
   if(allocated(bigarray7))deallocate(bigarray7)
   if(allocated(bigarray8))deallocate(bigarray8)

 end if

!--------------------------------------------------------------------

 if(option==2 .and. quit==1 )then

!  Estimation of the available memory
!  
!  A quarter of Mbyte is 256*1024/8 real(dp) numbers,
!  that is 32*1024 dp numbers.
!  One begins with the allocation of 4 Mbytes. If successful,
!  one increases that number, until the allocation is not successfull
!  any more. Unfortunately, on a P6 with the pghpf compiler, the
!  allocate instruction generate a core dump, instead of returning
!  an error code, so that this part of code has been made optional.

   nquarter_mbytes=16
   nmbytes=nquarter_mbytes/4.0_dp

!  With an increase ratio of 1.25_dp (see below), ii=5 leads to 9 MB,
!  ii=10 leads to 28 MB, ii=15 leads to 85 MB, ii=18 leads to 165 MB,
!  ii=30 is over 2 GB
   do ii=1,30
     allocate(bigarray(32*1024,nquarter_mbytes),STAT=ier)
     if(ier/=0)then
       write(message,'(a,i9,a)') &
&       ' memana : failed to allocate ',nmbytes,' Mbytes'
       call wrtout(std_out,message,'PERS')
       exit
     end if
     write(message,'(a,i9,a)')&
&     ' memana : succeeded to allocate ',nmbytes,' Mbytes'
     call wrtout(std_out,message,'PERS')
!    Here really test the space
!    do kk=1,nquarter_mbytes
!    do jj=1,32*1024,37
!    bigarray(jj,kk)=0.0_dp
!    end do
!    write(6,*)' memana : wrote ',kk,' quarter of mbytes'
!    end do
     deallocate(bigarray)
     nquarter_mbytes=dble(nquarter_mbytes)*1.25_dp
     nmbytes=nquarter_mbytes/4.0_dp
   end do
   if(allocated(bigarray))deallocate(bigarray)
   call leave_new('PERS')
   if(mpi_enreg%paral_compil_kpt==1)then
!    BEGIN TF_CHANGES
     call leave_test(mpi_enreg)
!    END TF_CHANGES
   end if

!  End the test of the available memory
 end if

!--------------------------------------------------------------------

 deallocate(cdpfftf,cdpfft,cdpmpw,cintfftf,cintfft,cintmpw,cdpadd,cintadd)
 deallocate(mbdpadd,mbdpfftf,mbdpfft,mbdpmpw,mbintadd,mbintfftf,mbintfft,mbintmpw)
 deallocate(mbother,mbtot)

!DEBUG
!write(6,*)' memana : exit '
!ENDDEBUG

end subroutine memana
!!***
