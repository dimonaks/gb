!{\src2tex{textfont=tt}}
!!****f* ABINIT/outqmc
!! NAME
!! outqmc
!!
!! FUNCTION
!! Write the wave function to a file in 'pwfn.data' format. This file can be
!! read by the Cambridge quantum Monte Carlo program 'CASINO' and used as
!! trial wave function input for a variational or diffusion Monte Carlo
!! calculation.
!!
!! See www.tcm.phy.cam.ac.uk/~mdt26/casino.html for more details.
!!
!! M.D.Towler (mdt26 at cam.ac.uk)
!! November 2003
!! N.D.M.Hine (nicholas.hine at imperial.ac.uk)
!! November 2004
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wavefunction coefficients
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mkmem=maximum number of k-points in core memory
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum number of plane waves
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nspinor=number of spinorial components of the wavefunctions(unsupported)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k
!!  psps <type(pseudopotential_type)>=all the information about psps
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation.
!!
!! OUTPUT
!! Writes the file pwfn.data
!!
!! PARENTS
!!      gstate,gstateimg
!!
!! CHILDREN
!!      leave_new,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outqmc(cg,dtset,eigen,gprimd,hdr,kg,mband,mkmem,mpi_enreg,&
                  &mpw,nkpt,npwarr,nspinor,nsppol,occ,psps,results_gs)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_57_iovars, except_this_one => outqmc
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer :: mband,mkmem,mpw,nkpt,nspinor,nsppol
 type(dataset_type) :: dtset
 type(hdr_type) :: hdr
 type(mpi_type) :: mpi_enreg
 type(pseudopotential_type) :: psps
 type(results_gs_type) :: results_gs
!arrays
 integer :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),eigen(mband*nkpt*nsppol)
 real(dp) :: gprimd(3,3),occ(mband*nkpt*nsppol)

!Local variables -------------------------
!scalars
 integer,parameter :: r2s_length=80
 integer :: iband,iband_kpt_shift,iband_sppol_shift,icg,icg_shift,icgfull
 integer :: icgfull_shift,ierr=0,ii,ikg,ikg_shift,ikgfull,ikpt,ikpt_shift,io
 integer :: iocc,isppol,jj,nband1,nband2,nelecs,nkgfull
 real(dp) :: norm
 logical :: already_open,am_master,foundkg
 character(50) :: pwfnfilename
 character(500) :: message
 character(80) :: dft_functional,pseudo_name,pseudo_type
 character(r2s_length) :: tmpr,tmpr2,tmpr3
!arrays
 integer :: kgfull(3,mpw*mkmem),kgmap(mpw*mkmem)
 real(dp) :: gcart_qmc(3),kptscart_qmc(3,nkpt)
 real(dp),allocatable :: cgfull(:,:),xcart_qmc(:,:)
! *********************************************************************

!Go away if I am not the master node.
 write(message,'(a,a)')ch10,' outqmc : enter '
 call wrtout(ab_out,message,'PERS')

 am_master=.true.
 if(mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then
   if(mpi_enreg%me/=0) am_master=.false.
   if(mpi_enreg%paral_level>1)then
     write(message,'(a,a,a,i5,a)')' outqmc : ERROR -',ch10,&
&     '  paral_level=',mpi_enreg%paral_level,&
&     ' : Parallelization over k-points is not currently supported'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
 end if
 if(.not.am_master)return

 do io=10,99
   inquire(unit=io,opened=already_open,iostat=ierr)
   if(ierr/=0.or.io>99)then
     write(message,'(a)')' outqmc : ERROR - problem finding free unit number.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
   if(.not.already_open)exit
 end do

 write(message,'(a)')' outqmc : will open CASINO pwfn.data file '
 call wrtout(ab_out,message,'PERS')

!DEBUG
!write(6,*)ch10,'outqmc: DEBUG: dtset%ndtset = ',dtset%ndtset,ch10
!ENDDEBUG
!Open CASINO pwfn.data file
 if (dtset%ndtset<2)then
   pwfnfilename='pwfn.data'
 else
   pwfnfilename='pwfn'//trim(i2s(dtset%jdtset))//'.data'
 end if
 open(io,file=pwfnfilename,form='formatted',recl=300,status='unknown', &
& iostat=ierr)
 if(ierr/=0)then
   write(message,'(a)')' outqmc : ERROR - Unable to open the pwfn.data file.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!Check if the full set of k vectors has been used in this calculation
 if (dtset%kptopt==1) then
   close(io,status='delete')
   write(message,'(a,a,a)')' outqmc : ERROR - kptopt=1 so k-points have been ',&
&   'generated in the irreducible Brillouin Zone only. ',&
&   'Set kptopt=2 to obtain full set of k-points.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!Check if the full set of G vectors has been used in this calculation
 do ikpt=1,nkpt
   if (dtset%istwfk(ikpt)/=1) then
     close(io,status='delete')
     write(message,'(a,a,a,i5,a,i2,a,a,a,a,a)')' outqmc : ERROR -',ch10,&
&     '  istwfk(',ikpt,')=',dtset%istwfk(ikpt),' (ie /= 1) so some ',&
&     'G-vectors may not be present.',ch10,'  Set istwfk=1 for each ',&
&     'k-point to obtain full set.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
 end do !ikpt

 write(message,'(a)')' outqmc : QMC trial wave function file for CASINO generated by ABINIT'
 call wrtout(ab_out,message,'PERS')

!Write the required quantities to pwfn.data
 write(io,"('QMC trial wave function file for CASINO generated by ABINIT (www.abinit.org).')")

 write(io,fmt="(/'BASIC INFO'/'----------')")
 write(io,fmt="('Generated by:')")
 write(io,fmt="(' ABINIT ',a)")trim(hdr%codvsn)
 write(io,fmt="('Method:'/' DFT')")

 write(io,fmt="('DFT Functional:')")
 select case (dtset%ixc)
   case(0)
     dft_functional='No exchange-correlation.'
   case(1)
     dft_functional='L(S)DA (Teter/Pade parametrization)'
   case(2)
     dft_functional='LDA (Perdew-Zunger-Ceperley-Alder parametrization)'
   case(3)
     dft_functional='LDA (old Teter rational polynomial parametrization)'
   case(4)
     dft_functional='LDA (Wigner)'
   case(5)
     dft_functional='LDA (Hedin-Lundqvist)'
   case(6)
     dft_functional='LDA (X-alpha)'
   case(7)
     dft_functional='L(S)DA (Perdew-Wang 92)'
   case(8)
     dft_functional='L(S)DA (Perdew-Wang 92, exchange-only)'
   case(9)
     dft_functional='L(S)DA (Perdew-Wang 92, exchange- and RPA-correlation)'
   case(10)
     dft_functional='Diff. between ixc=7 and 9; use with accurate RPA corr. energy'
   case(11)
     dft_functional='GGA (Perdew-Burke-Ernzerhof)'
   case(12)
     dft_functional='GGA (Perdew-Burke-Ernzerhof, exchange-only)'
   case(13)
     dft_functional='GGA (potential: van Leeuwen-Baerends ; energy: Perdew-Wang 92)'
   case(14)
     dft_functional='GGA (RPBE of Zhang and Yang)'
   case(15)
     dft_functional='GGA (RPBE of Hammer, Hansen and Norskov)'
   case(16)
     dft_functional='GGA (HTCH)'
   case(17)
     dft_functional='Not defined (as of 11/2003).'
   case(18)
     dft_functional='Not defined (as of 11/2003).'
   case(19)
     dft_functional='Not defined (as of 11/2003).'
   case(20)
     dft_functional='Fermi-Amaldi xc for TDDFT'
   case(21)
     dft_functional='Fermi-Amaldi xc for TDDFT with LDA xc kernel'
   case(22)
     dft_functional='Fermi-Amaldi xc for TDDFT with Burke-Petersilka-Gross hybrid xc kernel'
     case default
     dft_functional='Unknown type.'
 end select
 write(io,"(' ABINIT ixc = ',i2,' : ',a)")dtset%ixc,trim(dft_functional)

 write(io,"('Pseudopotential (of first atom type)')")
 select case(psps%pspcod(1))
   case(1)
     pseudo_type='ABINIT type 1' ; pseudo_name='Troullier-Martins'
   case(2)
     pseudo_type='ABINIT type 2' ; pseudo_name='Goedecker-Teter-Hutter (GTH)'
   case(3)
     pseudo_type='ABINIT type 3' ; pseudo_name='Hartwigsen-Goedecker-Hutter'
   case(4)
     pseudo_type='ABINIT type 4'
     pseudo_name='Teter pseudo generated using the ATOM code'
   case(5)
     pseudo_type='ABINIT type 5'
     pseudo_name='"Phoney" pseudo built on a Hamman grid'
   case(6)
     pseudo_type='ABINIT type 6'
     pseudo_name='Fritz-Haber Institut (Troullier Martins)'
     case default
     pseudo_type='Unknown pseudopotential type (as of 11/2003).' ; pseudo_name=''
 end select
 if(dtset%ixc<10)then
   write(io,"(1x,a,2x,': ',a)")trim(pseudo_type),trim(pseudo_name)
 else
   write(io,"(1x,a,3x,': ',a)")trim(pseudo_type),trim(pseudo_name)
 end if

 write(io,"('Plane wave cutoff (au)')")
 tmpr=r2s(hdr%ecut_eff,'(f12.3)')
 write(io,'(1x,a)')trim(tmpr)

 write(io,"('Spin polarized:')")
 select case(dtset%nspden)
   case(1)
     write(io,"(' .false.')")
   case(2)
     write(io,"(' .true.')")
   case(4)
     close(io,status='delete')
     write(message,'(a)')' outqmc : ERROR - nspden=4 but CASINO cannot yet deal with non-collinear spins.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
     case default
     close(io,status='delete')
     write(message,'(a)')' outqmc : ERROR - Unrecognized value of nspden.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
 end select

 write(io,"('Total energy (au per primitive cell)')")
 tmpr=r2s(results_gs%etotal,'(f24.14)')
 write(io,'(1x,a)')trim(tmpr)
 write(io,"('Kinetic energy')")
 tmpr=r2s(results_gs%energies%e_kinetic,'(f24.14)')
 write(io,'(1x,a)')trim(tmpr)
 write(io,"('Local potential energy (Local pseudopotential energy eei + pseudopotential core-core energy eii)')")
 tmpr=r2s((results_gs%energies%e_localpsp+results_gs%energies%e_corepsp),'(f24.14)')
 write(io,'(1x,a)')trim(tmpr)
 write(io,"('Non-local potential energy')")
 tmpr=r2s(results_gs%energies%e_nonlocalpsp,'(f24.14)')
 write(io,'(1x,a)')trim(tmpr)
 write(io,"('Electron-electron energy (Hartree Energy + Exchange-Correlation Energy)')")
 tmpr=r2s((results_gs%energies%e_hartree+results_gs%energies%e_xc),'(f24.14)')
 write(io,'(1x,a)')trim(tmpr)
 write(io,"('Ion-ion energy')")
 tmpr=r2s(results_gs%energies%e_ewald,'(f24.14)')
 write(io,'(1x,a)')trim(tmpr)
 write(io,"('Number of electrons per primitive cell')")

 nelecs=0
 do ii=1,dtset%natom
   nelecs=nelecs+psps%ziontypat(dtset%typat(ii))
 end do

 write(io,'(1x,i3)')nelecs
 write(io,*)
 write(io,"('GEOMETRY'/'--------')")
 write(io,"('Number of atoms per primitive cell')")
 write(io,'(1x,i3)')dtset%natom
 write(io,"('Atomic numbers and positions of atoms (au)')")

 allocate(xcart_qmc(3,dtset%natom))
 call xredxcart(dtset%natom,1,hdr%rprimd,xcart_qmc,hdr%xred)
 do ii=1,dtset%natom
   tmpr=r2s(xcart_qmc(1,ii),'(f24.14)')
   tmpr2=r2s(xcart_qmc(2,ii),'(f24.14)')
   tmpr3=r2s(xcart_qmc(3,ii),'(f24.14)')
   jj=psps%znucltypat(dtset%typat(ii))
   write(io,'(1x,i2,3(1x,a))')jj,trim(tmpr),trim(tmpr2),trim(tmpr3)
 end do
 deallocate(xcart_qmc)

 write(io,"('Primitive lattice vectors (au)')")
 do ii=1,3
   tmpr=r2s(hdr%rprimd(1,ii),'(f24.14)')
   tmpr2=r2s(hdr%rprimd(2,ii),'(f24.14)')
   tmpr3=r2s(hdr%rprimd(3,ii),'(f24.14)')
   write(io,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
 end do

!Copy the G vectors for the first k point into kgfull
 ikgfull=0
 do ikg=1,npwarr(1)
   ikgfull=ikgfull+1
   kgfull(1:3,ikgfull) = kg(1:3,ikg)
   kgmap(ikg)=ikgfull
 end do
 ikg_shift = npwarr(1)
!Go through the other k points and look for any G vectors that haven't
!already been found and add them to the end of kgfull
 do ikpt=2,nkpt
   do ikg=ikg_shift,ikg_shift+npwarr(ikpt)
     foundkg = .false.
     do ii=1,ikgfull
       if(kg(1,ikg)==kgfull(1,ii).and.kg(2,ikg)==kgfull(2,ii) &
&       .and.kg(3,ikg)==kgfull(3,ii)) then
         foundkg=.true.
         kgmap(ikg)=ii
         exit
       end if
     end do
     if(.not.foundkg)then
       ikgfull=ikgfull+1
       kgfull(1:3,ikgfull)=kg(1:3,ikg)
       kgmap(ikg)=ikgfull
     end if
   end do
   ikg_shift=ikg_shift+npwarr(ikpt)
 end do
 nkgfull=ikgfull

 write(io,*)
 write(io,"('G VECTORS'/'---------')")
 write(io,"('Number of G-vectors')")
 write(io,'(1x,i7)')nkgfull
 write(io,"('Gx Gy Gz (au)')")

 do ikgfull=1,nkgfull
   gcart_qmc=2*pi*(kgfull(1,ikgfull)*gprimd(1:3,1)&
&   +kgfull(2,ikgfull)*gprimd(1:3,2)+kgfull(3,ikgfull)*gprimd(1:3,3))
   write(io,*)gcart_qmc(1:3) ! '(3e26.16)'
 end do

!Populate the cgfull array, using the mapping in kgmap between the
!coefficients for kg in the per-kpoint list and the ones in the full list
!The number of xxx_shift's might seem excessive but the re-ordering of the
!list from (spin, kpt, band, kg) to (kpt, spin, band, kgfull) is quite
!complicated
 allocate(cgfull(2,nkgfull*nspinor*nsppol*mband*nkpt))
 cgfull(1:2,1:nkgfull*nspinor*nsppol*mband*nkpt)=0
 icg_shift=1
 do isppol=1,nsppol
   ikg_shift=1
   if(isppol==2)then
!    Go back to the beginning of cgfull but skip the first set of isppol=1 bands
     icgfull_shift=nkgfull*dtset%nband(1)
     ikpt_shift=nkpt
   else
     icgfull_shift=0 ! Start at the beginning of cgfull
     ikpt_shift=0
   end if
   do ikpt=1,nkpt
     do iband=1,dtset%nband(ikpt+ikpt_shift)
       ikg=ikg_shift
!      Find the index in Abinit's coefficient list
       do icg=icg_shift,icg_shift+npwarr(ikpt)-1
!        Map it to an index in the full CASINO list with the mapping recorded
!        when kgfull was read in
         icgfull = kgmap(ikg)+icgfull_shift
         cgfull(1:2,icgfull)=cg(1:2,icg)
         ikg=ikg+1
       end do !icg
       icg_shift=icg_shift+npwarr(ikpt)
       icgfull_shift=icgfull_shift+nkgfull
     end do !iband
     if(isppol==2)then
!      Skip the isppol==1 bands
       icgfull_shift=icgfull_shift+nkgfull*dtset%nband(ikpt)
     else
       if(nsppol==2)then
!        Skip the isppol==2 bands
         icgfull_shift=icgfull_shift+nkgfull*dtset%nband(ikpt+nkpt)
       else
         icgfull_shift=icgfull_shift ! Nothing to be skipped
       end if
     end if
     ikg_shift=ikg_shift+npwarr(ikpt)
   end do !ikpt
 end do !isppol

!See if each orbital is normalised and check for integer occupancy of orbitals.
!These are checked by CASINO and it will complain if they are not as expected.
 icgfull_shift=1
 ii=0
 iocc=1
 call wrtout(std_out,message,'COMMENT')
 do ikpt=1,nkpt
   do isppol=1,nsppol
     do iband=1,dtset%nband(ikpt+(isppol-1)*nkpt)
       if(occ(iocc)/=int(occ(iocc)))then
         write(message,'(a,a,a,i5,a,i1,a,i5,a,f10.8,a)')' outqmc : WARNING -',&
&         ch10,'Non-integer occupation number for kpt ',ikpt,', sppol ',&
&         isppol,', band ',iband,': occ=',occ(iocc),'.'
         call wrtout(std_out,message,'WARNING')
       end if
       iocc=iocc+1
       norm=0
       do icgfull=icgfull_shift,icgfull_shift+nkgfull-1
         norm=norm+cgfull(1,icgfull)**2+cgfull(2,icgfull)**2
       end do !icgfull
       icgfull_shift=icgfull_shift+nkgfull
       if((norm<0.999).or.(norm>1.001))then
         write(message,'(a,a,a,i5,a,i1,a,i5,a,f10.8,a)')' outqmc : WARNING -',&
&         ch10,'Incorrectly normalised orbital for kpt ',ikpt,', sppol ',&
&         isppol,', band ',iband,': norm=',norm,'.'
         call wrtout(std_out,message,'WARNING')
       end if
     end do !iband
   end do !isppol
 end do !ikpt

 write(io,*)
 write(io,"('WAVE FUNCTION'/'-------------')")
 write(io,"('Number of k-points')")
 write(io,'(1x,i5)')nkpt

 do ikpt=1,nkpt
   kptscart_qmc(1:3,ikpt)=2*pi*(dtset%kptns(1,ikpt)*gprimd(1:3,1)&
&   +dtset%kptns(2,ikpt)*gprimd(1:3,2)+dtset%kptns(3,ikpt)*gprimd(1:3,3))
 end do !ikpt

 iband_kpt_shift=0
 icg_shift=0
 icgfull_shift=0
 do ikpt=1,nkpt
   write(io,"('k-point # ; # of bands (up spin/down spin) ; k-point coords (au)')")
   if(nsppol==2)then
     nband1=dtset%nband(ikpt)
     nband2=dtset%nband(ikpt+nkpt)
   else
     nband1=dtset%nband(ikpt)
     nband2=0
   end if
   tmpr=r2s(kptscart_qmc(1,ikpt),'(f24.14)')
   tmpr2=r2s(kptscart_qmc(2,ikpt),'(f24.14)')
   tmpr3=r2s(kptscart_qmc(3,ikpt),'(f24.14)')
   write(io,'(3(1x,i5),3(1x,a))')ikpt,nband1,nband2,trim(tmpr),trim(tmpr2), &
&   trim(tmpr3)
   do isppol=1,nsppol
     if (isppol==2) then
       iband_sppol_shift=sum(dtset%nband(1:nkpt))
       iband_kpt_shift=sum(dtset%nband((nkpt+1):(nkpt+ikpt-1)))
     else
       iband_sppol_shift=0
       iband_kpt_shift=sum(dtset%nband(1:(ikpt-1)))
     end if
     do iband=1,dtset%nband(ikpt)
       write(io,"('Band, spin, eigenvalue (au)')")
       tmpr=r2s(eigen(iband_kpt_shift+iband_sppol_shift+iband),'(f24.14)')
       write(io,'(2(1x,i5),1x,a)')iband,isppol,tmpr
       write(io,"('Eigenvector coefficients')")
       do icgfull=1,nkgfull
         write(io,"(1x,'(',e23.16,',',e23.16,')')")cgfull(1:2,icgfull+icgfull_shift)
       end do !icgfull
       icgfull_shift=icgfull_shift+nkgfull
     end do !iband
   end do !isppol
 end do !ikpt

 write(message,'(a,a)')' outqmc : done with writing of QMC trial wave function file for CASINO',ch10
 call wrtout(ab_out,message,'PERS')

end subroutine outqmc
!!***

!!****f* ABINIT/i2s
!! NAME
!! i2s
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

 function i2s(n)
!----------------------------------------------------------------------- !
!I2S                                                                   !
!===                                                                   !
!Convert integers to left justified strings that can be printed in the !
!middle of a sentence without introducing large amounts of white space. !
!!
!Calling routine is intended to include something like:                !
!integer i                                                             !
!i=12                                                                  !
!write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!!
!MDT 2002                                                              !
!----------------------------------------------------------------------- !

 implicit none

!Arguments ----------------------
 integer, intent(in) :: n
 character(len=20) :: i2s

!Local variables ----------------
 integer :: i,j
 character :: tmp,sign

! *********************************************************************

 if(n==0)then
   i2s='0' ; return
 end if
 sign=' ' ; if(n<0)sign='-'

 do i=1,len(i2s)
   i2s(i:i)=' '
 end do

 i=abs(n)
 do j=1,len(i2s)
   if(i==0)exit
   i2s(j:j)=achar(ichar('0')+mod(i,10))
   i=i/10
 end do

 i=1 ; j=len_trim(i2s)
 do
   if(i>=j)exit
   tmp=i2s(j:j)
   i2s(j:j)=i2s(i:i)
   i2s(i:i)=tmp
   i=i+1
   j=j-1
 end do

 i2s=trim(sign)//i2s

end function i2s
!!***

!!****f* ABINIT/r2s
!! NAME
!! r2s
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE
 function r2s(r,real_format)
!------------------------------------------------------------------------- !
!Converts real variable with arbitrary format to string that can be      !
!trimmed and printed in the middle of a sentence without introducing     !
!large amounts of white space, as you would if you did                   !
!write(6,'(f12.6)')12.0 or similar. Note you need to pass through the    !
!format string e.g. f12.6 .                                              !
!!
!Calling routine is intended to include something like:                  !
!USE utilities                                                           !
!REAL(dp) r                                                              !
!r=12._dp                                                                 !
!tmpr=r2s(r,'(f12.6)')                                                   !
!write(6,*)'Real number ',trim(tmpr),' with words at the end.'           !
!!
!Note : DON'T USE R2S IN A WRITE STATEMENT SINCE THIS IS ILLEGAL         !
!IN FORTRAN90 (ALTHOUGH NOT IN FORTRAN200X). IF ANYONE HAS TIME, FEEL    !
!FREE TO WRITE A VERSION OF THIS WHICH ISN'T ILLEGAL - SIMILAR TO        !
!I2S ABOVE - SO THAT PEOPLE WHO HAVEN'T READ THIS NOTE DON'T FEEL        !
!TEMPTED TO CALL R2S IN A WRITE STATEMENT.                               !
!!
!MDT 2002                                                                !
!------------------------------------------------------------------------- !

 use defs_basis

 implicit none

!Arguments -------------------------
 real(dp),intent(in) :: r
 character(len=*),intent(in) :: real_format
 character(len=80) :: r2s

!Local variables -------------------

! *********************************************************************


 if(len(r2s)>0)then
   write(r2s,real_format)r
   r2s=adjustl(r2s)
 end if

end function r2s
!!***
