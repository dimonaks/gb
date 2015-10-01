!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_memory
!! NAME
!! wvl_memory
!!
!! FUNCTION
!! Estimation of the memory needed for waelet based computation job.
!! According to the value of the option variable,
!! might also try to allocate this amount of memory, and if it fails,
!! might estimate the available memory.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset=<type datafiles_type>contains all input variables.
!!  idtset=number of the current dataset
!!  mpi_enreg=informations about MPI parallelization
!!  npsp=number of pseudopotentials
!!  option : if 0 , no test of available memory
!!           if 1 , the routine tries to allocate the estimated memory, for testing
!!                    purposes, and if a failure occurs, the routine stops.
!!           if 2 , like 1, but before stopping, the routine will provide
!!                    an estimation of the available memory.
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! The estimator is the one provided by BigDFT.
!!
!! PARENTS
!!      invars2m
!!
!! CHILDREN
!!      dtsetcopy,dtsetfree,eleconf,leave_new,memoryestimator,mkradim,wrtout
!!      wvl_setboxgeometry,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_memory(dtset, idtset, mpi_enreg, npsp, option, pspheads)
  
  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API, only: eleconf, MemoryEstimator
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_53_abiutil
 use interfaces_57_iovars, except_this_one => wvl_memory
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idtset, npsp, option
  type(dataset_type),intent(in) :: dtset
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
  !scalars
  integer :: ityp, i, iasctype, mu, nstates, mxpl, mxchg
  character(len=2) :: symbol
  character(len=500) :: message
  real(dp) :: rcov, rprb, ehomo, radfine, peakmem
  integer, parameter :: nmax=6, lmax=3
  type(dataset_type) :: tmpDt
  !arrays
  integer :: neleconf(nmax, 0:lmax)
  real(dp) :: acell(3), rprimd(3,3), rprim(3,3)
  real(dp), allocatable :: radii_cf(:,:)
  real(dp), allocatable :: xred(:,:), xcart(:,:)
  character(len = 20), allocatable :: atomnames(:)

! **************************************************************************

 if(option<0 .or. option>2)then
   write(message, '(A,A,A,A,I0,A)') ch10,&
&   ' wvl_memory : BUG -',ch10,&
&   '  option=',option,' while the only allowed values are 0, 1, or 2.'
   call wrtout(std_out,message,'COLL')
 end if

 write(message,*)' wvl_memory : analysis of memory needs '
 call wrtout(std_out,message,'COLL')

 if(idtset/=0)then
   write(message,'(80a,a,a,i3,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need for DATASET', idtset,&
&   ' (WVL).'
 else
   write(message,'(80a,a,a,a)')('=',mu=1,80),ch10,&
&   ' Values of the parameters that define the memory need of the present run',&
&   ' (WVL).'
 end if
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'( a,f7.3,a,i7,2(a,F7.3),a,a,f7.3,a,i7,2(a,F7.3) )' ) &
& '  wvl_hgrid =', dtset%wvl_hgrid , '   nwfshist =', dtset%nwfshist, &
& ' wvl_crmult =', dtset%wvl_crmult, ' wvl_frmult =', dtset%wvl_frmult, ch10,&
& '  tl_radius =', dtset%tl_radius , '  tl_nprccg =', dtset%tl_nprccg, &
& ' wvl_cpmult =', dtset%wvl_cpmult, ' wvl_fpmult =', dtset%wvl_fpmult
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if (dtset%nsppol == 2) then
   nstates = dtset%nelect
 else
   nstates = dtset%mband
 end if
 write(message,'(4(a,i7))')&
& '      natom =', dtset%natom, '     ntypat =', dtset%ntypat, &
& '    nstates =', nstates,     '     nsppol =', dtset%nsppol
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'(80a)') ('=',mu=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

#if defined HAVE_BIGDFT
!First, use eleconf to get radii_cf().
 allocate(atomnames(dtset%ntypat))
 allocate(radii_cf(npsp, 2))
 do ityp = 1, npsp, 1
   call eleconf(int(pspheads(ityp)%znuclpsp), int(pspheads(ityp)%zionpsp), &
&   symbol, rcov, rprb, ehomo, neleconf, iasctype, mxpl, mxchg)
   write(atomnames(ityp), "(A)") symbol

!  new method for assigning the radii
   radii_cf(ityp, 1) = one / sqrt(abs(two * ehomo))
   radfine = 100.d0
   do i = 0, 4, 1
     if (pspheads(ityp)%GTHradii(i) /= zero) then
       radfine = min(radfine, pspheads(ityp)%GTHradii(i))
     end if
   end do
   radii_cf(ityp,2) = radfine
 end do

!Compute the shifted positions and acell
 acell = dtset%acell_orig(1:3,1)
 allocate(xred(3, dtset%natom))
 xred = dtset%xred_orig(:,:,1)
 rprimd = dtset%rprimd_orig(1:3,1:3,1)
 call dtsetcopy(tmpDt, dtset)
 tmpDt%prtvol = 1
 call wvl_setBoxGeometry(tmpDt, mpi_enreg%me, radii_cf, rprimd, xred)
!Compute acell and rprim from rprimd
 call mkradim(acell,rprim,rprimd)
 allocate(xcart(3, tmpDt%natom))
 call xredxcart(tmpDt%natom, 1, rprimd, xcart, xred)
 call MemoryEstimator(tmpDt%wvl%atoms%geocode, mpi_enreg%nproc, tmpDt%nwfshist, &
& tmpDt%wvl%n(1), tmpDt%wvl%n(2), tmpDt%wvl%n(3), &
& acell(1), acell(2), acell(3), &
& tmpDt%wvl%h(1), tmpDt%wvl%h(2), tmpDt%wvl%h(3), tmpDt%natom, tmpDt%ntypat, &
& tmpDt%typat, xcart, radii_cf, tmpDt%wvl_crmult, &
& tmpDt%wvl_frmult, tmpDt%mband, 0, atomnames, 0, &
& tmpDt%nsppol, peakmem)

 call dtsetfree(tmpDt)
 deallocate(radii_cf)
 deallocate(xred)
 deallocate(xcart)
 deallocate(atomnames)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_memory : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif

 write(message,'(80a,a)') ('=',mu=1,80), ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

end subroutine wvl_memory
!!***
