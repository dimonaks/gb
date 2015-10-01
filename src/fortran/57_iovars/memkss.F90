!{\src2tex{textfont=tt}}
!!****f* ABINIT/memkss
!! NAME
!! memkss
!!
!! FUNCTION
!! This routine evaluates the additional amount of memory required
!! by routine 'outkss'.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mproj=maximum dimension for number of projection operators for each
!!   angular momentum for nonlocal pseudopotential
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points.
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms in unit cell.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine is not available for paw calculations
!!
!! PARENTS
!!      outkss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine memkss(mband,mgfft,mkmem,mpi_enreg,mproj,mpsang,mpw,natom,&
&                 ngfft,nkpt,nspinor,nsym,ntypat)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mkmem,mproj,mpsang,mpw,natom,nkpt,nspinor
 integer,intent(in) :: nsym,ntypat
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer(i8b) :: isize,memsize
 character(len=500) :: msg

! *********************************************************************
!
 isize=580+fnlen+4*(81+nkpt+9*nsym)+8*15    !non allocatable var.
 if(mpi_enreg%paral_compil_kpt==1)then
   isize=isize+4*4                           !kpt_distrb
 end if
 memsize=isize
 isize=isize+4*nkpt+12*mpw+20*nkpt*mpw      !nbasek,gbasek,cnormk,gcurr
 memsize=max(memsize,isize)
 if(mpi_enreg%paral_compil_kpt==1)then
   isize=isize+12*mpw*nkpt                   !ibuf1,ibuf2,rbuf1
   memsize=max(memsize,isize)
   isize=isize-12*mpw*nkpt                   !ibuf1,ibuf2,rbuf1
 end if
 isize=isize+40*mpw                         !gbase,cnorm
 memsize=max(memsize,isize)
 isize=isize-4*nkpt-20*mpw*nkpt             !nbasek,gbasek,cnormk
 isize=isize+4*mpw                          !insort
 memsize=max(memsize,isize)
 isize=isize-16*mpw                         !cnorm
 isize=isize+28*mpw+24*nsym                 !gbig,nshell,gshell
 memsize=max(memsize,isize)
 isize=isize+4*mpw                          !shlim
 memsize=max(memsize,isize)
 isize=isize-44*mpw-24*nsym                 !gcurr,gbase,gshell,insort,nshell
 isize=isize-4*mpw                          !shlim
 isize=isize+8*mpw*nspinor&
& +16*mpw*nspinor*(mpw*nspinor+1)       !eigval,eigvec
 memsize=max(memsize,isize)
 isize=isize+8*mpw+8*ngfft(4)&
& *ngfft(5)*ngfft(6)      !ts,vlocal
 memsize=max(memsize,isize)
 isize=isize+8*mgfft+4+28*mpw               !gbound,indpw_k,kg_k
 memsize=max(memsize,isize)
 if (mkmem==0) then
   isize=isize+4*(mpw+2*mgfft+4)             !indpw_disk
   memsize=max(memsize,isize)
   isize=isize-4*(mpw+2*mgfft+4)             !indpw_disk
 end if
 isize=isize+8*natom&
& +24*mpw*ntypat*mpsang*mproj      !phkxred,ffnl,kinpw
 memsize=max(memsize,isize)
 isize=isize+16*mpw*natom                   !ph3d
 memsize=max(memsize,isize)
 isize=isize+48*mpw*nspinor&
& +8*mpw*nspinor*(mpw*nspinor+1)        !pwave,subghg,gvnlg
 if (nspinor==2)&
& isize=isize+40*mpw*nspinor                !pwave_so,subghg_so
 memsize=max(memsize,isize)
 isize=isize+8*mpw*nspinor*(mpw*nspinor+1)  !ghg
 memsize=max(memsize,isize)
 isize=isize+8*ngfft(4)*ngfft(5)*ngfft(6)   !work
 memsize=max(memsize,isize)
 isize=isize-8*ngfft(4)*ngfft(5)*ngfft(6)   !work
 isize=isize-8*mgfft+4+28*mpw&              !gbound,indpw_k,kg_k
&-8*natom-24*mpw*ntypat*mpsang*mproj&  !phkxred,ffnl,kinpw
&-16*mpw*natom                        !ph3d
 isize=isize-48*mpw*nspinor&
& -8*mpw*nspinor*(mpw*nspinor+1)        !pwave,subghg,gvnlg
 if (nspinor==2)&
& isize=isize-40*mpw*nspinor                !pwave_so,subghg_so

 isize=isize+56*mpw*nspinor                 !cwork,rwork
 memsize=max(memsize,isize)
 isize=isize-56*mpw*nspinor                 !cwork,rwork
 isize=isize+112*mpw*nspinor                !cwork,rwork,iwork,ifail
 memsize=max(memsize,isize)
 isize=isize-112*mpw*nspinor                !cwork,rwork,iwork,ifail
 isize=isize-8*mpw*nspinor*(mpw*nspinor+1)  !ghg
 isize=isize+8*mband                        !occ_k
 memsize=max(memsize,isize)
 isize=isize-8*mband                        !occ_k
 isize=isize-8*mpw*nspinor&
& -16*mpw*nspinor*(mpw*nspinor+1)       !eigval,eigvec
 isize=isize-32*mpw-8*ngfft(4)&
& *ngfft(5)*ngfft(6)     !gbig,ts,vlocal
 if(mpi_enreg%paral_compil_kpt==1)then
   isize=isize-4*4                           !kpt_distrb
 end if
 isize=isize-580-fnlen-4*(81+nkpt+9*nsym)-8*15   !non allocatable var.
!
 write(msg,'(2a,f8.2,a)')ch10,&
& ' Additional amount of memory required by "outkss" routine=',memsize*b2Mb,' Mbytes.'
 call wrtout(std_out,msg,'COLL')
!
end subroutine memkss
!!***
