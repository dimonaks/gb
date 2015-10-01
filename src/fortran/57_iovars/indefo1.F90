!{\src2tex{textfont=tt}}
!!****f* ABINIT/indefo1
!! NAME
!! indefo1
!!
!! FUNCTION
!! Initialisation phase : defaults values for a first batch of input variables
!! (especially dimensions, needed to allocate other parts of dtsets, as well
!!  as other input variables whose existence is needed for other initialisations to proceed).
!!
!! COPYRIGHT
!! Copyright (C) 2009-2010 ABINIT group (XG,MM,FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!  dtset=<type datafiles_type>contains all input variables for one dataset,
!!   some of which are given a default value here.
!!
!! PARENTS
!!      invars1m,newsp
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine indefo1(dtset)

 use defs_basis
 use defs_abitypes

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(dataset_type),intent(out) :: dtset

!Local variables -------------------------------
!scalars

!******************************************************************
!
!Set up default values. All variables to be output in outvars.f
!should have a default, even if a nonsensible one can be
!chosen to garantee print in that routine.

!DEBUG
!write(6,*)' indefo1 : enter '
!ENDDEBUG

!Use alphabetic order

!A
 dtset%amu(:)=-one
 dtset%acell_orig(:,:)=zero
!B
 dtset%bandpp=1
 dtset%berryopt=0
 dtset%bfield(:)=zero
!C
!D
 dtset%densty(:,:)=zero
!E
 dtset%efield(:)=zero
!F
!G
 dtset%gw_nqlwl=0
!I
 dtset%iatfix(:,:)=0
 dtset%icoulomb=0
!J
 dtset%jellslab=0
!K
 dtset%kptopt=0
!L
 dtset%lexexch(:)=-1
 dtset%lpawu(:)=-1
!M
 dtset%mkmem=-1
 dtset%mkqmem=-1
 dtset%mk1mem=-1
!N
 dtset%nnos=0
 dtset%natpawu=0
 dtset%natsph=0
 dtset%natvshift=0
 dtset%nbdblock=1
 dtset%nconeq=0
 dtset%nkpt=-1
 dtset%nkptgw=0
 dtset%normpawu(:)=0
 dtset%npband=1
 dtset%npfft=1
 dtset%npimage=1
 dtset%npkpt=1
 dtset%nqptdm=0
 dtset%nspden=1
 dtset%nspinor=1
 dtset%nsppol=1
 dtset%nsym=0     ! Actually, this default value is not used : it is to be reimposed before each call to ingeo in invars1
!O
 dtset%optdriver=0
!P
 dtset%pawspnorb=0
!Q
!R
 dtset%rprim_orig(:,:,:)=zero
 dtset%rprim_orig(1,1,:)=one
 dtset%rprim_orig(2,2,:)=one
 dtset%rprim_orig(3,3,:)=one
!S
 dtset%slabzbeg=zero
 dtset%slabzend=zero
 dtset%so_psp(:)=1
 dtset%spinat(:,:)=zero
 dtset%symmorphi=1
!T
 dtset%typat(:)=0
!U
 dtset%usedmatpu=0
 dtset%useexexch=0
 dtset%usepawu=0
!V
 dtset%vel_orig(:,:,:)=zero
!W
 dtset%wfoptalg=0
!X
 dtset%xred_orig(:,:,:)=zero
!Y
!Z
 dtset%zeemanfield(:)=zero
!DEBUG
!write(6,*)' indefo1 : exit '
!ENDDEBUG

end subroutine indefo1
!!***
