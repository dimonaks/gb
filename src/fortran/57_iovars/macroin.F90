!{\src2tex{textfont=tt}}
!!****f* ABINIT/macroin
!! NAME
!! macroin
!!
!! FUNCTION
!! Treat "macro" input variables, that can :
!! - initialize several other input variables for one given dataset
!! - initialize several other input variables for a set of datasets.
!! Note that the treatment of these different types of macro input variables is different.
!! Documentation of such input variables is very important, including the
!! proper echo, in the output file, of what such input variables have done. 
!!
!! Important information : all the "macro" input variables should be properly
!! identifiable to be so, and it is proposed to make them start with the string "macro".
!!
!! COPYRIGHT
!! Copyright (C) 2009-2010 ABINIT group (XG)
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
!!   some of which are given a value here.
!!   The dataset with number 0 should NOT be modified in the present routine.
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine macroin(dtsets,ndtset_alloc)

 use defs_basis
 use defs_datatypes
 use defs_abitypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc
!arrays
 type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------------
!scalars
 integer :: idtset

!******************************************************************
!
!DEBUG
!write(6,*)' macroin1 : enter '
!ENDDEBUG

 do idtset=1,ndtset_alloc
   
   if (dtsets(idtset)%macro_uj>0) then
     dtsets(idtset)%irdwfk   = 1        ! preconverged wave function compulsory 
!    dtsets(idtset)%nline    = maxval((/ int(dtsets(idtset)%natom/2) , 6 /))   ! using default value: \DeltaU< 1%
!    dtsets(idtset)%nnsclo   = 4        ! using default value: \DeltaU< 1% 
     dtsets(idtset)%tolvrs   = 10d-8    ! convergence on the potential; 10d-8^= 10d-5 on occupation 
     dtsets(idtset)%diemix   = 0.45_dp  ! fastest convergence: dn= E^(-istep * 0.229 )
     dtsets(idtset)%dmatpuopt= 3        ! normalization of the occupation operator 
!    dtsets(idtset)%nstep    = 255      ! expected convergence after 10 \pm 3, 30 as in default normally suficient
!    dtsets(idtset)%iscf     = 17       ! mixing on potential, 17: default for PAW
   end if ! macro_uj

 end do
 
!DEBUG
!write(6,*)' macroin1 : exit '
!ENDDEBUG

end subroutine macroin
!!***
