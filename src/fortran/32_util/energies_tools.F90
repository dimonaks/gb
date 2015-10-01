!{\src2tex{textfont=tt}}
!!****f* ABINIT/energies_tools
!! This module contains functions used to manipulate
!! type(energies_type) objects (initialization, copy,...)
!!
!!***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!!****f* ABINIT/energies_init
!!
!! NAME
!! energies_init
!!
!! FUNCTION
!! Set zero in all values of a type(energies_type) object.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!   energies <type(energies_type)>=values to initialise
!!
!! PARENTS
!!      driver,gstate,scfcv,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine energies_init(energies)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(out) :: energies

!Local variables-------------------------------

! *************************************************************************

!@energies_type

 energies%entropy       = zero
 energies%e_entropy     = zero
 energies%e_fermie      = zero
 energies%e_paw         = zero
 energies%e_pawdc       = zero
 energies%e_kinetic     = zero
 energies%e_localpsp    = zero
 energies%e_nonlocalpsp = zero
 energies%e_eigenvalues = zero
 energies%e_hartree     = zero
 energies%e_ewald       = zero
 energies%e_xc          = zero
 energies%e_vxc         = zero
 energies%e_xcdc        = zero
 energies%e_elecfield   = zero
 energies%e_corepsp     = zero
 energies%e_electronpositron   = zero
 energies%edc_electronpositron = zero
 energies%e0_electronpositron  = zero

end subroutine energies_init
!!***

!!****f* ABINIT/energies_copy
!!
!! NAME
!! energies_copy
!!
!! FUNCTION
!! Copy two type(energies_type) objects.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   energies_in <type(energies_type)>=input values (to copy)
!!
!! OUTPUT
!!   energies_out <type(energies_type)>=output values
!!
!! PARENTS
!!      setup_positron
!!
!! CHILDREN
!!
!! SOURCE

 subroutine energies_copy(energies_in,energies_out)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(in)  :: energies_in
 type(energies_type),intent(out) :: energies_out

!Local variables-------------------------------

!*************************************************************************

!@energies_type

 energies_out%entropy              = energies_in%entropy
 energies_out%e_entropy            = energies_in%e_entropy
 energies_out%e_fermie             = energies_in%e_fermie
 energies_out%e_paw                = energies_in%e_paw
 energies_out%e_pawdc              = energies_in%e_pawdc
 energies_out%e_kinetic            = energies_in%e_kinetic
 energies_out%e_localpsp           = energies_in%e_localpsp
 energies_out%e_nonlocalpsp        = energies_in%e_nonlocalpsp
 energies_out%e_eigenvalues        = energies_in%e_eigenvalues
 energies_out%e_hartree            = energies_in%e_hartree
 energies_out%e_ewald              = energies_in%e_ewald
 energies_out%e_xc                 = energies_in%e_xc
 energies_out%e_vxc                = energies_in%e_vxc
 energies_out%e_xcdc               = energies_in%e_xcdc
 energies_out%e_elecfield          = energies_in%e_elecfield
 energies_out%e_corepsp            = energies_in%e_corepsp
 energies_out%e_electronpositron   = energies_in%e_electronpositron
 energies_out%edc_electronpositron = energies_in%edc_electronpositron
 energies_out%e0_electronpositron  = energies_in%e0_electronpositron

end subroutine energies_copy
!!***
!!***
