!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_radmesh
!! NAME
!! m_radmesh
!!
!! FUNCTION
!! Module containing all the functions related to the radial meshes
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (MT,FJ,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  * Routines tagged with "@type_name" are strongly connected to the definition of the data type. 
!!    Strongly connected means that the proper functioning of the implementation relies on the 
!!    assumption that the tagged procedure is consistent with the type declaration.
!!    Every time a developer changes the structure "type_name" adding new entries, he/she has to make sure 
!!    that all the strongly connected routines are changed accordingly to accommodate the modification of the data type. 
!!    Typical examples of strongly connected routines are creation, destruction or reset methods.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_radmesh

 use defs_basis
 use m_errors

 use defs_datatypes, only : pawrad_type

 implicit none

 private

 public :: init_radmesh         ! Main creation method
 public :: nullify_radmesh      ! Nullify the object
 public :: destroy_radmesh      ! Frees the allocated memory
 public :: print_radmesh        ! Printout of the basic info
 public :: isame_radmesh        ! Checks whether two meshes are equivalent or have the same equation.
 public :: calc_slatradl        ! Calculates the radial part of Slater integrals.

 public :: compmesh            ! Computes all points of the radial mesh.
 public :: copymesh            ! Returns a copy of the mesh.
 public :: simp_gen            ! Performs integral on a given (generalized) grid using Simpson rule.
 public :: deducer0            ! Extrapolate r=0 value of a function from values near r=0.
 public :: nderiv_gen          ! Do corrected first (and -if requested- second) derivation on a given (generalized) grid.
 public :: poisson             ! Solves Poisson equation for angularly dependent charge distribution of angular momentum l.
 public :: ifromr              ! Retrieve the Index FROM a given R value in a radial grid.

 ! TODO: Might use bit flags, but all radmesh stuff should be encapsulated here!
 integer,private,parameter :: RMESH_LINEAR = 1
 integer,private,parameter :: RMESH_LOG1   = 2
 integer,private,parameter :: RMESH_LOG2   = 3
 integer,private,parameter :: RMESH_LOG3   = 4

CONTAINS  !===========================================================
!!***

!!****f* m_radmesh/init_radmesh
!! NAME
!!  init_radmesh
!!
!! FUNCTION
!!  Creation method for radial meshes.
!!
!! INPUTS
!!  mesh_size=Dimension of the radial mesh.
!!  mesh_type=Type of mesh
!!     1=regular grid: r(i)=(i-1)*AA
!!     2=logarithmic grid: r(i)=AA*(exp[BB*(i-1)]-1)
!!     3=logarithmic grid: r(i>1)=AA*exp[BB*(i-1)] and r(1)=0
!!     4=logarithmic grid: r(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!  rstep=Radial step of the mesh (AA parameter above)
!!  lstep=Exponential step of the mesh (BB parameter above)
!!    Needed only if mesh type is logarithmic.
!!  r_for_intg=Mesh size used in integrals computation
!!    Integrals will be computed up to rr(r_for_intg).
!!
!! OUTPUT
!!  Rmesh<pawrad_type>=The object completely initialized.
!!
!! PARENTS
!!      m_atom
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine init_radmesh(Rmesh,mesh_size,mesh_type,rstep,lstep,r_for_intg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mesh_size,mesh_type
 real(dp),intent(in) :: rstep,lstep
 real(dp),intent(in) :: r_for_intg
 type(pawrad_type),intent(out) :: Rmesh
!arrays
 
! *************************************************************************

 DBG_ENTER("COLL")

 !@pawrad_type
 call nullify_radmesh(Rmesh)

 Rmesh%mesh_size = mesh_size
 Rmesh%mesh_type = mesh_type
 Rmesh%rstep     = rstep
 Rmesh%lstep     = lstep

 allocate(Rmesh%rad    (Rmesh%mesh_size))
 allocate(Rmesh%radfact(Rmesh%mesh_size))
 allocate(Rmesh%simfact(Rmesh%mesh_size))

 call compmesh(Rmesh,r_for_intg)

 DBG_EXIT("COLL")

end subroutine init_radmesh
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/nullify_radmesh
!! NAME
!!  nullify_radmesh 
!!
!! FUNCTION
!!  Nullify all pointers in the object
!!
!! PARENTS
!!      m_atom,m_radmesh
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine nullify_radmesh(Rmesh)


 implicit none

!Arguments ------------------------------------
 type(pawrad_type),intent(inout) :: Rmesh
 
! *************************************************************************

 !@pawrad_type
 nullify(Rmesh%rad    )
 nullify(Rmesh%radfact)
 nullify(Rmesh%simfact)

end subroutine nullify_radmesh
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/destroy_radmesh
!! NAME
!!  destroy_radmesh
!!
!! FUNCTION
!!  Frees all memory allocated
!!
!! PARENTS
!!      m_atom
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine destroy_radmesh(Rmesh)


 implicit none

!Arguments ------------------------------------
 type(pawrad_type),intent(inout) :: Rmesh

!Local variables-------------------------------
 
! *************************************************************************

 DBG_ENTER("COLL")

 !@pawrad_type
 if (associated(Rmesh%rad    )) deallocate(Rmesh%rad    )
 if (associated(Rmesh%radfact)) deallocate(Rmesh%radfact)
 if (associated(Rmesh%simfact)) deallocate(Rmesh%simfact)

 DBG_EXIT("COLL")

end subroutine destroy_radmesh
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/print_radmesh
!! NAME
!!  print_radmesh
!!
!! FUNCTION
!!  Reports basic info on the object.
!!
!! INPUTS
!!  Rmesh<pawrad_type>=Object defining the radial mesh
!!  header=String for the header provided by the user.
!!  [unit]=Unit number for output, defaults to std_out
!!  [prtvol]=Verbosity level, minimal if not specified.
!!  [mode_paral]=Either "COLL" or "PERS". Passed to wrtout. Defaults to "COLL"
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_atom
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine print_radmesh(Rmesh,header,unit,prtvol,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(pawrad_type),intent(in) :: Rmesh

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg
 
! *************************************************************************

 !@pawrad_type
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=ch10//' ==== Info on the Radial Mesh ==== '
 if (PRESENT(header)) msg=ch10//' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 SELECT CASE (Rmesh%mesh_type)

 CASE (RMESH_LINEAR)
   write(msg,'(a,i4,a,g12.5)')&
&    ' - Linear mesh: r(i)=step*(i-1), size=',Rmesh%mesh_size,', step=',Rmesh%rstep

 CASE (RMESH_LOG1)
   write(msg,'(a,i4,2(a,g12.5))')&
&    ' - Logarithimc mesh: r(i)=AA*[exp(BB*(i-1))-1], size=',Rmesh%mesh_size,', AA=',Rmesh%rstep,' BB=',Rmesh%lstep

 CASE (RMESH_LOG2)
   write(msg,'(a,i4,2(a,g12.5))')&
&    ' - Logarithimc mesh: r(i)=AA*exp(BB*(i-2)), size=',Rmesh%mesh_size,', AA=',Rmesh%rstep,' BB=',Rmesh%lstep

 CASE (RMESH_LOG3)
   write(msg,'(a,i1,a,i4,a,g12.5)')&
&    ' - Logarithimc mesh: r(i)=-AA*ln(1-(i-1)/n), n=size=',Rmesh%mesh_size,', AA=',Rmesh%rstep

 CASE DEFAULT 
   msg = ' Unknown mesh type! Action : check your pseudopotential or input file.'
   MSG_ERROR(msg)
 END SELECT

 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol>1) then
   write(msg,'(a,i4)')' Mesh size for integrals = ',Rmesh%int_meshsz
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,g12.5)')' rmax=rad(mesh_size)     = ',Rmesh%rmax
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,g12.5)')' Value of stepint        = ',Rmesh%stepint
   call wrtout(my_unt,msg,my_mode)
 end if

end subroutine print_radmesh
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/isame_radmesh
!! NAME
!!  isame_radmesh
!!
!! FUNCTION
!!  Check two radial meshes, returns a logical flag defining
!!  whether the meshes have the same equation and an integer 
!!  flag
!!
!! INPUTS
!!  Rmesh1,Rmesh2<pawrad_type>=The two radial meshes.
!!
!! OUTPUT
!!  hasameq=.true. if the two meshes are defined by the same equation.
!!  whichdenser=
!!    * 0 if meshes are not compatible
!!    * 1 if Rmesh1 is denser than Rmesh2
!!    * 2 if Rmesh2 is denser than Rmesh1
!!
!! PARENTS
!!      m_atom,m_paw_slater
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine isame_radmesh(Rmesh1,Rmesh2,hasameq,whichdenser)


 implicit none

!Arguments ------------------------------------
 integer,intent(out) :: whichdenser 
 logical,intent(out) :: hasameq
 type(pawrad_type),intent(in) :: Rmesh1
 type(pawrad_type),intent(in) :: Rmesh2

!Local variables-------------------------------
 
! *************************************************************************

 !@pawrad_type
 whichdenser =0 ; hasameq=.FALSE.

 if (Rmesh1%mesh_type /= Rmesh2%mesh_type) RETURN
 
 SELECT CASE (Rmesh1%mesh_type)

 CASE (RMESH_LINEAR) !check for linear meshes
   hasameq = (Rmesh1%rstep == Rmesh2%rstep)
 
 CASE (RMESH_LOG1,&  !check for logarithmic meshes
&      RMESH_LOG2,&
&      RMESH_LOG3 )

   hasameq = (    Rmesh1%rstep == Rmesh2%rstep &
&            .and.Rmesh1%lstep == Rmesh2%lstep )

 CASE DEFAULT
   MSG_ERROR("Unknown mesh type")
 END SELECT

 ! === If meshes have same equation, check whether they are equal ===
 ! * Note that also int_meshsz must be equal
 if (hasameq) then 
   whichdenser= 1
   if (Rmesh2%mesh_size  > Rmesh1%mesh_size ) whichdenser = 2
   !if (Rmesh1%int_meshsz == Rmesh2%int_meshsz) whichdenser = 2
 end if

end subroutine isame_radmesh
!!***

!----------------------------------------------------------------------

!!****m* m_radmesh/calc_slatradl
!! NAME
!!  calc_slatradl
!!
!! FUNCTION
!!  Calculate the radial part of Slater integrals. See below.
!!
!! INPUTS
!!  ll= l quantum number in the expansion of the Coulomb term.
!!  mesh_size=Number of points on the radial mesh.
!!  ff1(radmesh), ff2(radmesh)= The two functions to be integrated.
!!  Pawrad <type(pawrad_type)>=Structure containing radial grid information.
!!
!! OUTPUT
!!  integral = 
!!   $ \dfrac{4\pi}{2L+1} \int ff1(r1) \dfrac{r_<^L}{r_>^{L+1}} ff2(r2) dr1 dr2 $
!!  where $r_< = min(r1,r2)$ and $r_> = Max(r1,r2)$.
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine calc_slatradl(ll,mesh_size,ff1,ff2,Pawrad,integral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!scalars
 integer,intent(in) :: mesh_size,ll
 real(dp),intent(out) :: integral
!arrays
 real(dp),intent(in) :: ff1(mesh_size),ff2(mesh_size)
 type(pawrad_type),intent(in) :: Pawrad

!Local variables ---------------------------------------
!scalars
 integer :: int_meshsz
 real(dp) :: qq
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: hh(:),gg(:)
 
! *************************************************************************

 if (mesh_size /= Pawrad%mesh_size) then 
  MSG_BUG("mesh_size /= Pawrad%mesh_size")
 end if

 allocate(hh(mesh_size),gg(mesh_size))
 !hh = zero
 !gg = zero

 int_meshsz=Pawrad%int_meshsz
 !$int_meshsz=Pawrad%mesh_size

 ! the line below requires hh as work array. 
 hh = ff2  

 ! TODO find where int_meshsz is calculated and if it can affects the results.
 if (int_meshsz<mesh_size) hh(int_meshsz+1:mesh_size)=zero

 call poisson(hh,ll,qq,Pawrad,gg)

 gg(2:mesh_size) = gg(2:mesh_size)/Pawrad%rad(2:mesh_size)
                                                                  
 hh(1)          = zero
 hh(2:mesh_size)= ff1(2:mesh_size) * gg(2:mesh_size)
 deallocate(gg)

 call simp_gen(integral,hh,Pawrad)
 integral = four_pi * integral

 deallocate(hh)

end subroutine calc_slatradl
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/compmesh
!! NAME
!! compmesh
!!
!! FUNCTION
!! Compute all points of a radial mesh
!! Grid can be regular or logarithimc
!!
!! SIDE EFFECTS
!!  mesh <type(pawrad_type)>=data containing radial grid information
!!
!! PARENTS
!!      m_atom,m_paw_pwij,m_radmesh,psp7in
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! NOTES
!!  Possible mesh types (mesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! INPUTS
!! mesh=data containing radial grid information
!! r_for_intg=upper boundary for (future) integration over the radial grid
!!            (can be negative for an integration over the whole grid)
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS 
!!  The following quantities are calculated inside the routine:
!!    %stepint = Radial step used to convert any function from the
!!               present grid onto a regular grid in order to integrate it using trapeze method
!!    %rad(mesh_size) = Coordinates of all the points of the mesh.
!!    %radfact(mesh_size) = Factors used to compute radial integrals.
!!    %int_meshsz = Integrals will be computed up to r(int_meshsz)
!!    %simfact(mesh_size) = Factor used to compute radial integrals by the a Simpson scheme
!!                          Integral[f] = Sum_i [simfact(i)*f(i)]
!!    %rmax = Max. value of r = rad(mesh_size)
!!
!! SOURCE

subroutine compmesh(mesh,r_for_intg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: r_for_intg
 type(pawrad_type),intent(inout) :: mesh

!Local variables-------------------------------
!scalars
 integer :: ir,ir_last,isim
 real(dp) :: hh
 character(len=500) :: msg

!************************************************************************

 if (mesh%mesh_type==1) then
   isim=3
   mesh%stepint=mesh%rstep
   mesh%rad(1)=zero;mesh%radfact(1)=one
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =mesh%rstep*dble(ir-1)
     mesh%radfact(ir)=one
   end do
 else if (mesh%mesh_type==2) then
   isim=3
   mesh%stepint=mesh%lstep
   mesh%rad(1)=zero;mesh%radfact(1)=mesh%rstep
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =mesh%rstep*(exp(mesh%lstep*dble(ir-1))-one)
     mesh%radfact(ir)=mesh%rad(ir)+mesh%rstep
   end do
 else if (mesh%mesh_type==3) then
   isim=4
   mesh%stepint=mesh%lstep
   mesh%rad(1)=zero;mesh%radfact(1)=zero
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =mesh%rstep*exp(mesh%lstep*dble(ir-2))
     mesh%radfact(ir)=mesh%rad(ir)
   end do
 else if (mesh%mesh_type==4) then
   isim=3
   mesh%lstep=one/dble(mesh%mesh_size)
   mesh%stepint=mesh%lstep
   mesh%rad(1)=zero;mesh%radfact(1)=mesh%rstep
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =-mesh%rstep*log(one-mesh%lstep*dble(ir-1))
     mesh%radfact(ir)=mesh%rstep/(one-mesh%lstep*dble(ir-1))
   end do
 else !  Other values of mesh_type are not allowed (see psp7in.F90)
  write(msg,'(a,i0)')" Unknow value of mesh_type: ",mesh%mesh_type
  MSG_ERROR(msg)
 end if

 mesh%int_meshsz=mesh%mesh_size
 if (r_for_intg>0.d0) then
   ir=min(ifromr(mesh,r_for_intg),mesh%mesh_size)
   if (ir<mesh%mesh_size) then
     if (abs(mesh%rad(ir+1)-r_for_intg)<abs(mesh%rad(ir)-r_for_intg)) ir=ir+1
   end if
   if (ir>1) then
     if (abs(mesh%rad(ir-1)-r_for_intg)<abs(mesh%rad(ir)-r_for_intg)) ir=ir-1
   end if
   mesh%int_meshsz=ir
 end if

 hh=mesh%stepint/3.d0
 mesh%simfact(mesh%int_meshsz)=hh*mesh%radfact(mesh%int_meshsz)
 mesh%simfact(1:isim-2)=zero
 ir_last=1
 do ir=mesh%int_meshsz,isim,-2
   mesh%simfact(ir-1)=4.d0*hh*mesh%radfact(ir-1)
   mesh%simfact(ir-2)=2.d0*hh*mesh%radfact(ir-2)
   ir_last=ir-2
 end do
 mesh%simfact(ir_last)=half*mesh%simfact(ir_last)
 if (mesh%int_meshsz<mesh%mesh_size) mesh%simfact(mesh%int_meshsz+1:mesh%mesh_size)=zero

 mesh%rmax=mesh%rad(mesh%mesh_size)

end subroutine compmesh
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/copymesh
!! NAME
!! copymesh
!!
!! FUNCTION
!! Copy one radial mesh (in a generalized format) to another
!!
!! INPUTS
!!  mesh1 <type(pawrad_type)>=data containing radial grid information of input mesh
!!
!! OUTPUT
!!  mesh2 <type(pawrad_type)>=data containing radial grid information of output mesh
!!
!! PARENTS
!!      m_atom,m_paw_pwij,psp7in,psp7nl
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! NOTES
!!  Possible mesh types (mesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

subroutine copymesh(mesh1,mesh2)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 type(pawrad_type),intent(in) :: mesh1
 type(pawrad_type),intent(out) :: mesh2

!Local variables-------------------------------
!scalars
 integer :: ir

! *************************************************************************
 mesh2%mesh_type =mesh1%mesh_type
 mesh2%mesh_size =mesh1%mesh_size
 mesh2%int_meshsz=mesh1%int_meshsz
 mesh2%lstep     =mesh1%lstep
 mesh2%rstep     =mesh1%rstep
 mesh2%stepint   =mesh1%stepint
 mesh2%rmax      =mesh1%rmax
!If you need the following lines, please put ierr as a local variable (integer)
!those lines are not useful and the behavior is undefined. => problems with g95 (PMA)
!if (associated(mesh2%rad)) deallocate(mesh2%rad,STAT=ierr)
!if (associated(mesh2%radfact)) deallocate(mesh2%radfact,STAT=ierr)
!if (associated(mesh2%simfact)) deallocate(mesh2%simfact,STAT=ierr)
 allocate(mesh2%rad(mesh1%mesh_size))
 allocate(mesh2%radfact(mesh1%mesh_size))
 allocate(mesh2%simfact(mesh1%mesh_size))
 do ir=1,mesh1%mesh_size
   mesh2%rad(ir)    =mesh1%rad(ir)
   mesh2%radfact(ir)=mesh1%radfact(ir)
   mesh2%simfact(ir)=mesh1%simfact(ir)
 end do
end subroutine copymesh
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/simp_gen
!! NAME
!! simp_gen
!!
!! FUNCTION
!! Do integral on a given (generalized) grid using Simpson rule
!!
!! INPUTS
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  func(:)=integrand values
!!
!! OUTPUT
!!  intg=resulting integral by Simpson rule
!!
!! PARENTS
!!      Lij,m_atom,m_paw_commutator,m_paw_pwij,m_paw_slater,m_radmesh
!!      make_cs_dia,make_efg_paw,mlwfovlp_projpaw,optics_paw,optics_paw_core
!!      partial_dos_fractions_paw,pawdensities,pawdij,pawdij0,pawdijso
!!      pawfrnhat,pawinit,pawkij,pawnabla_init,pawpuxinit,pawshpfun,pawxc
!!      pawxc3,pawxcm,pawxcmpositron,pawxcpositron,poslifetime,psp7cg,psp7in
!!      psp7lo,psp7nl,qijb,qijb_bk,qijb_kk,smatrix_pawinit
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

subroutine simp_gen(intg,func,radmesh)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: intg
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: func(radmesh%int_meshsz)

!Local variables-------------------------------
!scalars
 integer :: ii,nn
 real(dp) :: resid,simp

! *************************************************************************

 nn=radmesh%int_meshsz

 simp=zero
 do ii=1,nn
   simp=simp+func(ii)*radmesh%simfact(ii)
 end do

 resid=zero
 if (radmesh%mesh_type==3) then
   resid=half*(func(2)+func(1))*(radmesh%rad(2)-radmesh%rad(1))
   if (mod(nn,2)==1) resid=resid+radmesh%stepint/3.d0*(1.25d0*func(2)*radmesh%radfact(2) &
&   +2.d0*func(3)*radmesh%radfact(3)-0.25d0*func(4)*radmesh%radfact(4))
 else if (mod(nn,2)==0) then
   resid=radmesh%stepint/3.d0*(1.25d0*func(1)*radmesh%radfact(1)+2.d0*func(2)*radmesh%radfact(2) &
&   -0.25d0*func(3)*radmesh%radfact(3))
 end if

 intg=simp+resid

end subroutine simp_gen
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/deducer0
!! NAME
!! deducer0
!!
!! FUNCTION
!! Extrapolate r=0 value of a function from values near r=0
!! using a 3 points formula
!!
!! INPUTS
!!  funcsz=size of array func
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! SIDE EFFECTS
!!  func(funcsz)=array containing values of function to extrapolate
!!
!! PARENTS
!!      Lij,denfgr,m_paw_slater,m_radmesh,make_cs_dia,make_efg_paw,optics_paw
!!      optics_paw_core,pawdenpot,pawdensities,pawdij0,pawdijso,pawfrnhat
!!      pawkij,pawnabla_init,pawxc,pawxcsph,pawxcsphpositron,psp7in
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine deducer0(func,funcsz,radmesh)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: funcsz
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(inout) :: func(funcsz)

!Local variables-------------------------------

! *************************************************************************

 if (radmesh%mesh_type==1.or.radmesh%mesh_type==2.or.radmesh%mesh_type==4) then
   func(1)=func(4)+3*(func(2)-func(3))
 else if (radmesh%mesh_type==3) then
   func(1)=func(4)+exp(two*radmesh%lstep)/(exp(radmesh%lstep)-one)*(func(2)-func(3))
 end if

end subroutine deducer0
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/nderiv_gen
!! NAME
!! nderiv_gen
!!
!! FUNCTION
!! Do corrected first (and -if requested- second) derivation on a given (generalized) grid.
!! This routine interfaces nderiv (derivation on a regular grid).
!!
!! INPUTS
!!  func(:)=input function
!!  nder= order of the derivation (1 or 2)
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  der(:,nder)=resulting derived function
!!
!! PARENTS
!!      optics_paw,optics_paw_core,pawdijso,pawinit,pawnabla_init,pawxc
!!      pawxcsph,pawxcsphpositron,poslifetime,psp7cc,spline_paw_fncs
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

subroutine nderiv_gen(der,func,nder,radmesh)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nder
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: func(radmesh%mesh_size)
 real(dp),intent(out) :: der(radmesh%mesh_size,nder)

!Local variables-------------------------------
!scalars
 integer :: msz
 character(len=500) :: msg

! *************************************************************************

 if(nder/=1 .and. nder/=2)then
   write(msg,'(a,i0)')' first or second derivatives are allowed while the argument nder=',nder
   MSG_BUG(msg)
 end if

 msz=radmesh%mesh_size

 if (radmesh%mesh_type==1) then

   call nderiv(radmesh%rstep,func,der(1:msz,1),radmesh%mesh_size,1)
   if (nder==2) call nderiv(radmesh%rstep,func,der(1:msz,2),radmesh%mesh_size,2)

 else if (radmesh%mesh_type==2) then

   call nderiv(radmesh%lstep,func,der(1:msz,1),radmesh%mesh_size,1)
   der(1:msz,1)=der(1:msz,1)/radmesh%radfact(1:msz)
   if (nder==2)then
     call nderiv(radmesh%lstep,func,der(1:msz,2),radmesh%mesh_size,2)
     der(1:msz,2)=(der(1:msz,2)/radmesh%radfact(1:msz)-der(1:msz,1))/radmesh%radfact(1:msz)
   end if

 else if (radmesh%mesh_type==3) then

   call nderiv(radmesh%lstep,func(2:msz),der(2:msz,1),msz-1,1)
   der(2:msz,1)=der(2:msz,1)/radmesh%radfact(2:msz)
   call deducer0(der(:,1),msz,radmesh)
   if (nder==2)then
     call nderiv(radmesh%lstep,func(2:msz),der(2:msz,2),msz-1,2)
     der(2:msz,2)=(der(2:msz,2)/radmesh%radfact(2:msz)-der(2:msz,1))/radmesh%radfact(2:msz)
     call deducer0(der(:,2),msz,radmesh)
   end if

 else if (radmesh%mesh_type==4) then

   call nderiv(radmesh%lstep,func,der(1:msz,1),radmesh%mesh_size,1)
   der(1:msz,1)=der(1:msz,1)/radmesh%radfact(1:msz)
   if (nder==2)then
     call nderiv(radmesh%lstep,func,der(1:msz,2),radmesh%mesh_size,2)
     der(1:msz,2)=der(1:msz,2)/radmesh%radfact(1:msz)**2-der(1:msz,1)/radmesh%rstep
   end if

 end if

end subroutine nderiv_gen
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/poisson
!! NAME
!! poisson
!!
!! FUNCTION
!!  Solve poisson equation for angularly dependent charge
!!  distribution of angular momentum l
!!  Densities and potentials are given on a (generalized) radial grid
!!
!! INPUTS
!!  den(radmesh%mesh_size)= electron density * (4*pi*r**2) appropriate for l
!!  ll= l quantum number
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  qq= lth moment of the charge
!!  rv(radmesh%mesh_size)= electrostatic potential * r in (Hartree*Bohr) units
!!          where v(r)=\frac{1}{2l+1}(\frac{int[(r''^(l+2))g(r'')dr'']} {r^(l+1)}
!!                                   +(r^l) int[r''^(1-l)g(r'')dr''])
!!
!! PARENTS
!!      m_radmesh,pawdenpot,pawdij0,pawinit,pawkij,pawpuxinit,psp7in
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

subroutine poisson(den,ll,qq,radmesh,rv)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll
 real(dp),intent(out) :: qq
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: den(radmesh%mesh_size)
 real(dp),intent(out) :: rv(radmesh%mesh_size)

!Local variables ---------------------------------------
!scalars
 integer :: ir,jr,mm,nn
 real(dp) :: angm,hh,rr
!arrays
 real(dp) :: aa(radmesh%mesh_size),bb(radmesh%mesh_size),cc(radmesh%mesh_size)
 real(dp) :: dd(radmesh%mesh_size+1),ee(4)
 real(dp),allocatable :: radl(:),radl1(:)

! ***************************************************************************

!==============================================
!UNIFORM GRID - NUMEROV ALGORITHM
!==============================================
 if (radmesh%mesh_type==1) then
   hh=radmesh%rstep
   nn=radmesh%int_meshsz-1
   do ir=1,nn
     aa(ir)=two*hh*den(ir+1)/(ir)
     bb(ir)=den(ir+1)*((ir*hh)**ll)
   end do
   dd(1)=zero
   dd(2:nn+1)=bb(1:nn)
   call simp_gen(qq,dd,radmesh)
   qq=qq/(2*ll+1)
   rv(1)=aa(1)+0.1_dp*aa(2)
   do ir=2,nn-1
     rv(ir)=aa(ir)+0.1_dp*(aa(ir+1)+aa(ir-1))
   end do
   rv(nn)=aa(nn)+0.1_dp*aa(nn-1)
   angm=dble(ll*(ll+1))
   rr=(nn+1)*hh
   rv(nn)=rv(nn)+(2.4_dp-0.2_dp*angm/((nn+1)**2))*qq/(rr**ll)
   do ir=1,nn
     aa(ir)=angm/(ir*ir)
     bb(ir)=2.4_dp+aa(ir)
   end do
   do ir=1,nn-1
     cc(ir)=-1.2_dp+0.1_dp*aa(ir+1)
   end do
   do ir=nn,2,-1
     aa(ir)=-1.2_dp+0.1_dp*aa(ir-1)
   end do
   if (nn.eq.1) then
     rv(2)=rv(1)/bb(1)
     rv(1)=zero
   else
     do ir=2,nn
       bb(ir)=bb(ir)-aa(ir)*cc(ir-1)/bb(ir-1)
     end do
     rv(1)=rv(1)/bb(1)
     do ir=2,nn
       rv(ir)=(rv(ir)-aa(ir)*rv(ir-1))/bb(ir)
     end do
     do ir=nn-1,1,-1
       rv(ir)=rv(ir)-cc(ir)*rv(ir+1)/bb(ir)
     end do
     do ir=nn+1,2,-1
       rv(ir)=rv(ir-1)
     end do
     rv(1)=zero
     if (nn+1<radmesh%mesh_size) rv(nn+2:radmesh%mesh_size)=zero
   end if
   rv(:)=half*rv(:)

!  ==============================================
!  ANY OTHER GRID - SIMPSON ALGORITHM
!  ==============================================
 else
   nn=radmesh%mesh_size;hh=third*radmesh%stepint
   do while (abs(den(nn))<tol16.and.nn>radmesh%int_meshsz)
     nn=nn-1
   end do
   mm=nn;if (radmesh%mesh_type==3) mm=mm-1
   allocate(radl(nn),radl1(nn))
   do jr=nn,2,-1
     ir=nn-jr+1
     radl(jr) =radmesh%rad(jr)**ll
     radl1(jr)=radmesh%rad(jr)*radl(jr)
     aa(ir)=den(jr)*radmesh%radfact(jr)*radl(jr)
     bb(ir)=den(jr)*radmesh%radfact(jr)/radl1(jr)
   end do
   radl(1)=zero;radl1(1)=zero
   ee(2)=aa(nn-1);ee(3)=aa(nn-2);ee(4)=aa(nn-3)
   call deducer0(ee,4,radmesh)
   aa(nn)=ee(1)
   ee(2)=bb(nn-1);ee(3)=bb(nn-2);ee(4)=bb(nn-3)
   call deducer0(ee,4,radmesh)
   bb(nn)=ee(1)
   cc(1)=zero;dd(1)=zero
   do ir=3,mm,2
     cc(ir)  =cc(ir-2)+hh*(aa(ir-2)+four*aa(ir-1)+aa(ir))
     cc(ir-1)=cc(ir-2)+hh*(1.25_dp*aa(ir-2)+two*aa(ir-1)-quarter*aa(ir))
     dd(ir)  =dd(ir-2)+hh*(bb(ir-2)+four*bb(ir-1)+bb(ir))
     dd(ir-1)=dd(ir-2)+hh*(1.25_dp*bb(ir-2)+two*bb(ir-1)-quarter*bb(ir))
   end do
   if (mod(mm,2)==0) then
     cc(mm)=cc(mm-2)+hh*(aa(mm-2)+four*aa(mm-1)+aa(mm))
     dd(mm)=dd(mm-2)+hh*(bb(mm-2)+four*bb(mm-1)+bb(mm))
   end if
   if (mm<nn) then
     cc(nn)=cc(mm)+half*(aa(mm)+aa(nn))*(radmesh%rad(1+nn-mm)-radmesh%rad(1))
     dd(nn)=dd(mm)+half*(bb(mm)+bb(nn))*(radmesh%rad(1+nn-mm)-radmesh%rad(1))
   end if
   rv(1)=zero
   do ir=2,nn
     jr=nn-ir+1
     rv(ir)=(dd(jr)*radl1(ir)+(cc(nn)-cc(jr))/radl(ir))/(two*ll+one)
   end do
   if (nn<radmesh%mesh_size) rv(nn+1:radmesh%mesh_size)=rv(nn)
   qq=cc(nn)/(two*ll+one)
   deallocate(radl,radl1)
 end if

end subroutine poisson
!!***

!----------------------------------------------------------------------

!!****f* m_radmesh/ifromr
!! NAME
!! ifromr
!!
!! FUNCTION
!! Retreive Index FROM a given R value in a radial grid
!! Grid can be regular or logarithimc
!!
!! INPUTS
!!  rr=given input r value
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  ifromr=index of rr in radial grid
!!
!! PARENTS
!!      compmesh,psp7in
!!
!! CHILDREN
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

function ifromr(radmesh,rr)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: ifromr
 real(dp),intent(in) :: rr
 type(pawrad_type),intent(in) :: radmesh

!Local variables-------------------------------

! *************************************************************************

 if (radmesh%mesh_type==1) then
   ifromr=int(tol8+rr/radmesh%rstep)+1
 else if (radmesh%mesh_type==2) then
   ifromr=int(tol8+log(1.d0+rr/radmesh%rstep)/radmesh%lstep)+1
 else if (radmesh%mesh_type==3) then
   if (rr<radmesh%rstep) then
     ifromr=1
   else
     ifromr=int(tol8+log(rr/radmesh%rstep)/radmesh%lstep)+2
   end if
 else if (radmesh%mesh_type==4) then
   ifromr=int(tol8+(1.d0-exp(-rr/radmesh%rstep))/radmesh%lstep)+1
 else 
!  Other values of mesh_type are not allowed (see psp7in.F90)
   stop
 end if

end function ifromr
!!***
  
!----------------------------------------------------------------------

END MODULE m_radmesh
!!***
