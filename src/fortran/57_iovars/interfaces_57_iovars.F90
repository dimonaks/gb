!!****m* ABINIT/interfaces_57_iovars
!! NAME
!! interfaces_57_iovars
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/57_iovars
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_57_iovars

 implicit none

interface
 subroutine append_xyz(dtset_char,lenstr,string,xyz_fname,strln)
  use defs_basis
  implicit none
  integer,intent(inout) :: lenstr
  integer,intent(in) :: strln
  character(len=2),intent(in) :: dtset_char
  character(len=strln),intent(inout) :: string
  character(len=fnlen),intent(in) :: xyz_fname
 end subroutine append_xyz
end interface

interface
 subroutine chkdpr(advice_change_cond,cond_number,cond_string,cond_values,&  
  &  ierr,input_name,input_value,minimal_flag,minimal_value,unit)
  use defs_basis
  implicit none
  integer,intent(in) :: advice_change_cond
  integer,intent(in) :: cond_number
  integer,intent(inout) :: ierr
  integer,intent(in) :: minimal_flag
  integer,intent(in) :: unit
  character(len=*),intent(in) :: input_name
  real(dp),intent(in) :: input_value
  real(dp),intent(in) :: minimal_value
  character(len=*),intent(in) :: cond_string(3)
  integer,intent(in) :: cond_values(3)
 end subroutine chkdpr
end interface

interface
 subroutine chkinp(dtsets,iout,mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  type(mpi_type),intent(in) :: mpi_enreg
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine chkinp
end interface

interface
 subroutine chkint(advice_change_cond,cond_number,cond_string,cond_values,&  
  &  ierr,input_name,input_value,&  
  &  list_number,list_values,minmax_flag,minmax_value,unit)
  implicit none
  integer,intent(in) :: advice_change_cond
  integer,intent(in) :: cond_number
  integer,intent(inout) :: ierr
  integer,intent(in) :: input_value
  integer,intent(in) :: list_number
  integer,intent(in) :: minmax_flag
  integer,intent(in) :: minmax_value
  integer,intent(in) :: unit
  character(len=*),intent(in) :: input_name
  character(len=*),intent(in) :: cond_string(3)
  integer,intent(in) :: cond_values(3)
  integer,intent(in) :: list_values(list_number)
 end subroutine chkint
end interface

interface
 subroutine chkint_eq(advice_change_cond,cond_number,cond_string,cond_values,&  
  &  ierr,input_name,input_value,list_number,list_values,unit)
  implicit none
  integer,intent(in) :: advice_change_cond
  integer,intent(in) :: cond_number
  integer,intent(inout) :: ierr
  integer,intent(in) :: input_value
  integer,intent(in) :: list_number
  integer,intent(in) :: unit
  character(len=*),intent(in) :: input_name
  character(len=*),intent(in) :: cond_string(3)
  integer,intent(in) :: cond_values(3)
  integer,intent(in) :: list_values(list_number)
 end subroutine chkint_eq
end interface

interface
 subroutine chkint_ge(advice_change_cond,cond_number,cond_string,cond_values,&  
  &  ierr,input_name,input_value,minmax_value,unit)
  implicit none
  integer,intent(in) :: advice_change_cond
  integer,intent(in) :: cond_number
  integer,intent(inout) :: ierr
  integer,intent(in) :: input_value
  integer,intent(in) :: minmax_value
  integer,intent(in) :: unit
  character(len=*),intent(in) :: input_name
  character(len=*),intent(in) :: cond_string(3)
  integer,intent(in) :: cond_values(3)
 end subroutine chkint_ge
end interface

interface
 subroutine chkint_le(advice_change_cond,cond_number,cond_string,cond_values,&  
  &  ierr,input_name,input_value,&  
  &  minmax_value,unit)
  implicit none
  integer,intent(in) :: advice_change_cond
  integer,intent(in) :: cond_number
  integer,intent(inout) :: ierr
  integer,intent(in) :: input_value
  integer,intent(in) :: minmax_value
  integer,intent(in) :: unit
  character(len=*),intent(in) :: input_name
  character(len=*),intent(in) :: cond_string(3)
  integer,intent(in) :: cond_values(3)
 end subroutine chkint_le
end interface

interface
 subroutine chkint_ne(advice_change_cond,cond_number,cond_string,cond_values,&  
  &  ierr,input_name,input_value,&  
  &  list_number,list_values,unit)
  implicit none
  integer,intent(in) :: advice_change_cond
  integer,intent(in) :: cond_number
  integer,intent(inout) :: ierr
  integer,intent(in) :: input_value
  integer,intent(in) :: list_number
  integer,intent(in) :: unit
  character(len=*),intent(in) :: input_name
  character(len=*),intent(in) :: cond_string(3)
  integer,intent(in) :: cond_values(3)
  integer,intent(in) :: list_values(list_number)
 end subroutine chkint_ne
end interface

interface
 subroutine chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&  
  &  ierr,input_name,input_value,&  
  &  list_number,list_values,minmax_flag,minmax_value,unit)
  implicit none
  integer,intent(in) :: advice_change_cond
  integer,intent(in) :: cond_number
  integer,intent(inout) :: ierr
  integer,intent(in) :: input_value
  integer,intent(in) :: list_number
  integer,intent(in) :: minmax_flag
  integer,intent(in) :: minmax_value
  integer,intent(in) :: unit
  character(len=*),intent(in) :: input_name
  character(len=*),intent(in) :: cond_string(3)
  integer,intent(in) :: cond_values(3)
  integer,intent(in) :: list_values(list_number)
 end subroutine chkint_prt
end interface

interface
 subroutine chkneu(charge,dtset,occopt)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: occopt
  real(dp),intent(in) :: charge
  type(dataset_type),intent(inout) :: dtset
 end subroutine chkneu
end interface

interface
 subroutine chkvars (string)
  implicit none
  character(len=*),intent(in) :: string
 end subroutine chkvars
end interface

interface
 subroutine dsksta(ishm,usepaw,nbandkss,mpsang,natom,ntypat,npwkss,nkpt,nspinor,nsppol,nsym2,dimlmn)
  implicit none
  integer,intent(in) :: ishm
  integer,intent(in) :: mpsang
  integer,intent(in) :: natom
  integer,intent(in) :: nbandkss
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwkss
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym2
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  integer,intent(in) :: dimlmn(natom*usepaw)
 end subroutine dsksta
end interface

interface
 subroutine finddistrproc(dtset,mband,mpi_enreg)
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine finddistrproc
end interface

interface
 subroutine getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,mixalch,npsp,npspalch,&  
  &  ntypat,ntypalch,pspheads)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: lmnmax
  integer,intent(out) :: lmnmaxso
  integer,intent(out) :: lnmax
  integer,intent(out) :: lnmaxso
  integer,intent(in) :: npsp
  integer,intent(in) :: npspalch
  integer,intent(in) :: ntypalch
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: mixalch(npspalch,ntypalch)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine getdim_nloc
end interface

interface
 subroutine importxyz (lenstr,string_raw,string_upper,strln)
  implicit none
  integer,intent(inout) :: lenstr
  integer,intent(in) :: strln
  character(len=*),intent(in) :: string_raw
  character(len=*),intent(inout) :: string_upper
 end subroutine importxyz
end interface

interface
 subroutine indefo(dtsets,ndtset_alloc)
  use defs_abitypes
  implicit none
  integer,intent(in) :: ndtset_alloc
  type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)
 end subroutine indefo
end interface

interface
 subroutine indefo1(dtset)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(out) :: dtset
 end subroutine indefo1
end interface

interface
 subroutine ingeo (acell,berryopt,bravais,efield_xcart,&  
  &  genafm,iatfix,icoulomb,iimage,iout,jdtset,jellslab,lenstr,&  
  &  msym,natom,nimage,nspden,nsppol,nsym,ntypat,pawspnorb,&  
  &  ptgroupma,rprim,slabzbeg,slabzend,spgroup,spinat,string,symafm,&  
  &  symmorphi,symrel,tnons,tolsym,typat,vel,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(out) :: icoulomb
  integer,intent(in) :: iimage
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(out) :: jellslab
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(inout) :: natom
  integer,intent(in) :: nimage
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(out) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: pawspnorb
  integer,intent(out) :: ptgroupma
  integer,intent(out) :: spgroup
  integer,intent(inout) :: symmorphi
  real(dp),intent(out) :: slabzbeg
  real(dp),intent(out) :: slabzend
  character(len=*),intent(in) :: string
  real(dp),intent(out) :: tolsym
  integer,intent(out) :: bravais(11)
  real(dp),intent(out) :: acell(3)
  real(dp),intent(inout) :: efield_xcart(3)
  real(dp),intent(out) :: genafm(3)
  integer,intent(out) :: iatfix(3,natom)
  real(dp),intent(out) :: rprim(3,3)
  real(dp),intent(inout) :: spinat(3,natom)
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
  integer,intent(out) :: typat(natom)
  real(dp),intent(out) :: vel(3,natom)
  real(dp),intent(out) :: xred(3,natom)
 end subroutine ingeo
end interface

interface
 subroutine ingeobld (iout,jdtset,lenstr,natrd,natom,&  
  &  nobj,string,typat,typat_read,xcart,xcart_read)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(in) :: natom
  integer,intent(in) :: natrd
  integer,intent(in) :: nobj
  character(len=*),intent(in) :: string
  integer,intent(out) :: typat(natom)
  integer,intent(in) :: typat_read(natrd)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(in) :: xcart_read(3,natrd)
 end subroutine ingeobld
end interface

interface
 subroutine inkpts(bravais,chksymbreak,iout,iscf,istwfk,jdtset,&  
  &  kpt,kptopt,kptnrm,kptrlatt,&  
  &  kptrlen,lenstr,msym,nkpt,nqpt,nshiftk,nsym,&  
  &  occopt,qpt,qptnrm,response,rprimd,&  
  &  shiftk,string,symafm,symrel,vacuum,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: chksymbreak
  integer,intent(in) :: iout
  integer,intent(in) :: iscf
  integer,intent(in) :: jdtset
  integer,intent(in) :: kptopt
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(inout) :: nkpt
  integer,intent(in) :: nqpt
  integer,intent(out) :: nshiftk
  integer,intent(in) :: nsym
  integer,intent(in) :: occopt
  integer,intent(in) :: response
  real(dp),intent(out) :: kptnrm
  real(dp),intent(out) :: kptrlen
  real(dp),intent(out) :: qptnrm
  character(len=*),intent(in) :: string
  integer,intent(in) :: bravais(11)
  integer,intent(out) :: kptrlatt(3,3)
  integer,intent(in) :: vacuum(3)
  integer,intent(out) :: istwfk(nkpt)
  real(dp),intent(out) :: kpt(3,nkpt)
  real(dp),intent(out) :: qpt(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: shiftk(3,8)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(out) :: wtk(nkpt)
 end subroutine inkpts
end interface

interface
 subroutine invacuum(jdtset,lenstr,natom,rprimd,string,vacuum,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(in) :: natom
  character(len=*),intent(in) :: string
  integer,intent(out) :: vacuum(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine invacuum
end interface

interface
 subroutine invars0(dtsets,istatr,istatshft,lenstr,&  
  &  msym,mxnatom,mxnimage,mxntypat,ndtset,ndtset_alloc,npsp,papiopt,timopt,string)
  use defs_abitypes
  implicit none
  integer,intent(out) :: istatr
  integer,intent(out) :: istatshft
  integer,intent(in) :: lenstr
  integer,intent(out) :: msym
  integer,intent(out) :: mxnatom
  integer,intent(out) :: mxnimage
  integer,intent(out) :: mxntypat
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(out) :: npsp
  integer,intent(out) :: papiopt
  integer,intent(inout) :: timopt
  character(len=*),intent(in) :: string
  type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)
 end subroutine invars0
end interface

interface
 subroutine invars1(bravais,dtset,iout,jdtset,lenstr,mband_upper,mpi_enreg,msym,&  
  &  string,symafm,symrel,tnons,zion_max)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(out) :: mband_upper
  integer,intent(in) :: msym
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  character(len=*),intent(inout) :: string
  real(dp),intent(in) :: zion_max
  integer,intent(inout) :: bravais(11)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine invars1
end interface

interface
 subroutine invars1m(bravais_,dmatpuflag,dtsets,iout,lenstr,mband_upper_,mpi_enreg,&  
  &  msym,mxgw_nqlwl,mxlpawu,mxmband_upper,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxnconeq,&  
  &  mxnkpt,mxnkptgw,mxnimage,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,&  
  &  ndtset,ndtset_alloc,string,zion_max)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: dmatpuflag
  integer,intent(in) :: iout
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(out) :: mxgw_nqlwl
  integer,intent(out) :: mxlpawu
  integer,intent(out) :: mxmband_upper
  integer,intent(in) :: mxnatom
  integer,intent(out) :: mxnatpawu
  integer,intent(out) :: mxnatsph
  integer,intent(out) :: mxnatvshift
  integer,intent(out) :: mxnconeq
  integer,intent(in) :: mxnimage
  integer,intent(out) :: mxnkpt
  integer,intent(out) :: mxnkptgw
  integer,intent(out) :: mxnnos
  integer,intent(out) :: mxnqptdm
  integer,intent(out) :: mxnspinor
  integer,intent(out) :: mxnsppol
  integer,intent(out) :: mxnsym
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  type(mpi_type),intent(inout) :: mpi_enreg
  character(len=*),intent(inout) :: string
  real(dp),intent(in) :: zion_max
  integer,intent(out) :: bravais_(11,0:ndtset_alloc)
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  integer,intent(out) :: mband_upper_(0:ndtset_alloc)
 end subroutine invars1m
end interface

interface
 subroutine invars2(bravais,dtset,iout,jdtset,lenstr,&  
  &  mband,msym,npsp,pspheads,string,usepaw,zionpsp)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(in) :: mband
  integer,intent(in) :: msym
  integer,intent(in) :: npsp
  integer,intent(in) :: usepaw
  type(dataset_type),intent(inout) :: dtset
  character(len=*),intent(in) :: string
  integer,intent(in) :: bravais(11)
  type(pspheader_type),intent(in) :: pspheads(npsp)
  real(dp),intent(in) :: zionpsp(npsp)
 end subroutine invars2
end interface

interface
 subroutine invars2m(bravais_,dtsets,iout,lenstr,&  
  &  mband_upper_,mpi_enreg,msym,ndtset,ndtset_alloc,npsp,pspheads,string)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  type(mpi_type),intent(inout) :: mpi_enreg
  character(len=*),intent(in) :: string
  integer,intent(in) :: bravais_(11,0:ndtset_alloc)
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: mband_upper_(0:ndtset_alloc)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine invars2m
end interface

interface
 subroutine iofn2(npsp,filnam)
  use defs_basis
  implicit none
  integer,intent(in) :: npsp
  character(len=fnlen), intent(out) :: filnam(npsp)
 end subroutine iofn2
end interface

interface
 function is_input_variable(invar)
  implicit none
  character(len=*),intent(in) :: invar
  logical :: is_input_variable
 end function is_input_variable
end interface

interface
 subroutine macroin(dtsets,ndtset_alloc)
  use defs_abitypes
  implicit none
  integer,intent(in) :: ndtset_alloc
  type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)
 end subroutine macroin
end interface

interface
 subroutine macroin2(dtsets,ndtset_alloc)
  use defs_abitypes
  implicit none
  integer,intent(in) :: ndtset_alloc
  type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)
 end subroutine macroin2
end interface

interface
 subroutine memana(cadd,cfft,cfftf,chain,cmpw,dttyp,iout,iprcel,iscf,&  
  &  marrays,mbcg,mbdiskpd,mbdiskwf,mbf_fftgr,mbgylm,mffmem,&  
  &  mkmem,mpi_enreg,mpw,natom,nchain,nfft,nfftf,occopt,option,prtvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: iprcel
  integer,intent(in) :: iscf
  integer,intent(in) :: marrays
  integer,intent(in) :: mffmem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nchain
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: mbcg
  real(dp),intent(in) :: mbdiskpd
  real(dp),intent(in) :: mbdiskwf
  real(dp),intent(in) :: mbf_fftgr
  real(dp),intent(in) :: mbgylm
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: cadd(marrays)
  integer,intent(in) :: cfft(marrays)
  integer,intent(in) :: cfftf(marrays)
  logical,intent(in) :: chain(marrays,nchain)
  integer,intent(in) :: cmpw(marrays)
  integer,intent(in) :: dttyp(marrays)
 end subroutine memana
end interface

interface
 subroutine memkss(mband,mgfft,mkmem,mpi_enreg,mproj,mpsang,mpw,natom,&  
  &  ngfft,nkpt,nspinor,nsym,ntypat)
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mproj
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
 end subroutine memkss
end interface

interface
 subroutine memorf(cplex,n1xccc,getcell,idtset,intxc,iout,iprcel,&  
  &  iscf,jdtset,lmnmax,lnmax,mband,mffmem,mgfft,&  
  &  mkmems,mpi_enreg,mpsang,mpssoang,mpw,mqgrid,&  
  &  natom,nband,nfft,ngfft,&  
  &  nkpt,nloalg,nspden,nspinor,nsppol,nsym,ntypat,&  
  &  occopt,optddk,optphon,option,optstrs,prtvol,useylm,xclevel)
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: getcell
  integer,intent(in) :: idtset
  integer,intent(in) :: intxc
  integer,intent(in) :: iout
  integer,intent(in) :: iprcel
  integer,intent(in) :: iscf
  integer,intent(in) :: jdtset
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mband
  integer,intent(in) :: mffmem
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: mqgrid
  integer,intent(in) :: n1xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: optddk
  integer,intent(in) :: option
  integer,intent(in) :: optphon
  integer,intent(in) :: optstrs
  integer,intent(in) :: prtvol
  integer,intent(in) :: useylm
  integer,intent(in) :: xclevel
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: mkmems(3)
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine memorf
end interface

interface
 subroutine memory(n1xccc,getcell,idtset,icoulomb,intxc,ionmov,iout,iprcch,iprcel,&  
  &  iscf,jdtset,lmnmax,lnmax,&  
  &  mband,mffmem,mgfft,mgfftdiel,mgfftf,mkmem,mpi_enreg,mpsang,mpssoang,mpw,mqgrid_ff,mqgrid_vl,&  
  &  natom,nband,nfft,nfftdiel,nfftf,ngfft,ngfftdiel,ngfftf,&  
  &  nkpt,nloalg,npsp,npulayit,npwdiel,nspden,nspinor,nsppol,nsym,ntypat,&  
  &  occopt,optforces,option,optstress,pawcpxocc,pawmixdg,pawnhatxc,pawspnorb,pawstgylm,&  
  &  prtvol,pspheads,typat,ucvol,usepaw,useylm,xclevel)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: getcell
  integer,intent(in) :: icoulomb
  integer,intent(in) :: idtset
  integer,intent(in) :: intxc
  integer,intent(in) :: ionmov
  integer,intent(in) :: iout
  integer,intent(in) :: iprcch
  integer,intent(in) :: iprcel
  integer,intent(in) :: iscf
  integer,intent(in) :: jdtset
  integer,intent(in) :: lmnmax
  integer,intent(in) :: lnmax
  integer,intent(in) :: mband
  integer,intent(in) :: mffmem
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mgfftf
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: mpw
  integer,intent(in) :: mqgrid_ff
  integer,intent(in) :: mqgrid_vl
  integer,intent(in) :: n1xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkpt
  integer,intent(in) :: npsp
  integer,intent(in) :: npulayit
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: optforces
  integer,intent(in) :: option
  integer,intent(in) :: optstress
  integer,intent(in) :: pawcpxocc
  integer,intent(in) :: pawmixdg
  integer,intent(in) :: pawnhatxc
  integer,intent(in) :: pawspnorb
  integer,intent(in) :: pawstgylm
  integer,intent(in) :: prtvol
  integer,intent(in) :: usepaw
  integer,intent(in) :: useylm
  integer,intent(in) :: xclevel
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: nband(nkpt*nsppol)
  type(pspheader_type) :: pspheads(npsp)
  integer,intent(in) :: typat(natom)
 end subroutine memory
end interface

interface
 subroutine out_acknowl(dtsets,iout,ndtset_alloc,npsp,pspheads) 
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine out_acknowl
end interface

interface
 subroutine outqmc(cg,dtset,eigen,gprimd,hdr,kg,mband,mkmem,mpi_enreg,&  
  &  mpw,nkpt,npwarr,nspinor,nsppol,occ,psps,results_gs)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: mband
  integer :: mkmem
  integer :: mpw
  integer :: nkpt
  integer :: nspinor
  integer :: nsppol
  type(dataset_type) :: dtset
  type(hdr_type) :: hdr
  type(mpi_type) :: mpi_enreg
  type(pseudopotential_type) :: psps
  type(results_gs_type) :: results_gs
  real(dp) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp) :: eigen(mband*nkpt*nsppol)
  real(dp) :: gprimd(3,3)
  integer :: kg(3,mpw*mkmem)
  integer :: npwarr(nkpt)
  real(dp) :: occ(mband*nkpt*nsppol)
 end subroutine outqmc
end interface

interface
 function i2s(n)
  implicit none
  integer, intent(in) :: n
  character(len=20) :: i2s
 end function i2s
end interface

interface
 function r2s(r,real_format)
  use defs_basis
  implicit none
  real(dp),intent(in) :: r
  character(len=80) :: r2s
  character(len=*),intent(in) :: real_format
 end function r2s
end interface

interface
 subroutine outvar1 (choice,dtsets,iout,istatr,istatshft,&  
  &  jdtset_,mxgw_nqlwl,mxmband,mxnatom,mxnatsph,mxnatvshift,mxnimage,mxnkptgw,mxnkpt,mxnqptdm,mxnsppol,mxnsym,mxntypat,&  
  &  mxnatpawu,ndtset,ndtset_alloc,npsp,prtvol_glob,response,response_,results_out,usepaw)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: iout
  integer,intent(in) :: istatr
  integer,intent(in) :: istatshft
  integer,intent(in) :: mxgw_nqlwl
  integer,intent(in) :: mxmband
  integer,intent(in) :: mxnatom
  integer,intent(in) :: mxnatpawu
  integer,intent(in) :: mxnatsph
  integer,intent(in) :: mxnatvshift
  integer,intent(in) :: mxnimage
  integer,intent(in) :: mxnkpt
  integer,intent(in) :: mxnkptgw
  integer,intent(in) :: mxnqptdm
  integer,intent(in) :: mxnsppol
  integer,intent(in) :: mxnsym
  integer,intent(in) :: mxntypat
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  integer,intent(in) :: prtvol_glob
  integer,intent(in) :: response
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  integer,intent(in) :: response_(ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 end subroutine outvar1
end interface

interface
 subroutine outvars (choice,dmatpuflag,dtsets,iout,istatr,istatshft,mpi_enreg,&  
  &  mxgw_nqlwl,mxlpawu,mxmband,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxnconeq,mxnimage,mxnkptgw,mxnkpt,&  
  &  mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,mxntypat,&  
  &  ndtset,ndtset_alloc,npsp,results_out,timopt)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: dmatpuflag
  integer,intent(in) :: iout
  integer,intent(in) :: istatr
  integer,intent(in) :: istatshft
  integer,intent(in) :: mxgw_nqlwl
  integer,intent(in) :: mxlpawu
  integer,intent(in) :: mxmband
  integer,intent(in) :: mxnatom
  integer,intent(in) :: mxnatpawu
  integer,intent(in) :: mxnatsph
  integer,intent(in) :: mxnatvshift
  integer,intent(in) :: mxnconeq
  integer,intent(in) :: mxnimage
  integer,intent(in) :: mxnkpt
  integer,intent(in) :: mxnkptgw
  integer,intent(in) :: mxnnos
  integer,intent(in) :: mxnqptdm
  integer,intent(in) :: mxnspinor
  integer,intent(in) :: mxnsppol
  integer,intent(in) :: mxnsym
  integer,intent(in) :: mxntypat
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  integer,intent(in) :: timopt
  type(mpi_type),intent(in) :: mpi_enreg
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 end subroutine outvars
end interface

interface
 subroutine outxml_open(filename)
  implicit none
  character(len = *), intent(in) :: filename
 end subroutine outxml_open
end interface

interface
 subroutine outxml_finalise(tsec, values)
  use defs_basis
  implicit none
  integer, intent(in) :: values(8)
  real(dp), intent(in) :: tsec(2)
 end subroutine outxml_finalise
end interface

interface
 subroutine out_resultsgs_XML(dtset, level, results_gs, usepaw)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: level
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  type(results_gs_type),intent(inout) :: results_gs
 end subroutine out_resultsgs_XML
end interface

interface
 subroutine out_geometry_XML(dtset, level, natom, rprimd, xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: level
  integer,intent(in) :: natom
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine out_geometry_XML
end interface

interface
 subroutine parsefile(filnamin, lenstr, ndtset, string)
  use defs_basis
  implicit none
  integer, intent(out) :: lenstr
  integer, intent(out) :: ndtset
  character(len = *), intent(in) :: filnamin
  character(len = strlen), intent(out) :: string
 end subroutine parsefile
end interface

interface
 subroutine prtocc(dtsets,iout,jdtset_,&  
  &  ndtset_alloc,prtvol_glob,results_out)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: prtvol_glob
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 end subroutine prtocc
end interface

interface
 subroutine prttagm(dprarr,intarr,iout,jdtset_,length,marr,narr,ndtset_alloc,token,typevarphys)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: length
  integer,intent(in) :: marr
  integer,intent(in) :: narr
  integer,intent(in) :: ndtset_alloc
  character(len=*),intent(in) :: token
  character(len=3),intent(in) :: typevarphys
  real(dp),intent(in) :: dprarr(marr,0:ndtset_alloc)
  integer,intent(in) :: intarr(marr,0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
 end subroutine prttagm
end interface

interface
 subroutine wvl_memory(dtset, idtset, mpi_enreg, npsp, option, pspheads)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: idtset
  integer,intent(in) :: npsp
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine wvl_memory
end interface

interface
 subroutine wvl_setBoxGeometry(dtset, me, radii, rprimd, xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: me
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: radii(dtset%ntypat,2)
  real(dp),intent(inout) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine wvl_setBoxGeometry
end interface

end module interfaces_57_iovars
!!***
