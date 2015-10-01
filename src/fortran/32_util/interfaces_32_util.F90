!!****m* ABINIT/interfaces_32_util
!! NAME
!! interfaces_32_util
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/32_util
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

module interfaces_32_util

 implicit none

interface
 subroutine acrossb(a,b,c)
  use defs_basis
  implicit none
  real(dp),intent(in) :: a(3)
  real(dp),intent(in) :: b(3)
  real(dp),intent(out) :: c(3)
 end subroutine acrossb
end interface

interface
 subroutine appdig(integ,string,strinn)
  implicit none
  integer,intent(in) :: integ
  character(len=*),intent(in) :: string
  character(len=*),intent(out) :: strinn
 end subroutine appdig
end interface

interface
 subroutine atmdata(amu,rcov,symbol,znucl)
  use defs_basis
  implicit none
  real(dp),intent(out) :: amu
  real(dp),intent(out) :: rcov
  character(len=2),intent(out) :: symbol
  real(dp),intent(in) :: znucl
 end subroutine atmdata
end interface

interface
 subroutine symbol2znucl(amu,rcov,symbol,znucl)
  use defs_basis
  implicit none
  real(dp),intent(out) :: amu
  real(dp),intent(out) :: rcov
  character(len=2),intent(in) :: symbol
  real(dp),intent(out) :: znucl
 end subroutine symbol2znucl
end interface

interface
 subroutine atmlength(densty,length,zion,znucl)
  use defs_basis
  implicit none
  real(dp),intent(in) :: densty
  real(dp),intent(out) :: length
  real(dp),intent(in) :: zion
  real(dp),intent(in) :: znucl
 end subroutine atmlength
end interface

interface
 subroutine besjm(arg,besjx,cosx,nn,nx,sinx,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  integer,intent(in) :: nx
  real(dp),intent(in) :: arg
  real(dp),intent(out) :: besjx(nx)
  real(dp),intent(in) :: cosx(nx)
  real(dp),intent(in) :: sinx(nx)
  real(dp),intent(in) :: xx(nx)
 end subroutine besjm
end interface

interface
 subroutine blow_pawuj(mat,nj,matt)
  use defs_basis
  implicit none
  integer,intent(in) :: nj
  real(dp),intent(in) :: mat(nj,nj)
  real(dp),intent(out) :: matt(nj+1,nj+1)
 end subroutine blow_pawuj
end interface

interface
 subroutine bound_deriv(func,mesh,nn,yp1,ypn)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: nn
  type(pawrad_type),intent(in) :: mesh
  real(dp), intent(out) :: yp1
  real(dp), intent(out) :: ypn
  real(dp), intent(in) :: func(nn)
 end subroutine bound_deriv
end interface

interface
 subroutine chknm8(nmxpct,nmfond)
  implicit none
  character(len=9),intent(in) :: nmfond
  character(len=9),intent(in) :: nmxpct
 end subroutine chknm8
end interface

interface
 function clp(x)
  use defs_basis
  implicit none
  real(dp) :: clp
  real(dp),intent(in) :: x
 end function clp
end interface

interface
 subroutine ctrap(imax,ff,hh,ans)
  use defs_basis
  implicit none
  integer,intent(in) :: imax
  real(dp),intent(out) :: ans
  real(dp),intent(in) :: hh
  real(dp),intent(in) :: ff(imax)
 end subroutine ctrap
end interface

interface
 subroutine ctrap_gen(intg,func,radmesh)
  use defs_basis
  use defs_datatypes
  implicit none
  real(dp),intent(out) :: intg
  type(pawrad_type),intent(in) :: radmesh
  real(dp),intent(in) :: func(radmesh%int_meshsz)
 end subroutine ctrap_gen
end interface

interface
 subroutine derf(derf_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derf_yy
  real(dp),intent(in) :: yy
 end subroutine derf
end interface

interface
 subroutine derfc(derfc_yy,yy)
  use defs_basis
  implicit none
  real(dp),intent(out) :: derfc_yy
  real(dp),intent(in) :: yy
 end subroutine derfc
end interface

interface
 subroutine energies_init(energies)
  use defs_datatypes
  implicit none
  type(energies_type),intent(out) :: energies
 end subroutine energies_init
end interface

interface
 subroutine energies_copy(energies_in,energies_out)
  use defs_datatypes
  implicit none
  type(energies_type),intent(in) :: energies_in
  type(energies_type),intent(out) :: energies_out
 end subroutine energies_copy
end interface

interface
 function gaunt(ll,mm,l1,m1,l2,m2)
  use defs_basis
  implicit none
  integer,intent(in) :: l1
  integer,intent(in) :: l2
  integer,intent(in) :: ll
  integer,intent(in) :: m1
  integer,intent(in) :: m2
  integer,intent(in) :: mm
  real(dp) :: gaunt
 end function gaunt
end interface

interface
 subroutine hermit(chmin,chmout,ierr,ndim)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ndim
  real(dp),intent(inout) :: chmin(ndim*ndim+ndim)
  real(dp),intent(inout) :: chmout(ndim*ndim+ndim)
 end subroutine hermit
end interface

interface
 subroutine initylmr(mpsang,normchoice,npts,nrm,option,rr,ylmr,ylmr_gr)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  integer,intent(in) :: normchoice
  integer,intent(in) :: npts
  integer,intent(in) :: option
  real(dp),intent(in) :: nrm(npts)
  real(dp),intent(in) :: rr(3,npts)
  real(dp),intent(out) :: ylmr(mpsang*mpsang,npts)
  real(dp),optional,intent(out) :: ylmr_gr(3*(option/2)+6*(option/3),mpsang*mpsang,npts)
 end subroutine initylmr
end interface

interface
 subroutine interpol3d(r,nr1,nr2,nr3,denval,grid)
  use defs_basis
  implicit none
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  real(dp),intent(out) :: denval
  real(dp),intent(in) :: grid(nr1,nr2,nr3)
  real(dp),intent(in) :: r(3)
 end subroutine interpol3d
end interface

interface
 subroutine interpol3d_indices (r,nr1,nr2,nr3,&  
  &  ir1,ir2,ir3,  pr1,pr2,pr3)
  use defs_basis
  implicit none
  integer,intent(out) :: ir1
  integer,intent(out) :: ir2
  integer,intent(out) :: ir3
  integer,intent(in) :: nr1
  integer,intent(in) :: nr2
  integer,intent(in) :: nr3
  integer,intent(out) :: pr1
  integer,intent(out) :: pr2
  integer,intent(out) :: pr3
  real(dp),intent(in) :: r(3)
 end subroutine interpol3d_indices
end interface

interface
 subroutine inupper(string)
  implicit none
  character(len=*),intent(inout) :: string
 end subroutine inupper
end interface

interface
 subroutine isfile(filnam,status)
  use defs_basis
  implicit none
  character(len=fnlen),intent(inout) :: filnam
  character(len=3),intent(in) :: status
 end subroutine isfile
end interface

interface
 subroutine jbessel(bes,besp,bespp,ll,order,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: order
  real(dp),intent(out) :: bes
  real(dp),intent(out) :: besp
  real(dp),intent(out) :: bespp
  real(dp),intent(in) :: xx
 end subroutine jbessel
end interface

interface
 subroutine kramerskronig(nomega,omega,eps,method,only_check)
  use defs_basis
  implicit none
  integer,intent(in) :: method
  integer,intent(in) :: nomega
  integer,intent(in) :: only_check
  complex,intent(inout) :: eps(nomega)
  real(dp),intent(in) :: omega(nomega)
 end subroutine kramerskronig
end interface

interface
 subroutine lxyz(lp,mp,idir,ll,mm,lidir)
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ll
  integer,intent(in) :: lp
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  complex(dpc),intent(out) :: lidir
 end subroutine lxyz
end interface

interface
 subroutine mat_mlms2jmj(lcor,mat_mlms,mat_jmj,ndij,option,optspin,prtvol)
  use defs_basis
  implicit none
  integer,intent(in) :: lcor
  integer,intent(in) :: ndij
  integer,intent(in) :: option
  integer,intent(in) :: optspin
  integer,intent(in) :: prtvol
  complex(dpc),intent(inout) :: mat_jmj(2*(2*lcor+1),2*(2*lcor+1))
  complex(dpc),intent(inout) :: mat_mlms(2*lcor+1,2*lcor+1,ndij)
 end subroutine mat_mlms2jmj
end interface

interface
 subroutine mat_slm2ylm(lcor,mat_inp_c,mat_out_c,ndij,option,optspin,prtvol)
  use defs_basis
  implicit none
  integer,intent(in) :: lcor
  integer,intent(in) :: ndij
  integer,intent(in) :: option
  integer,intent(in) :: optspin
  integer,intent(in) :: prtvol
  complex(dpc) :: mat_inp_c(2*lcor+1,2*lcor+1,ndij)
  complex(dpc) :: mat_out_c(2*lcor+1,2*lcor+1,ndij)
 end subroutine mat_slm2ylm
end interface

interface
 subroutine matcginv(a,lda,n)
  use defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  complex(gwpc),intent(inout) :: a(lda,n)
 end subroutine matcginv
end interface

interface
 subroutine matcginv_dpc(a,lda,n)
  use defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  complex(dpc),intent(inout) :: a(lda,n)
 end subroutine matcginv_dpc
end interface

interface
 subroutine mati3inv(mm,mit)
  implicit none
  integer,intent(out) :: mit(3,3)
  integer,intent(in) :: mm(3,3)
 end subroutine mati3inv
end interface

interface
 subroutine matr3eigval(eigval,matr)
  use defs_basis
  implicit none
  real(dp),intent(out) :: eigval(3)
  real(dp),intent(in) :: matr(3,3)
 end subroutine matr3eigval
end interface

interface
 subroutine matr3inv(aa,ait)
  use defs_basis
  implicit none
  real(dp),intent(in) :: aa(3,3)
  real(dp),intent(out) :: ait(3,3)
 end subroutine matr3inv
end interface

interface
 subroutine matrginv(a,lda,n)
  use defs_basis
  implicit none
  integer,intent(in) :: lda
  integer,intent(in) :: n
  real(dp),intent(inout) :: a(lda,n)
 end subroutine matrginv
end interface

interface
 subroutine mkherm(array,ndim)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  real(dp),intent(inout) :: array(2,ndim,ndim)
 end subroutine mkherm
end interface

interface
 subroutine mknormpath(nbounds,bounds,gmet,ndiv_small,ndiv,npt_tot,path)
  use defs_basis
  implicit none
  integer,intent(in) :: nbounds
  integer,intent(in) :: ndiv_small
  integer,intent(inout) :: npt_tot
  real(dp),intent(in) :: bounds(3,nbounds)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(inout) :: ndiv(nbounds-1)
  real(dp),intent(out),optional :: path(3,npt_tot)
 end subroutine mknormpath
end interface

interface
 subroutine nderiv(hh,yy,zz,ndim,norder)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  integer,intent(in) :: norder
  real(dp),intent(in) :: hh
  real(dp),intent(in) :: yy(ndim)
  real(dp),intent(out) :: zz(ndim)
 end subroutine nderiv
end interface

interface
 subroutine normev(evec,ndim,num)
  use defs_basis
  implicit none
  integer,intent(in) :: ndim
  integer,intent(in) :: num
  real(dp),intent(inout) :: evec(2*ndim,num)
 end subroutine normev
end interface

interface
 subroutine pclock(itimpt,unit,mode_paral)
  implicit none
  integer,intent(in) :: itimpt
  integer,intent(in),optional :: unit
  character(len=4),intent(in),optional :: mode_paral
 end subroutine pclock
end interface

interface
 function permutations(nn,kk)
  use defs_basis
  implicit none
  integer,intent(in) :: kk
  integer,intent(in) :: nn
  real(dp) :: permutations
 end function permutations
end interface

interface
 function phim(costheta,sintheta,mm)
  use defs_basis
  implicit none
  integer,intent(in) :: mm
  real(dp),intent(in) :: costheta
  real(dp) :: phim
  real(dp),intent(in) :: sintheta
 end function phim
end interface

interface
 subroutine pl_deriv(mpsang,pl_d2,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: pl_d2(mpsang)
 end subroutine pl_deriv
end interface

interface
 subroutine plm_coeff(blm,mpsang,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: blm(5,mpsang*mpsang)
 end subroutine plm_coeff
end interface

interface
 subroutine plm_d2theta(mpsang,plm_d2t,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: mpsang
  real(dp),intent(in) :: xx
  real(dp),intent(out) :: plm_d2t(mpsang*mpsang)
 end subroutine plm_d2theta
end interface

interface
 function plm_dphi(ll,mm,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: mm
  real(dp) :: plm_dphi
  real(dp),intent(in) :: xx
 end function plm_dphi
end interface

interface
 function plm_dtheta(ll,mm,xx)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: mm
  real(dp) :: plm_dtheta
  real(dp),intent(in) :: xx
 end function plm_dtheta
end interface

interface
 subroutine print_ij(a_ij,adim,cplex,ndim,opt_io,opt_l,opt_l_index,opt_pack,opt_prtvol,pack2ij,test_value,unt,&  
  &  opt_sym,asym_ij)    !Optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: adim
  integer,intent(in) :: cplex
  integer,intent(in) :: ndim
  integer,intent(in) :: opt_io
  integer,intent(in) :: opt_l
  integer,intent(in) :: opt_pack
  integer,intent(in) :: opt_prtvol
  integer,intent(in),optional :: opt_sym
  integer,intent(in) :: unt
  real(dp),intent(in) :: test_value
  real(dp),intent(in) :: a_ij(cplex*adim)
  real(dp),intent(in),optional :: asym_ij(cplex*adim)
  integer,intent(in) :: opt_l_index(ndim*min(1+opt_l,1))
  integer,intent(in) :: pack2ij(adim*opt_pack)
 end subroutine print_ij
end interface

interface
 subroutine printxsf(n1,n2,n3,datagrid,basis,origin,natom, ntypat, typat, xcart, znucl, nunit,realrecip)
  use defs_basis
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  integer,intent(in) :: nunit
  integer,intent(in) :: realrecip
  real(dp),intent(in) :: basis(3,3)
  real(dp),intent(in) :: datagrid(n1*n2*n3)
  real(dp),intent(in) :: origin(3)
  integer,intent(in) :: typat(ntypat)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine printxsf
end interface

interface
 subroutine prmat (mat, ni, nj, mi, unitm)
  use defs_basis
  implicit none
  integer,intent(in) :: mi
  integer,intent(in) :: ni
  integer,intent(in) :: nj
  integer,intent(in), optional :: unitm
  real(dp),intent(in) :: mat(mi,nj)
 end subroutine prmat
end interface

interface
 subroutine ratint(npts,xin,xpt,yin,yerr,ypt)
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  real(dp),intent(in) :: xpt
  real(dp),intent(out) :: yerr
  real(dp),intent(out) :: ypt
  real(dp),intent(in) :: xin(npts)
  real(dp),intent(in) :: yin(npts)
 end subroutine ratint
end interface

interface
 subroutine realgaunt(l_max,ngnt,gntselect,realgnt)
  use defs_basis
  implicit none
  integer,intent(in) :: l_max
  integer,intent(out) :: ngnt
  integer,intent(out) :: gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)
  real(dp),intent(out) :: realgnt((2*l_max-1)**2*(l_max)**4)
 end subroutine realgaunt
end interface

interface
 subroutine rotmat(xaxis,zaxis,inversion_flag,umat)
  use defs_basis
  implicit none
  integer,intent(out) :: inversion_flag
  real(dp),intent(out) :: umat(3,3)
  real(dp),intent(in) :: xaxis(3)
  real(dp),intent(in) :: zaxis(3)
 end subroutine rotmat
end interface

interface
 function set_istwfk(kpoint) result(istwfk)
  use defs_basis
  implicit none
  integer :: istwfk
  real(dp),intent(in) :: kpoint(3)
 end function set_istwfk
end interface

interface
 subroutine simpson_int(nsimpson,simp_delta,simp_funct,simp_res)
  use defs_basis
  implicit none
  integer,intent(in) :: nsimpson
  real(dp),intent(in) :: simp_delta
  real(dp),intent(in) :: simp_funct(nsimpson)
  real(dp),intent(out) :: simp_res(nsimpson)
 end subroutine simpson_int
end interface

interface
 subroutine slxyzs(lp,mp,idir,ll,mm,sls_val)
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: ll
  integer,intent(in) :: lp
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  complex(dpc),intent(out) :: sls_val
 end subroutine slxyzs
end interface

interface
 subroutine solvbes(root,alpha,beta,ll,nq)
  use defs_basis
  implicit none
  integer :: ll
  integer :: nq
  real(dp) :: alpha
  real(dp) :: beta
  real(dp) :: root(nq)
 end subroutine solvbes
end interface

interface
 subroutine status(counter,filstat,istat,level,routine)
  use defs_basis
  implicit none
  integer,intent(in) :: counter
  integer,intent(in) :: istat
  integer,intent(in) :: level
  character(len=fnlen),intent(in) :: filstat
  character(len=*),intent(in) :: routine
 end subroutine status
end interface

interface
 subroutine symq3(nsym,qpt,symq,symrec,timrev,prtvol)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in),optional :: prtvol
  integer,intent(out) :: timrev
  real(dp),intent(in) :: qpt(3)
  integer,intent(out) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine symq3
end interface

interface
 subroutine wrap2_pmhalf(num,red,shift)
  use defs_basis
  implicit none
  real(dp),intent(in) :: num
  real(dp),intent(out) :: red
  real(dp),intent(out) :: shift
 end subroutine wrap2_pmhalf
end interface

interface
 subroutine wrap2_zero_one(num,red,shift)
  use defs_basis
  implicit none
  real(dp),intent(in) :: num
  real(dp),intent(out) :: red
  real(dp),intent(out) :: shift
 end subroutine wrap2_zero_one
end interface

interface
 subroutine ylm_cmplx(lx,ylm,xx,yy,zz)
  use defs_basis
  implicit none
  integer,intent(in) :: lx
  real(dp),intent(in) :: xx
  real(dp),intent(in) :: yy
  real(dp),intent(in) :: zz
  complex(dpc),intent(out) :: ylm((lx+1)*(lx+1))
 end subroutine ylm_cmplx
end interface

interface
 subroutine ys(lp,mp,ll,mm,ys_val)
  use defs_basis
  implicit none
  integer,intent(in) :: ll
  integer,intent(in) :: lp
  integer,intent(in) :: mm
  integer,intent(in) :: mp
  complex(dpc),intent(out) :: ys_val
 end subroutine ys
end interface

end module interfaces_32_util
!!***
