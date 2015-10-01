!{\src2tex{textfont=tt}}
!!****f* ABINIT/matcginv
!! NAME
!! matcginv
!!
!! FUNCTION
!! Invert a general matrix of complex elements.
!!  CGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
!!  CGETRI computes the inverse of a matrix using the LU factorization computed by CGETRF.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2010 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! lda=leading dimension of complex matrix a
!! n=size of complex matrix a
!! a=matrix of complex elements
!!
!! OUTPUT
!! a=inverse of a input matrix
!! SIDE EFFECTS
!! a(lda,n)= array of complex elements, input, inverted at output
!!
!! PARENTS
!!      calc_ffm,calc_rpa_functional,eps1_tc,make_epsm1_driver
!!
!! CHILDREN
!!      cbgmdi,cbgmlu,cgeicd,cgetrf,cgetri,wrtout,zgetrf,zgetri
!!
!! NOTES
!!  MG: This routine is obsolete, m_abilasi provides a better version supporting overloading and Scalapack
!!  We keep matcginv simply because m_abilasi is an high level module depending on MPI and Scalapack stuff.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine matcginv(a,lda,n)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lda,n
!arrays
 complex(gwpc),intent(inout) :: a(lda,n)

!Local variables-------------------------------
!scalars
 integer :: ierr,nwork
 character(len=500) :: message
!arrays
 integer :: ipvt(n)
 complex(gwpc),allocatable :: work(:)
!no_abirules
#if defined HAVE_LINALG_ASL
 real(dp) :: det
 complex(gwpc) :: cdet
#endif

! *************************************************************************

 nwork=n

 allocate(work(nwork))

#if defined HAVE_LINALG_ASL
 call cbgmlu(a,lda,n,ipvt,ierr)
 if(ierr/=0) then
   write(message,'(10a)')ch10,&
&   ' matcginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ASL routine cbgmlu failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if
 call cbgmdi(a,lda,n,ipvt,cdet,det,-1,work,ierr)
 if(ierr/=0) then
   write(message,'(10a)')ch10,&
&   ' matcginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ASL routine dbgmdi failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if
#else
#if defined HAVE_GW_DPC
 call zgetrf(n,n,a,lda,ipvt,ierr)
#else
 call cgetrf(n,n,a,lda,ipvt,ierr)
#endif
 if(ierr/=0) then
   write(message,'(10a)')ch10,&
&   ' matcginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine cgetrf failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if
#if defined HAVE_GW_DPC
 call zgetri(n,a,n,ipvt,work,n,ierr)
#else
 call cgetri(n,a,n,ipvt,work,n,ierr)
#endif
 if(ierr/=0) then
   write(message,'(10a)')ch10,&
&   ' matcginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine cgetri failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if
#endif

 deallocate(work)

end subroutine matcginv
!!***
