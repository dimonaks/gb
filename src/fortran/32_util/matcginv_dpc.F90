!{\src2tex{textfont=tt}}
!!****f* ABINIT/matcginv_dpc
!! NAME
!! matcginv_dpc
!!
!! FUNCTION
!! Invert a general matrix of complex elements.
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
!! OUTPUT
!! a=inverse of a input matrix
!! SIDE EFFECTS
!! a(lda,n)= array of complex elements, input, inverted at output
!!
!!
!! PARENTS
!!      m_matlu,m_matrix,m_oper
!!
!! CHILDREN
!!      leave_new,wrtout,zgetrf,zgetri
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine matcginv_dpc(a,lda,n)

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
 complex(dpc),intent(inout) :: a(lda,n)

!Local variables-------------------------------
!scalars
 integer :: ierr,nwork
 character(len=500) :: message
!arrays
 integer :: ipvt(n)
 complex(dpc),allocatable :: work(:)
!no_abirules

! *************************************************************************

 nwork=n

 allocate(work(nwork))

 call zgetrf(n,n,a,lda,ipvt,ierr)
 if(ierr/=0) then
   write(message,'(8a,i4,2a)')ch10,&
&   ' matcginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine zgetrf failed.',ierr,ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if
 call zgetri(n,a,n,ipvt,work,n,ierr)
 if(ierr/=0) then
   write(message,'(10a)')ch10,&
&   ' matcginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine zgetri failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

 deallocate(work)

end subroutine matcginv_dpc
!!***
