!{\src2tex{textfont=tt}}
!!****f* ABINIT/matrginv
!! NAME
!! matrginv
!!
!! FUNCTION
!! Invert a general matrix of real*8 elements.
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
!! a=matrix of real elements
!! OUTPUT
!! a=inverse of a input matrix
!!
!! SIDE EFFECTS
!! a(lda,n)= array of real elements, input, inverted at output
!!
!!
!! PARENTS
!!      elast9,linear_optics_paw,mlwfovlp_ylmfar,piezo9,q0fit,setup_bse,spectra
!!      strainsym,trashme
!!
!! CHILDREN
!!      dbgmdi,dbgmlu,dgeicd,dgetrf,dgetri,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine matrginv(a,lda,n)

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
 real(dp),intent(inout) :: a(lda,n)

!Local variables-------------------------------
!scalars
 integer :: ierr,nwork
#if defined HAVE_LINALG_ESSL
 real(dp) :: rcond
#endif
 character(len=500) :: message
!arrays
 integer,allocatable :: ipvt(:)
#if defined HAVE_LINALG_ESSL
 real(dp) :: det(2)
#elif defined HAVE_LINALG_ASL
 real(dp) :: det(2)
#endif
 real(dp),allocatable :: work(:)

! *************************************************************************

#if defined HAVE_LINALG_ESSL
 nwork=200*n
#else
 nwork=n
#endif

 allocate(work(nwork))
 allocate(ipvt(n))


#if defined HAVE_LINALG_ESSL

 call dgeicd(a,lda,n,0,rcond,det,work,nwork)
 if(abs(rcond)==zero) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ESSL routine dgeicd failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

#elif defined HAVE_LINALG_ASL

 call dbgmlu(a,lda,n,ipvt,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ASL routine dbgmlu failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 call dbgmdi(a,lda,n,ipvt,det,-1,work,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The ASL routine dbgmdi failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

#else

 call dgetrf(n,n,a,lda,ipvt,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine dgetrf failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 call dgetri(n,a,lda,ipvt,work,n,ierr)
 if(ierr /= 0) then
   write(message, '(10a)' ) ch10,&
&   ' matrginv : BUG -',ch10,&
&   '  The matrix that has been passed in argument of this subroutine',ch10,&
&   '  is probably either singular or nearly singular.',ch10,&
&   '  The LAPACK routine dgetri failed.',ch10,&
&   '  Action : Contact ABINIT group '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

#endif

 deallocate(work)
 deallocate(ipvt)         ! added by MM on Nov. 26

end subroutine matrginv
!!***
