!{\src2tex{textfont=tt}}
!!****f* ABINIT/dsksta
!! NAME
!! dsksta
!!
!! FUNCTION
!! This routine evaluates the amount of disk space required by the _KSS file.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (MT,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dimlmn(natom*usepaw)=Number of nlm partial waves for each atom.
!!  ishm=Number of G-shells to be saved in _KSS file.
!!  mpsang=Max angular momentum +1 for pseudos.
!!  natom=Number of atoms in the unit cell.
!!  nbandkss=Number of desired bands to be saved in _KSS file
!!  nkpt=Number of k points.
!!  npwkss=Number of desired G-vectors to be saved in _KSS file.
!!  nspinor=Number of spinorial components.
!!  nsppol=Number of independent spin polarizations.
!!  ntypat=Number of type of atoms.
!!  nsym2=Number of symmetries in space group, without INV
!!  usepaw=1 if PAW.
!!
!! OUTPUT
!!  Writes on standard output
!!
!! TODO 
!!  Should be called inside memory to give an estimation of the size of the file.
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

subroutine dsksta(ishm,usepaw,nbandkss,mpsang,natom,ntypat,npwkss,nkpt,nspinor,nsppol,nsym2,dimlmn)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: usepaw,ishm,nbandkss,mpsang,natom,ntypat,nkpt
 integer,intent(in) :: npwkss,nspinor,nsppol,nsym2
!arrays
 integer,intent(in) :: dimlmn(natom*usepaw)

!Local variables-------------------------------
!scalars
 integer :: bsize_tot,bsize_hdr,bsize_kb,bsize_wf,bsize_cprj
 character(len=500) :: msg

! *********************************************************************

!The Abinit header is not considered.
 bsize_hdr= 80*2 + & !title 
&5*4 + & !nsym2,nbandksseff,npwkss,ishm,mpsang
&nsym2*9*4 + & !symrel2
&nsym2*3*8 + & !tnons
&npwkss*3*4 + & !gbig
&ishm*4     !shlim

!NOTE: vkb does not depend on nsppol, however the elements are written for each spin.
 bsize_kb=0 
 if (usepaw==0) then
   bsize_kb= nsppol* &
&   (         mpsang*ntypat       *8 + & !vkbsign
&  2*(nkpt*mpsang*ntypat*npwkss*8)  & !vkbd,vkbd
&  )
 end if

 bsize_wf= nsppol* &
& ( nkpt*nbandkss                 *8 + & !energies
&nkpt*nbandkss*nspinor*npwkss*2*8   & !wfg
&)

!For PAW add space required by projectors.
 bsize_cprj=0
 if (usepaw==1) then 
   bsize_cprj=SUM(dimlmn(:))*(nsppol*nkpt*nspinor*nbandkss*2*8) 
 end if

 bsize_tot = bsize_hdr + bsize_kb + bsize_wf + bsize_cprj
 write(msg,'(2a,f8.2,4a,4(a,f8.2,2a))')ch10,&
& ' Total amount of disk space required by _KSS file = ',bsize_tot*b2Mb,' Mb.',ch10,&
& '  Subdivided into : ',ch10,&
& '  Header             = ',bsize_hdr *b2Mb,' Mb.',ch10,&
& '  KB elements        = ',bsize_kb  *b2Mb,' Mb.',ch10,&
& '  Wavefunctions (PW) = ',bsize_wf  *b2Mb,' Mb.',ch10,&
& '  PAW projectors     = ',bsize_cprj*b2Mb,' Mb.',ch10
 call wrtout(std_out,msg,'COLL')

end subroutine dsksta
!!***
