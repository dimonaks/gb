!{\src2tex{textfont=tt}}
!!****f* ABINIT/ingeo
!!
!! NAME
!! ingeo
!!
!! FUNCTION
!! Initialize geometry variables for the ABINIT code.
!! 1) set up unit cell : acell, rprim and rprimd ; deduce Bravais lattice
!! 2) (removed)
!! 3) Set up the number of atoms (natrd) in the primitive set, to be read.
!! 4) Read the type of each atom in the primitive set
!! 5) Read coordinates for each atom in the primitive set
!! 6) Eventually read the symmetries
!! 7) Checks whether the geometry builder must be used,
!!    and call it if needed. Call eventually the symmetry builder and analyser
!!    Make the adequate transfers if the geometry
!!    builder is not needed.
!! 8) Initialize the fixing of atoms,
!!    the initial velocities, and the initial atomic spin
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG, RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! berryopt == 4: electric field is on; berryopt /= 4: electric field is off
!! iimage= index of the current image
!! iout=unit number of output file
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! msym=default maximal number of symmetries
!! natom=number of atoms
!! nimage=number of images
!! nspden=number of spin-density components
!! nsppol=number of independent spin polarizations
!! ntypat=number of type of atoms
!! pawspnorb=1 when spin-orbit is activated within PAW 
!! string*(*)=character string containing all the input data, used
!!  only if choice=1 or 3. Initialized previously in instrng.
!!
!! OUTPUT
!! acell(3)=length of primitive vectors
!! bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!! genafm(3)=magnetic translation generator (in case of Shubnikov group type IV)
!! iatfix(3,natom)=indices for atoms fixed along some (or all) directions
!! jellslab=not zero if jellslab keyword is activated
!! slabzbeg, slabzend= the z coordinates of beginning / end of the jellium slab
!! nsym=actual number of symmetries
!! ptgroupma = magnetic point group number
!! rprim(3,3)=dimensionless real space primitive translations
!! spgroup=symmetry space group
!! spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!! symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!! symmorphi=if 0, only allows symmorphic symmetry operations
!! symrel(3,3,1:msym)=symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!! tolsym=tolerance for the symmetry operations
!! typat(natom)=type integer for each atom in cell
!! vel(3,natom)=initial velocity of atom in bohr/atomic time units
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! SIDE EFFECTS
!! efield_xcart=use when berryopt == 4, contains the electric field
!!              in cartesian coordinates (INOUT is due to xcartxred())
!!
!! NOTES
!! the parameters ntypat and natom have already been read in indims,
!! and were used to dimension the arrays needed here.
!!
!!
!! PARENTS
!!      invars1
!!
!! CHILDREN
!!      chkorthsy,fillcell,gensymshub,gensymshub4,gensymspgr,ingeobld,intagm
!!      leave_new,mati3inv,metric,mkrdim,symanal,symatm,symfind,symlatt
!!      symrelrot,symzat,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ingeo (acell,berryopt,bravais,efield_xcart,&
& genafm,iatfix,icoulomb,iimage,iout,jdtset,jellslab,lenstr,&
& msym,natom,nimage,nspden,nsppol,nsym,ntypat,pawspnorb,&
& ptgroupma,rprim,slabzbeg,slabzend,spgroup,spinat,string,symafm,&
& symmorphi,symrel,tnons,tolsym,typat,vel,xred)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_42_parser
 use interfaces_57_iovars, except_this_one => ingeo
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,iimage,iout,jdtset,lenstr,msym,nimage,nspden,nsppol,ntypat,pawspnorb
 integer,intent(inout) :: natom,symmorphi
 integer,intent(out) :: icoulomb,jellslab,nsym,ptgroupma,spgroup
 real(dp),intent(out) :: slabzbeg,slabzend,tolsym
 character(len=*),intent(in) :: string
!arrays
 integer,intent(out) :: bravais(11),iatfix(3,natom),symafm(msym)
 integer,intent(out) :: symrel(3,3,msym),typat(natom)
 real(dp),intent(inout) :: efield_xcart(3),spinat(3,natom)
 real(dp),intent(out) :: acell(3),genafm(3),rprim(3,3),tnons(3,msym)
 real(dp),intent(out) :: vel(3,natom),xred(3,natom)

!Local variables-------------------------------
 character(len=*), parameter :: format01110 ="(1x,a6,1x,(t9,8i8) )"
 character(len=*), parameter :: format01160 ="(1x,a6,1x,1p,(t9,3g18.10)) "
!scalars
 integer :: bckbrvltt,brvltt,chkprim,iatom,idir,ii,irreducible,isym
 integer :: jsym,marr,mu,natfix,natrd,nobj,noncoll
 integer :: nptsym,nsym_now,shubnikov,spgaxor,spgorig
 integer :: spgroupma,tacell,tangdeg,tgenafm,tnatrd,tread,trprim,tscalecart,tspgroupma
 integer :: txangst,txcart,txred,use_inversion
 real(dp) :: alpha,a2,aa,cc,cosang,ucvol
 character(len=30) :: token,strimg
 character(len=5) :: stringimage
 character(len=500) :: message
!arrays
 integer,allocatable :: ptsymrel(:,:,:),typat_read(:),symrec(:,:,:),indsym(:,:,:)
 real(dp) :: angdeg(3),efield_xred(3),gmet(3,3),gprimd(3,3),rmet(3,3),rcm(3)
 real(dp) :: rprimd(3,3),rprimd_new(3,3),scalecart(3)
!real(dp) :: tsec(2)
 real(dp),allocatable :: tnons_cart(:,:),xangst_read(:,:)
 real(dp),allocatable :: xcart(:,:),xcart_read(:,:),xred_read(:,:)
!no_abirules
!Dummy arguments for subroutine 'intagm' to parse input file
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

! *************************************************************************

!call timab(47,1,tsec)

!DEBUG
!write(6,*)' ingeo : enter '
!call flush(6)
!stop
!ENDDEBUG

 marr=max(12,3*natom,9*msym)
 allocate(intarr(marr),dprarr(marr))

 if(nimage/=1)then
   if(iimage<10)then
     write(stringimage,'(i1)')iimage
   else if(iimage<100)then
     write(stringimage,'(i2)')iimage
   else if(iimage<1000)then
     write(stringimage,'(i3)')iimage
   else if(iimage<10000)then
     write(stringimage,'(i4)')iimage
   else
     write(std_out,'(3a)')' ingeo : BUG -',ch10,&
&     ' Still to be coded '
     call leave_new('COLL')
   end if
   strimg='_'//trim(stringimage)//'img'
   alpha=real(iimage-1)/real(nimage-1)
 end if

!1) set up unit cell : acell, rprim and rprimd ---------------------

 acell(1:3)=one
 token = 'acell'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tacell,'LEN')
 if(tacell==1) acell(1:3)=dprarr(1:3)

 if(nimage/=1 .and. iimage/=1)then
   token = 'acell_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'LEN')
   if(tread==1)then
     acell(1:3)=acell(1:3)+alpha*(dprarr(1:3)-acell(1:3))
     tacell=1
   end if
   token = 'acell'//trim(strimg)
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'LEN')
   if(tread==1)then
     acell(1:3)=dprarr(1:3)
     tacell=1
   end if
 end if

 scalecart(1:3)=one
 token = 'scalecart'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tscalecart,'LEN')
 if(tscalecart==1) scalecart(1:3)=dprarr(1:3)

 if(nimage/=1 .and. iimage/=1)then
   token = 'scalecart_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'LEN')
   if(tread==1)then
     scalecart(1:3)=scalecart(1:3)+alpha*(dprarr(1:3)-scalecart(1:3))
     tscalecart=1
   end if
   token = 'scalecart'//trim(strimg)
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'LEN')
   if(tread==1)then
     scalecart(1:3)=dprarr(1:3)
     tscalecart=1
   end if
 end if


 if(nimage==1)then
!  Alternate SIESTA input variable for acell
   token = 'LatticeConstant'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
   if(tread==1)then
     if(tacell==1)then
       write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&       ' ingeo: ERROR -',ch10,&
&       '  acell and LatticeConstant cannot be specified simultaneously',ch10,&
&       '  for the same dataset.',ch10,&
&       '  Action : check the input file.'
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     else
!      Note that LatticeConstant is a scalar
       acell(1:3)=dprarr(1)
       tacell=1
     end if
   end if

!  Might initialize scalecart from CML file
   if(tacell==0)then
     token = '_acell'
     call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tacell,'LEN')
     if(tacell==1) acell(1:3)=dprarr(1:3)
   end if
 end if

!Check that input length scales acell(3) are > 0
 do mu=1,3
   if(acell(mu)<=0.0_dp) then
     write(message, '(a,a,a,a,i5,a, 1p,e14.6,a,a,a,a)' ) ch10,&
&     ' ingeo: ERROR -',ch10,&
&     '  Length scale',mu,' is input as acell=',acell(mu),ch10,&
&     '  However, length scales must be > 0 ==> stop',ch10,&
&     '  Action : correct acell in input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end do

!Initialize rprim, or read the angles
 tread=0
 token = 'rprim'
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),token,trprim,'DPR')
 if(trprim==1)rprim(:,:)=reshape( dprarr(1:9) , (/3,3/) )
 if(nimage/=1 .and. iimage/=1)then
   token = 'rprim_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     rprim(:,:)=rprim(:,:)+alpha*(reshape( dprarr(1:9) , (/3,3/) )-rprim(:,:))
     trprim=1
   end if
   token = 'rprim'//trim(strimg)
   call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     rprim(:,:)=reshape( dprarr(1:9) , (/3,3/) )
     trprim=1
   end if
 end if

!If none of the rprim were read ...
 if(trprim==0)then
   token = 'angdeg'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tangdeg,'DPR')
   angdeg(:)=dprarr(1:3)
   if(nimage/=1 .and. iimage/=1)then
     token = 'angdeg_lastimg'
     call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
     if(tread==1)then
       angdeg(:)=angdeg(:)+alpha*(dprarr(1:3)-angdeg(:))
       tangdeg=1
     end if
     token = 'angdeg'//trim(strimg)
     call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
     if(tangdeg==1)then
       angdeg(:)=dprarr(1:3)
       tangdeg=1
     end if
   end if

!  Might initialize angdeg from CML file
   if(tangdeg==0 .and. nimage==1)then
     token = '_angdeg'
     call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tangdeg,'DPR')
     angdeg(:)=dprarr(1:3)
   end if
   if(tangdeg==1)then
     write(message, '(a,a)' )ch10,&
&     ' ingeo : use angdeg to generate rprim.'
     call wrtout(std_out,message,'COLL')

!    Check that input angles are positive
     do mu=1,3
       if(angdeg(mu)<=0.0_dp) then
         write(message, '(a,a,a,a,i5,a,1p,e14.6,a,a,a,a)' ) ch10,&
&         ' ingeo: ERROR -',ch10,&
&         '  Angle number',mu,' is input as angdeg=',angdeg(mu),ch10,&
&         '  However, angles must be > 0 ==> stop',ch10,&
&         '  Action : correct angdeg in input file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do

!    Check that the sum of angles is smaller than 360 degrees
     if(angdeg(1)+angdeg(2)+angdeg(3)>=360.0_dp) then
       write(message, '(a,a,a,a,a,a,es14.4,a,a,a)' ) ch10,&
&       ' ingeo: ERROR -',ch10,&
&       '  The sum of input angles (angdeg(1:3)) must be lower than 360 degrees',ch10,&
&       '  while it is ',angdeg(1)+angdeg(2)+angdeg(3),'.',ch10,&
&       '  Action : correct angdeg in input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     if( abs(angdeg(1)-angdeg(2))<tol12 .and. &
&     abs(angdeg(2)-angdeg(3))<tol12 .and. &
&     abs(angdeg(1)-90._dp)+abs(angdeg(2)-90._dp)+abs(angdeg(3)-90._dp)>tol12 )then
!      Treat the case of equal angles (except all right angles) :
!      generates trigonal symmetry wrt third axis
       cosang=cos(pi*angdeg(1)/180.0_dp)
       a2=2.0_dp/3.0_dp*(1.0_dp-cosang)
       aa=sqrt(a2)
       cc=sqrt(1.0_dp-a2)
       rprim(1,1)=aa        ; rprim(2,1)=0.0_dp                 ; rprim(3,1)=cc
       rprim(1,2)=-0.5_dp*aa ; rprim(2,2)= sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,2)=cc
       rprim(1,3)=-0.5_dp*aa ; rprim(2,3)=-sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,3)=cc
!      DEBUG
!      write(6,*)' ingeo : angdeg=',angdeg(1:3)
!      write(6,*)' ingeo : aa,cc=',aa,cc
!      ENDDEBUG
     else
!      Treat all the other cases
       rprim(:,:)=0.0_dp
       rprim(1,1)=1.0_dp
       rprim(1,2)=cos(pi*angdeg(3)/180.0_dp)
       rprim(2,2)=sin(pi*angdeg(3)/180.0_dp)
       rprim(1,3)=cos(pi*angdeg(2)/180.0_dp)
       rprim(2,3)=(cos(pi*angdeg(1)/180.0_dp)-rprim(1,2)*rprim(1,3))/rprim(2,2)
       rprim(3,3)=sqrt(1.0_dp-rprim(1,3)**2-rprim(2,3)**2)
     end if

   end if
!  No problem if neither rprim nor angdeg are defined : use default rprim
 end if

!Rescale rprim using scalecart (and set scalecart to one)
 rprim(:,1)=scalecart(:)*rprim(:,1)
 rprim(:,2)=scalecart(:)*rprim(:,2)
 rprim(:,3)=scalecart(:)*rprim(:,3)
 scalecart(:)=one

!Compute different matrices in real and reciprocal space, also
!checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 tolsym=tol8
 token = 'tolsym'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) tolsym=dprarr(1)

!Find a tentative Bravais lattice and its point symmetries (might not use them)
!Note that the Bravais lattice might not be the correct one yet (because the
!actual atomic locations might lower the symattry obtained from the lattice parameters only)
 allocate(ptsymrel(3,3,msym))
 call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)

!Get efield in coordinates relative to lattice vectors
 if (berryopt==4) then
   call xredxcart(1,-1,rprimd,efield_xcart,efield_xred)
 end if

!3) Possibly, initialize a jellium slab
 jellslab=0
 token = 'jellslab'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) jellslab=intarr(1)

 if(jellslab/=0)then
   slabzbeg=zero
   token = 'slabzbeg'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) slabzbeg=dprarr(1)

   slabzend=zero
   token = 'slabzend'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) slabzend=dprarr(1)
 end if

!4) Set up the number of atoms in the primitive set, to be read.

!This is the default
 natrd=natom

 token = 'natrd'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tnatrd,'INT')
 if(tnatrd==1) natrd=intarr(1)

 if(natrd<1 .or. natrd>natom)then
   write(message, '(a,a,a,a,a,a,a,a,i5,a,i5,a,a,a)' ) ch10,&
&   ' ingeo : ERROR -',ch10,&
&   '  The number of atoms to be read (natrd) must be positive',&
&   ch10,'  and not bigger than natom.',ch10,&
&   '  This is not the case : natrd=',natrd,', natom=',natom,'.',ch10,&
&   '  Action : correct natrd or natom in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if


!5) Read the type and initial spin of each atom in the primitive set--------

!Check for the use of the old name of this variable
 token = 'type'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' ingeo : ERROR -',ch10,&
&   '  The use of the "type" input variable is forbidden since version 4.1 .',ch10,&
&   '  Action : replace "type" by "typat".'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 allocate(typat_read(natrd))
 typat_read(1)=1

 token = 'typat'
 call intagm(dprarr,intarr,jdtset,marr,natrd,string(1:lenstr),token,tread,'INT')

!If not read, try the CML data
 if(tread==0)then

!  DEBUG
!  write(6, '(a)' ) ' ingeo : before intagm _typat'
!  write(6, '(a)' ) string(1:lenstr)
!  ENDDEBUG
   token = '_typat'
   call intagm(dprarr,intarr,jdtset,marr,natrd,string(1:lenstr),token,tread,'INT')
 end if
 if(tread==1) typat_read(1:natrd)=intarr(1:natrd)

 do iatom=1,natrd
   if(typat_read(iatom)<1 .or. typat_read(iatom)>ntypat )then
     write(message, '(a,a,a,a,i4,a,i3,a,a,a,i3,a,a,a)' )ch10,&
&     ' ingeo : ERROR - ',ch10,&
&     '  The input type of atom number',iatom,&
&     ' is equal to ',typat_read(iatom),',',ch10,&
&     '  while it should be between 1 and ntypat=',ntypat,'.',ch10,&
&     '  Action : change either the variable typat or the variable ntypat.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end do

!6) Read coordinates for each atom in the primitive set--------

 allocate(xangst_read(3,natrd),xcart_read(3,natrd),xred_read(3,natrd))

 token = 'xred'
 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),&
& token,txred,'DPR')
 if(txred==1) xred_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 if(nimage/=1 .and. iimage/=1)then
   token = 'xred_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     xred_read(:,1:natrd)=xred_read(:,1:natrd)+&
&     alpha*(reshape( dprarr(1:3*natrd) , (/3,natrd/) )-xred_read(:,1:natrd))
     txred=1
   end if
   token = 'xred'//trim(strimg)
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     xred_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
     txred=1
   end if
 end if

!DEBUG
!write(6,*)' ingeo : xred_read ='
!write(6,*)xred_read(:,1:natrd)
!ENDDEBUG

 token = 'xangst'
 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),&
& token,txangst,'DPR')
 if(txangst==1) xangst_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 if(nimage/=1 .and. iimage/=1)then
   token = 'xangst_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     xangst_read(:,1:natrd)=xangst_read(:,1:natrd)+&
&     alpha*(reshape( dprarr(1:3*natrd) , (/3,natrd/) )-xangst_read(:,1:natrd))
     txangst=1
   end if
   token = 'xangst'//trim(strimg)
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     xangst_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
     txangst=1
   end if
 end if

 token = 'xcart'
 call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),&
& token,txcart,'LEN')
 if(txcart==1)xcart_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 if(nimage/=1 .and. iimage/=1)then
   token = 'xcart_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     xcart_read(:,1:natrd)=xcart_read(:,1:natrd)+&
&     alpha*(reshape( dprarr(1:3*natrd) , (/3,natrd/) )-xcart_read(:,1:natrd))
     txcart=1
   end if
   token = 'xcart'//trim(strimg)
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     xcart_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
     txcart=1
   end if
 end if

!Might initialize xred from CML or XYZ file
 if(txred+txcart+txangst==0)then
   token = '_xred'
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,txred,'DPR')
   if(txred==1) &
&   xred_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

   token = '_xangst'
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,txangst,'DPR')
   if(txangst==1) &
&   xangst_read(:,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )
 end if

 if (txred+txcart+txangst==0) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' ingeo : ERROR -',ch10,&
&   '  Neither xred nor xangst nor xcart are present in input file. ',ch10,&
&   '  Action : define one of these in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (txred==1)   write(message, '(a)' ) '  xred   is defined in input file'
 if (txangst==1) write(message, '(a)' ) '  xangst is defined in input file'
 if (txcart ==1) write(message, '(a)' ) '  xcart  is defined in input file'
 call wrtout(std_out,message,'COLL')

 if (txred+txcart+txangst>1)then
   write(message, '(6a)' )ch10,&
&   ' ingeo : ERROR -',ch10,&
&   '  Too many input channels for atomic positions are defined.',ch10,&
&   '  Action : choose to define only one of these.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(txred==1)then
   message = ' ingeo : takes atomic coordinates from input array xred '
   call wrtout(std_out,message,'COLL')
   call xredxcart(natrd,1,rprimd,xcart_read,xred_read)
 else
   if(txangst==1)then
     message = ' ingeo : takes atomic coordinates from input array xangst'
     call wrtout(std_out,message,'COLL')
     xcart_read(:,:)=xangst_read(:,:)/Bohr_Ang
   else
     message = ' ingeo : takes atomic coordinates from input array xcart'
     call wrtout(std_out,message,'COLL')
   end if
   txred=1
 end if
!At this stage, the cartesian coordinates are known, for the
!atoms whose coordinates where read.

!DEBUG
!write(6,*)' ingeo : xcart_read ='
!write(6,*)xcart_read
!ENDDEBUG

!Here, allocate the variable that will contain the completed
!sets of xcart, after the use of the geometry builder or the symmetry builder
 allocate(xcart(3,natom))

!7) Eventually read the symmetries

!Take care of the symmetries
 token = 'nsym'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
!Might initialize nsym from CML file
 if(tread==0)then
   token = '_nsym'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 end if

 if(tread==1) nsym=intarr(1)

!Check that nsym is not negative
 if (nsym<0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' ingeo: ERROR -',ch10,&
&   '  Input nsym must be positive or 0, but was ',nsym,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action: correct nsym in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
!Check that nsym is not bigger than msym
 if (nsym>msym) then
   write(message, '(a,a,a,a,i4,a,i4,a,a,a,a,a)' ) ch10,&
&   ' ingeo: ERROR -',ch10,&
&   '  Input nsym=',nsym,' exceeds msym=',msym,'.',ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action: correct nsym in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Read symmorphi
 token = 'symmorphi'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) symmorphi=intarr(1)

!Now, read the symmetry operations
 if(nsym>0)then
   token = 'symrel'
   call intagm(dprarr,intarr,jdtset,marr,9*nsym,string(1:lenstr),token,tread,'INT')
!  Might initialize symrel from CML file
   if(tread==0)then
     token = '_symrel'
     call intagm(dprarr,intarr,jdtset,marr,9*nsym,string(1:lenstr),token,tread,'INT')
   end if
   if(nsym>1 .and. tread==0)then
     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' ingeo : ERROR -',ch10,&
&     '  When nsym>1, symrel must be defined in the input file.',ch10,&
&     '  Action : either change nsym, or define symrel in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if(tread==1) symrel(:,:,1:nsym)=reshape( intarr(1:9*nsym) , (/3,3,nsym/) )

!  Take care of tnons
   tnons(:,1:nsym)=zero
   token = 'tnons'
   call intagm(dprarr,intarr,jdtset,marr,3*nsym,string(1:lenstr),token,tread,'DPR')
!  Might initialize tnons from CML file
   if(tread==0)then
     token = '_tnons'
     call intagm(dprarr,intarr,jdtset,marr,3*nsym,string(1:lenstr),token,tread,'DPR')
   end if
   if(tread==1) tnons(:,1:nsym)=reshape( dprarr(1:3*nsym) , (/3,nsym/) )

   if(symmorphi==0)then
     do isym=1,nsym
       if(sum(tnons(:,isym)**2)>tol6)then
         write(message, '(8a,i2,a,3f8.4,3a)' ) ch10,&
&         ' ingeo : ERROR -',ch10,&
&         '  When symmorph/=1, the vectors of translation (tnons)',ch10,&
&         '  a symmetry operation must vanish.',ch10,&
&         '  However, for the symmetry operation number ',isym,', tnons =',tnons(:,isym),'.',ch10,&
&         '  Action : either change your list of allowed symmetry operations, or use the symmetry finder (nsym=0).'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if

!  Take care of symafm
   token = 'symafm'
   call intagm(dprarr,intarr,jdtset,marr,nsym,string(1:lenstr),token,tread,'INT')
   if(tread==1) symafm(1:nsym)=intarr(1:nsym)

 end if


!8) Checks whether the geometry builder must be used,
!and call it if needed. Call the symmetry builder and analyzer if needed.

!At this stage, nsym might still contain the default 0
!msym contains the default 192.
!The cartesian coordinates of the atoms of the primitive set
!are contained in xcart_read.

 nobj=0
 token = 'nobj'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) nobj=intarr(1)

!If there are objects, chkprim will not be used immediately
!But, if there are no objects, but a space group, it will be used directly.
 chkprim=1
 token = 'chkprim'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) chkprim=intarr(1)

 if(nobj/=0)then

!  Spinat is read for each atom, from 1 to natom
   token = 'spinat'
   call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),token,tread,'DPR')
   if(tread==1) then
     spinat(1:3,1:natom) = reshape( dprarr(1:3*natom) , (/3,natom/) )
   else if (nspden==4.or.(nspden==2.and.nsppol==1)) then
     write(message, '(8a)' )ch10,&
&     ' ingeo: ERROR -',ch10,&
&     '  When nspden=4 or (nspden==2 and nsppol==1), the input variable spinat must be',ch10,&
&     '  defined in the input file, which is apparently not the case.',ch10,&
&     '  Action : define spinat or use nspden=1 in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Will use the geometry builder
   if(tnatrd/=1 .and. nobj/=0)then
     write(message, '(a,a,a,a,a,a,i8,a,a,a,a,a)' ) ch10,&
&     ' ingeo : ERROR -',ch10,&
&     '  The number of atoms to be read (natrd) must be initialized',ch10,&
&     '  in the input file, when nobj=',nobj,'.',ch10,&
&     '  This is not the case.',ch10,&
&     '  Action : initialize natrd in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   if(jellslab/=0)then
     write(message, '(4a,i8,3a)' ) ch10,&
&     ' ingeo : ERROR -',ch10,&
&     '  A jellium slab cannot be used when nobj=',nobj,'.',ch10,&
&     '  Action : change one of the input variables jellslab or nobj in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   call ingeobld (iout,jdtset,lenstr,natrd,natom,&
&   nobj,string,typat,typat_read,xcart,xcart_read)

!  Finalize the computation of coordinates : produce xred.
   call xredxcart(natom,-1,rprimd,xcart,xred)

 else ! nobj==0

!  Spinat is read for each irreducible atom, from 1 to natrd
   token = 'spinat'
   call intagm(dprarr,intarr,jdtset,marr,3*natrd,string(1:lenstr),token,tread,'DPR')
   if(tread==1)spinat(1:3,1:natrd) = reshape( dprarr(1:3*natrd) , (/3,natrd/) )

!  Get xred
   call xredxcart(natrd,-1,rprimd,xcart_read,xred)

   spgroup=0
   token = 'spgroup'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) spgroup=intarr(1)
!  Might initialize spgroup from CML file
   if(tread==0)then
     token = '_spgroup'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
     if(tread==1) spgroup=intarr(1)
   end if

   if(spgroup/=0 .or. nsym/=0)then

     if(jellslab/=0 .and. nsym/=1 .and. spgroup/=1)then
       write(message, '(4a,i8,3a)' ) ch10,&
&       ' ingeo : ERROR -',ch10,&
&       '  For the time being, a jellium slab can only be used either',ch10,&
&       '  either with the symmetry finder (nsym=0) or with the space group 1 (nsym=1)',ch10,&
&       '  Action : change one of the input variables jellslab or nsym or spgroup in your input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     typat(1:natrd)=typat_read(1:natrd)

     if(spgroup/=0 .and. nsym/=0)then
       write(message, '(a,a,a,a,i4,a,a,i4,a,a,a,a,a,a,a,a)' ) ch10,&
&       ' ingeo : ERROR -',ch10,&
&       '  The spatial group number spgroup=',spgroup,ch10,&
&       '  is specified, as well as the number of symmetries nsym=',nsym,ch10,&
&       '  This is not allowed, as you can define the symmetries',ch10,&
&       '  either using spgroup OR using nsym, but not both.',ch10,&
&       '  Action : modify your input file',ch10,&
&       '   (either set spgroup to 0, or nsym to 0)'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     brvltt=0

     if(spgroup/=0)then

!      Will generate the spatial group using spgroup
!      Assign default values
       spgaxor=1
       spgorig=1
       token = 'brvltt'
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
       if(tread==1) brvltt=intarr(1)
       token = 'spgaxor'
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
       if(tread==1) spgaxor=intarr(1)
       token = 'spgorig'
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
       if(tread==1) spgorig=intarr(1)

!      DEBUG
!      write(6,*)' ingeo : brvltt = ',brvltt
!      write(6,*)' ingeo : spgaxor = ',spgaxor
!      write(6,*)' ingeo : spgorig = ',spgorig
!      write(6,*)' ingeo : spgroup = ',spgroup
!      write(6,*) 'ingeo : before symmetry part, msym is msym = ',msym
!      ENDDEBUG

!      Treat the case of magnetic groups
       token = 'spgroupma'
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tspgroupma,'INT')
       if(tspgroupma==1) spgroupma=intarr(1)
       token = 'genafm'
       call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tgenafm,'DPR')
       if(tgenafm==1) genafm(1:3)=dprarr(1:3)
       if(tspgroupma/=0 .and. tgenafm/=0)then
         write(message, '(a,a,a,a,i4,a,a,3es8.2,a,a,a,a,a,a,a,a)' ) ch10,&
&         ' ingeo : ERROR -',ch10,&
&         '  The spatial group number spgroupma=',spgroupma,ch10,&
&         '  is specified, as well as the antiferromagnetic generator genafm=',genafm(1:3),ch10,&
&         '  This is not allowed, as you can define the magnetic space group',ch10,&
&         '  either using spgroupma OR using genafm, but not both.',ch10,&
&         '  Action : modify your input file',ch10,&
&         '   (either define spgroupma or genafm)'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if

!      TODO : all the symmetry generation operations should be in one big routine

!      If spgroupma is defined, check whether it is consistent
!      with spgroup, determine the Shubnikov type,
!      and, for type IV, find the corresponding genafm
       shubnikov=1
       if(tspgroupma==1)then
         call gensymshub(genafm,spgroup,spgroupma,shubnikov)
       else if(tgenafm==1)then
         shubnikov=4
       end if

!      Generate the spatial group of symmetries in a conventional cell
!      In case of Shubnikov space group type IV, only generate the
!      Fedorov (non-magnetic) group. For Shubnikov type III space group,
!      the magnetic part is generated here.
       bckbrvltt=brvltt
       if(brvltt==-1)brvltt=0
       call gensymspgr(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&       spgroupma,symafm,symrel,tnons)

!      For shubnikov type IV groups,
!      double the space group, using the antiferromagnetic translation generator
       if(shubnikov==4)then
         call gensymshub4(genafm,msym,nsym,symafm,symrel,tnons)
       end if

!      DEBUG
!      write(6,*)' after gensymshub4, nsym =',nsym
!      write(6,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!      do ii=1,nsym
!      write(6,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!      end do
!      ENDDEBUG

!      If brvltt was -1 at input, one should now change the conventional cell
!      to a primitive one, if brvltt/=1
       if(bckbrvltt==-1 .and. brvltt/=1)then
!        Will work with rprim only
         rprim(:,:)=rprimd(:,:)
         rprimd_new(:,:)=rprimd(:,:)
         acell(:)=1.0_dp
         select case(brvltt)
           case(5)
             rprimd_new(:,2)=(rprim(:,2)+rprim(:,3))*0.5_dp
             rprimd_new(:,3)=(rprim(:,3)-rprim(:,2))*0.5_dp
           case(6)
             rprimd_new(:,1)=(rprim(:,1)+rprim(:,3))*0.5_dp
             rprimd_new(:,3)=(rprim(:,3)-rprim(:,1))*0.5_dp
           case(4)
             rprimd_new(:,1)=(rprim(:,1)+rprim(:,2))*0.5_dp
             rprimd_new(:,2)=(rprim(:,2)-rprim(:,1))*0.5_dp
           case(3)
             rprimd_new(:,1)=(rprim(:,2)+rprim(:,3))*0.5_dp
             rprimd_new(:,2)=(rprim(:,1)+rprim(:,3))*0.5_dp
             rprimd_new(:,3)=(rprim(:,1)+rprim(:,2))*0.5_dp
           case(2)
             rprimd_new(:,1)=(-rprim(:,1)+rprim(:,2)+rprim(:,3))*0.5_dp
             rprimd_new(:,2)=( rprim(:,1)-rprim(:,2)+rprim(:,3))*0.5_dp
             rprimd_new(:,3)=( rprim(:,1)+rprim(:,2)-rprim(:,3))*0.5_dp
           case(7)
             rprimd_new(:,1)=( rprim(:,1)*2.0_dp+rprim(:,2)+rprim(:,3))/3.0_dp
             rprimd_new(:,2)=(-rprim(:,1)      +rprim(:,2)+rprim(:,3))/3.0_dp
             rprimd_new(:,3)=(-rprim(:,1)-rprim(:,2)*2.0_dp+rprim(:,3))/3.0_dp
         end select
         call symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)
!        Produce xred in the new system of coordinates
         call xredxcart(natrd, 1,rprimd,xcart,xred)
         call xredxcart(natrd,-1,rprimd_new,xcart,xred)
!        Produce tnons in the new system of coordinates
         allocate(tnons_cart(3,nsym))
         call xredxcart(nsym, 1,rprimd,tnons_cart,tnons)
         call xredxcart(nsym,-1,rprimd_new,tnons_cart,tnons)
         deallocate(tnons_cart)

!        DEBUG
!        write(6,*)' after change of coordinates, nsym =',nsym
!        write(6,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!        do ii=1,nsym
!        write(6,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!        end do
!        ENDDEBUG

!        Prune the symmetry operations : suppress those with
!        exactly the same point and magnetic part
         nsym_now=1
         do isym=2,nsym
           irreducible=1
           do jsym=1,nsym_now
             if(sum(abs(symrel(:,:,isym)-symrel(:,:,jsym)))==0 .and. &
&             symafm(isym)==symafm(jsym)                          )then
               irreducible=0
               exit
             end if
           end do
           if(irreducible==1)then
             nsym_now=nsym_now+1
             symrel(:,:,nsym_now)=symrel(:,:,isym)
             tnons(:,nsym_now)=tnons(:,isym)
             symafm(nsym_now)=symafm(isym)
           end if
         end do
         nsym=nsym_now
!        Translate tnons in the ]-0.5,0.5] interval
         tnons(:,1:nsym)=tnons(:,1:nsym)-nint(tnons(:,1:nsym)-1.0d-8)

!        DEBUG
!        write(6,*)' after reduction, nsym =',nsym
!        write(6,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!        do ii=1,nsym
!        write(6,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!        end do
!        ENDDEBUG

!        Now that symrel, tnons and xred are expressed in the primitive
!        axis system, update the geometric quantities
         rprimd(:,:)=rprimd_new(:,:)
         rprim(:,:)=rprimd_new(:,:)
         call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
         call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
       end if

     end if

     if(natom/=natrd)then
!      Generate the full set of atoms from its knowledge in the irreducible part.
       call fillcell(natom,natrd,nsym,spinat,symafm,symrel,tnons,tolsym,typat,xred)
     end if

!    Check whether the symmetry operations are consistent with the lattice vectors
     call chkorthsy(gprimd,iout,nsym,rmet,rprimd,symrel)

   else ! spgroup==0 and nsym==0

!    Here, spgroup==0 as well as nsym==0, so must generate
!    the spatial group of symmetry. However, all the atom
!    positions must be known, so the number
!    of atoms to be read must equal the total number of atoms.
     if(natrd/=natom)then

       write(message, '(a,a,a,a,i4,a,a,i4,a,a,a,a,a,a,a,a,a)' ) ch10,&
&       ' ingeo : ERROR -',ch10,&
&       '  The number of atoms to be read (natrd)=',natrd,ch10,&
&       '  differs from the total number of atoms (natom)=',natom,ch10,&
&       '  while spgroup=0 and nsym=0.',&
&       '  This is not allowed, since the information needed to',ch10,&
&       '  generate the missing atomic coordinates is not available.',ch10,&
&       '  Action : modify your input file',ch10,&
&       '   (either natrd, or natom, or spgroup, or nsym)'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')

     else

       typat(:)=typat_read(:)
!      Find the symmetry operations : nsym, symafm, symrel and tnons.
!      Use nptsym and ptsymrel, as determined by symlatt
       noncoll=0;if (nspden==4) noncoll=1
       use_inversion=1;if (nspden==4.or.pawspnorb>0) use_inversion=0

!      DEBUG
!      write(6,*)' ingeo : tolsym=',tolsym
!      ENDDEBUG

       call symfind(berryopt,efield_xred,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&
&       ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred)

!      If the tolerance on symmetries is bigger than 1.e-8, symmetrize the atomic positions
       if(tolsym>1.00001e-8)then
         allocate(indsym(4,natom,nsym),symrec(3,3,nsym))
         do isym=1,nsym
           call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
         end do
         call symatm(indsym,natom,nsym,symrec,tnons,tolsym,typat,xred)
         call symzat(indsym,natom,nsym,symrel,tnons,xred)
         deallocate(indsym,symrec)

         write(message, '(a,a,a,a,es14.6,10a)' ) ch10,&
&         ' ingeo : WARNING -',ch10,&
&         '  The tolerance on symmetries =',tolsym,ch10,&
&         '  is bigger than the usual tolerance, i.e. 1.0e-8 .',ch10,&
&         '  In order to avoid spurious effect, the atomic coordinates have been',ch10,&
&         '  symmetrized before storing them in the dataset internal variable.',ch10,&
&         '  So, do not be surprised by the fact that your input variables (xcart, xrec, ...)',ch10,&
&         '  do not correspond the the ones echoes by ABINIT, the latter being used to do the calculations.'
         call wrtout(std_out,message,'COLL')

       end if

     end if

   end if

!  Finalize the computation of coordinates : produce xcart
   call xredxcart(natom,1,rprimd,xcart,xred)

 end if ! check of existence of an object

 deallocate(ptsymrel)
 deallocate(xangst_read,xcart_read,xcart,xred_read)
 deallocate(typat_read)

!Correct the default nsym value, if a symmetry group has not been generated.
 if(nsym==0)nsym=1

!--------------------------------------------------------------------------------------------------------

 icoulomb=0
 token = 'icoulomb'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1)icoulomb=intarr(1)

!calculate the center of the atomic system such as to put the
!atoms in the middle of the simulation box for the free BC case.
 if (icoulomb == 1) then
   rcm(:)=zero
   do iatom=1,natom
     rcm(:)=rcm(:)+xred(:,iatom)
   end do 
   rcm(:)=rcm(:)/real(natom,dp)-half
   do iatom=1,natom
     xred(:,iatom)=xred(:,iatom)-rcm(:)
   end do
!  Also modify the tnons
   do isym=1,nsym
     tnons(:,isym)=matmul(symrel(:,:,isym),rcm(:))-rcm(:)+tnons(:,isym)
   end do

   write(message,'(3a)')' ingeo : WARNING -',ch10,&
&   ' Because icoulomb is 1, the average center of coordinates of the system has been translated to (0.5,0.5,0.5) '
   call wrtout(std_out,message,'COLL')
 end if

!========================================================================================================
!
!At this stage, the cell parameters and atomic coordinates are known, as well as the symmetry operations
!There has been a preliminary analysis of the holohedry (not definitive, though ...)
!
!========================================================================================================

!Here, determine correctly the Bravais lattice and other space group or shubnikov group characteristics
 call symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

!--------------------------------------------------------------------------------------

!DEBUG
!write(6,*)' ingeo : before symmorphi filter '
!do isym=1,nsym
!write(6,'(i2,9i3,3f8.4,i3)' )isym,symrel(:,:,isym),tnons(:,isym),symafm(isym)
!end do
!symmorphi=0
!ENDDEBUG

!Finally prune the set of symmetry in case non-symmorphic operations must be excluded
 if(symmorphi==0)then
   jsym=0
   do isym=1,nsym
     if(sum(tnons(:,isym)**2)<tol6)then
       jsym=jsym+1
!      This symmetry operation is non-symmorphic, and can be kept
       if(isym/=jsym)then
         symrel(:,:,jsym)=symrel(:,:,isym)
         tnons(:,jsym)=tnons(:,isym)
         symafm(jsym)=symafm(isym)
       end if
     end if
   end do
   nsym=jsym
 end if

!DEBUG
!write(6,*)' ingeo : after symmorphi filter '
!do isym=1,nsym
!write(6,'(i2,9i3,3f8.4,i3)' )isym,symrel(:,:,isym),tnons(:,isym),symafm(isym)
!end do
!ENDDEBUG

!DEBUG
!call symmultsg(nsym,symafm,symrel,tnons)
!ENDDEBUG

!9) initialize the list of fixed atoms, and initial velocities -----------------
!Note : these inputs do not influence the previous generation of
!symmetry operations. This might be changed in the future

!idir=0 is for iatfix , idir=1 is for iatfixx,
!idir=2 is for iatfixy, idir=3 is for iatfixz
 iatfix(:,:)=0

 do idir=0,3

   if(idir==0)then
     token = 'natfix'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   else if(idir==1)then
     token = 'natfixx'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   else if(idir==2)then
     token = 'natfixy'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   else if(idir==3)then
     token = 'natfixz'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   end if

!  Use natfix also for natfixx,natfixy,natfixz
   natfix=0
   if(tread==1) natfix=intarr(1)


!  Checks the validity of natfix
   if (natfix<0 .or. natfix>natom) then
     write(message, '(a,a,a,a,a,a,i4,a,i4,a,a,a)' ) ch10,&
&     ' ingeo: ERROR -',ch10,&
&     '  The input variables natfix, natfixx, natfixy and natfixz must be',ch10,&
&     '  between 0 and natom (=',natom,'), while one of them is',natfix,'.',ch10,&
&     '  Action: correct that occurence in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Read iatfix 
   if(idir==0)then
     token = 'iatfix'
     call intagm(dprarr,intarr,jdtset,marr,natfix,&
&     string(1:lenstr),token,tread,'INT')
   else if(idir==1)then
     token = 'iatfixx'
     call intagm(dprarr,intarr,jdtset,marr,natfix,&
&     string(1:lenstr),token,tread,'INT')
   else if(idir==2)then
     token = 'iatfixy'
     call intagm(dprarr,intarr,jdtset,marr,natfix,&
&     string(1:lenstr),token,tread,'INT')
   else if(idir==3)then
     token = 'iatfixz'
     call intagm(dprarr,intarr,jdtset,marr,natfix,&
&     string(1:lenstr),token,tread,'INT')
   end if

!  If some iatfix was read, natfix must vanish
   if (natfix==0 .and. tread==1)then
     write(message, '(4a,i1,5a)' ) ch10,&
&     ' ingeo: ERROR -',ch10,&
&     '  For direction ',idir,' the corresponding natfix is zero,',ch10,&
&     '  while iatfix specifies some atoms to be fixed.',ch10,&
&     '  Action: either specify a non-zero natfix(x,y,z) or suppress iatfix(x,y,z).'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  If natfix is non-zero, iatfix must be defined
   if (natfix>0 .and. tread==0)then
     write(message, '(4a,i1,3a,i4,3a)' ) ch10,&
&     ' ingeo: ERROR -',ch10,&
&     '  For direction ',idir,' no iatfix has been specified,',ch10,&
&     '  while natfix specifies that some atoms to be fixed, natfix=',natfix,'.',ch10,&
&     '  Action: either set natfix(x,y,z) to zero or define iatfix(x,y,z).'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   if(tread==1)then
     do ii=1,natfix
!      Checks the validity of the input iatfix
       if (intarr(ii)<1 .or. intarr(ii)>natom) then
         write(message, '(a,a,a,a,a,a,i4,a,a,a)' ) ch10,&
&         ' ingeo: ERROR -',ch10,&
&         '  The input variables iatfix, iatfixx, iatfixy and iatfixz must be',&
&         ch10,&
&         '  between 1 and natom, while one of them is',intarr(ii),'.',ch10,&
&         '  Action: correct that occurence in your input file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
!      Finally set the value of the internal iatfix array
       do iatom=1,natom
         if(intarr(ii)==iatom)then
           if(idir==0)iatfix(1:3,iatom)=1
           if(idir/=0)iatfix(idir,iatom)=1
         end if
       end do
     end do
   end if

!  DEBUG
!  write(6,*)' ingeo  : idir ',idir,'tread ',tread,'natfix',natfix,'iatfix ',iatfix
!  if(idir==3)stop
!  ENDDEBUG

 end do

 token = 'vel'
 vel(:,:)=zero
 call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),token,tread,'DPR')
 if(tread==1)vel(:,:)=reshape( dprarr(1:3*natom) , (/3,natom/) )

 if(nimage/=1 .and. iimage/=1)then
   token = 'vel_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),token,tread,'DPR')
   if(tread==1)vel(:,:)=vel(:,:)&
&   +alpha*(reshape( dprarr(1:3*natom) , (/3,natom/) )-vel(:,:))
   token = 'vel'//trim(strimg)
   call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),token,tread,'DPR')
   if(tread==1)vel(:,:) = reshape( dprarr(1:3*natom) , (/3,natom/) )
 end if

 deallocate(intarr,dprarr)

!DEBUG
!write(6,*)' ingeo : exit '
!do isym=1,nsym
!write(6,'(i2,9i3,3f8.4,i3)' )isym,symrel(:,:,isym),tnons(:,isym),symafm(isym)
!end do
!ENDDEBUG

!DEBUG
!write(6,*)' ingeo : end of subroutine, nsym= ',nsym
!write(6,*)' ingeo : ptgroupma=',ptgroupma
!write(6,*)' ingeo : natfix=',natfix
!write(6,*)' ingeo : iatfix=',iatfix
!if(.true.)stop
!ENDDEBUG

!call timab(47,2,tsec)

end subroutine ingeo
!!***
