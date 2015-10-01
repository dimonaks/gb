!Simple Fortran subroutine to read input tokens
 subroutine readin(filnamin, &
 &acell, rprim, xred, typat, znucl, &
 &periodcell,firotatez,&
 &hkl1, uvtw1, uvtw2, uvtw3, uvw1, uvw2, uvw3, mul, test,&
 &slicepms, shiftatoms, deldist, shiftcell, shiftgrain, a_c_conv, sixvalues)

 use defs_basis
 implicit none
 
 
 
!for parser
 integer,parameter :: max_at=11000
 character(len = fnlen), intent(inout) :: filnamin

 character(len = strlen) :: string
 integer :: lenstr
 integer :: option, marr, tread,narr,trprim,txred
!  character(len = strlen) :: string_raw
 character(len = 30) :: token

 
 integer :: jdtset,ds_input
 integer intarr(max_at*3)
 real(dp) dprarr(max_at*3)
 real(dp) :: dprarr1(1)
 integer :: intarr1(1)
!For input
 real(dp),intent(out) :: acell(3),firotatez,deldist,shiftcell(2), shiftgrain(2), a_c_conv(4), sixvalues(6)

 real(dp),intent(out) :: xred(3,max_at),rprim(3,3), mul(3), slicepms(3), shiftatoms(4),  znucl(100)
! real(dp) :: xred(3,max_at)
 integer,intent(out) :: typat(max_at), periodcell(6)
 integer,intent(out) :: hkl1(3), uvtw1(4), uvtw2(4), uvtw3(4), uvw1(3), uvw2(3), uvw3(3)
 integer,intent(out) :: test
 integer :: natom, ntypat

 jdtset=0
 option=1
 !filnamin = trim(filnamin)

 call instrng (filnamin,lenstr,option,strlen,string)
! write(std_out,'(3a)') ch10,' Here=',ch10
!To make case-insensitive, map characters of string to upper case:

 call inupper(string(1:lenstr))


!! INPUTS
!!  jdtset=see the notes section
!!  marr=dimension of the intarr and dprarr arrays, as declared in the
!!   calling subroutine.
!!  narr=actual size of array to be read in.
!!  string=character string containing 'tags' and data.
!!  typevarphys= variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR'=>real(dp) (no special treatment)
!!   'LEN'=>real(dp) (expect a "length", identify bohr, au or angstrom,
!!       and return in au -atomic units=bohr- )
!!   'ENE'=>real(dp) (expect a "energy", identify Ha, hartree, eV, Ry, Rydberg)
!!   'LOG'=>integer, but read logical variable T,F,.true., or .false.
!!   'KEY'=>character, returned in token


 token = 'natom'
 marr=1; narr=1
 call intagm(dprarr1,intarr1,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread.ne.1)  stop ' readin : no natom!'
 natom = intarr1(1)

 marr = natom*3
 !allocate(dprarr(marr),intarr(marr))

! acell(1:3)=one
 token = 'acell'; narr=3
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'LEN',ds_input)
 if(tread==1) acell(1:narr) = dprarr(1:narr)
 

 
 
! rprim(:,:)=0.0_dp; rprim(1,1)=1.0_dp; rprim(2,2)=1.0_dp; rprim(3,3)=1.0_dp
 token = 'rprim'
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),token,trprim,'DPR',ds_input)
 if(trprim==1)rprim(:,:)=reshape( dprarr(1:9) , (/3,3/) ) 
 
 token = 'typat'; narr=natom
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) typat(1:narr) = intarr(1:narr)
 
 ntypat = maxval(typat)


! allocate(xred(3,natom))
 token = 'xred'; 
 call intagm(dprarr,intarr,jdtset,marr,3*natom,string(1:lenstr),token,txred,'DPR',ds_input)
 if(txred==1) xred(:,1:natom) = reshape( dprarr(1:3*natom) , (/3,natom/) )

 marr = 100
 token = 'znucl'; narr = ntypat
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) znucl(1:narr) = dprarr(1:narr)


 marr = 6
 token = 'periodcell'; narr=6
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) periodcell(1:narr) = intarr(1:narr)

 token = 'firotatez'; narr=1
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
! if(tread.ne.1) stop 'Error There are no periodcell variable in input file ' 
 if(tread==1)firotatez = dprarr(1)

 
 !Boundary plane and Three directions in 4-index notation

 token = 'hkl1'; narr=3
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) hkl1(1:narr) = intarr(1:narr)
 
 uvw1 = 0
 uvw2 = 0
 uvw3 = 0
 uvtw1 = 0
 uvtw2 = 0
 uvtw3 = 0
 
 token = 'uvtw1'; narr=4
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) uvtw1(1:narr) = intarr(1:narr)
 
 token = 'uvtw2'; narr=4
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) uvtw2(1:narr) = intarr(1:narr)
 
 token = 'uvtw3'; narr=4
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) uvtw3(1:narr) = intarr(1:narr)
 
 token = 'uvw1'; narr=3
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) uvw1(1:narr) = intarr(1:narr)
 
 token = 'uvw2'; narr=3
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) uvw2(1:narr) = intarr(1:narr)
 
 token = 'uvw3'; narr=3
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) uvw3(1:narr) = intarr(1:narr)
 




!The multiplicators, by which new vectors parrallel to the directions will be multiplied.
 token = 'mul'; narr=3
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) mul(1:narr) = dprarr(1:narr)

! slice parameters 
  token = 'slicepms'; narr=3
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) slicepms(1:narr) = dprarr(1:narr)
 

! shift of atoms parameters
 shiftatoms = 0
 token = 'shiftatoms'; narr=4
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) shiftatoms(1:narr) = dprarr(1:narr)

 
 token = 'test'; narr=1
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread==1) test = intarr(1)

 deldist = 3
 token = 'deldist'; narr=1
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) deldist = dprarr(1)
 
 shiftcell(1) = 1
 shiftcell(2) = 0
 token = 'shiftcell'; narr=2
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) shiftcell(1:narr) = dprarr(1:narr)
 
 shiftgrain = 0
 token = 'shiftgrain'; narr=2
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) shiftgrain(1:narr) = dprarr(1:narr)
!write(6,*) shiftgrain(1:narr)


!Parameters for generation of cells for manual search of lattice constants
 a_c_conv = 0
 token = 'a_c_conv'; narr=4
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) a_c_conv(1:narr) = dprarr(1:narr)
 !write(6,*) a_c_conv(1:narr)


!parameter for different use
 sixvalues = 0
 token = 'sixvalues'; narr=6
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'DPR',ds_input)
 if(tread==1) sixvalues(1:narr) = dprarr(1:narr)




 !Always last
 token = 'name'; narr=1
 call intagm(dprarr,intarr,jdtset,marr,narr,string(1:lenstr),token,tread,'KEY',ds_input)
 filnamin = token
 filnamin = trim(filnamin)
 write(6,*) token
! deallocate(dprarr,intarr)
! deallocate(xred)
 end subroutine readin







 subroutine readinnatom(filnamin,&
&natom)
 use defs_basis
 implicit none
 
!for parser
 character(len = fnlen),intent(in):: filnamin

 character(len = strlen) :: string
 integer :: lenstr,ds_input
 integer :: option, marr, tread,narr,trprim
!  character(len = strlen) :: string_raw
 character(len = 30) :: token

 
 integer :: jdtset
 real(dp) :: dprarr1(1)
 integer :: intarr1(1)
!For input

 integer,intent(out) :: natom
 jdtset=0
 option=1

 call instrng (filnamin,lenstr,option,strlen,string)

 call inupper(string(1:lenstr))


 token = 'natom'
 marr=1; narr=1
 call intagm(dprarr1,intarr1,jdtset,marr,narr,string(1:lenstr),token,tread,'INT',ds_input)
 if(tread.ne.1)  stop ' readin : no natom!'
 natom = intarr1(1)
!write(6,*)' readinnatom : natom_read ='
 end subroutine readinnatom






