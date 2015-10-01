!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_ab6_invars_f90
!! NAME
!! m_ab6_invars_f90
!!
!! FUNCTION
!! driver for the parser
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_ab6_invars

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use interfaces_57_iovars

  implicit none

  private

  ! We store here a list of dtset arrays to be able to
  ! parse several ABINIT files without freeing it.
  ! The simplest portable way to do it, is to create
  ! a list of dtsets arrays and to use the list index
  ! as an identifier that can be given to the other languages.
  type, private :: dtsets_list
     integer                     :: id
     type(dtsets_list),  pointer :: next => null()
     type(dataset_type), pointer :: dtsets(:)
     type(pspheader_type), pointer :: pspheads(:)
     integer :: mxgw_nqlwl, mxnatom, mxntypat, mxlpawu, mxmband_upper, mxnatpawu, &
         & mxnatsph,  mxnconeq, mxnimage, mxnkptgw, mxnatvshift, &
         & mxnkpt,  mxnnos, mxnqptdm, mxnsppol, mxnsym, mxnspinor, mxmband
     integer :: istatr, istatshft, dmatpuflag, papiopt, timopt
  end type dtsets_list
  type(dtsets_list), pointer :: my_dtsets => null()
  integer :: nb_dtsets = 0

  ! These flags should be .true. inside ABINIT.
  ! Use ab6_invars_set_flags() to change them.
  logical, private :: call_status = .false.
  character(len = fnlen), private :: opt_status_file
  logical, private :: call_timab = .false.
  real(dp), pointer, private :: opt_timab_tsec(:)

  ! These pointers are used only inside ABINIT and refers to their
  ! equivalent. The C binding don't support them.
  ! Use ab6_invars_get_abinit_vars() to get them.
  type(MPI_type), pointer :: mpi_enreg

  logical, private, parameter :: AB_DBG = .false.

  include "ab6_invars_f90.inc"

  ! The following group is used for Fortran bindings only,
  ! and specifically its usage inside ABINIT. They have no
  ! C or Python equivalent.
  public :: ab6_invars_set_flags
  public :: ab6_invars_set_mpi
  public :: ab6_invars_get_abinit_vars
  public :: ab6_invars_load

  ! The following routines are the main creation routines, having also
  ! an equivalent in C or Python.
  public :: ab6_invars_new_from_file
  public :: ab6_invars_new_from_string
  public :: ab6_invars_free

  ! The following routines are the main getter functions, also available
  ! in C or Python.
  public :: ab6_invars_get_ndtset
  public :: ab6_invars_get_integer
  public :: ab6_invars_get_real
  public :: ab6_invars_get_shape
  public :: ab6_invars_get_integer_array
  public :: ab6_invars_get_real_array

contains

  subroutine ab6_invars_set_flags(status, timab, status_file, timab_tsec)


    logical, intent(in) :: status, timab
    character(len = fnlen), intent(in), optional :: status_file
    real(dp), intent(in), target, optional :: timab_tsec(:)

    call_status = status
    if (present(status_file)) then
       write(opt_status_file, "(A)") status_file
    else
       write(opt_status_file, "(A)") "status"
    end if
    call_timab  = timab
    if (present(timab_tsec)) then
       opt_timab_tsec => timab_tsec
    else
       allocate(opt_timab_tsec(2))
    end if
  end subroutine ab6_invars_set_flags
  
  subroutine ab6_invars_set_mpi(mpi_enreg_)


    type(MPI_type), intent(in), target :: mpi_enreg_

    mpi_enreg => mpi_enreg_
  end subroutine ab6_invars_set_mpi

  subroutine ab6_invars_get_abinit_vars(dtsetsId, dtsets, pspheads, &
       & mxgw_nqlwl, mxlpawu, mxmband, mxmband_upper, mxnatom, mxnatpawu, &
       & mxnatsph,  mxnatvshift, mxnconeq, mxnimage, mxnkpt, mxnkptgw, &
       & mxnnos, mxnqptdm, mxnsppol, mxnspinor, mxnsym, mxntypat, &
       & istatr, istatshft, papiopt, timopt, dmatpuflag)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: dtsetsId
    type(dataset_type), pointer :: dtsets(:)
    type(pspheader_type), pointer :: pspheads(:)
    integer, intent(out) :: mxgw_nqlwl, mxnatom, mxntypat, mxlpawu, &
         & mxmband_upper, mxnatpawu, mxnatsph,  mxnconeq, mxnimage, mxnkptgw, mxnatvshift, &
         & mxnkpt,  mxnnos, mxnqptdm, mxnsppol, mxnsym, mxnspinor,mxmband, &
         & istatr, istatshft, papiopt, timopt, dmatpuflag

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)

    if (associated(token)) then
       dtsets => token%dtsets
       pspheads => token%pspheads
       mxgw_nqlwl = token%mxgw_nqlwl
       mxnatom = token%mxnatom
       mxnimage = token%mxnimage
       mxntypat = token%mxntypat
       mxlpawu = token%mxlpawu
       mxmband_upper = token%mxmband_upper
       mxnatpawu = token%mxnatpawu
       mxnatsph = token%mxnatsph
       mxnconeq = token%mxnconeq
       mxnkptgw = token%mxnkptgw
       mxnatvshift = token%mxnatvshift
       mxnkpt = token%mxnkpt
       mxnnos = token%mxnnos
       mxnqptdm = token%mxnqptdm
       mxnsppol = token%mxnsppol
       mxnsym = token%mxnsym
       mxnspinor = token%mxnspinor
       mxmband = token%mxmband
       istatr = token%istatr
       istatshft = token%istatshft
       papiopt = token%papiopt
       timopt = token%timopt
       dmatpuflag = token%dmatpuflag
    else
       nullify(dtsets)
       nullify(pspheads)
    end if
  end subroutine ab6_invars_get_abinit_vars
  
  subroutine new_token(token)


    type(dtsets_list), pointer :: token

    ! We allocate a new list token and prepend it.
    if (AB_DBG) write(0,*) "AB module: create a new token."
    nb_dtsets = nb_dtsets + 1

    allocate(token)
    token%id = nb_dtsets
    nullify(token%dtsets)
    token%next => my_dtsets

    my_dtsets => token
    if (AB_DBG) write(0,*) "AB module: creation OK with id ", token%id
  end subroutine new_token

  subroutine free_token(token)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_53_abiutil
!End of the abilint section

    type(dtsets_list), pointer :: token

    integer :: idtset

    if (.not. associated(token)) then
       return
    end if

    ! We free a token list.
    if (AB_DBG) write(0,*) "AB module: free request on dataset array ", token%id
    if (associated(token%dtsets)) then
       if (AB_DBG) write(0,*) " | ", size(token%dtsets), "dtsets found."
       do idtset = 0, size(token%dtsets) - 1, 1
          if (AB_DBG) write(0,*) " | free dtset ", idtset
          call dtsetFree(token%dtsets(idtset))
          if (AB_DBG) write(0,*) " | free OK"
       end do
       deallocate(token%dtsets)
       nullify(token%dtsets)
       deallocate(token%pspheads)
       nullify(token%pspheads)
       if (AB_DBG) write(0,*) " | general free OK"

       ! For now, we let the list grow.
!!$       ! We remove token from the list.
!!$       call _dtset_get_from_id(prev, dtsetsId - 1)
!!$       if (associated(prev)) then
!!$          prev%next => token%next
!!$       else
!!$          my_dtsets => token%next
!!$       end if
!!$       deallocate(token)
!!$       nb_dtsets = max(nb_dtsets - 1, 0)
    end if
    if (AB_DBG) write(0,*) "AB module: free done"
  end subroutine free_token

  subroutine get_token(token, id)


    type(dtsets_list), pointer :: token
    integer, intent(in) :: id

    type(dtsets_list), pointer :: tmpLst

    if (AB_DBG) write(0,*) "AB module: request list element ", id
    nullify(token)
    ! List element are prepended so element id is at (nb - id) position.
    tmpLst => my_dtsets
    do
       if (.not. associated(tmpLst)) then
          exit
       end if
       if (tmpLst%id == id .and. associated(tmpLst%dtsets)) then
          token => tmpLst
          return
       end if
       tmpLst => tmpLst%next
    end do
  end subroutine get_token

  subroutine ab6_invars_new_from_string(dtsetsId, instr, len)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_47_xml
 use interfaces_57_iovars
!End of the abilint section

    integer, intent(out) :: dtsetsId
    integer, intent(in) :: len
    character(len = len), intent(in) :: instr

    character(len = strlen) :: string
    integer :: lenstr, ndtset
    integer :: marr, tread
    character(len = 30) :: token
    integer :: intarr(1)
    real(dp) :: dprarr(1)
    character(len=500) :: message

    dtsetsId = 0

    if (len > strlen) then
       return
    end if

    write(string,*) instr

    !To make case-insensitive, map characters of string to upper case:
    call inupper(string(1:len))

    !Might import data from CML file(s) into string
    !Need string_raw to deal properly with CML filenames
    lenstr = len
    call importcml(lenstr, instr, string, len)
    call importxyz(lenstr, instr, string, len)

    !6) Take ndtset from the input string, then allocate
    !the arrays whose dimensions depends only on ndtset and msym.

    ndtset=0 ; marr=1
    token = 'ndtset'
    call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
    if(tread==1) ndtset=intarr(1)
    !Check that ndtset is not negative
    if (ndtset<0 .or. ndtset>99) then
       write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
            &  ' abinit : ERROR -',ch10,&
            &  '  Input ndtset must be non-negative and < 100, but was ',ndtset,ch10,&
            &  '  This is not allowed.  ',ch10,&
            &  '  Action : modify ndtset in the input file.'
       call wrtout(06,  message,'COLL')
       call leave_new('COLL')
    end if

    allocate(mpi_enreg)
    mpi_enreg%nproc = 1
    mpi_enreg%me = 0
    mpi_enreg%paral_compil_mpio = 0
    call ab6_invars_load(dtsetsId, string, lenstr, ndtset, .false., .false.)
    deallocate(mpi_enreg)
  end subroutine ab6_invars_new_from_string

  subroutine ab6_invars_new_from_file(dtsetsId, filename, n)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_57_iovars
!End of the abilint section

    integer, intent(out) :: dtsetsId
    integer, intent(in) :: n
    character(len = n), intent(in) :: filename

    character(len = strlen) :: string
    integer :: lenstr, ndtset

    dtsetsId = 0

    if (call_status) call status(0,opt_status_file,99,1,'call parsefile')

    if (AB_DBG) write(0,*) "AB module: read '", trim(filename), "' to string."
    call parsefile(filename, lenstr, ndtset, string)
    if (AB_DBG) write(0,*) "AB module: read OK, string length ", lenstr

    allocate(mpi_enreg)
    mpi_enreg%nproc = 1
    mpi_enreg%me = 0
    mpi_enreg%paral_compil = 0
    mpi_enreg%paral_compil_mpio = 0
    call ab6_invars_load(dtsetsId, string, lenstr, ndtset, .false., .false.)
    deallocate(mpi_enreg)
  end subroutine ab6_invars_new_from_file

  subroutine ab6_invars_load(dtsetsId, string, lenstr, ndtset, &
       & with_psp, with_mem, pspfilnam)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_45_psp_parser
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_57_iovars
!End of the abilint section

    integer, intent(out) :: dtsetsId
    character(len = strlen), intent(inout) :: string
    integer, intent(in) :: lenstr, ndtset
    logical, intent(in) :: with_psp, with_mem
    character(len = fnlen), intent(in), optional :: pspfilnam(:)

    type(dtsets_list), pointer :: token
    integer :: jdtset
    integer :: ndtset_alloc
    integer :: npsp, ii, idtset, msym, usepaw
    integer,allocatable :: bravais_(:,:),mband_upper_(:)
    real(dp) :: zion_max
    character(len = fnlen), allocatable :: pspfilnam_(:)

    ! We allocate a new list token and prepend it.
    if (AB_DBG) write(0,*) "AB module: allocate a new object."
    call new_token(token)
    dtsetsId = token%id

    ndtset_alloc=ndtset ; if(ndtset==0)ndtset_alloc=1
    allocate(token%dtsets(0:ndtset_alloc))
    do idtset = 0, size(token%dtsets) - 1, 1
       call dtset_nullify(token%dtsets(idtset))
    end do
    if (AB_DBG) write(0,*) "AB module: allocation OK at ", dtsetsId

    if (AB_DBG) write(0,*) "AB module: call invars0()."
    token%timopt = 1
    if(mpi_enreg%paral_compil==1) token%timopt = 0
    !7) Continue to analyze the input string, get upper dimensions,
    !and allocate the remaining arrays.
    if (call_status) call status(0,opt_status_file,99,1,'call invars0  ')
    call invars0(token%dtsets,token%istatr,token%istatshft,lenstr,&
         & msym,token%mxnatom,token%mxnimage,token%mxntypat,ndtset,ndtset_alloc,npsp,&
         & token%papiopt, token%timopt, string)

    !Be careful : at these fourth and fifth calls of status, istatr and istatshft taken
    !from the input variables will be saved definitively.
    if (call_status) call status(0,opt_status_file,token%istatr,1,'init istatr   ')
    if (call_status) call status(0,opt_status_file,token%istatshft,1,'init istatshft')

    if (call_timab) call timab(41,2,opt_timab_tsec)
    if (call_timab) call timab(token%timopt,5,opt_timab_tsec)

    !8) Finish to read the "file" file completely, as npsp is known,
    !and also initialize pspheads, that contains the important information
    !from the pseudopotential headers, as well as the psp filename

    if (call_status) call status(0,opt_status_file,99,1,'call iofn2    ')

    if (call_timab) call timab(42,1,opt_timab_tsec)

    usepaw=0
    allocate(token%pspheads(npsp))
    if (with_psp) then
       if (mpi_enreg%me == 0) then
          if (.not. present(pspfilnam)) then
             allocate(pspfilnam_(npsp))
             call iofn2(npsp, pspfilnam_)
             call inpspheads(pspfilnam_,npsp,token%pspheads)

!DEBUG
!      write(6,*)' ab6_invars_f90 : token%pspheads(1)%nproj(0:3)=',token%pspheads(1)%nproj(0:3)
!ENDDEBUG

             deallocate(pspfilnam_)
          else
             call inpspheads(pspfilnam,npsp,token%pspheads)
          end if
          if(minval(abs(token%pspheads(1:npsp)%pspcod-7))==0) usepaw=1
       end if
       !Communicate pspheads to all processors
       call pspheads_comm(npsp,mpi_enreg,token%pspheads,usepaw)
    else
       ! No psp files are given, we put default values into pspheads.
       token%pspheads(:)%zionpsp = 1
       token%pspheads(:)%pspxc   = token%dtsets(1)%ixc
       token%pspheads(:)%pspso   = 0
       token%pspheads(:)%xccc    = 0
    end if
    !If (all) pspcod are 7 then this is a PAW calculation
    do idtset=0,ndtset_alloc
       token%dtsets(idtset)%usepaw=usepaw
    end do

    !Take care of other dimensions, and part of the content of dtsets
    !that is or might be needed early.
    !zion_max=maxval(pspheads(1:npsp)%zionpsp) ! This might not work properly with HP compiler
    zion_max=token%pspheads(1)%zionpsp
    do ii=1,npsp
       zion_max=max(token%pspheads(ii)%zionpsp,zion_max)
    end do
    if (AB_DBG) write(0,*) "AB module: OK."

    allocate(bravais_(11,0:ndtset_alloc))
    allocate(mband_upper_ (  0:ndtset_alloc))

    if (AB_DBG) write(0,*) "AB module: call invars1m()."

!DEBUG
!      write(6,*)' ab6_invars_f90 , before invars1m : token%pspheads(1)%nproj(0:3)=',token%pspheads(1)%nproj(0:3)
!ENDDEBUG
    call invars1m(bravais_,token%dmatpuflag,token%dtsets,ab_out,lenstr,mband_upper_,&
         & mpi_enreg,msym,token%mxgw_nqlwl,token%mxlpawu,token%mxmband_upper,&
         & token%mxnatom,token%mxnatpawu,token%mxnatsph,token%mxnatvshift,&
         & token%mxnconeq,token%mxnkpt,token%mxnkptgw,token%mxnimage,token%mxnnos,token%mxnqptdm,&
         & token%mxnspinor,token%mxnsppol,token%mxnsym,ndtset,ndtset_alloc,string,&
         & zion_max)

!DEBUG
!      write(6,*)' ab6_invars_f90 , after : token%pspheads(1)%nproj(0:3)=',token%pspheads(1)%nproj(0:3)
!      stop
!ENDDEBUG

    if (AB_DBG) write(0,*) "AB module: OK."

    !9) Provide defaults for the variables that have not yet been initialized.
    if (AB_DBG) write(0,*) "AB module: call indefo()."
    if (call_status) call status(0,opt_status_file,99,1,'call indefo   ')

    call indefo(token%dtsets,ndtset_alloc)
    if (AB_DBG) write(0,*) "AB module: OK."
   
    if (call_status) call status(0,opt_status_file,99,1,'call macroin  ')

    call macroin(token%dtsets,ndtset_alloc) 

    !10) Perform some global initialization, depending on the value of
    ! pseudopotentials, parallelism variables, or macro input variables
    
    !If all the pseudopotentials have the same pspxc, override the default
    !value for dtsets 1 to ndtset
    if(with_psp .and. minval(abs((token%pspheads(1:npsp)%pspxc-token%pspheads(1)%pspxc)))==0)then
       token%dtsets(1:ndtset_alloc)%ixc=token%pspheads(1)%pspxc
    end if

    !Bands-FFT parallelism is activated, override the default for particular datasets, imposing MPI-IO
    if (mpi_enreg%paral_compil_mpio==1) then
       do idtset=1,ndtset_alloc
          if (token%dtsets(idtset)%paral_kgb==1) token%dtsets(idtset)%accesswff=IO_MODE_MPI
       end do
    end if

    if (call_timab) call timab(42,2,opt_timab_tsec)

    !11) Call the main input routine.
    if (AB_DBG) write(0,*) "AB module: call invars2()."

    if (call_status) call status(0,opt_status_file,99,1,'call invars2m ')
    if (call_timab) call timab(43,1,opt_timab_tsec)

    if (with_mem) then

!DEBUG
!      write(6,*)' ab6_invars_f90 : token%pspheads(1)%nproj(0:3)=',token%pspheads(1)%nproj(0:3)
!ENDDEBUG

       call invars2m(bravais_,token%dtsets,ab_out,lenstr,&
            & mband_upper_,mpi_enreg,msym,ndtset,ndtset_alloc,npsp,token%pspheads,string)
    else
       do idtset = 1, ndtset_alloc, 1
          jdtset=token%dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
          call invars2(bravais_(:, idtset),token%dtsets(idtset),ab_out,jdtset,lenstr,&
               & mband_upper_(idtset),msym,npsp,token%pspheads,string,usepaw,&
               & token%pspheads(1:npsp)%zionpsp)
       end do
    end if
    deallocate(bravais_, mband_upper_)
    if (AB_DBG) write(0,*) "AB module: OK."

    if (call_status) call status(0,opt_status_file,99,1,'call macroin2  ')

    call macroin2(token%dtsets,ndtset_alloc)

    !mxmband=maxval(dtsets(1:ndtset_alloc)%mband) ! This might not work with the HP compiler
    token%mxmband=token%dtsets(1)%mband
    do ii=1,ndtset_alloc
       token%mxmband=max(token%dtsets(ii)%mband,token%mxmband)
    end do

    if (call_timab) call timab(43,2,opt_timab_tsec)

  end subroutine ab6_invars_load

  subroutine ab6_invars_free(dtsetsId)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: dtsetsId

    type(dtsets_list), pointer :: token

    nullify(token)
    call get_token(token, dtsetsId)
    call free_token(token)
  end subroutine ab6_invars_free

  subroutine ab6_invars_get_ndtset(dtsetsId, value, errno)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

    integer, intent(in) :: dtsetsId
    integer, intent(out) :: value
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)
    if (associated(token)) then
       value = size(token%dtsets) - 1
       errno = AB6_NO_ERROR
    else
       errno = AB6_ERROR_OBJ
    end if
  end subroutine ab6_invars_get_ndtset

  include "ab6_invars_f90_get.f90"
end module m_ab6_invars
!!***
