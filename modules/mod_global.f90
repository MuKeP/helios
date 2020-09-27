#   if defined(__INTEL_COMPILER)
#       define __COMPILER 1
#   elif defined(__GFORTRAN__)
#       define __COMPILER 2
#       define __unix 1
#   elif defined(__SUNPRO_F90)
#       define __COMPILER 3
#   elif defined(__SUNPRO_F95)
#       define __COMPILER 3
#   else
#       define __COMPILER 4
#   endif

#   if defined(__unix)
#       define __OS 2
#   else
#       define __OS 1
#   endif

#   if defined(_OPENMP)
#       define __opnmp 1
#   else
#       define __opnmp 0
#   endif

    module glob

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

#   if (__COMPILER==1)
        use ifport,  only: signal,system,getpid,getuid,splitpathqq,hostnam
        use ifport,  only: SIGABRT,SIGINT,SIGTERM
        use ifposix, only: PXFSTRUCTCREATE
#       if (__OS==2)
            use ifposix, only: PXFSIGADDSET,PXFCONST,PXFSYSCONF,PXFGETPPID
            use ifposix, only: PXFGETPWUID,PXFSTRGET
#       endif
#   elif (__COMPILER==2)
        !
#   elif (__COMPILER==3)
        !
#   elif (__COMPILER==4)
        use DFPort, only: signal,system,getpid
        use DFLib,  only: splitpathqq
        use DFPort, only: SIGABRT,SIGINT,SIGTERM
#   endif

#   if(__opnmp==1)
        include "omp_lib.h"
!        use omp_lib
#   endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: glVersion='4.541'
    character (len=*), parameter :: glDate   ='2018.12.17'
    character (len=*), parameter :: glAuthor ='Anton B. Zakharov'

    integer*4, parameter :: r16kind=16, r8kind=8, r4kind=4
    integer*4, parameter :: i16kind=16, i8kind=8, i4kind=4, i2kind=2, i1kind=1
    integer*4, parameter ::             l8kind=8, l4kind=4, l2kind=2, l1kind=1

    public :: r16kind, r8kind, r4kind
    public :: i16kind, i8kind, i4kind, i2kind, i1kind
    public ::          l8kind, l4kind, l2kind, l1kind

    ! *glu => for global use  (storage and non-accuracy-demanding procedures)
    ! *spu => for special use (accuracy-demanding procedures)
    integer(kind=i4kind), parameter :: rglu=r8kind, rspu=r16kind
    integer(kind=i4kind), parameter :: iglu=i4kind, ispu=i8kind
    integer(kind=i4kind), parameter :: lglu=l1kind

    !   ~~~~ Global data settings ~~~~ !
#   if(__OS==1)
        character (len=*), parameter :: os='win',osSeparator=char(92),osMove='move',osCopy='copy'
        logical*1        , parameter :: isPosix=.FALSE.
#   else
        character (len=*), parameter :: os='nix',osSeparator=char(47),osMove='mv'  ,osCopy='cp'
        logical*1        , parameter :: isPosix=.TRUE.
#   endif

    real(kind=rglu), parameter :: gluZero=0,gluUnity=1
    real(kind=rspu), parameter :: spuZero=0,spuUnity=1

    !tolerance in comparison of two real variables
    real(kind=rglu), parameter :: gluCompare=10._rglu**floor(log10(epsilon(gluUnity))+1._rglu)
    real(kind=rspu), parameter :: spuCompare=10._rspu**floor(log10(epsilon(spuUnity))+1._rspu)

    logical(kind=lglu), parameter :: true=.true., false=.false.
    real   (kind=rglu), parameter :: NaN=0./0. !keep carefully since some compilers may show different behavior
    character  (len=*), parameter :: months(12)=['Jan','Feb','Mar','Apr','May','Jun',&
                                                 'Jul','Aug','Sep','Oct','Nov','Dec']
    integer(kind=iglu), parameter :: monthDays(12)=[31,28,31,30,31,30,31,31,30,31,30,31]

    character (len=*) , parameter :: days(7)   =['Mon','Tue','Wed','Thu','Fri','Sat','Sun']

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    abstract interface
        subroutine outOfMemory(label, req, left)
        character (len=*), intent(in) :: label
        integer(kind=8)  , intent(in) :: req,left
        end subroutine outOfMemory
    end interface

    abstract interface
        subroutine changeMemory(label, size)
        character (len=*), intent(in) :: label
        integer(kind=8)  , intent(in) :: size
        end subroutine changeMemory
    end interface

    type uch
        character (len=1), allocatable :: ch(:)
        integer(kind=iglu)             :: ln

        contains

        procedure :: del=>uchDel
        procedure :: get=>uchGet
    end type uch

    type ivarVector
        integer(kind=iglu), allocatable :: v(:)
        integer(kind=iglu)              :: ln

        contains

        procedure :: del=>ivarDel
        procedure :: get=>ivarGet
    end type ivarVector

    type rvarVector
        real   (kind=rglu), allocatable :: v(:)
        integer(kind=iglu)              :: ln

        contains

        procedure :: del=>rvarDel
        procedure :: get=>rvarGet
    end type rvarVector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    procedure(outOfMemory),  pointer :: oomSub => null()
    procedure(changeMemory), pointer :: cmemSub => null()

    integer(kind=i8kind), protected :: glSharedMemory
    real   (kind=rglu)   :: lastTimerPoint=0


    character (len=200) :: fmtstr
    integer(kind=iglu)  :: defaultIO=6
    real   (kind=rglu)  :: pi,dtr
    character           :: glPrintString*500
    integer(kind=iglu)  :: iounit=6

    integer(kind=iglu)  :: void
    real   (kind=rglu)  :: voidr
    character           :: voidc*1
    logical(kind=lglu)  :: voidl

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    interface find
        module procedure find_c,find_uc,find_r16,find_r8,find_r4,find_i8,find_i4,find_i2,find_i1
    end interface find

    interface isEven
        module procedure isEven_r16,isEven_r8,isEven_r4,isEven_i8,isEven_i4,isEven_i2,isEven_i1
    end interface isEven

    interface same
        module procedure same_i8,same_i4,same_i2,same_i1
    end interface same

    interface mid
        module procedure mid_i8,mid_i4,mid_i2,mid_i1,mid_r16,&
                         mid_r8,mid_r4,mid_l8,mid_l4,mid_l2 ,&
                         mid_l1,mid_ch,mid_uc
    end interface mid

    interface lenTrimArray
        module procedure lenTrimArray_i8,lenTrimArray_i4 ,lenTrimArray_i2,&
                         lenTrimArray_i1,lenTrimArray_r16,lenTrimArray_r8,&
                         lenTrimArray_r4
    end interface lenTrimArray

    interface collectArray
        module procedure collectArray_i8,collectArray_i4 ,collectArray_i2,&
                         collectArray_i1,collectArray_r16,collectArray_r8,&
                         collectArray_r4
    end interface collectArray

    interface rangen
        module procedure random_generate_i8,random_generate_i4,&
                         random_generate_r8,random_generate_r4
    end interface rangen

    interface arrayAnalize
        module procedure arrayAnalize1,arrayAnalize2,arrayAnalize3,&
                         arrayAnalize4,arrayAnalize5,arrayAnalize6
    end interface arrayAnalize

    interface compareReal
        module procedure compare_r16,compare_r8,compare_r4
    end interface compareReal

    interface assignment (=)
        module procedure assign_uchSet ,assign_uchGet, &
                         assign_ivarSet,assign_ivarGet,&
                         assign_rvarSet,assign_rvarGet,&
                         assign_uchExchange
    end interface

    interface vec_set
        module procedure rvar_set,ivar_set
    end interface

    interface call_wrapper
        module procedure call_wrapper_i8,call_wrapper_i4,call_wrapper_i2,call_wrapper_i1,&
                         call_wrapper_l8,call_wrapper_l4,call_wrapper_l2,call_wrapper_l1,&
                         call_wrapper_r16,call_wrapper_r8,call_wrapper_ch,call_wrapper_uc
    end interface call_wrapper

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    ! module attributes
    public :: glVersion,glDate,glAuthor
    public :: glSetIOunit,glGetIOunit
    public :: glFinalize

    ! data type presets
    public :: iglu,ispu,lglu,NaN
    public :: rglu,gluZero,gluUnity,gluCompare
    public :: rspu,spuZero,spuUnity,spuCompare

    ! to be compatible with GNU Fortran. FUCK GNU FORTRAN.
    !public :: fmtstr

    ! types
    public :: uch,ivarVector,rvarVector
    ! type "methods"
    !public :: uchGet,uchSet,uchDel,delVector,setVector !,uch_set,vec_set

    public :: assignment (=)

    ! platform constants
    public :: os,osSeparator,osMove,osCopy
    ! platform functions
    public :: getPath,isPosix,getPOSIXinfo

    ! signals
    public :: signal,system,getpid
    public :: SIGABRT,SIGINT,SIGTERM,PXFSTRUCTCREATE,PXFSIGADDSET

    ! OpenMP functions
    public :: setThreadsNumber,getThreadNumber

    ! variables
    public :: defaultIO,glPrintString
    public :: true,false,void,voidr,voidc,voidl,pi,dtr

    ! tools
    public :: find,isEven,same,mid,lenTrimArray,collectArray,rangen
    public :: random_generate_array_r8,random_generate_array_r4
    public :: random_generate_array_i8,random_generate_array_i4
    public :: compareStrings,arrayAnalize,compareReal,purifyValues,definePi

    ! almost useless
    public :: infiniteVoidLoop,nullSub,screenProgress

    ! time and date
    public :: convertTime,timeControl,dayOfWeek,isLeapYear,timeStamp,date_time,timeControlCheckpoint

    ! memory
    public :: glControlMemory,glShareMemory,glMemoryLeft,glFromBytes,glSharedMemory

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MEMORY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function glControlMemory(bytesUsed,label,action) result(ret)
    implicit none

    integer(kind=i8kind), intent(in)           :: bytesUsed
    character (len=*)   , intent(in)           :: label
    character (len=*)   , intent(in), optional :: action
    integer(kind=i8kind)                       :: ubytes


    ubytes=bytesUsed
    if (present(action)) then
        if (action.EQ.'free') then
            ubytes=-abs(bytesUsed)
        elseif(action.EQ.'use') then
            ubytes= abs(bytesUsed)
        endif
    endif

    glSharedMemory=glSharedMemory-ubytes
    call cmemSub(label,ubytes)

    if (glSharedMemory.LT.0) call oomSub(label, ubytes, glSharedMemory+ubytes)

    ret=0; return
    end function glControlMemory

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function glShareMemory(bytes,violsub,chsub) result(ret)
    implicit none

    integer(kind=i8kind), intent(in) :: bytes
    external                         :: violsub,chsub


    oomSub => violsub; cmemSub => chsub
    glSharedMemory=bytes

    ret=0; return
    end function glShareMemory

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=r16kind) function glMemoryLeft(units) result(ret)
    implicit none

    character (len=*), optional    :: units
    character (len=:), allocatable :: dunits


    if (present(units)) then
        dunits=units
    else
        dunits='none'
    endif

    ret=glFromBytes(glSharedMemory, dunits)

    return
    end function glMemoryLeft

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure real(kind=r16kind) function glFromBytes(size, units) result(ret)
    implicit none

    integer(kind=i8kind), intent(in) :: size
    character (len=*),    intent(in) :: units
    real(kind=r16kind)               :: divisor


    divisor=1024._r16kind

    select case (units)
        case ('B','b','BYTES','bytes','BYTE','byte')
            ret=real(size, kind=r16kind)
        case ('KB','kb','Kb','KILOBYTES','kilobytes','KILOBYTE','kilobyte')
            ret=real(size, kind=r16kind)/divisor
        case ('MB','mb','Mb','MEGABYTES','megabytes','MEGABYTE','megabyte')
            ret=real(size, kind=r16kind)/(divisor*divisor)
        case ('GB','gb','Gb','GIGABYTES','gigabytes','GIGABYTE','gigabyte')
            ret=real(size, kind=r16kind)/(divisor*divisor*divisor)
        case default
            ret=real(size, kind=r16kind)
    end select

    return
    end function glFromBytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ UCH METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    elemental subroutine assign_uchSet(this,str)
    implicit none

    type(uch)         , intent(out) :: this
    character(len=*)  , intent(in)  :: str
    integer(kind=iglu)              :: i


    this%ln=len(str); allocate (this%ch(len(str)))
    do i = 1,len(str)
        this%ch(i)=str(i:i)
    enddo

    return
    end subroutine assign_uchSet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function uch_set(str) result(ret)
    implicit none

    type(uch)                      :: ret
    character (len=*), intent(in)  :: str
    integer(kind=iglu)             :: i


    ret%ln=len(str); allocate (ret%ch(len(str)))
    do i = 1,len(str)
        ret%ch(i)=str(i:i)
    enddo

    return
    end function uch_set

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    elemental subroutine assign_uchExchange(this,other)
    implicit none

    type(uch), intent(out) :: this
    type(uch), intent(in)  :: other
    integer(kind=iglu)     :: i


    this%ln=other%ln; allocate (this%ch(other%ln))
    do i = 1,other%ln
        this%ch(i)=other%ch(i)
    enddo

    return
    end subroutine assign_uchExchange

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine assign_uchGet(str,this)
    implicit none

    character (len=:), allocatable, intent(out) :: str
    type(uch)                     , intent(in)  :: this
    integer(kind=iglu)                          :: i


    allocate (character (len=this%ln) :: str)

    do i = 1,this%ln
        str(i:i)=this%ch(i)
    enddo

    return
    end subroutine assign_uchGet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    elemental subroutine uchDel(this)
    implicit none
    class(uch), intent(out) :: this


    if (allocated(this%ch)) then
        this%ln=0; deallocate (this%ch)
    endif
    return
    end subroutine uchDel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function uchGet(this,ustart,ustop,ustep) result(ret)
    implicit none
    class(uch)        , intent(in)             :: this
    integer(kind=iglu), intent(in), optional   :: ustart,ustop,ustep
    character (len=:) , allocatable            :: ret
    integer(kind=iglu)                         :: dstart,dstop,dstep,i,ln,cnt


    dstart=1     ; if (present(ustart)) dstart=ustart
    dstop=this%ln; if (present(ustop) ) dstop =ustop
    dstep=1      ; if (present(ustep) ) dstep =ustep

    ln=(dstop-dstart)/dstep+1; allocate (character (len=ln) :: ret)

    cnt=0
    do i = dstart,dstop,dstep
        cnt=cnt+1; ret(cnt:cnt)=this%ch(i)
    enddo

    return
    end function uchGet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function uchLen(this) result(ret)
    implicit none
    class(uch), intent(in) :: this


    ret=this%ln; return
    end function uchLen

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VAR VECTOR METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine assign_ivarGet(arr,this)
    implicit none
    integer(kind=iglu), allocatable, intent(out) :: arr(:)
    type(ivarVector)               , intent(in)  :: this

    if (allocated(arr)) deallocate (arr)
    allocate (arr(this%ln)); arr=this%v

    return
    end subroutine assign_ivarGet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine assign_ivarSet(this,arr)
    implicit none
    integer(kind=iglu), intent(in)  :: arr(:)
    type(ivarVector)  , intent(out) :: this

    if (allocated(this%v)) deallocate (this%v)
    this%ln=UBound(arr,1)
    allocate (this%v(UBound(arr,1))); this%v=arr

    return
    end subroutine assign_ivarset

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function ivar_set(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    type(ivarVector)               :: ret


    ret%ln=UBound(arr,1)
    allocate (ret%v(ret%ln)); ret%v=arr

    return
    end function ivar_set

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function rvar_set(arr) result(ret)
    implicit none

    real(kind=rglu), intent(in) :: arr(:)
    type(rvarVector)            :: ret


    ret%ln=UBound(arr,1)
    allocate (ret%v(ret%ln)); ret%v=arr

    return
    end function rvar_set

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    elemental subroutine ivarDel(this)
    implicit none
    class(ivarVector), intent(out) :: this


    if (allocated(this%v)) then
        this%ln=0; deallocate (this%v)
    endif

    return
    end subroutine ivarDel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function ivarGet(this,ustart,ustop,ustep) result(ret)
    implicit none
    class(ivarVector) , intent(in)           :: this
    integer(kind=iglu), intent(in), optional :: ustart,ustop,ustep
    integer(kind=iglu), allocatable          :: ret(:)
    integer(kind=iglu)                       :: dstart,dstop,dstep,ln,i


    dstart=1     ; if (present(ustart)) dstart=ustart
    dstop=this%ln; if (present(ustop))  dstop =ustop
    dstep=1      ; if (present(ustep))  dstep =ustep

    ln=(dstop-dstart+1)/dstep; allocate (ret(ln))

    do i = dstart,dstop,dstep
        ret(i:i)=this%v(i)
    enddo

    return
    end function ivarGet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function ivarLen(this) result(ret)
    implicit none
    class(ivarVector), intent(in) :: this

    ret=this%ln; return
    end function ivarlen

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine assign_rvarGet(arr,this)
    implicit none
    real(kind=rglu) , allocatable, intent(out) :: arr(:)
    type(rvarVector)             , intent(in)  :: this

    if (allocated(arr)) deallocate (arr)
    allocate (arr(this%ln)); arr=this%v

    return
    end subroutine assign_rvarGet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine assign_rvarSet(this,arr)
    implicit none
    real(kind=rglu) , intent(in)  :: arr(:)
    type(rvarVector), intent(out) :: this

    if (allocated(this%v)) deallocate (this%v)
    this%ln=UBound(arr,1)
    allocate (this%v(UBound(arr,1))); this%v=arr

    return
    end subroutine assign_rvarset

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    elemental subroutine rvarDel(this)
    implicit none
    class(rvarVector), intent(out) :: this


    if (allocated(this%v)) then
        this%ln=0; deallocate (this%v)
    endif

    return
    end subroutine rvarDel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function rvarGet(this,ustart,ustop,ustep) result(ret)
    implicit none
    class(rvarVector) , intent(in)           :: this
    integer(kind=iglu), intent(in), optional :: ustart,ustop,ustep
    real(kind=rglu)   , allocatable          :: ret(:)
    integer(kind=iglu)                       :: dstart,dstop,dstep,ln,i


    dstart=1     ; if (present(ustart)) dstart=ustart
    dstop=this%ln; if (present(ustop))  dstop =ustop
    dstep=1      ; if (present(ustep))  dstep =ustep

    ln=(dstop-dstart+1)/dstep; allocate (ret(ln))

    do i = dstart,dstop,dstep
        ret(i:i)=this%v(i)
    enddo

    return
    end function rvarGet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function rvarLen(this) result(ret)
    implicit none
    class(rvarVector), intent(in) :: this

    ret=this%ln; return
    end function rvarlen

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VALUE LEN FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_i8(myInt)  result(ret) ! my integer's digits (old name)
    implicit none
    integer(kind=i8kind), intent(in) :: myInt


    if (myInt.EQ.0) then
        ret=1
        return
    endif
    if (myInt.LT.0) then
        ret=int( log10(abs(dble(myInt))) + 1 )+1
        return
    endif

    ret=int( log10(dble(myInt)) + 1 ); return
    end function mid_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_i4(myInt) result(ret)
    implicit none
    integer(kind=i4kind), intent(in) :: myInt


    if (myInt.EQ.0) then
        ret=1
        return
    endif
    if (myInt.LT.0) then
        ret=int( log10(abs(dble(myInt))) + 1 )+1
        return
    endif

    ret=int( log10(dble(myInt)) + 1 ); return
    end function mid_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_i2(myInt) result(ret)
    implicit none
    integer(kind=i2kind), intent(in) :: myInt


    if (myInt.EQ.0) then
        ret=1
        return
    endif
    if (myInt.LT.0) then
        ret=int( log10(abs(dble(myInt))) + 1 )+1
        return
    endif

    ret=int( log10(dble(myInt)) + 1 ); return
    end function mid_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_i1(myInt) result(ret)
    implicit none
    integer(kind=i1kind), intent(in) :: myInt


    if (myInt.EQ.0) then
        ret=1
        return
    endif
    if (myInt.LT.0) then
        ret=int( log10(abs(dble(myInt))) + 1 )+1
        return
    endif

    ret=int( log10(dble(myInt)) + 1 ); return
    end function mid_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_r16(myReal) result(ret)
    implicit none
    real(kind=r16kind), intent(in) :: myReal


    if (abs(myReal).LT.1._r16kind) then
        if (myReal.LT.0._r16kind) then
            ret=2
        else
            ret=1
        endif
        return
    endif
    if (myReal.LT.0._r16kind) then
        ret=int( log10(abs(myReal))+1 )+1
        return
    endif
    ret=int( log10(myReal)+1 ); return
    end function mid_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_r8(myReal) result(ret)
    implicit none
    real(kind=r8kind), intent(in) :: myReal


    if (abs(myReal).LT.1._r8kind) then
        if (myReal.LT.0._r8kind) then
            ret=2
        else
            ret=1
        endif
        return
    endif
    if (myReal.LT.0._r8kind) then
        ret=int( log10(abs(myReal))+1 )+1
        return
    endif
    ret=int( log10(myReal)+1 ); return
    end function mid_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_r4(myReal) result(ret)
    implicit none
    real(kind=r4kind), intent(in) :: myReal


    if (abs(myReal).LT.1._r4kind) then
        if (myReal.LT.0._r4kind) then
            ret=2
        else
            ret=1
        endif
        return
    endif
    if (myReal.LT.0._r4kind) then
        ret=int( log10(abs(myReal))+1 )+1
        return
    endif

    ret=int( log10(myReal)+1 ); return
    end function mid_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_l8(logi) result(ret)
    implicit none
    logical(kind=l8kind), intent(in) :: logi


    ret=5; if (logi) ret=4; return
    end function mid_l8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_l4(logi) result(ret)
    implicit none
    logical(kind=l4kind), intent(in) :: logi


    ret=5; if (logi) ret=4; return
    end function mid_l4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_l2(logi) result(ret)
    implicit none
    logical(kind=l2kind), intent(in) :: logi


    ret=5; if (logi) ret=4; return
    end function mid_l2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_l1(logi) result(ret)
    implicit none
    logical(kind=l1kind), intent(in) :: logi


    ret=5; if (logi) ret=4; return
    end function mid_l1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_ch(str) result(ret)
    implicit none
    character (len=*), intent(in) :: str


    ret=len(str); return
    end function mid_ch

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function mid_uc(str) result(ret)
    implicit none
    type(uch), intent(in) :: str


    ret=str%ln; return
    end function mid_uc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIM ARRAY FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function lenTrimArray_i8(arr) result(ret)
    implicit none
    integer(kind=i8kind), intent(in) :: arr(:)
    integer(kind=iglu)               :: ln,i


    ln=UBound(arr,1)

    ret=0
    do i = ln,1,-1
        if (arr(i).NE.0) then
            ret=i
            exit
        endif
    enddo

    return  ! return the position of last non-zero element
    end function lenTrimArray_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function lenTrimArray_i4(arr) result(ret)
    implicit none
    integer(kind=i4kind), intent(in) :: arr(:)
    integer(kind=iglu)               :: ln,i


    ln=UBound(arr,1)

    ret=0
    do i = ln,1,-1
        if (arr(i).NE.0) then
            ret=i
            exit
        endif
    enddo

    return  ! return the position of last non-zero element
    end function lenTrimArray_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function lenTrimArray_i2(arr) result(ret)
    implicit none
    integer(kind=i2kind), intent(in) :: arr(:)
    integer(kind=iglu)               :: ln,i


    ln=UBound(arr,1)

    ret=0
    do i = ln,1,-1
        if (arr(i).NE.0) then
            ret=i
            exit
        endif
    enddo

    return  ! return the position of last non-zero element
    end function lenTrimArray_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function lenTrimArray_i1(arr) result(ret)
    implicit none
    integer(kind=i1kind), intent(in) :: arr(:)
    integer(kind=iglu)               :: ln,i


    ln=UBound(arr,1)

    ret=0
    do i = ln,1,-1
        if (arr(i).NE.0) then
            ret=i
            exit
        endif
    enddo

    return  ! return the position of last non-zero element
    end function lenTrimArray_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function lenTrimArray_r16(arr,tol) result(ret)
    implicit none

    real   (kind=r16kind), intent(in)           :: arr(:)
    real   (kind=r16kind), intent(in), optional :: tol
    real   (kind=r16kind)                       :: deftol
    integer(kind=iglu)                          :: ln,i


    ln=UBound(arr,1)

    deftol=epsilon(deftol)*10._r16kind; if (present(tol)) deftol=tol

    ret=0
    do i = ln,1,-1
        if (arr(i).GT.deftol) then
            ret=i
            exit
        endif
    enddo

    return ! return the position of last non-zero element
    end function lenTrimArray_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function lenTrimArray_r8(arr,tol) result(ret)
    implicit none

    real   (kind=r8kind), intent(in)           :: arr(:)
    real   (kind=r8kind), intent(in), optional :: tol
    real   (kind=r8kind)                       :: deftol
    integer(kind=iglu)                         :: ln,i


    ln=UBound(arr,1)

    deftol=epsilon(deftol)*10._r8kind; if (present(tol)) deftol=tol

    ret=0
    do i = ln,1,-1
        if (arr(i).GT.deftol) then
            ret=i
            exit
        endif
    enddo

    return  ! return the position of last non-zero element
    end function lenTrimArray_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function lenTrimArray_r4(arr,tol) result(ret)
    implicit none

    real   (kind=r4kind), intent(in)           :: arr(:)
    real   (kind=r4kind), intent(in), optional :: tol
    real   (kind=r4kind)                       :: deftol
    integer(kind=iglu)                         :: ln,i


    ln=UBound(arr,1)

    deftol=epsilon(deftol)*10._r4kind; if (present(tol)) deftol=tol

    ret=0
    do i = ln,1,-1
        if (arr(i).GT.deftol) then
            ret=i
            exit
        endif
    enddo

    return  ! return the position of last non-zero element
    end function lenTrimArray_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COLLECT ARRAY FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function collectArray_i8(arr) result(ret)
    implicit none
    integer(kind=i8kind) :: arr(:)
    integer(kind=iglu)   :: ln,i,lst,els
    integer(kind=iglu)   :: fzero,fnzero


    ln=UBound(arr,1) ! define size of arr

    els=0
    do i = 1,ln ! define number of non-zero elements
        if (arr(i).NE.0) els=els+1
    enddo

    if (els.EQ.0) then
        ret=0; return ! all elements equal zero
    endif

    do i = 1,ln ! find first non-zero element
        if (arr(i).NE.0) then
            fnzero=i
            exit
        endif
    enddo

    if (fnzero.NE.1) then ! put this element to the first position
        arr(1)=arr(fnzero); arr(fnzero)=0
    endif

    do
        lst=lenTrimArray(arr) ! check "trim" of arr
        if (lst.EQ.els) exit ! done.

        do i = 1,lst
            if (arr(i).EQ.0) then ! first zero element
                fzero=i
                exit
            endif
        enddo

        do i = fzero+1,lst
            if (arr(i).NE.0) then ! first non-zero element
                fnzero=i
                exit
            endif
        enddo
        arr(fzero)=arr(fnzero); arr(fnzero)=0 ! swap.
    enddo

    ret=lenTrimArray(arr) ! length of trim for "trimed" arr.
    return
    end function collectArray_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function collectArray_i4(arr) result(ret)
    implicit none
    integer(kind=i4kind) :: arr(:)
    integer(kind=iglu)   :: ln,i,lst,els
    integer(kind=iglu)   :: fzero,fnzero


    ln=UBound(arr,1) ! define size of arr

    els=0
    do i = 1,ln ! define number of non-zero elements
        if (arr(i).NE.0) els=els+1
    enddo

    if (els.EQ.0) then
        ret=0; return ! all elements equal zero
    endif

    do i = 1,ln ! find first non-zero element
        if (arr(i).NE.0) then
            fnzero=i
            exit
        endif
    enddo

    if (fnzero.NE.1) then ! put this element to the first position
        arr(1)=arr(fnzero); arr(fnzero)=0
    endif

    do
        lst=lenTrimArray(arr) ! check "trim" of arr
        if (lst.EQ.els) exit ! done.

        do i = 1,lst
            if (arr(i).EQ.0) then ! first zero element
                fzero=i
                exit
            endif
        enddo

        do i = fzero+1,lst
            if (arr(i).NE.0) then ! first non-zero element
                fnzero=i
                exit
            endif
        enddo
        arr(fzero)=arr(fnzero); arr(fnzero)=0 ! swap.
    enddo

    ret=lenTrimArray(arr) ! length of trim for "trimed" arr.
    return
    end function collectArray_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function collectArray_i2(arr) result(ret)
    implicit none
    integer(kind=i2kind) :: arr(:)
    integer(kind=iglu)   :: ln,i,lst,els
    integer(kind=iglu)   :: fzero,fnzero


    ln=UBound(arr,1) ! define size of arr

    els=0
    do i = 1,ln ! define number of non-zero elements
        if (arr(i).NE.0) els=els+1
    enddo

    if (els.EQ.0) then
        ret=0; return ! all elements equal zero
    endif

    do i = 1,ln ! find first non-zero element
        if (arr(i).NE.0) then
            fnzero=i
            exit
        endif
    enddo

    if (fnzero.NE.1) then ! put this element to the first position
        arr(1)=arr(fnzero); arr(fnzero)=0
    endif

    do
        lst=lenTrimArray(arr) ! check "trim" of arr
        if (lst.EQ.els) exit ! done.

        do i = 1,lst
            if (arr(i).EQ.0) then ! first zero element
                fzero=i
                exit
            endif
        enddo

        do i = fzero+1,lst
            if (arr(i).NE.0) then ! first non-zero element
                fnzero=i
                exit
            endif
        enddo
        arr(fzero)=arr(fnzero); arr(fnzero)=0 ! swap.
    enddo

    ret=lenTrimArray(arr) ! length of trim for "trimed" arr.
    return
    end function collectArray_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function collectArray_i1(arr) result(ret)
    implicit none
    integer(kind=i1kind) :: arr(:)
    integer(kind=iglu)   :: ln,i,lst,els
    integer(kind=iglu)   :: fzero,fnzero


    ln=UBound(arr,1) ! define size of arr

    els=0
    do i = 1,ln ! define number of non-zero elements
        if (arr(i).NE.0) els=els+1
    enddo

    if (els.EQ.0) then
        ret=0; return ! all elements equal zero
    endif

    do i = 1,ln ! find first non-zero element
        if (arr(i).NE.0) then
            fnzero=i
            exit
        endif
    enddo

    if (fnzero.NE.1) then ! put this element to the first position
        arr(1)=arr(fnzero); arr(fnzero)=0
    endif

    do
        lst=lenTrimArray(arr) ! check "trim" of arr
        if (lst.EQ.els) exit ! done.

        do i = 1,lst
            if (arr(i).EQ.0) then ! first zero element
                fzero=i
                exit
            endif
        enddo

        do i = fzero+1,lst
            if (arr(i).NE.0) then ! first non-zero element
                fnzero=i
                exit
            endif
        enddo
        arr(fzero)=arr(fnzero); arr(fnzero)=0 ! swap.
    enddo

    ret=lenTrimArray(arr) ! length of trim for "trimed" arr.
    return
    end function collectArray_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function collectArray_r16(arr,tol) result(ret)
    implicit none
    real(kind=r16kind)           :: deftol,arr(:),hold
    real(kind=r16kind), optional :: tol
    integer(kind=iglu)           :: ln,i,lst,els
    integer(kind=iglu)           :: fzero,fnzero


    ln=UBound(arr,1) ! define size of arr
    deftol=epsilon(deftol)*10._r16kind; if (present(tol)) deftol=tol

    els=0
    do i = 1,ln ! define number of non-zero elements
        if (arr(i).GT.deftol) els=els+1
    enddo

    if (els.EQ.0) then
        ret=0; return ! all elements equal zero
    endif

    do i = 1,ln ! find first non-zero element
        if (arr(i).GT.deftol) then
            fnzero=i
            exit
        endif
    enddo

    if (fnzero.NE.1) then ! put this element to the first position
        hold=arr(1); arr(1)=arr(fnzero); arr(fnzero)=hold
    endif

    do
        lst=lenTrimArray(arr) ! check "trim" of arr
        if (lst.EQ.els) exit ! done.

        do i = 1,lst
            if (arr(i).LE.deftol) then ! first zero element
                fzero=i
                exit
            endif
        enddo

        do i = fzero+1,lst
            if (arr(i).GT.deftol) then ! first non-zero element
                fnzero=i
                exit
            endif
        enddo
        hold=arr(fzero); arr(fzero)=arr(fnzero); arr(fnzero)=hold ! swap.
    enddo

    ret=lenTrimArray(arr) ! length of trim for "trimed" arr.
    return
    end function collectArray_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function collectArray_r8(arr,tol) result(ret)
    implicit none
    real(kind=r8kind)           :: deftol,arr(:),hold
    real(kind=r8kind), optional :: tol
    integer(kind=iglu)          :: ln,i,lst,els
    integer(kind=iglu)          :: fzero,fnzero


    ln=UBound(arr,1) ! define size of arr
    deftol=epsilon(deftol)*10._r8kind; if (present(tol)) deftol=tol

    els=0
    do i = 1,ln ! define number of non-zero elements
        if (arr(i).GT.deftol) els=els+1
    enddo

    if (els.EQ.0) then
        ret=0; return ! all elements equal zero
    endif

    do i = 1,ln ! find first non-zero element
        if (arr(i).GT.deftol) then
            fnzero=i
            exit
        endif
    enddo

    if (fnzero.NE.1) then ! put this element to the first position
        hold=arr(1); arr(1)=arr(fnzero); arr(fnzero)=hold
    endif

    do
        lst=lenTrimArray(arr) ! check "trim" of arr
        if (lst.EQ.els) exit ! done.

        do i = 1,lst
            if (arr(i).LE.deftol) then ! first zero element
                fzero=i
                exit
            endif
        enddo

        do i = fzero+1,lst
            if (arr(i).GT.deftol) then ! first non-zero element
                fnzero=i
                exit
            endif
        enddo
        hold=arr(fzero); arr(fzero)=arr(fnzero); arr(fnzero)=hold ! swap.
    enddo

    ret=lenTrimArray(arr) ! length of trim for "trimed" arr.
    return
    end function collectArray_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function collectArray_r4(arr,tol) result(ret)
    implicit none
    real(kind=r4kind)           :: deftol,arr(:),hold
    real(kind=r4kind), optional :: tol
    integer(kind=iglu)          :: ln,i,lst,els
    integer(kind=iglu)          :: fzero,fnzero


    ln=UBound(arr,1) ! define size of arr
    deftol=epsilon(deftol)*10; if (present(tol)) deftol=tol

    els=0
    do i = 1,ln ! define number of non-zero elements
        if (arr(i).GT.deftol) els=els+1
    enddo

    if (els.EQ.0) then
        ret=0; return ! all elements equal zero
    endif

    do i = 1,ln ! find first non-zero element
        if (arr(i).GT.deftol) then
            fnzero=i
            exit
        endif
    enddo

    if (fnzero.NE.1) then ! put this element to the first position
        hold=arr(1); arr(1)=arr(fnzero); arr(fnzero)=hold
    endif

    do
        lst=lenTrimArray(arr) ! check "trim" of arr
        if (lst.EQ.els) exit ! done.

        do i = 1,lst
            if (arr(i).LE.deftol) then ! first zero element
                fzero=i
                exit
            endif
        enddo

        do i = fzero+1,lst
            if (arr(i).GT.deftol) then ! first non-zero element
                fnzero=i
                exit
            endif
        enddo
        hold=arr(fzero); arr(fzero)=arr(fnzero); arr(fnzero)=hold ! swap.
    enddo

    ret=lenTrimArray(arr) ! length of trim for "trimed" arr.
    return
    end function collectArray_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAME FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function same_i8(x,y) result(ret)
    implicit none
    integer(kind=i8kind), intent(in) :: x,y


    ret=.NOT.btest(x+y,0); return
    end function same_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function same_i4(x,y) result(ret)
    implicit none
    integer(kind=i4kind), intent(in) :: x,y


    ret=.NOT.btest(x+y,0); return
    end function same_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function same_i2(x,y) result(ret)
    implicit none
    integer(kind=i2kind), intent(in) :: x,y


    ret=.NOT.btest(x+y,0); return
    end function same_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function same_i1(x,y) result(ret)
    implicit none
    integer(kind=i1kind), intent(in) :: x,y


    ret=.NOT.btest(x+y,0); return
    end function same_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ IS EVEN FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function isEven_i8(x) result(ret)
    implicit none
    integer(kind=i8kind), intent(in) :: x


    ret=.NOT.btest(x,0); return
    end function isEven_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function isEven_i4(x) result(ret)
    implicit none
    integer(kind=i4kind), intent(in) :: x


    ret=.NOT.btest(x,0); return
    end function isEven_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function isEven_i2(x) result(ret)
    implicit none
    integer(kind=i2kind), intent(in) :: x


    ret=.NOT.btest(x,0); return
    end function isEven_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function isEven_i1(x) result(ret)
    implicit none
    integer(kind=i1kind), intent(in) :: x


    ret=.NOT.btest(x,0); return
    end function isEven_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function isEven_r16(x) result(ret)
    implicit none
    real(kind=r16kind), intent(in) :: x


    ret=isEven_i8(int(x,kind=i8kind)); return
    end function isEven_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function isEven_r8(x) result(ret)
    implicit none
    real(kind=r8kind), intent(in) :: x


    ret=isEven_i8(int(x,kind=i8kind)); return
    end function isEven_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function isEven_r4(x) result(ret)
    implicit none
    real(kind=r4kind), intent(in) :: x


    ret=isEven_i8(int(x,kind=i8kind)); return
    end function isEven_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PSEUDO RANDOM VALUE FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=i8kind) function random_generate_i8(insta,insto) result(ret)
    implicit none
    integer(kind=i8kind), intent(in) :: insta,insto
    integer(kind=i8kind)             :: setG,sigsw,inum,multip
    real   (kind=r8kind)             :: rnum


    if (insta.GT.insto) then
        ret=0; return
    endif

    if (insto-insta.EQ.0) then
        ret=insta; return
    endif

    setG=max(abs(insta),abs(insto))
    setG=int(log10(real(setG,r8kind)),i8kind)+1

    if ((insta.LT.0).AND.(insto.LE.0)) sigsw=1
    if ((insta.LT.0).AND.(insto.GT.0)) sigsw=2
    if ((insta.GE.0).AND.(insto.GT.0)) sigsw=3

    multip=int(10**setG,i8kind)

    do
        call random_number(rnum)
        inum=int(rnum*multip,i8kind)
        if ((inum.GE.insta).AND.(inum.LE.insto)) exit
    enddo

    select case (sigsw)
        case (1); inum=-inum
        case (2)
            call random_number(rnum)
            if (rnum.GT.0.5_rglu) then
                continue
            else
                inum=-inum
            endif
            if (inum.LE.insta) inum=-inum
            if (inum.GE.insto) inum=-inum
        case (3); continue
    end select

    ret=inum; return
    end function random_generate_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=i4kind) function random_generate_i4(insta,insto) result(ret)
    implicit none
    integer(kind=i4kind), intent(in) :: insta,insto
    integer(kind=i4kind)             :: setG,sigsw,inum,multip
    real   (kind=r8kind)             :: rnum


    if (insta.GT.insto) then
        ret=0; return
    endif

    if (insto-insta.EQ.0) then
        ret=insta; return
    endif

    setG=max(abs(insta),abs(insto))
    setG=int(log10(real(setG,r8kind)),i4kind)+1

    if ((insta.LT.0).AND.(insto.LE.0)) sigsw=1
    if ((insta.LT.0).AND.(insto.GT.0)) sigsw=2
    if ((insta.GE.0).AND.(insto.GT.0)) sigsw=3

    multip=int(10**setG,i4kind)

    do
        call random_number(rnum)
        inum=int(rnum*multip,i4kind)
        if ((inum.GE.insta).AND.(inum.LE.insto)) exit
    enddo

    select case (sigsw)
        case (1); inum=-inum
        case (2)
            call random_number(rnum)
            if (rnum.GT.0.5_rglu) then
                continue
            else
                inum=-inum
            endif
            if (inum.LE.insta) inum=-inum
            if (inum.GE.insto) inum=-inum
        case (3); continue
    end select

    ret=inum; return
    end function random_generate_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=r8kind) function random_generate_r8(insta,insto) result(ret)
    implicit none
    real   (kind=r8kind), intent(in) :: insta,insto
    integer(kind=iglu)               :: sigsw
    real   (kind=r8kind)             :: setG,rnum,renum,multip


    if (insta.GT.insto) then; ret=0; return; endif
    if (insto-insta.LE.gluCompare) then; ret=insta; return; endif

    setG=max(abs(insta),abs(insto))
    setG=int(log10(setG))+1

    if ((insta.LT.0).AND.(insto.LE.0)) sigsw=1
    if ((insta.LT.0).AND.(insto.GT.0)) sigsw=2
    if ((insta.GE.0).AND.(insto.GT.0)) sigsw=3

    multip=10.**setG

    do
        call random_number(rnum)
        rnum=rnum*multip
        if ((rnum.GE.insta).AND.(rnum.LE.insto)) then
            renum=rnum; exit
        endif
    enddo

    select case (sigsw)
        case (1); renum=-renum
        case (2)
            call random_number(rnum)
            if (rnum.GT.0.5_rglu) then
                continue
            else
                renum=-renum
            endif
            if (renum.LE.insta) renum=-renum
            if (renum.GE.insto) renum=-renum
        case (3); continue
    end select

    ret=renum; return
    end function random_generate_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=r4kind) function random_generate_r4(insta,insto) result(ret)
    implicit none
    real   (kind=r4kind), intent(in) :: insta,insto
    integer(kind=iglu)               :: sigsw
    real   (kind=r4kind)             :: setG,rnum,renum,multip


    if (insta.GT.insto) then; ret=0; return; endif
    if (insto-insta.LE.gluCompare) then; ret=insta; return; endif

    setG=max(abs(insta),abs(insto))
    setG=int(log10(setG))+1

    if ((insta.LT.0).AND.(insto.LE.0)) sigsw=1
    if ((insta.LT.0).AND.(insto.GT.0)) sigsw=2
    if ((insta.GE.0).AND.(insto.GT.0)) sigsw=3

    multip=10.**setG

    do
        call random_number(rnum)
        rnum=rnum*multip
        if ((rnum.GE.insta).AND.(rnum.LE.insto)) then
            renum=rnum; exit
        endif
    enddo

    select case (sigsw)
        case (1); renum=-renum
        case (2)
            call random_number(rnum)
            if (rnum.GT.0.5_rglu) then
                continue
            else
                renum=-renum
            endif
            if (renum.LE.insta) renum=-renum
            if (renum.GE.insto) renum=-renum
        case (3); continue
    end select

    ret=renum; return
    end function random_generate_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine random_generate_array_r8(size,array,insta,insto)
    implicit none
    integer(kind=iglu)  , intent(in) :: size
    real   (kind=r8kind)             :: array(size)

    real   (kind=r8kind), intent(in) :: insta,insto
    integer(kind=iglu)               :: k


    if (insta.GT.insto) then; array=0; return; endif
    if (insto-insta.LE.gluCompare) then; array=insta; return; endif

    do k = 1,size
        array(k)=random_generate_r8(insta,insto)
    enddo

    return
    end subroutine random_generate_array_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine random_generate_array_r4(size,array,insta,insto)
    implicit none
    integer(kind=iglu)  , intent(in) :: size
    real   (kind=r4kind)             :: array(size)

    real   (kind=r4kind), intent(in) :: insta,insto
    integer(kind=iglu)               :: k


    if (insta.GT.insto) then; array=0; return; endif
    if (insto-insta.LE.gluCompare) then; array=insta; return; endif

    do k = 1,size
        array(k)=random_generate_r4(insta,insto)
    enddo

    return
    end subroutine random_generate_array_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine random_generate_array_i8(size,array,insta,insto)
    implicit none
    integer(kind=iglu)  , intent(in) :: size
    integer(kind=i8kind), intent(in) :: insta,insto
    integer(kind=i8kind)             :: array(size)
    integer(kind=iglu)               :: k


    if (insta.GT.insto) then
        array=0; return
    endif

    if (insto-insta.EQ.0) then
        array=insta; return
    endif

    do k = 1,size
        array(k)=random_generate_i8(insta,insto)
    enddo

    return
    end subroutine random_generate_array_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine random_generate_array_i4(size,array,insta,insto)
    implicit none
    integer(kind=iglu)  , intent(in) :: size
    integer(kind=i4kind), intent(in) :: insta,insto
    integer(kind=i4kind)             :: array(size)
    integer(kind=iglu)               :: k


    if (insta.GT.insto) then
        array=0; return
    endif

    if (insto-insta.EQ.0) then
        array=insta; return
    endif

    do k = 1,size
        array(k)=random_generate_i4(insta,insto)
    enddo

    return
    end subroutine random_generate_array_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIND FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_uc(array,val) result(ret)
    implicit none
    type(uch)        , intent(in) :: array(:)
    character (len=*), intent(in) :: val
    integer(kind=iglu)            :: i


    do i = lbound(array,1),ubound(array,1)
        if ( trim(val) .EQ. trim(uchGet(array(i))) ) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_uc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_c(array,val) result(ret)
    implicit none
    character (len=*), intent(in) :: array(:)
    character (len=*), intent(in) :: val
    integer(kind=iglu)            :: i


    do i = lbound(array,1),ubound(array,1)
        if ( trim(val) .EQ. trim(array(i)) ) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_r16(array,val) result(ret)
    implicit none
    real   (kind=r16kind), intent(in) :: val,array(:)
    integer(kind=iglu)                :: i


    do i = lbound(array,1),ubound(array,1)
        if (abs(array(i)-val).LE.epsilon(val)*10._r16kind) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_r8(array,val) result(ret)
    implicit none
    real   (kind=r8kind), intent(in) :: val,array(:)
    integer(kind=iglu)               :: i


    do i = lbound(array,1),ubound(array,1)
        if (abs(array(i)-val).LE.epsilon(val)*10) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_r4(array,val) result(ret)
    implicit none
    real   (kind=r4kind), intent(in) :: val,array(:)
    integer(kind=iglu)               :: i


    do i = lbound(array,1),ubound(array,1)
        if (abs(array(i)-val).LE.epsilon(val)*10) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_i8(array,val) result(ret)
    implicit none
    integer(kind=i8kind), intent(in) :: val,array(:)
    integer(kind=iglu)               :: i


    do i = lbound(array,1),ubound(array,1)
        if (array(i).EQ.val) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_i4(array,val) result(ret)
    implicit none
    integer(kind=i4kind), intent(in) :: val,array(:)
    integer(kind=iglu)               :: i


    do i = lbound(array,1),ubound(array,1)
        if (array(i).EQ.val) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_i2(array,val) result(ret)
    implicit none
    integer(kind=i2kind), intent(in) :: val,array(:)
    integer(kind=iglu)               :: i


    do i = lbound(array,1),ubound(array,1)
        if (array(i).EQ.val) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function find_i1(array,val) result(ret)
    implicit none
    integer(kind=i1kind), intent(in) :: val,array(:)
    integer(kind=iglu)               :: i


    do i = lbound(array,1),ubound(array,1)
        if (array(i).EQ.val) then
            ret=i; return
        endif
    enddo

    ret=-1; return
    end function find_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ POSIX FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getPOSIXinfo(hostname,login,memory,ppid)
    implicit none

    type(uch),            intent(out) :: hostname, login
    integer(kind=i4kind), intent(out) :: memory,ppid

    character(len=64)                 :: chostname, clogin
    integer(kind=i4kind)              :: ierr,pgsz,pgcnt,varnum,loglen,uid
    integer(kind=i4kind)              :: upasswd_structure



    chostname=repeat(' ', len(chostname)); clogin=repeat(' ', len(clogin))
    hostname=''; login=''
    memory=0; ppid=0
#   if(__OS==2)
        ! hostname
        ierr=hostnam(chostname)
        hostname=trim(chostname)

        ! login
        ! more safe way to get login, since PXFGETLOGIN sometimes leads to segmentation fault
        call PXFSTRUCTCREATE('passwd',upasswd_structure,ierr)
        call PXFGETPWUID(GETUID(),upasswd_structure,ierr)
        call PXFSTRGET(upasswd_structure,'pw_name',clogin,loglen,ierr)

        login=clogin(1:loglen)

        ! memory
        call PXFCONST('_SC_PAGESIZE',varnum,ierr)
        call PXFSYSCONF(varnum,pgsz,ierr)
        call PXFCONST('_SC_PHYS_PAGES',varnum,ierr)
        call PXFSYSCONF(varnum,pgcnt,ierr)

        memory=int(real(pgsz, kind=r8kind)*real(pgcnt, kind=r8kind)/1024.**2)

        ! parental process PID
        call PXFGETPPID(ppid,ierr)
#   endif


    return
    end subroutine getPOSIXinfo

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TIME AND DATE FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function timeControl(cputime) result(ret)
    implicit none

    real(kind=rglu), optional :: cputime


#   if(__opnmp==1)
        ret=omp_get_wtime()
#   else
        call cpu_time(ret)
#   endif
    if (present(cputime)) call cpu_time(cputime)

    return
    end function timeControl

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function timeControlCheckpoint(msg,save,drop,raw,rcall) result(ret)
    implicit none

    character (len=*),  optional, intent(in) :: msg
    logical(kind=lglu), optional, intent(in) :: save,drop,raw
    integer(kind=iglu), optional, intent(in) :: rcall
    logical(kind=lglu)                       :: usave,udrop,uraw,prnt


    prnt=present(msg)
    udrop=false; if (present(drop)) udrop=drop
    usave=true ; if (present(save)) usave=save
    uraw =false; if (present(raw))  uraw =raw

    if (udrop) lastTimerPoint=0

    if (lastTimerPoint.EQ.0) then
        lastTimerPoint=timecontrol()
        ret=lastTimerPoint
        if (prnt) then
            if (uraw) then
                write (iounit,"(' [ ',F13.3,' ] ',A)") 0.,msg !fmt on more than 30 years.
            else
                write (iounit,'(A,1X,A)') msg,timeStamp()
            endif
        endif
        return
    endif

    ret=timecontrol()-lastTimerPoint
    if (.NOT.usave) lastTimerPoint=timecontrol()
    if (prnt) then
        if (uraw) then
            write (iounit,"(' [ ',F13.3,' ] ',A)") ret,msg !fmt on more than 30 years.
        else
            write (iounit,'(A,1X,A)') msg,convertTime(ret)
        endif
    endif

    return
    end function timeControlCheckpoint

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function isLeapYear(year) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: year


    if ( ((mod(year,4).EQ.0).AND.(mod(year,100).NE.0)) .OR. (mod(year,400).EQ.0)) then
        ret=1
    else
        ret=0
    endif

    return
    end function isLeapYear

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function dayOfWeek(year,month,day,past) result(ret)
    implicit none

    integer(kind=iglu), parameter  :: offsetYear=1900,offsetMonth    =1,&
                                          offsetDay =1   ,offsetDayOfWeek=1 ! 1=Monday

    integer(kind=iglu), intent(in)            :: year,month,day
    integer(kind=iglu), intent(out), optional :: past
    integer(kind=iglu)                        :: i,dyear,dmonth,dday,daysPast


    dyear =year -offsetYear
    dmonth=month-offsetMonth
    dday  =day  -offsetDay

    daysPast=offsetDayOfWeek

    daysPast=daysPast+dyear*365
    do i = 1,dmonth
        daysPast=daysPast+monthDays(i)
    enddo
    daysPast=daysPast+dday

    do i = offsetYear,year
        daysPast=daysPast+isLeapYear(i)
    enddo

    if ( (isLeapYear(year).EQ.1) .AND. ( (month.LT.2) .OR. ((month.EQ.2).AND.(day.LE.29)) ) ) then
        daysPast=daysPast-1
    endif

    if (present(past)) past=daysPast
    ret=mod(daysPast,7); if (ret.EQ.0) ret=7

    return
    end function dayOfWeek

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function convertTime(secs,spec) result(ret)
    implicit none
    character (len=25)                        :: ret
    real   (kind=rglu), intent(in)            :: secs
    logical(kind=lglu), optional, intent(in)  :: spec

    integer(kind=iglu)                        :: days,hours,minutes,seconds,mseconds,inLen
    character (len=4)                         :: dayid
    logical(kind=lglu)                        :: uspec


    ret=repeat(' ',len(ret))
    uspec=false
    if (present(spec)) uspec=spec

    days    =int( secs/86400. )
    hours   =int( secs/3600. - days*24. )
    minutes =int( secs/60. - hours*60. - days*1440. )
    seconds =int( secs - minutes*60. - hours*3600. - days*86400. )
    mseconds=int( 100.*(secs-int(secs)) )

    if (uspec) then
        hours=hours+days*24
        inLen=mid(hours); if (hours.LT.2) inLen=2
        write (ret,100) hours,minutes,seconds,mseconds
        100 format (i<inLen>.2,':',i2.2,':',i2.2,'.',i2.2)
    else
        dayid='day '
        if ((days.GT.1).OR.(days.EQ.0)) dayid='days'
        write (ret,101) days,dayid,hours,minutes,seconds,mseconds
        101 format (i<mid(days)>,1X,A4,1X,i2.2,':',i2.2,':',i2.2,'.',i2.2)
    endif

    return
    end function convertTime

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function date_time() result(ret)
    implicit none

    character (len=25)   :: ret
    integer(kind=i4kind) :: values(8)


    values=0; call date_and_time(values=values)

    write (ret,100) days(dayOfWeek(values(1),values(2),values(3))),&
                    values(1),months(values(2)),values(3),values(5:7)

100 format (A3,1X,i4,'-',A3,'-',i2.2,2X,i2.2,':',i2.2,':',i2.2)

    return
    end function date_time

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function timeStamp() result(ret)
    implicit none
    character (len=19)   :: ret
    integer(kind=i4kind) :: values(8)


    values=0; call date_and_time(values=values)

    !VALUE(1): The year
    !VALUE(2): The month
    !VALUE(3): The day of the month
    !VALUE(4): Time difference with UTC in minutes
    !VALUE(5): The hour of the day
    !VALUE(6): The minutes of the hour
    !VALUE(7): The seconds of the minute
    !VALUE(8): The milliseconds of the second

    write (ret,100) values(1:3),values(5:7)

100 format (i4,'.',i2.2,'.',i2.2,1X,i2.2,':',i2.2,':',i2.2)

    return
    end function timeStamp

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OPENMP FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function setThreadsNumber(n_threads) result(ret)
    implicit none

    integer(kind=iglu) :: n_threads


#   if(__opnmp==1)
        call omp_set_num_threads(n_threads)
#   endif

    ret=0; return
    end function setThreadsNumber

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function getThreadNumber() result(ret)
    implicit none


    ret=1
#   if(__opnmp==1)
        ret=omp_get_thread_num()+1
#   endif

    return
    end function getThreadNumber

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAY ANALIZE SUBROUTINE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine arrayAnalize1(arr,distributionWidth,unit)
    implicit none

    integer(kind=iglu), parameter   :: bnd=1
    integer(kind=iglu), allocatable :: pntDistribution(:)

    real   (kind=rglu)   :: arr(:)
    real   (kind=rglu)   :: ul,dl
    integer(kind=iglu)   :: distributionWidth,unit
    integer(kind=iglu)   :: bounds(2,bnd)
    integer(kind=iglu)   :: a ,i,j,k,l, m,n,o,p
    integer(kind=i8kind) :: combinations


    do k = 1,bnd
        bounds(1,k)=LBound(arr,k); bounds(2,k)=UBound(arr,k)
    enddo

    allocate (pntDistribution(0:distributionWidth))
    do a = 0,distributionWidth
        pntDistribution(a)=0; ul=real(10,rglu)**real(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=-1

        do i = bounds(1,1),bounds(2,1)
            if ( (abs(arr(i)).GT.dl) .AND. (abs(arr(i)).LE.ul) ) then
                pntDistribution(a)=pntDistribution(a)+1
            endif
        enddo
    enddo

    combinations=1
    do a = 1,bnd
        combinations=combinations*(bounds(2,a)-bounds(1,a)+1)
    enddo

    write (unit,100) combinations
    do a = 0,distributionWidth
        ul=real(10,rglu)**(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=0
        write (unit,200) ul,dl,pntDistribution(a)
    enddo
    write (unit,*)

    deallocate (pntDistribution)

100 format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200 format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

    return
    end subroutine arrayAnalize1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine arrayAnalize2(arr,distributionWidth,unit)
    implicit none

    integer(kind=iglu), parameter   :: bnd=2
    integer(kind=iglu), allocatable :: pntDistribution(:)

    real   (kind=rglu)   :: arr(:,:)
    real   (kind=rglu)   :: ul,dl
    integer(kind=iglu)   :: distributionWidth,unit
    integer(kind=iglu)   :: bounds(2,bnd)
    integer(kind=iglu)   :: a ,i,j,k,l, m,n,o,p
    integer(kind=i8kind) :: combinations


    do k = 1,bnd
        bounds(1,k)=LBound(arr,k); bounds(2,k)=UBound(arr,k)
    enddo

    allocate (pntDistribution(0:distributionWidth))
    do a = 0,distributionWidth
        pntDistribution(a)=0; ul=real(10,rglu)**real(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=-1

        do i = bounds(1,1),bounds(2,1)
        do j = bounds(1,2),bounds(2,2)
            if ( (abs(arr(i,j)).GT.dl) .AND. (abs(arr(i,j)).LE.ul) ) then
                pntDistribution(a)=pntDistribution(a)+1
            endif
        enddo
        enddo
    enddo

    combinations=1
    do a = 1,bnd
        combinations=combinations*(bounds(2,a)-bounds(1,a)+1)
    enddo

    write (unit,100) combinations
    do a = 0,distributionWidth
        ul=real(10,rglu)**(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=0
        write (unit,200) ul,dl,pntDistribution(a)
    enddo
    write (unit,*)

    deallocate (pntDistribution)

100 format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200 format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

    return
    end subroutine arrayAnalize2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine arrayAnalize3(arr,distributionWidth,unit)
    implicit none

    integer(kind=iglu), parameter   :: bnd=3
    integer(kind=iglu), allocatable :: pntDistribution(:)

    real   (kind=rglu)   :: arr(:,:,:)
    real   (kind=rglu)   :: ul,dl
    integer(kind=iglu)   :: distributionWidth,unit
    integer(kind=iglu)   :: bounds(2,bnd)
    integer(kind=iglu)   :: a ,i,j,k,l, m,n,o,p
    integer(kind=i8kind) :: combinations


    do k = 1,bnd
        bounds(1,k)=LBound(arr,k); bounds(2,k)=UBound(arr,k)
    enddo

    allocate (pntDistribution(0:distributionWidth))
    do a = 0,distributionWidth
        pntDistribution(a)=0; ul=real(10,rglu)**real(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=-1

        do i = bounds(1,1),bounds(2,1)
        do j = bounds(1,2),bounds(2,2)
        do k = bounds(1,3),bounds(2,3)
            if ( (abs(arr(i,j,k)).GT.dl) .AND. (abs(arr(i,j,k)).LE.ul) ) then
                pntDistribution(a)=pntDistribution(a)+1
            endif
        enddo
        enddo
        enddo
    enddo

    combinations=1
    do a = 1,bnd
        combinations=combinations*(bounds(2,a)-bounds(1,a)+1)
    enddo

    write (unit,100) combinations
    do a = 0,distributionWidth
        ul=real(10,rglu)**(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=0
        write (unit,200) ul,dl,pntDistribution(a)
    enddo
    write (unit,*)

    deallocate (pntDistribution)

100 format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200 format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

    return
    end subroutine arrayAnalize3

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine arrayAnalize4(arr,distributionWidth,unit)
    implicit none

    integer(kind=iglu), parameter   :: bnd=4
    integer(kind=iglu), allocatable :: pntDistribution(:)

    real   (kind=rglu)   :: arr(:,:,:,:)
    real   (kind=rglu)   :: ul,dl
    integer(kind=iglu)   :: distributionWidth,unit
    integer(kind=iglu)   :: bounds(2,bnd)
    integer(kind=iglu)   :: a ,i,j,k,l, m,n,o,p
    integer(kind=i8kind) :: combinations


    do k = 1,bnd
        bounds(1,k)=LBound(arr,k); bounds(2,k)=UBound(arr,k)
    enddo

    allocate (pntDistribution(0:distributionWidth))
    do a = 0,distributionWidth
        pntDistribution(a)=0; ul=real(10,rglu)**real(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=-1

        do i = bounds(1,1),bounds(2,1)
        do j = bounds(1,2),bounds(2,2)
        do k = bounds(1,3),bounds(2,3)
        do l = bounds(1,4),bounds(2,4)
            if ( (abs(arr(i,j,k,l)).GT.dl) .AND. (abs(arr(i,j,k,l)).LE.ul) ) then
                pntDistribution(a)=pntDistribution(a)+1
            endif
        enddo
        enddo
        enddo
        enddo
    enddo

    combinations=1
    do a = 1,bnd
        combinations=combinations*(bounds(2,a)-bounds(1,a)+1)
    enddo

    write (unit,100) combinations
    do a = 0,distributionWidth
        ul=real(10,rglu)**(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=0
        write (unit,200) ul,dl,pntDistribution(a)
    enddo
    write (unit,*)

    deallocate (pntDistribution)

100 format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200 format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

    return
    end subroutine arrayAnalize4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine arrayAnalize5(arr,distributionWidth,unit)
    implicit none

    integer(kind=iglu), parameter   :: bnd=5
    integer(kind=iglu), allocatable :: pntDistribution(:)

    real   (kind=rglu)   :: arr(:,:,:,:,:)
    real   (kind=rglu)   :: ul,dl
    integer(kind=iglu)   :: distributionWidth,unit
    integer(kind=iglu)   :: bounds(2,bnd)
    integer(kind=iglu)   :: a ,i,j,k,l, m,n,o,p
    integer(kind=i8kind) :: combinations


    do k = 1,bnd
        bounds(1,k)=LBound(arr,k); bounds(2,k)=UBound(arr,k)
    enddo

    allocate (pntDistribution(0:distributionWidth))
    do a = 0,distributionWidth
        pntDistribution(a)=0; ul=real(10,rglu)**real(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=-1

        do i = bounds(1,1),bounds(2,1)
        do j = bounds(1,2),bounds(2,2)
        do k = bounds(1,3),bounds(2,3)
        do l = bounds(1,4),bounds(2,4)
        do m = bounds(1,5),bounds(2,5)
            if ( (abs(arr(i,j,k,l,m)).GT.dl) .AND. (abs(arr(i,j,k,l,m)).LE.ul) ) then
                pntDistribution(a)=pntDistribution(a)+1
            endif
        enddo
        enddo
        enddo
        enddo
        enddo
    enddo

    combinations=1
    do a = 1,bnd
        combinations=combinations*(bounds(2,a)-bounds(1,a)+1)
    enddo

    write (unit,100) combinations
    do a = 0,distributionWidth
        ul=real(10,rglu)**(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=0
        write (unit,200) ul,dl,pntDistribution(a)
    enddo
    write (unit,*)

    deallocate (pntDistribution)

100 format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200 format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

    return
    end subroutine arrayAnalize5

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine arrayAnalize6(arr,distributionWidth,unit)
    implicit none

    integer(kind=iglu), parameter   :: bnd=6
    integer(kind=iglu), allocatable :: pntDistribution(:)

    real   (kind=rglu)   :: arr(:,:,:,:,:,:)
    real   (kind=rglu)   :: ul,dl
    integer(kind=iglu)   :: distributionWidth,unit
    integer(kind=iglu)   :: bounds(2,bnd)
    integer(kind=iglu)   :: a ,i,j,k,l, m,n,o,p
    integer(kind=i8kind) :: combinations


    do k = 1,bnd
        bounds(1,k)=LBound(arr,k); bounds(2,k)=UBound(arr,k)
    enddo

    allocate (pntDistribution(0:distributionWidth))
    do a = 0,distributionWidth
        pntDistribution(a)=0; ul=real(10,rglu)**real(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=-1

        do i = bounds(1,1),bounds(2,1)
        do j = bounds(1,2),bounds(2,2)
        do k = bounds(1,3),bounds(2,3)
        do l = bounds(1,4),bounds(2,4)
        do m = bounds(1,5),bounds(2,5)
        do n = bounds(1,6),bounds(2,6)
            if ( (abs(arr(i,j,k,l,m,n)).GT.dl) .AND. (abs(arr(i,j,k,l,m,n)).LE.ul) ) then
                pntDistribution(a)=pntDistribution(a)+1
            endif
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
    enddo

    combinations=1
    do a = 1,bnd
        combinations=combinations*(bounds(2,a)-bounds(1,a)+1)
    enddo

    write (unit,100) combinations
    do a = 0,distributionWidth
        ul=real(10,rglu)**(-a); dl=real(10,rglu)**(-a-1)
        if (a.EQ.distributionWidth) dl=0
        write (unit,200) ul,dl,pntDistribution(a)
    enddo
    write (unit,*)

    deallocate (pntDistribution)

100 format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200 format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

    return
    end subroutine arrayAnalize6

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function compare_r16(a,b,add) result(ret)
    implicit none
    real   (kind=r16kind), intent(in)           :: a,b
    integer(kind=iglu)   , intent(in), optional :: add
    real   (kind=rglu)                          :: wa,wb,d
    integer(kind=iglu)                          :: power


    wa=abs(a); wb=abs(b)
    if (present(add)) then
        power=int(log10(float(int(max(wa,wb)))))+1+add
    else
        power=int(log10(float(int(max(wa,wb)))))+1
    endif
    d=epsilon(a)*10**(power-1)

    ret=abs(a-b).LT.d; return
    end function compare_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function compare_r8(a,b,add) result(ret)
    implicit none
    real   (kind=r8kind), intent(in)           :: a,b
    integer(kind=iglu)  , intent(in), optional :: add
    real   (kind=rglu)                         :: wa,wb,d
    integer(kind=iglu)                         :: power


    wa=abs(a); wb=abs(b)
    if (present(add)) then
        power=int(log10(float(int(max(wa,wb)))))+1+add
    else
        power=int(log10(float(int(max(wa,wb)))))+1
    endif
    d=epsilon(a)*10**(power-1)

    ret=abs(a-b).LT.d; return
    end function compare_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical(kind=lglu) function compare_r4(a,b,add) result(ret)
    implicit none
    real   (kind=r4kind), intent(in)           :: a,b
    integer(kind=iglu)  , intent(in), optional :: add
    real   (kind=rglu)                         :: wa,wb,d
    integer(kind=iglu)                         :: power


    wa=abs(a); wb=abs(b)
    if (present(add)) then
        power=int(log10(float(int(max(wa,wb)))))+1+add
    else
        power=int(log10(float(int(max(wa,wb)))))+1
    endif
    d=epsilon(a)*10**(power-1)

    ret=abs(a-b).LT.d; return
    end function compare_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer(kind=iglu) function compareStrings(fstr,sstr) result(ret)
    implicit none

    character (len=*)  , intent(in) :: fstr,sstr
    integer(kind=iglu)              :: k


    do k = 1,min(len(fstr),len(sstr))

        if     ( iachar(fstr(k:k)) .GT. iachar(sstr(k:k)) ) then
            ret=1
            exit
        elseif ( iachar(fstr(k:k)) .EQ. iachar(sstr(k:k)) ) then
            ret=0
            cycle
        else
            ret=-1
            exit
        endif

    enddo

    if     ( (ret.EQ.0) .AND. (len(fstr).GT.len(sstr)) ) then
        ret=1
    elseif ( (ret.EQ.0) .AND. (len(fstr).EQ.len(sstr)) ) then
        ret=0
    elseif ( (ret.EQ.0) .AND. (len(fstr).LT.len(sstr)) ) then
        ret=-1
    endif

    return
    end function compareStrings

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function screenProgress(string) result(ret)
    implicit none

    character (len=*) :: string


    write (*, '(A\)') char(13)//string

    ret=0; return
    end function screenProgress

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine infiniteVoidLoop()
    implicit none


    do
        continue
    enddo

    end subroutine infiniteVoidLoop

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine nullSub
    implicit none


    continue

    return
    end subroutine nullSub

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_i8(var) result(ret)
    implicit none
    integer(kind=i8kind) :: var


    ret=0; return
    end function call_wrapper_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_i4(var) result(ret)
    implicit none
    integer(kind=i4kind) :: var


    ret=0; return
    end function call_wrapper_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_i2(var) result(ret)
    implicit none
    integer(kind=i2kind) :: var


    ret=0; return
    end function call_wrapper_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_i1(var) result(ret)
    implicit none
    integer(kind=i1kind) :: var


    ret=0; return
    end function call_wrapper_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_r16(var) result(ret)
    implicit none
    real(kind=r16kind) :: var


    ret=0; return
    end function call_wrapper_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_r8(var) result(ret)
    implicit none
    real(kind=r8kind) :: var


    ret=0; return
    end function call_wrapper_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_r4(var) result(ret)
    implicit none
    real(kind=r4kind) :: var


    ret=0; return
    end function call_wrapper_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_l8(var) result(ret)
    implicit none
    logical(kind=l8kind) :: var


    ret=0; return
    end function call_wrapper_l8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_l4(var) result(ret)
    implicit none
    logical(kind=l4kind) :: var


    ret=0; return
    end function call_wrapper_l4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_l2(var) result(ret)
    implicit none
    logical(kind=l2kind) :: var


    ret=0; return
    end function call_wrapper_l2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_l1(var) result(ret)
    implicit none
    logical(kind=l1kind) :: var


    ret=0; return
    end function call_wrapper_l1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_ch(var) result(ret)
    implicit none
    character (len=*) :: var


    ret=0; return
    end function call_wrapper_ch

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function call_wrapper_uc(var) result(ret)
    implicit none
    type(uch) :: var


    ret=0; return
    end function call_wrapper_uc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function definePi() ! define pi and degr->radian koef.
    implicit none


    pi=2*acos(gluZero); dtr=pi/180_r8kind

    definePi=0; return
    end function definePi

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function purifyValues(threshold,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) result(ret)
    implicit none

    real(kind=rglu), optional :: threshold,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
    real(kind=rglu)           :: dthreshold


    dthreshold=1E-10_rglu; if (present(threshold)) dthreshold=threshold

    if (present(a1 )) then; if (abs(a1 ).LT.dthreshold) a1 =0; endif
    if (present(a2 )) then; if (abs(a2 ).LT.dthreshold) a2 =0; endif
    if (present(a3 )) then; if (abs(a3 ).LT.dthreshold) a3 =0; endif
    if (present(a4 )) then; if (abs(a4 ).LT.dthreshold) a4 =0; endif
    if (present(a5 )) then; if (abs(a5 ).LT.dthreshold) a5 =0; endif
    if (present(a6 )) then; if (abs(a6 ).LT.dthreshold) a6 =0; endif
    if (present(a7 )) then; if (abs(a7 ).LT.dthreshold) a7 =0; endif
    if (present(a8 )) then; if (abs(a8 ).LT.dthreshold) a8 =0; endif
    if (present(a9 )) then; if (abs(a9 ).LT.dthreshold) a9 =0; endif
    if (present(a10)) then; if (abs(a10).LT.dthreshold) a10=0; endif

    ret=0; return
    end function purifyValues

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function getPath() result(ret)
    implicit none

    type(uch)           :: ret
    character (len=512) :: fpath,rpath,drive
    character (len=1)   :: cnull


    call getarg(0,fpath)
    void=splitpathqq(fpath,drive,rpath,cnull,cnull)

    ret=trim(drive)//trim(rpath)

    return
    end function getPath

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function glGetIOunit() result(ret)
    implicit none


    ret=iounit; return
    end function glGetIOunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine glSetIOunit(iunt)
    implicit none

    integer(kind=iglu), intent(in) :: iunt
    integer(kind=iglu)             :: unt


    unt=iunt
    if (unt.EQ.5) unt=6 !stdin  (forbidden)
    if (unt.EQ.0) unt=6 !stdout (for compatibility)
    iounit=unt

    return
    end subroutine glSetIOunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine glFinalize
    implicit none


    iounit=6

    return
    end subroutine glFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module glob
