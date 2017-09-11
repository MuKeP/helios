	!DEC$if defined(__unix)
		!DEC$define OS=2
	!DEC$else
		!DEC$define OS=1
	!DEC$endif

	!DEC$if defined(__INTEL_COMPILER)
		!DEC$define compiler=1
	!DEC$elseif defined(__GFORTRAN__)
		!DEC$define compiler=2
	!DEC$elseif defined(__SUNPRO_F90)
		!DEC$define compiler=3
	!DEC$elseif defined(__SUNPRO_F95)
		!DEC$define compiler=3
	!DEC$else
		!DEC$define compiler=4
	!DEC$endif
	
	!DEC$if defined(_OPENMP)
		!DEC$define opnmp=1
	!DEC$else
		!DEC$define opnmp=0
	!DEC$endif

	module glob

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	!DEC$if (compiler.EQ.1) ! Intel Fortran Compiler
		use ifport , only: signal,system,getpid,splitpathqq
		use ifport , only: SRT$REAL8,SRT$REAL4
		use ifport , only: SRT$INTEGER8,SRT$INTEGER4,SRT$INTEGER2,SRT$INTEGER1
		use ifport , only: SIGABRT,SIGINT,SIGTERM
		use ifposix, only: PXFSTRUCTCREATE
		!DEC$if (os.EQ.2)
			use ifposix, only: PXFSIGADDSET
		!DEC$endif
	!DEC$elseif (compiler.EQ.2)
		!
	!DEC$elseif (compiler.EQ.3)
		!
	!DEC$elseif (compiler.EQ.4) ! Compaq© Visual Fortran
		use DFPort, only: signal,system,getpid
		use DFLib , only: splitpathqq
		use DFLib , only: SRT$REAL8,SRT$REAL4
		use DFLib , only: SRT$INTEGER4,SRT$INTEGER2,SRT$INTEGER1
		use DFPort, only: SIGABRT,SIGINT,SIGTERM
	!DEC$endif

	!MS$if(opnmp.EQ.1)
		include "omp_lib.h"
	!MS$endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: glVersion='4.240'
	character (len=*), parameter :: glDate   ='2017.09.03'
	character (len=*), parameter :: glAuthor ='Anton B. Zakharov'

	integer*4, parameter :: r16kind=16, r8kind=8, r4kind=4
	integer*4, parameter :: i16kind=16, i8kind=8, i4kind=4, i2kind=2, i1kind=1
	integer*4, parameter ::             l8kind=8, l4kind=4, l2kind=2, l1kind=1

	public :: r16kind, r8kind, r4kind
	public :: i16kind, i8kind, i4kind, i2kind, i1kind
	public ::          l8kind, l4kind, l2kind, l1kind

	! glu = real for global use  (storage and non-accuracy-demanding procedures)
	! spu = real for special use (accuracy-demanding procedures)
	integer(kind=i4kind), parameter :: rglu=r8kind, rspu=r8kind
	integer(kind=i4kind), parameter :: iglu=i4kind, ispu=i8kind
	integer(kind=i4kind), parameter :: lglu=l1kind

	!   ~~~~ Global data settings ~~~~ !
	!MS$if(OS.EQ.1)
		character (len=*), parameter :: os='win',osSeparator='\',osMove='move',osCopy='copy'
	!MS$else
		character (len=*), parameter :: os='nix',osSeparator='/',osMove='mv'  ,osCopy='cp'
	!MS$endif

	real(kind=rglu), parameter :: gluZero=0,gluUnity=1
	real(kind=rspu), parameter :: spuZero=0,spuUnity=1

	!tolerance in comparison of two real variables
	real(kind=rglu), parameter :: gluCompare=real(10,kind=rglu)**floor(log10(epsilon(gluUnity))+real(1,rglu))
	real(kind=rspu), parameter :: spuCompare=real(10,kind=rspu)**floor(log10(epsilon(spuUnity))+real(1,rspu))

	logical(kind=lglu), parameter :: true=.true., false=.false.

	character (len=*) , parameter :: months(12)=['Jan','Feb','Mar','Apr','May','Jun',&
	                                             'Jul','Aug','Sep','Oct','Nov','Dec']
	integer(kind=iglu), parameter :: monthDays(12)=[31,28,31,30,31,30,31,31,30,31,30,31]

	character (len=*) , parameter :: days(7)   =['Mon','Tue','Wed','Thu','Fri','Sat','Sun']

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type uch
		character (len=1) , allocatable :: ch(:)
		integer(kind=iglu)              :: ln
	end type uch

	type ivarVector
		integer(kind=iglu), allocatable :: v(:)
		integer(kind=iglu)              :: ln
	end type ivarVector

	type rvarVector
		real   (kind=rglu), allocatable :: v(:)
		integer(kind=iglu)              :: ln
	end type rvarVector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) :: defaultIO=6
	real   (kind=rglu) :: pi,dtr
	character          :: glPrintString*500
	integer(kind=iglu) :: iounit=6

	integer(kind=iglu) :: void
	real   (kind=rglu) :: voidr
	character          :: voidc*1
	logical(kind=lglu) :: voidl

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	interface uchGet
		module procedure uchGet_def,uchGet_opt
	end interface uchGet

	interface find
		module procedure find_c,find_uc,find_r8,find_r4,find_i8,find_i4,find_i2,find_i1
	end interface find

	interface isEven
		module procedure isEven_i8,isEven_i4,isEven_i2,isEven_i1
	end interface isEven

	interface same
		module procedure same_i8,same_i4,same_i2,same_i1
	end interface same

	interface mid
		module procedure mid_i8,mid_i4,mid_i2,mid_i1,mid_r8,&
		                 mid_r4,mid_l8,mid_l4,mid_l2,mid_l1,&
		                 mid_c
	end interface mid

	interface lenTrimArray
		module procedure lenTrimArray_i8,lenTrimArray_i4,lenTrimArray_i2,&
		                 lenTrimArray_i1,lenTrimArray_r8,lenTrimArray_r4
	end interface lenTrimArray

	interface collectArray
		module procedure collectArray_i8,collectArray_i4,collectArray_i2,&
		                 collectArray_i1,collectArray_r8,collectArray_r4
	end interface collectArray

	interface rangen
		module procedure random_generate_i8,random_generate_i4,&
		                 random_generate_r8,random_generate_r4
	end interface rangen

	interface sort
		module procedure sort_uch,sort_r8,sort_r4,sort_i4,sort_i2,sort_i1
	end interface sort

	interface arrayAnalize
		module procedure arrayAnalize1,arrayAnalize2,arrayAnalize3,&
		                 arrayAnalize4,arrayAnalize5,arrayAnalize6
	end interface arrayAnalize

	interface delVector
		module procedure delVector_r,delVector_i
	end interface delVector

	interface setVector
		module procedure setVector_r,setVector_i
	end interface setVector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private

	! module attributes
	public :: glVersion,glDate,glAuthor
	public :: glSetIOunit
	public :: glFinalize

	! data type presets
	public :: iglu,ispu,lglu
	public :: rglu,gluZero,gluUnity,gluCompare
	public :: rspu,spuZero,spuUnity,spuCompare

	! types
	public :: uch,ivarVector,rvarVector
	! type "methods"
	public :: uchGet,uchSet,uchDel,delVector,setVector

	! platform constants
	public :: os,osSeparator,osMove,osCopy
	! platform functions
	public :: getPath

	! signals
	public :: signal,system,getpid
	public :: SIGABRT,SIGINT,SIGTERM,PXFSTRUCTCREATE,PXFSIGADDSET

	! OpenMP functions
	public :: setThreadsNumber,getThreadNumber

	! sort function types
	public :: SRT$REAL8,SRT$REAL4
	public :: SRT$INTEGER4,SRT$INTEGER2,SRT$INTEGER1

	! variables
	public :: defaultIO,glPrintString
	public :: true,false,void,voidr,voidc,voidl,pi,dtr

	! tools
	public :: find,isEven,same,mid,lenTrimArray,collectArray,rangen,sort
	public :: compareStrings,arrayAnalize,definePi
	! almost useless
	public :: infiniteVoidLoop,nullSub

	! time and date
	public :: convertTime,timeControl,dayOfWeek,isLeapYear,timeStamp,date_time

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ UCH METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	pure function uchGet_def(uchar) result(ret)
	implicit none

	type(uch), intent(in)    :: uchar
	character (len=uchar%ln) :: ret
	integer(kind=iglu)       :: i


	if (uchar%ln.EQ.0) return
	ret=uchar%ch(1)
	do i = 1,uchar%ln-1
		ret=ret(:i)//uchar%ch(i+1)
	enddo

	return
	end function uchGet_def

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function uchGet_opt(uchar,ustart,uend) result(ret)
	implicit none

	type(uch)         , intent(in) :: uchar
	integer(kind=iglu), intent(in) :: ustart,uend
	character (len=uend-ustart+1)  :: ret
	integer(kind=iglu)             :: i


	if (uchar%ln.EQ.0) return
	ret=uchar%ch(ustart)
	do i = ustart,uend-1
		ret=ret(:i-ustart+1)//uchar%ch(i+1)
	enddo

	return
	end function uchGet_opt

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function uchSet(str,ch) result(ret)
	implicit none

	character (len=*), intent(in)    :: str
	type(uch), optional              :: ch
	type(uch)                        :: ret
	integer(kind=iglu)               :: i


	if (present(ch)) deallocate (ch%ch, stat=i)
	allocate ( ret%ch(len(str)) )
	ret%ln=len(str)

	do i = 1,ret%ln
		ret%ch(i)=str(i:i)
	enddo

	return
	end function uchSet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function uchDel(ch) result(ret)
	implicit none

	type(uch)          :: ch
	integer(kind=iglu) :: ret


	ret=0; deallocate (ch%ch, stat=ret); ch%ln=0

	return
	end function uchDel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VAR VECTOR METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function setVector_r(vec) result(ret)
	type(rvarVector) :: ret
	real(kind=rglu)  :: vec(:)


	ret%ln=UBound(vec,1)

	allocate (ret%v(ret%ln))
	ret%v=vec

	return
	end function setVector_r

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function setVector_i(vec) result(ret)
	type(ivarVector)   :: ret
	integer(kind=iglu) :: vec(:)


	ret%ln=UBound(vec,1)

	allocate (ret%v(ret%ln))
	ret%v=vec

	return
	end function setVector_i

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	elemental subroutine delVector_r(vector)
	implicit none

	type(rvarVector), intent(out) :: vector
	integer(kind=iglu)            :: err

		
	deallocate (vector%v, stat=err)

	return
	end subroutine delVector_r

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	elemental subroutine delVector_i(vector)
	implicit none

	type(ivarVector), intent(inout) :: vector
	integer(kind=iglu)              :: err

		
	deallocate (vector%v, stat=err)

	return
	end subroutine delVector_i

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

	pure integer(kind=iglu) function mid_r8(myReal) result(ret)
	implicit none
	real(kind=r8kind), intent(in) :: myReal


	if (abs(myReal).LT.1d0) then
		if (myReal.LT.0.d0) then
			ret=2
		else
			ret=1
		endif
		return
	endif
	if (myReal.LT.0) then
		ret=int( log10(abs(myReal))+1 )+1
		return
	endif
	ret=int( log10(myReal)+1 ); return
	end function mid_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	pure integer(kind=iglu) function mid_r4(myReal) result(ret)
	implicit none
	real(kind=r4kind), intent(in) :: myReal


	if (abs(myReal).LT.1e0) then
		if (myReal.LT.0e0) then
			ret=2
		else
			ret=1
		endif
		return
	endif
	if (myReal.LT.0) then
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

	pure integer(kind=iglu) function mid_c(str) result(ret)
	implicit none
	character (len=*), intent(in) :: str


	ret=len(str); return
	end function mid_c

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

	pure integer(kind=iglu) function lenTrimArray_r8(arr,tol) result(ret)
	implicit none

	real   (kind=r8kind), intent(in)           :: arr(:)
	real   (kind=r8kind), intent(in), optional :: tol
	real   (kind=r8kind)                       :: deftol
	integer(kind=iglu)                         :: ln,i


	ln=UBound(arr,1)

	deftol=epsilon(deftol)*10; if (present(tol)) deftol=tol

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

	deftol=epsilon(deftol)*10; if (present(tol)) deftol=tol

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

	integer(kind=iglu) function collectArray_r8(arr,tol) result(ret)
	implicit none
	real(kind=r8kind)           :: deftol,arr(:),hold
	real(kind=r8kind), optional :: tol
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
			if (rnum.GT.0.5d0) then
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
			if (rnum.GT.0.5d0) then
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
			if (rnum.GT.0.5d0) then
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
			if (rnum.GT.0.5d0) then
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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SORT FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function sort_r8(vector,rev) result(ret)
	implicit none
	real   (kind=r8kind)           :: vector(:),swap
	logical(kind=iglu)  , optional :: rev
	logical(kind=lglu)             :: urev
	integer(kind=iglu)             :: send,k,vlen


	urev=false; if (present(rev)) urev=rev

	send=ubound(vector,1)
	call sortqq(loc(vector),send,SRT$REAL8)
	ret=send

	if (urev) then
		vlen=ubound(vector,1)
		do k = 1,vlen/2
			swap=vector(k); vector(k)=vector(vlen-k+1); vector(vlen-k+1)=swap
		enddo		
	endif

	return
	end function sort_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function sort_r4(vector,rev) result(ret)
	implicit none
	real   (kind=r4kind)           :: vector(:),swap
	logical(kind=iglu)  , optional :: rev
	logical(kind=lglu)             :: urev
	integer(kind=iglu)             :: send,k,vlen


	urev=false; if (present(rev)) urev=rev

	send=UBound(vector,1)
	call sortqq(loc(vector),send,SRT$REAL4)
	ret=send

	if (urev) then
		vlen=UBound(vector,1)
		do k = 1,vlen/2
			swap=vector(k); vector(k)=vector(vlen-k+1); vector(vlen-k+1)=swap
		enddo		
	endif

	return
	end function sort_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function sort_i4(vector,rev) result(ret)
	implicit none
	integer(kind=i4kind)           :: vector(:),swap
	logical(kind=lglu)  , optional :: rev
	logical(kind=lglu)             :: urev
	integer(kind=iglu)             :: send,k,vlen


	urev=false; if (present(rev)) urev=rev

	send=UBound(vector,1)
	call sortqq(loc(vector),send,SRT$INTEGER4)
	ret=send

	if (urev) then
		vlen=UBound(vector,1)
		do k = 1,vlen/2
			swap=vector(k)
			vector(k)=vector(vlen-k+1)
			vector(vlen-k+1)=swap
		enddo		
	endif

	return
	end function sort_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function sort_i2(vector,rev) result(ret)
	implicit none
	integer(kind=i2kind)           :: vector(:),swap
	logical(kind=lglu)  , optional :: rev
	logical(kind=lglu)             :: urev
	integer(kind=iglu)             :: send,k,vlen


	urev=false; if (present(rev)) urev=rev

	send=UBound(vector,1)
	call sortqq(loc(vector),send,SRT$INTEGER2)
	ret=send

	if (urev) then
		vlen=UBound(vector,1)
		do k = 1,vlen/2
			swap=vector(k)
			vector(k)=vector(vlen-k+1)
			vector(vlen-k+1)=swap
		enddo		
	endif

	return
	end function sort_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function sort_i1(vector,rev) result(ret)
	implicit none
	integer(kind=i1kind)           :: vector(:),swap
	logical(kind=lglu)  , optional :: rev
	logical(kind=lglu)             :: urev
	integer(kind=iglu)             :: send,k,vlen


	urev=false; if (present(rev)) urev=rev

	send=UBound(vector,1)
	call sortqq(loc(vector),send,SRT$INTEGER1)
	ret=send

	if (urev) then
		vlen=UBound(vector,1)
		do k = 1,vlen/2
			swap=vector(k)
			vector(k)=vector(vlen-k+1)
			vector(vlen-k+1)=swap
		enddo		
	endif

	return
	end function sort_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function sort_uch(vector,rev,attribute) result(ret)
	implicit none
	type(uch)                       :: vector(:),swap
	character (len=*), optional     :: attribute
	character (len=len(attribute))  :: uattribute
	logical(kind=lglu), optional    :: rev
	logical(kind=lglu)              :: urev,cnd
	integer(kind=iglu)              :: N,k,cnt


	urev=false; if (present(rev)) urev=rev
	uattribute=repeat(' ',len(uattribute)); uattribute='length'
	if (present(attribute)) uattribute=attribute
	N=UBound(vector,1)

	if ( (uattribute.NE.'alphabet').AND.(uattribute.NE.'length') ) then
		ret=-1
		return
	endif

	ret=0
	select case (uattribute) ! assume that vector is small, bubble sort.
		case ('length')
			do
				cnt=0
				do k = 1,N-1

					if (urev) then
						cnd=(vector(k)%ln.LT.vector(k+1)%ln)
					else
						cnd=(vector(k)%ln.GT.vector(k+1)%ln)
					endif

					if (cnd) then
						swap=uchSet( uchGet(vector(k)) )
						vector(k)=uchSet( uchGet(vector(k+1)) )
						vector(k+1)=uchSet( uchGet(swap) )
						cnt=cnt+1
					endif
				enddo
				ret=ret+cnt
				if (cnt.EQ.0) exit
			enddo

		case ('alphabet')
			do
				cnt=0
				do k = 1,N-1

					if (urev) then
						cnd=compareStrings( uchGet(vector(k)), uchGet(vector(k+1)) ) .LT. 0
					else
						cnd=compareStrings( uchGet(vector(k)), uchGet(vector(k+1)) ) .GT. 0
					endif

					if (cnd) then
						swap=uchSet( uchGet(vector(k)) )
						vector(k)=uchSet( uchGet(vector(k+1)) )
						vector(k+1)=uchSet( uchGet(swap) )
						cnt=cnt+1
					endif
				enddo
				ret=ret+cnt
				if (cnt.EQ.0) exit
			enddo

	end select

	return
	end function sort_uch

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TIME AND DATE FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu) function timeControl(cputime) result(ret)
	implicit none
		
	real(kind=rglu), optional :: cputime
		
		
	!MS$if(opnmp.EQ.1)
		ret=omp_get_wtime()
	!MS$else
		call cpu_time(ret)
	!MS$endif				
	if (present(cputime)) call cpu_time(cputime)
		
	return
	end function timeControl

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
	character (len=20)                        :: ret
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

100	format (A3,1X,i4,'-',A3,'-',i2.2,2X,i2.2,':',i2.2,':',i2.2)

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

100	format (i4,'.',i2.2,'.',i2.2,1X,i2.2,':',i2.2,':',i2.2)

	return
	end function timeStamp

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OPENMP FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function setThreadsNumber(n_threads) result(ret)
	implicit none

	integer(kind=iglu) :: n_threads


	!MS$if(opnmp.EQ.1)
		call omp_set_num_threads(n_threads)
	!MS$endif
		
	ret=0; return
	end function setThreadsNumber

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function getThreadNumber() result(ret)
	implicit none


	ret=1
	!MS$if(opnmp.EQ.1)
		ret=omp_get_thread_num()+1
	!MS$endif
		
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

100	format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200	format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

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

100	format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200	format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

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

100	format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200	format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

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

100	format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200	format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

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

100	format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200	format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

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

100	format (/4X,'Elements ',i<mid(combinations)>,'. Distribution:')
200	format ( 4X,ES9.3,' <--> ',ES9.3,i<mid(combinations)+1>,' elements')

	return
	end subroutine arrayAnalize6

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

	integer(kind=iglu) function definePi() ! define pi and degr->radian koef.
	implicit none


	pi=2*acos(gluZero); dtr=pi/180

	definePi=0; return
	end function definePi

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function getPath() result(ret)
	implicit none

	type(uch)           :: ret
	character (len=512) :: fpath,rpath,drive
	character (len=1)   :: cnull


	call getarg(0,fpath)
	void=splitpathqq(fpath,drive,rpath,cnull,cnull)

	ret=uchSet(trim(drive)//trim(rpath))

	return
	end function getPath

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
