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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	!DEC$if (compiler.EQ.1) ! Intel Fortran Compiler
		use ifport, only: signal,system,getpid
		use ifport, only: SRT$REAL8,SRT$REAL4
		use ifport, only: SRT$INTEGER8,SRT$INTEGER4,SRT$INTEGER2,SRT$INTEGER1
		use ifport, only: SIGABRT,SIGINT,SIGTERM
		!use ifcore
		!use ifposix
	!DEC$elseif (compiler.EQ.2)
		!
	!DEC$elseif (compiler.EQ.3)
		!
	!DEC$elseif (compiler.EQ.4) ! Compaq© Visual Fortran
		use DFLib   
	!DEC$endif

	!MS$if(opnmp.EQ.1)
		include "omp_lib.h"
	!MS$endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: glVersion='4.210'
	character (len=*), parameter :: glDate   ='2017.07.07'
	character (len=*), parameter :: glAuthor ='Anton B. Zakharov'

	integer*4, parameter :: r16kind=16, r8kind=8, r4kind=4
	integer*4, parameter :: i16kind=16, i8kind=8, i4kind=4, i2kind=2, i1kind=1
	integer*4, parameter ::             l8kind=8, l4kind=4, l2kind=2, l1kind=1

	public :: r16kind, r8kind, r4kind
	public :: i16kind, i8kind, i4kind, i2kind, i1kind
	public ::          l8kind, l4kind, l2kind, l1kind

	! glu = real for global use  (storage and non-accuracy-demanding procedures).
	! spu = real for special use (accuracy-demanding procedures).
	integer*4, parameter :: rglu=r8kind, rspu=r8kind
	integer*4, parameter :: iglu=i4kind, ispu=i8kind
	integer*4, parameter :: lglu=l1kind

	!   ~~~~ Global data settings ~~~~ !
	!MS$if(OS.EQ.1)
		character (len=*), parameter :: os='win',osSeparator='\',osMove='move',osCopy='copy'
	!MS$else
		character (len=*), parameter :: os='nix',osSeparator='/',osMove='mv  ',osCopy='cp  '
	!MS$endif

	real(kind=rglu), parameter :: gluUnity=real(1,rglu),gluZero=real(0,rglu)
	real(kind=rspu), parameter :: spuUnity=real(1,rspu),spuZero=real(0,rspu)

	logical(kind=lglu), parameter :: true=.true., false=.false.

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type uch
		character (len=1), allocatable :: ch(:)
		integer*4                      :: ln
	end type uch

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) :: defaultIO=6
	real   (kind=rglu) :: pi,dtr
	character          :: glPrintString*500
	integer(kind=iglu) :: iounit=6

	integer(kind=iglu) :: void
	real   (kind=rglu) :: voidr
	character          :: voidc*1
	logical(kind=lglu) :: voidl

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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
		module procedure igen, rgen
	end interface rangen

	interface sort
		module procedure sort_uch,sort_r8,sort_r4,sort_i4,sort_i2,sort_i1
	end interface sort

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private

	!DEC$if (compiler.EQ.1)
		public :: signal,system,getpid
		public :: SRT$REAL8,SRT$REAL4
		public :: SRT$INTEGER8,SRT$INTEGER4,SRT$INTEGER2,SRT$INTEGER1
		public :: SIGABRT,SIGINT,SIGTERM
	!DEC$endif

	public :: rglu,rspu,iglu,ispu,lglu,gluUnity,gluZero,spuUnity,spuZero

	public :: glVersion,glDate,glAuthor

	! types
	public :: uch
	public :: uchGet,uchSet,uchDel

	! variables
	public :: defaultIO
	public :: pi,dtr,void,voidr,voidc,voidl,glPrintString,true,false

	! polymorphic
	public :: find,isEven,same,mid,lenTrimArray,collectArray,rangen,sort

	! simple
	public :: definePi,convertTime,glSetIOunit,glFinalize
	public :: timeControl,compareStrings,dayOfWeek,isLeapYear,timeStamp

	! platform constants
	public :: os,osSeparator,osMove,osCopy,osCall

	! OpenMP
	public :: setThreadsNumber,getThreadNumber

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~uch functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure function uchGet_def(uchar) result(ret)
		implicit none

		type(uch), intent(in)    :: uchar
		character (len=uchar%ln) :: ret
		integer*4                :: i


		if (uchar%ln.EQ.0) return
		ret=uchar%ch(1)
		do i = 1,uchar%ln-1
			ret=ret(:i)//uchar%ch(i+1)
		enddo

		return
		end function uchGet_def

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		function uchGet_opt(uchar,ustart,uend) result(ret)
		implicit none

		type(uch), intent(in)         :: uchar
		integer*4, intent(in)         :: ustart,uend
		character (len=uend-ustart+1) :: ret
		integer*4                     :: i


		if (uchar%ln.EQ.0) return
		ret=uchar%ch(ustart)
		do i = ustart,uend-1
			ret=ret(:i-ustart+1)//uchar%ch(i+1)
		enddo

		return
		end function uchGet_opt

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		function uchSet(str,ch) result(ret)
		implicit none

		character (len=*), intent(in)    :: str
		type(uch), optional              :: ch
		type(uch)                        :: ret
		integer*4                        :: i


		if (present(ch)) deallocate (ch%ch, stat=i)
		allocate ( ret%ch(len(str)) )
		ret%ln=len(str)

		do i = 1,ret%ln
			ret%ch(i)=str(i:i)
		enddo

		return
		end function uchSet

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		function uchDel(ch) result(ret)
		implicit none

		type(uch) :: ch
		integer*4 :: ret


		ret=0; deallocate (ch%ch, stat=ret); ch%ln=0

		return
		end function uchDel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~mid functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_i8(myInt)  result(ret) ! my integer's digits (old name)
		implicit none
		integer*8, intent(in) :: myInt


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_i4(myInt) result(ret)
		implicit none
		integer*4, intent(in) :: myInt


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_i2(myInt) result(ret)
		implicit none
		integer*2, intent(in) :: myInt


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_i1(myInt) result(ret)
		implicit none
		integer*1, intent(in) :: myInt


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_r8(myReal) result(ret)
		implicit none
		real*8, intent(in) :: myReal


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_r4(myReal) result(ret)
		implicit none
		real*4, intent(in) :: myReal


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_l8(logi) result(ret)
		implicit none
		logical*8, intent(in) :: logi


		ret=5; if (logi) ret=4; return
		end function mid_l8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_l4(logi) result(ret)
		implicit none
		logical*4, intent(in) :: logi


		ret=5; if (logi) ret=4; return
		end function mid_l4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_l2(logi) result(ret)
		implicit none
		logical*2, intent(in) :: logi


		ret=5; if (logi) ret=4; return
		end function mid_l2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_l1(logi) result(ret)
		implicit none
		logical*1, intent(in) :: logi


		ret=5; if (logi) ret=4; return
		end function mid_l1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function mid_c(str) result(ret)
		implicit none
		character (len=*), intent(in) :: str


		ret=len(str); return
		end function mid_c

!   ~~~~~~~~~~~~~~~~~~~~~~~trim array functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function lenTrimArray_i8(arr) result(ret)
		implicit none
		integer*8, intent(in) :: arr(:)
		integer*4             :: ln,i


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function lenTrimArray_i4(arr) result(ret)
		implicit none
		integer*4, intent(in) :: arr(:)
		integer*4             :: ln,i


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function lenTrimArray_i2(arr) result(ret)
		implicit none
		integer*2, intent(in) :: arr(:)
		integer*4             :: ln,i


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function lenTrimArray_i1(arr) result(ret)
		implicit none
		integer*1, intent(in) :: arr(:)
		integer*4             :: ln,i


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function lenTrimArray_r8(arr,tol) result(ret)
		implicit none

		real*8, intent(in)           :: arr(:)
		real*8, intent(in), optional :: tol
		real*8                       :: deftol
		integer*4                    :: ln,i


		ln=UBound(arr,1)

		deftol=1D-10; if (present(tol)) deftol=tol

		ret=0
		do i = ln,1,-1
			if (arr(i).GT.deftol) then
				ret=i
				exit
			endif
		enddo

		return  ! return the position of last non-zero element
		end function lenTrimArray_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function lenTrimArray_r4(arr,tol) result(ret)
		implicit none

		real*4, intent(in)           :: arr(:)
		real*4, intent(in), optional :: tol
		real*4                       :: deftol
		integer*4                    :: ln,i


		ln=UBound(arr,1)

		deftol=1E-10; if (present(tol)) deftol=tol

		ret=0
		do i = ln,1,-1
			if (arr(i).GT.deftol) then
				ret=i
				exit
			endif
		enddo

		return  ! return the position of last non-zero element
		end function lenTrimArray_r4

!   ~~~~~~~~~~~~~~~~~~~~~~collect array functions~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function collectArray_i8(arr) result(ret)
		implicit none
		integer*8 :: arr(:)
		integer*4 :: ln,i,lst,els
		integer*4 :: fzero,fnzero


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function collectArray_i4(arr) result(ret)
		implicit none
		integer*4 :: arr(:)
		integer*4 :: ln,i,lst,els
		integer*4 :: fzero,fnzero


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function collectArray_i2(arr) result(ret)
		implicit none
		integer*2 :: arr(:)
		integer*4 :: ln,i,lst,els
		integer*4 :: fzero,fnzero


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function collectArray_i1(arr) result(ret)
		implicit none
		integer*1 :: arr(:)
		integer*4 :: ln,i,lst,els
		integer*4 :: fzero,fnzero


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function collectArray_r8(arr,tol) result(ret)
		implicit none
		real*8           :: deftol,arr(:),hold
		real*8, optional :: tol
		integer*4        :: ln,i,lst,els
		integer*4        :: fzero,fnzero


		ln=UBound(arr,1) ! define size of arr
		deftol=1D-10; if (present(tol)) deftol=tol

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function collectArray_r4(arr,tol) result(ret)
		implicit none
		real*4           :: deftol,arr(:),hold
		real*4, optional :: tol
		integer*4        :: ln,i,lst,els
		integer*4        :: fzero,fnzero


		ln=UBound(arr,1) ! define size of arr
		deftol=1E-10; if (present(tol)) deftol=tol

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure function strByLogical_l1(logi) result(ret)
		implicit none
		logical*1, intent(in)     :: logi
		character (len=mid(logi)) :: ret


		if (logi) then; ret='true'
		else;           ret='false'; endif; return
		end function strByLogical_l1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function definePi() ! define pi and degr->radian koef.
		implicit none


		pi=2.d0*acos(0.d0); dtr=pi/180.d0

		definePi=0; return
		end function definePi

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure function convertTime(secs,spec)
		implicit none
		character (len=20)               :: convertTime
		real*8, intent(in)               :: secs
		logical*1, optional, intent(in)  :: spec

		integer*4                        :: days,hours,minutes,seconds,mseconds,inLen
		character (len=4)                :: dayid
		logical                          :: uspec


		uspec=false
		if (present(spec)) uspec=spec

		days    =int( secs/86400.d0 )
		hours   =int( secs/3600.d0 - days*24.d0 )
		minutes =int( secs/60.d0 - hours*60.d0 - days*1440.d0 )
		seconds =int( secs - minutes*60.d0 - hours*3600.d0 - days*86400.d0 )
		mseconds=int( 100.d0*(secs-int(secs)) )

		if (uspec) then
			hours=hours+days*24
			inLen=mid(hours); if (hours.LT.2) inLen=2
			write (convertTime,100) hours,minutes,seconds,mseconds
			100 format (i<inLen>.2,':',i2.2,':',i2.2,'.',i2.2)
		else
			dayid='day '
			if ((days.GT.1).OR.(days.EQ.0)) dayid='days'
			write (convertTime,101) days,dayid,hours,minutes,seconds,mseconds
			101 format (i<mid(days)>,1X,A4,1X,i2.2,':',i2.2,':',i2.2,'.',i2.2)
		endif

		return
		end function convertTime

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function same_i8(x,y) result(ret)
		implicit none
		integer*8, intent(in) :: x,y


		ret=.NOT.btest(x+y,0); return
		end function same_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function same_i4(x,y) result(ret)
		implicit none
		integer*4, intent(in) :: x,y


		ret=.NOT.btest(x+y,0); return
		end function same_i4

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function same_i2(x,y) result(ret)
		implicit none
		integer*2, intent(in) :: x,y


		ret=.NOT.btest(x+y,0); return
		end function same_i2

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function same_i1(x,y) result(ret)
		implicit none
		integer*1, intent(in) :: x,y


		ret=.NOT.btest(x+y,0); return
		end function same_i1

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function isEven_i8(x) result(ret)
		implicit none
		integer*8, intent(in) :: x


		ret=.NOT.btest(x,0); return
		end function isEven_i8

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function isEven_i4(x) result(ret)
		implicit none
		integer*4, intent(in) :: x


		ret=.NOT.btest(x,0); return
		end function isEven_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function isEven_i2(x) result(ret)
		implicit none
		integer*2, intent(in) :: x


		ret=.NOT.btest(x,0); return
		end function isEven_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure logical*1 function isEven_i1(x) result(ret)
		implicit none
		integer*1, intent(in) :: x


		ret=.NOT.btest(x,0); return
		end function isEven_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function igen(insta,insto)
		implicit none
		integer*4, intent(in) :: insta,insto
		integer*4             :: setG,sigsw,inum,multip
		real*8                :: rnum


		if (insta.GT.insto) then
			igen=0; return
		endif

		if (insto-insta.EQ.0) then
			igen=insta; return
		endif

		setG=max(abs(insta),abs(insto))
		setG=int(log10(dble(setG)))+1

		if ((insta.LT.0).AND.(insto.LE.0)) sigsw=1
		if ((insta.LT.0).AND.(insto.GT.0)) sigsw=2
		if ((insta.GE.0).AND.(insto.GT.0)) sigsw=3

		multip=int(10**setG)

		do
			call random_number(rnum)
			inum=int(rnum*multip)
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

		igen=inum; return
		end function igen

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		real*8 function rgen(insta,insto)
		implicit none
		real*8, intent(in) :: insta,insto
		integer*4          :: sigsw
		real*8             :: setG,rnum,renum,multip


		if (insta.GT.insto) then; rgen=0; return; endif
		if (insto-insta.LE.1D-10) then; rgen=insta; return; endif

		setG=max(abs(insta),abs(insto))
		setG=int(log10(setG))+1

		if ((insta.LT.0).AND.(insto.LE.0)) sigsw=1
		if ((insta.LT.0).AND.(insto.GT.0)) sigsw=2
		if ((insta.GE.0).AND.(insto.GT.0)) sigsw=3

		multip=10.d0**setG

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

		rgen=renum; return
		end function rgen

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_uc(array,val) result(ret)
		implicit none
		type(uch)        , intent(in) :: array(:)
		character (len=*), intent(in) :: val
		integer*4                     :: i


		do i = lbound(array,1),ubound(array,1)
			if ( trim(val) .EQ. trim(uchGet(array(i))) ) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_uc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_c(array,val) result(ret)
		implicit none
		character (len=*), intent(in) :: array(:)
		character (len=*), intent(in) :: val
		integer*4                     :: i


		do i = lbound(array,1),ubound(array,1)
			if ( trim(val) .EQ. trim(array(i)) ) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_r8(array,val) result(ret)
		implicit none
		real*8, intent(in) :: val,array(:)
		real*8, parameter  :: rTol=1D-15
		integer*4          :: i


		do i = lbound(array,1),ubound(array,1)
			if (abs(array(i)-val).LE.rTol) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_r4(array,val) result(ret)
		implicit none
		real*4, intent(in) :: val,array(:)
		real*4, parameter  :: rTol=1D-8
		integer*4          :: i


		do i = lbound(array,1),ubound(array,1)
			if (abs(array(i)-val).LE.rTol) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_i8(array,val) result(ret)
		implicit none
		integer*8, intent(in) :: val,array(:)
		integer*4             :: i


		do i = lbound(array,1),ubound(array,1)
			if (array(i).EQ.val) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_i4(array,val) result(ret)
		implicit none
		integer*4, intent(in) :: val,array(:)
		integer*4             :: i


		do i = lbound(array,1),ubound(array,1)
			if (array(i).EQ.val) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_i2(array,val) result(ret)
		implicit none
		integer*2, intent(in) :: val,array(:)
		integer*4             :: i


		do i = lbound(array,1),ubound(array,1)
			if (array(i).EQ.val) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function find_i1(array,val) result(ret)
		implicit none
		integer*1, intent(in) :: val,array(:)
		integer*4             :: i


		do i = lbound(array,1),ubound(array,1)
			if (array(i).EQ.val) then
				ret=i; return
			endif
		enddo

		ret=-1; return
		end function find_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function sort_r8(vector,rev) result(ret)
		implicit none
		real*8              :: vector(:),swap
		logical*1, optional :: rev
		logical*1           :: urev
		integer*4           :: send,k,vlen


		urev=false; if (present(rev)) urev=rev

		send=UBound(vector,1)
		call sortqq(loc(vector),send,SRT$REAL8)
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
		end function sort_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function sort_r4(vector,rev) result(ret)
		implicit none
		real*4              :: vector(:),swap
		logical*1, optional :: rev
		logical*1           :: urev
		integer*4           :: send,k,vlen


		urev=false; if (present(rev)) urev=rev

		send=UBound(vector,1)
		call sortqq(loc(vector),send,SRT$REAL4)
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
		end function sort_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function sort_i4(vector,rev) result(ret)
		implicit none
		integer*4           :: vector(:),swap
		logical*1, optional :: rev
		logical*1           :: urev
		integer*4           :: send,k,vlen


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function sort_i2(vector,rev) result(ret)
		implicit none
		integer*2           :: vector(:),swap
		logical*1, optional :: rev
		logical*1           :: urev
		integer*4           :: send,k,vlen


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function sort_i1(vector,rev) result(ret)
		implicit none
		integer*1           :: vector(:),swap
		logical*1, optional :: rev
		logical*1           :: urev
		integer*4           :: send,k,vlen


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function sort_uch(vector,rev,attribute) result(ret)
		implicit none
		type(uch)                      :: vector(:),swap
		character (len=*), optional    :: attribute
		character (len=len(attribute)) :: uattribute
		logical*1        , optional    :: rev
		logical*1                      :: urev,cnd
		integer*4                      :: N,k,cnt


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		pure integer*4 function compareStrings(fstr,sstr) result(ret)
		implicit none

		character (len=*)  , intent(in) :: fstr,sstr
		integer*4                       :: fnd,k


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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		real(kind=rglu) function timeControl(cputime) result(ret)
		implicit none
		
		real(kind=rglu), optional :: cputime
		
		
		!!MS$if(OS.EQ.1)
			call cpu_time(ret)
			if (present(cputime)) cputime=ret
		!!MS$else
			!MS$if(opnmp.EQ.1)
				ret=omp_get_wtime()
				if (present(cputime)) call cpu_time(cputime)
			!MS$else
				call cpu_time(ret)
				if (present(cputime)) cputime=ret
			!MS$endif				
		!!MS$endif
		
		return
		end function timeControl

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		
		pure integer*4 function isLeapYear(year) result(ret)
		implicit none

		integer*4, intent(in) :: year


		if ( ((mod(year,4).EQ.0).AND.(mod(year,100).NE.0)) .OR. (mod(year,400).EQ.0)) then
			ret=1
		else
			ret=0
		endif

		return
		end function isLeapYear

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function dayOfWeek(year,month,day,past) result(ret)
		implicit none

		integer*4, parameter  :: months(12)=[31,28,31,30,31,30,31,31,30,31,30,31]
		integer*4, parameter  :: offsetYear=1900,offsetMonth    =1,&
		                         offsetDay =1   ,offsetDayOfWeek=1 ! 1=Monday

		integer*4, intent(in)            :: year,month,day
		integer*4, intent(out), optional :: past

		integer*4                        :: i,dyear,dmonth,dday,daysPast


		dyear =year -offsetYear
		dmonth=month-offsetMonth
		dday  =day  -offsetDay

		daysPast=offsetDayOfWeek

		daysPast=daysPast+dyear*365
		do i = 1,dmonth
			daysPast=daysPast+months(i)
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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		function timeStamp() result(ret)
		character (len=19) :: ret
		integer*4          :: values(8)


		values=0; call date_and_time(values=values)

		!VALUE(1): The year
		!VALUE(2): The month
		!VALUE(3): The day of the month
		!VALUE(4): Time difference with UTC in minutes
		!VALUE(5): The hour of the day
		!VALUE(6): The minutes of the hour
		!VALUE(7): The seconds of the minute
		!VALUE(8): The milliseconds of the second

		write (ret,1) values(1:3),values(5:7)

		1 format (i4,'.',i2.2,'.',i2.2,1X,i2.2,':',i2.2,':',i2.2)

		return
		end function timeStamp

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function setThreadsNumber(n_threads) result(ret)
		implicit none

		integer(kind=iglu) :: n_threads


		!!MS$if(OS.EQ.2)
			!MS$if(opnmp.EQ.1)
				call omp_set_num_threads(n_threads)
			!MS$endif
		!!MS$endif
		
		ret=0; return
		end function setThreadsNumber

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function getThreadNumber() result(ret)
		implicit none


		ret=1
		!!MS$if(OS.EQ.2)
			!MS$if(opnmp.EQ.1)
				ret=omp_get_thread_num()+1
			!MS$endif
		!!MS$endif
		
		return
		end function getThreadNumber

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function osCall(msg) result(ret)
		implicit none
		character (len=*) :: msg
		integer*4         :: sstat


		!MS$if(OS.EQ.1)
			call  system(msg)
		!MS$else
			sStat=system(msg)
		!MS$endif
		ret=0; return
		end function osCall

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine glSetIOunit(iunt)
		implicit none

		integer*4, intent(in) :: iunt
		integer*4             :: unt


		unt=iunt
		if (unt.EQ.5) unt=6 !stdin  (forbidden)
		if (unt.EQ.0) unt=6 !stdout (for compatibility)
		iounit=unt

		return
		end subroutine glSetIOunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine glFinalize
		implicit none


		!continue

		return
		end subroutine glFinalize


!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module glob
