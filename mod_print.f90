	module printmod

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	use glob, only: mid,sort,void,true,false,collectArray,lenTrimArray
	use fcontrol
	use txtParser

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: prVersion='2.212'
	character (len=*), parameter :: prDate   ='2017.07.10'
	character (len=*), parameter :: prAuthor ='Anton B. Zakharov'

	character (len=*), parameter :: prMatrixPrefix=' '

	integer*4, parameter      :: prBufferStrLen=1000

	integer*4, parameter      :: maxfmtlen=128,maxVariableNameLen=128,defiDigits=7,&
	                             defrDigits=10,defrPresc=5,defmaxWidth=150,        &
								 maxBufferLen=100,fmtPrintLen=32
	real*8   , parameter      :: defsptol=1D-8

	integer*4, parameter      :: fixedDistance=2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type variable
		character (len=maxfmtlen)          :: fmt
		character (len=maxVariableNameLen) :: label
		integer*4                          :: vtype
		real*8                             :: rvalue
		integer*4                          :: ivalue
		integer*4                          :: labelLen,dotShift
	end type variable

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type(variable)                 :: varBuffer(maxBufferLen)
	character (len=prBufferStrLen) :: prReadStr

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=maxfmtlen) :: ufmt,ueigenfmt
	integer*4                 :: arrType,vLen,internalIOunit,iMid,rMid,sigL,dotShift
	integer*4                 :: carriageControl,bufferLen=0,prBuffer=0

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	interface prMatrix
		module procedure r8printer,r4printer,i4printer,i2printer,i1printer
	end interface prMatrix

	interface prVector
		module procedure r8vector,r4vector,i4vector,i2vector,i1vector
	end interface prVector

	interface prAccumulateValues
		module procedure prAccumulateValuesSingler,prAccumulateValuesArrayr,&
		                 prAccumulateValuesSinglei,prAccumulateValuesArrayi
	end interface prAccumulateValues

	interface prFlushBuffer
		module procedure prFlushBufferSingle,prFlushBufferArray
	end interface prFlushBuffer

	interface prStrByVal
		module procedure nomStrByInt_i8,nomStrByInt_i4,nomStrByInt_i2,nomStrByInt_i1,    &
		                 defStrByInt_i8,defStrByInt_i4,defStrByInt_i2,defStrByInt_i1,    &
		                 strByRealf_r8,strByReale_r8,strByRealf_r4,strByReale_r4,        &
		                 strByLogical_l8,strByLogical_l4,strByLogical_l2,strByLogical_l1,&
						 strByRealfmt_r8,strByRealfmt_r4,strByIntfmt_i8,strByIntfmt_i4,  &
						 strByIntfmt_i2,strByIntfmt_i1
	end interface prStrByVal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private

	public :: prVersion,prDate,prAuthor,prFinalize
	public :: prUndoOutput
	public :: prMatrix,prVector,prEigenProblem
	public :: prAccumulateValues,prFlushValues
	public :: prLongText,prTable,prEchoFile,prStrByVal
	public :: prBuffer,prOpenBuffer,prCloseBuffer,prFlushBuffer,prAppendFile

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function prOpenBuffer() result(ret)
		implicit none


		if (prBuffer.NE.0) void=fcNullID(prBuffer)
		prBuffer=fcNewID(); open(prBuffer,status='scratch')

		ret=prBuffer; return
		end function prOpenBuffer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prCloseBuffer()
		implicit none


		if (prBuffer.EQ.0) return
		void=fcNullID(prBuffer)

		return
		end subroutine prCloseBuffer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prFlushBufferSingle(iounit)
		implicit none

		integer*4 :: iounit

		
		if (prBuffer.EQ.0) return

		rewind(prBuffer)
		do
			if (eof(prBuffer)) exit
			read  (prBuffer,'(A)') prReadStr
			write (iounit,'(A)') trim(prReadStr)
		enddo
		rewind(prBuffer)

		return
		end subroutine prFlushBufferSingle

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prFlushBufferArray(iounit)
		implicit none

		integer*4 :: iounit(:)
		integer*4 :: k

		
		if (prBuffer.EQ.0) return

		rewind(prBuffer)
		do
			if (eof(prBuffer)) exit
			read  (prBuffer,'(A)') prReadStr
			do k = 1,UBound(iounit,1)
				write (iounit(k),'(A)') trim(prReadStr)
			enddo
		enddo
		rewind(prBuffer)

		return
		end subroutine prFlushBufferArray

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function prAccumulateValuesSingler(label,value,fmt) result(ret)
		implicit none

		character (len=*)           :: label
		real*8                      :: value
		character (len=*), optional :: fmt


		rMid=mid(abs(value)); arrType=1
		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif

		bufferLen=bufferLen+1
		if (bufferLen.GT.maxBufferLen) then; ret=-1; return; endif

		varBuffer(bufferLen)%vtype=arrType
		varBuffer(bufferLen)%rvalue=value
		varBuffer(bufferLen)%label=trim(label)
		varBuffer(bufferLen)%labelLen=len_trim(label)
		varBuffer(bufferLen)%fmt=trim(ufmt)
		varBuffer(bufferLen)%dotShift=dotShift

		ret=bufferLen; return
		end function prAccumulateValuesSingler

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function prAccumulateValuesSinglei(label,value,fmt) result(ret)
		implicit none

		character (len=*)           :: label
		integer*4                   :: value
		character (len=*), optional :: fmt


		iMid=mid(abs(value)); arrType=2
		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif

		bufferLen=bufferLen+1
		if (bufferLen.GT.maxBufferLen) then; ret=-1; return; endif

		varBuffer(bufferLen)%vtype=arrType
		varBuffer(bufferLen)%ivalue=value
		varBuffer(bufferLen)%label=trim(label)
		varBuffer(bufferLen)%labelLen=len_trim(label)
		varBuffer(bufferLen)%fmt=trim(ufmt)
		varBuffer(bufferLen)%dotShift=dotShift

		ret=bufferLen; return
		end function prAccumulateValuesSinglei

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function prAccumulateValuesArrayr(label,value,fmt) result(ret)
		implicit none

		character (len=*)           :: label(:)
		real*8                      :: value(:)
		character (len=*), optional :: fmt
		integer*4                   :: i


		arrType=1
		do i = 1,min(UBound(label,1),UBound(value,1))
			rMid=mid(abs(value(i)))

			if ( present(fmt) ) then
				void=parsfmt(fmt)
			else
				void=parsfmt()
			endif

			bufferLen=bufferLen+1
			if (bufferLen.GT.maxBufferLen) then; ret=-1; return; endif

			varBuffer(bufferLen)%vtype=arrType
			varBuffer(bufferLen)%rvalue=value(i)
			varBuffer(bufferLen)%label=trim(label(i))
			varBuffer(bufferLen)%labelLen=len_trim(label(i))
			varBuffer(bufferLen)%fmt=trim(ufmt)
			varBuffer(bufferLen)%dotShift=dotShift
		enddo

		ret=bufferLen; return
		end function prAccumulateValuesArrayr

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function prAccumulateValuesArrayi(label,value,fmt) result(ret)
		implicit none

		character (len=*)           :: label(:)
		integer*4                   :: value(:)
		character (len=*), optional :: fmt
		integer*4                   :: i


		arrType=2
		do i = 1,min(UBound(label,1),UBound(value,1))
			iMid=mid(abs(value(i)))

			if ( present(fmt) ) then
				void=parsfmt(fmt)
			else
				void=parsfmt()
			endif

			bufferLen=bufferLen+1
			if (bufferLen.GT.maxBufferLen) then; ret=-1; return; endif

			varBuffer(bufferLen)%vtype=arrType
			varBuffer(bufferLen)%ivalue=value(i)
			varBuffer(bufferLen)%label=trim(label(i))
			varBuffer(bufferLen)%labelLen=len_trim(label(i))
			varBuffer(bufferLen)%fmt=trim(ufmt)
			varBuffer(bufferLen)%dotShift=dotShift
		enddo

		ret=bufferLen; return
		end function prAccumulateValuesArrayi

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prFlushValues(iounit,indent,inverse)
		implicit none

		integer*4           :: iounit
		integer*4, optional :: indent
		logical*1, optional :: inverse
		logical*1           :: dinverse
		character (len=100) :: ifmt

		integer*4           :: labelShift,valueShift,rDiff,&
		                       i,sta,sto,ste,dindent


		carriageControl=0; internalIOunit=iounit
		dinverse=false; if (present(inverse)) dinverse=inverse
		dindent=1; if (present(indent)) dindent=indent

		if (dinverse) then
			sta=bufferLen; sto=1; ste=-1
		else
			sta=1; sto=bufferLen; ste=1
		endif

		labelShift=0; valueShift=0
		do i = 1,bufferLen
			labelShift=max(labelShift,varBuffer(i)%labelLen)
			valueShift=max(valueShift,varBuffer(i)%dotShift)
		enddo

		do i = sta,sto,ste
			rDiff=(labelShift-varBuffer(i)%labelLen)+(valueShift-varBuffer(i)%dotShift)+1
			ifmt=tpFill(ifmt)
			ifmt='('//prStrByVal(dindent)//'X,A,'//prStrByVal(rDiff)//'X,'//trim(varBuffer(i)%fmt)//')'

			call advanceCarriage
			select case (varBuffer(i)%vtype)
				case (1)
					write (iounit,ifmt) trim(varBuffer(i)%label),varBuffer(i)%rvalue

				case (2)
					write (iounit,ifmt) trim(varBuffer(i)%label),varBuffer(i)%ivalue
			end select
		enddo
		bufferLen=0

		return
		end subroutine prFlushValues

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prAppendFile(src,trgt,prefix)
		implicit none

		character (len=*)           :: src
		integer*4                   :: trgt,iounit,err
		character (len=*), optional :: prefix
		character (len=64)          :: dprefix
		character (len=300)         :: readline


		dprefix=tpFill(dprefix); dprefix='INPUT FILE'
		if (present(prefix)) dprefix=trim(prefix)

		iounit=fcNewID()
		open (iounit,file=src,status='old',iostat=err)

		if (err.NE.0) then
			return
		endif

		do
			if (eof(iounit)) exit
			readline=tpFill(readline); read (iounit,'(A)') readline
			write (trgt,'(A,1X,A)') trim(dprefix),trim(readline)
		enddo
		close (iounit); void=fcNullID(iounit)

		return
		end subroutine prAppendFile

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine r8vector(vector,iounit,label,fmt,orient,maxwidth,maxcols,err)
		implicit none

		real*8                      :: vector(:)
		integer*4                   :: iounit
		character (len=*)           :: label

		character (len=*), optional :: fmt,orient
		integer*4        , optional :: maxwidth,maxcols,err

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=10)          :: uorient
		integer*4                   :: umaxwidth,ucolumns

		integer*4                   :: N,indexLen,columnLen,Nrows
		
		integer*4                   :: sta,sto,k,i

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(vector,1); arrType=1; rMid=mid(maxval(abs(vector))); indexLen=mid(N)

		umaxwidth=defmaxWidth  ; if ( present(maxwidth) ) umaxwidth=maxwidth
		uorient=tpFill(uorient); uorient='horizontal'; if ( present(orient) ) uorient=trim(orient)
		carriageControl=0; internalIOunit=iounit

		if (.NOT.tpIsStrInList(trim(uorient),['horizontal','vertical']) )   then
			if(present(err)) err=-1
			return
		endif

		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif		

		columnLen=1+indexLen+1+fixedDistance-1+vLen

		if (present(maxcols)) then
			ucolumns=maxcols
			if ( present(maxwidth) .AND. (ucolumns.GT.int( (umaxwidth-1)/( columnLen ) ))) then
				ucolumns=int( (umaxwidth-1)/( columnLen ) )
			endif
		else
			ucolumns=int( (umaxwidth-1)/( columnLen ) )
		endif
		Nrows=N/ucolumns; if (mod(N,ucolumns).NE.0) Nrows=Nrows+1

		ufmt='(1X,'//prStrByVal(ucolumns)//'(1X,i'//prStrByVal(indexLen)//',")",'//prStrByVal(fixedDistance-1)//'X,'//trim(ufmt)//'))'

		write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(3)
		select case (trim(uorient))
			case ('horizontal')
				do k = 1,nrows
					sta=ucolumns*(k-1)+1; sto=k*ucolumns; if (sto.GT.N) sto=N
					write (iounit,ufmt) (i,vector(i), i=sta,sto); call advanceCarriage
				enddo

			case ('vertical')
				do k = 1,nrows
					write (iounit,ufmt) (i,vector(i), i=k,N,Nrows); call advanceCarriage
				enddo
		end select
		write (iounit,100); call advanceCarriage(3)

 99 format (/ 2X,A,1X,A/)
100	format (//)

		return
		end subroutine r8vector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine r4vector(vector,iounit,label,fmt,orient,maxwidth,maxcols,err)
		implicit none

		real*4                      :: vector(:)
		integer*4                   :: iounit
		character (len=*)           :: label

		character (len=*), optional :: fmt,orient
		integer*4        , optional :: maxwidth,maxcols,err

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=10)          :: uorient
		integer*4                   :: umaxwidth,ucolumns

		integer*4                   :: N,indexLen,columnLen,Nrows
		
		integer*4                   :: sta,sto,k,i

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(vector,1); arrType=1; rMid=mid(maxval(abs(vector))); indexLen=mid(N)

		umaxwidth=defmaxWidth  ; if ( present(maxwidth) ) umaxwidth=maxwidth
		uorient=tpFill(uorient); uorient='horizontal'; if ( present(orient) ) uorient=trim(orient)
		carriageControl=0; internalIOunit=iounit

		if (.NOT.tpIsStrInList(trim(uorient),['horizontal','vertical']) )   then
			if(present(err)) err=-1
			return
		endif

		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif		

		columnLen=1+indexLen+1+fixedDistance-1+vLen

		if (present(maxcols)) then
			ucolumns=maxcols
			if ( present(maxwidth) .AND. (ucolumns.GT.int( (umaxwidth-1)/( columnLen ) ))) then
				ucolumns=int( (umaxwidth-1)/( columnLen ) )
			endif
		else
			ucolumns=int( (umaxwidth-1)/( columnLen ) )
		endif
		Nrows=N/ucolumns; if (mod(N,ucolumns).NE.0) Nrows=Nrows+1

		ufmt='(1X,'//prStrByVal(ucolumns)//'(1X,i'//prStrByVal(indexLen)//',")",'//prStrByVal(fixedDistance-1)//'X,'//trim(ufmt)//'))'

		write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(3)
		select case (trim(uorient))
			case ('horizontal')
				do k = 1,nrows
					sta=ucolumns*(k-1)+1; sto=k*ucolumns; if (sto.GT.N) sto=N
					write (iounit,ufmt) (i,vector(i), i=sta,sto); call advanceCarriage
				enddo

			case ('vertical')
				do k = 1,nrows
					write (iounit,ufmt) (i,vector(i), i=k,N,Nrows); call advanceCarriage
				enddo
		end select
		write (iounit,100); call advanceCarriage(3)

 99 format (/ 2X,A,1X,A/)
100	format (//)

		return
		end subroutine r4vector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine i4vector(vector,iounit,label,fmt,orient,maxwidth,maxcols,err)
		implicit none

		integer*4                   :: vector(:)
		integer*4                   :: iounit
		character (len=*)           :: label

		character (len=*), optional :: fmt,orient
		integer*4        , optional :: maxwidth,maxcols,err

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=10)          :: uorient
		integer*4                   :: umaxwidth,ucolumns

		integer*4                   :: N,indexLen,columnLen,Nrows
		
		integer*4                   :: sta,sto,k,i

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(vector,1); arrType=2; iMid=mid(maxval(abs(vector))); indexLen=mid(N)

		umaxwidth=defmaxWidth  ; if ( present(maxwidth) ) umaxwidth=maxwidth
		uorient=tpFill(uorient); uorient='horizontal'; if ( present(orient) ) uorient=trim(orient)
		carriageControl=0; internalIOunit=iounit

		if (.NOT.tpIsStrInList(trim(uorient),['horizontal','vertical']) )   then
			if(present(err)) err=-1
			return
		endif

		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif		

		columnLen=1+indexLen+1+fixedDistance-1+vLen

		if (present(maxcols)) then
			ucolumns=maxcols
			if ( present(maxwidth) .AND. (ucolumns.GT.int( (umaxwidth-1)/( columnLen ) ))) then
				ucolumns=int( (umaxwidth-1)/( columnLen ) )
			endif
		else
			ucolumns=int( (umaxwidth-1)/( columnLen ) )
		endif
		Nrows=N/ucolumns; if (mod(N,ucolumns).NE.0) Nrows=Nrows+1

		ufmt='(1X,'//prStrByVal(ucolumns)//'(1X,i'//prStrByVal(indexLen)//',")",'//prStrByVal(fixedDistance-1)//'X,'//trim(ufmt)//'))'

		write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(3)
		select case (trim(uorient))
			case ('horizontal')
				do k = 1,nrows
					sta=ucolumns*(k-1)+1; sto=k*ucolumns; if (sto.GT.N) sto=N
					write (iounit,ufmt) (i,vector(i), i=sta,sto); call advanceCarriage
				enddo

			case ('vertical')
				do k = 1,nrows
					write (iounit,ufmt) (i,vector(i), i=k,N,Nrows); call advanceCarriage
				enddo
		end select
		write (iounit,100); call advanceCarriage(3)

 99 format (/ 2X,A,1X,A/)
100	format (//)

		return
		end subroutine i4vector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine i2vector(vector,iounit,label,fmt,orient,maxwidth,maxcols,err)
		implicit none

		integer*2                   :: vector(:)
		integer*4                   :: iounit
		character (len=*)           :: label

		character (len=*), optional :: fmt,orient
		integer*4        , optional :: maxwidth,maxcols,err

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=10)          :: uorient
		integer*4                   :: umaxwidth,ucolumns

		integer*4                   :: N,indexLen,columnLen,Nrows
		
		integer*4                   :: sta,sto,k,i

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(vector,1); arrType=2; iMid=mid(maxval(abs(vector))); indexLen=mid(N)

		umaxwidth=defmaxWidth  ; if ( present(maxwidth) ) umaxwidth=maxwidth
		uorient=tpFill(uorient); uorient='horizontal'; if ( present(orient) ) uorient=trim(orient)
		carriageControl=0; internalIOunit=iounit

		if (.NOT.tpIsStrInList(trim(uorient),['horizontal','vertical']) )   then
			if(present(err)) err=-1
			return
		endif

		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif		

		columnLen=1+indexLen+1+fixedDistance-1+vLen

		if (present(maxcols)) then
			ucolumns=maxcols
			if ( present(maxwidth) .AND. (ucolumns.GT.int( (umaxwidth-1)/( columnLen ) ))) then
				ucolumns=int( (umaxwidth-1)/( columnLen ) )
			endif
		else
			ucolumns=int( (umaxwidth-1)/( columnLen ) )
		endif
		Nrows=N/ucolumns; if (mod(N,ucolumns).NE.0) Nrows=Nrows+1

		ufmt='(1X,'//prStrByVal(ucolumns)//'(1X,i'//prStrByVal(indexLen)//',")",'//prStrByVal(fixedDistance-1)//'X,'//trim(ufmt)//'))'

		write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(3)
		select case (trim(uorient))
			case ('horizontal')
				do k = 1,nrows
					sta=ucolumns*(k-1)+1; sto=k*ucolumns; if (sto.GT.N) sto=N
					write (iounit,ufmt) (i,vector(i), i=sta,sto); call advanceCarriage
				enddo

			case ('vertical')
				do k = 1,nrows
					write (iounit,ufmt) (i,vector(i), i=k,N,Nrows); call advanceCarriage
				enddo
		end select
		write (iounit,100); call advanceCarriage(3)

 99 format (/ 2X,A,1X,A/)
100	format (//)

		return
		end subroutine i2vector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine i1vector(vector,iounit,label,fmt,orient,maxwidth,maxcols,err)
		implicit none

		integer*1                   :: vector(:)
		integer*4                   :: iounit
		character (len=*)           :: label

		character (len=*), optional :: fmt,orient
		integer*4        , optional :: maxwidth,maxcols,err

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=10)          :: uorient
		integer*4                   :: umaxwidth,ucolumns

		integer*4                   :: N,indexLen,columnLen,Nrows
		
		integer*4                   :: sta,sto,k,i

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(vector,1); arrType=2; iMid=mid(maxval(abs(vector))); indexLen=mid(N)

		umaxwidth=defmaxWidth  ; if ( present(maxwidth) ) umaxwidth=maxwidth
		uorient=tpFill(uorient); uorient='horizontal'; if ( present(orient) ) uorient=trim(orient)
		carriageControl=0; internalIOunit=iounit

		if (.NOT.tpIsStrInList(trim(uorient),['horizontal','vertical']) )   then
			if(present(err)) err=-1
			return
		endif

		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif		

		columnLen=1+indexLen+1+fixedDistance-1+vLen

		if (present(maxcols)) then
			ucolumns=maxcols
			if ( present(maxwidth) .AND. (ucolumns.GT.int( (umaxwidth-1)/( columnLen ) ))) then
				ucolumns=int( (umaxwidth-1)/( columnLen ) )
			endif
		else
			ucolumns=int( (umaxwidth-1)/( columnLen ) )
		endif
		Nrows=N/ucolumns; if (mod(N,ucolumns).NE.0) Nrows=Nrows+1

		ufmt='(1X,'//prStrByVal(ucolumns)//'(1X,i'//prStrByVal(indexLen)//',")",'//prStrByVal(fixedDistance-1)//'X,'//trim(ufmt)//'))'

		write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(3)
		select case (trim(uorient))
			case ('horizontal')
				do k = 1,nrows
					sta=ucolumns*(k-1)+1; sto=k*ucolumns; if (sto.GT.N) sto=N
					write (iounit,ufmt) (i,vector(i), i=sta,sto); call advanceCarriage
				enddo

			case ('vertical')
				do k = 1,nrows
					write (iounit,ufmt) (i,vector(i), i=k,N,Nrows); call advanceCarriage
				enddo
		end select
		write (iounit,100); call advanceCarriage(3)

 99 format (/ 2X,A,1X,A/)
100	format (//)

		return
		end subroutine i1vector

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prEigenProblem(vecs,vals,iounit,label,fmt,maxwidth,maxcols,sparse,sptol,transpose,iostat)

		implicit none

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4        , intent(in)           :: iounit
		real*8           , intent(in)           :: vecs(:,:),vals(:)
		character (len=*), intent(in)           :: label
		real*8           , intent(in), optional :: sptol
		integer*4        , intent(in), optional :: maxwidth,maxcols
		character (len=*), intent(in), optional :: fmt
		logical*1        , intent(in), optional :: transpose,sparse
		integer*4        ,             optional :: iostat

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4                     :: N,M,rMid2
		real*8                        :: usptol
		integer*4                     :: i,j,k,a,umaxwidth,ucolumns,lblLen,uN,uM
		logical*1                     :: dsparse,diostat,dtranspose
		integer*4                     :: nmbrRowLen,nmbrColLen,fmtCols
		integer*4                     :: firstShift,indShift,valShift

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(vecs,1); M=N
		carriageControl=0; arrType=1
		internalIOunit=iounit
		rMid =mid(maxval(abs(vecs)))
		rMid2=mid(maxval(abs(vals)))

		if (rMid2.GT.rMid) rMid=rMid2

		usptol=defsptol;       if ( present(sptol)    ) usptol=sptol
		dsparse=false;         if ( present(sparse)   ) dsparse=sparse;       dsparse=false
		dtranspose=false;      if ( present(transpose)) dtranspose=transpose; dtranspose=false
		umaxwidth=defmaxWidth; if ( present(maxwidth) ) umaxwidth=maxwidth

		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif
		diostat=present(iostat)
		lblLen=len_trim(label)

		uN=N; uM=M
		nmbrRowLen=mid(N); nmbrColLen=mid(M)

		if (present(maxcols)) then
			ucolumns=maxcols
		else
			ucolumns=int( (umaxwidth-nmbrRowLen-1-4)/(vLen+fixedDistance) )
		endif
		firstShift=1+nmbrRowLen+4

		if (vLen-nmbrColLen.LT.0) then
			indShift=fixedDistance
			valShift=nmbrColLen-vLen+fixedDistance
		else
			indShift=vLen-nmbrColLen+fixedDistance
			valShift=fixedDistance
		endif
			
		write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(4)

		ueigenfmt='(1X,' //prStrByVal(nmbrRowLen)//'X,4X,'//prStrByVal(ucolumns)//'('//prStrByVal(valShift)//'X,'//trim(ufmt)//'))'
		ufmt     ='(1X,i'//prStrByVal(nmbrRowLen)// ',4X,'//prStrByVal(ucolumns)//'('//prStrByVal(valShift)//'X,'//trim(ufmt)//'))'

		do i = 1,uM,ucolumns
			if (uM-i.LT.ucolumns) then
				k=uM
			else
				k=i+ucolumns-1
			end if
			fmtCols=k-i+1
			write (iounit,ueigenfmt) (vals(a), a=i,k); call advanceCarriage
			write (iounit,100)       (a,       a=i,k); call advanceCarriage
			do j = 1,uN
				write (iounit,ufmt) j,(vecs(j,a), a=i,k); call advanceCarriage
			enddo
			write (iounit,101); call advanceCarriage(3)
		enddo

 99 format (// 2X,A,1X,A<lblLen>/)
199 format (// 2X,A,1X,A<lblLen>,1X,'(transposed)'/)
299 format (// 2X,A,1X,A<lblLen>,1X,'(unrarefied, print tol=',ES7.0,')'/)

100 format (<firstShift>X,<fmtCols>(<indShift>X,i<nmbrColLen>))
101	format (//)

105 format (1X,i<nmbrRowLen>,4X,<fmtCols>(<valShift>X,F9.4))
600 format (4X,'Error: ',i3,' occured during array output procedure.')

		return
		end subroutine prEigenProblem

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine r8printer(array,iounit,label,fmt,maxwidth,maxcols,sparse,sptol,transpose,iostat)

		implicit none

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4        , intent(in)           :: iounit
		real*8           , intent(in)           :: array(:,:)
		character (len=*), intent(in)           :: label
		real*8           , intent(in), optional :: sptol
		integer*4        , intent(in), optional :: maxwidth,maxcols
		character (len=*), intent(in), optional :: fmt
		logical*1        , intent(in), optional :: transpose,sparse
		integer*4        ,             optional :: iostat

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4                     :: N,M
		real*8                        :: usptol
		integer*4                     :: i,j,k,a,umaxwidth,ucolumns,lblLen,uN,uM
		logical*1                     :: dsparse,diostat,dtranspose
		integer*4                     :: nmbrRowLen,nmbrColLen,fmtCols
		integer*4                     :: firstShift,indShift,valShift

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(array,1); M=UBound(array,2)
		carriageControl=0; arrType=1
		internalIOunit=iounit
		rMid=mid(maxval(abs(array)))

		usptol=defsptol;       if ( present(sptol)    ) usptol=sptol
		dsparse=false;         if ( present(sparse)   ) dsparse=sparse
		dtranspose=false;      if ( present(transpose)) dtranspose=transpose
		umaxwidth=defmaxWidth; if ( present(maxwidth) ) umaxwidth=maxwidth
		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif
		diostat=present(iostat)
		lblLen=len_trim(label)

		if (dsparse) then
			write (iounit,299) prMatrixPrefix,trim(label),usptol; call advanceCarriage(3)
			nmbrRowLen=mid(N); nmbrColLen=mid(M)
			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',1X,i'//prStrByVal(nmbrColLen)//','//prStrByVal(fixedDistance)//'X,'//trim(ufmt)//')'
			do i = 1,N
			do j = 1,M
				if (abs(array(i,j)).GT.usptol) then
					call advanceCarriage
					write (iounit,ufmt) i,j,array(i,j)
				endif
			enddo
			enddo
		else
			if (.NOT.dtranspose) then
				uN=N; uM=M
				nmbrRowLen=mid(N); nmbrColLen=mid(M)
			else
				uN=M; uM=N
				nmbrRowLen=mid(M); nmbrColLen=mid(N)
			endif

			if (present(maxcols)) then
				ucolumns=maxcols
			else
				ucolumns=int( (umaxwidth-nmbrRowLen-1-4)/(vLen+fixedDistance) )
			endif
			firstShift=1+nmbrRowLen+4

			if (vLen-nmbrColLen.LT.0) then
				indShift=fixedDistance
				valShift=nmbrColLen-vLen+fixedDistance
			else
				indShift=vLen-nmbrColLen+fixedDistance
				valShift=fixedDistance
			endif
			
			if (.NOT.dtranspose) then
				write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(4)
			else
				write (iounit,199) prMatrixPrefix,trim(label); call advanceCarriage(4)
			endif

			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',4X,'//prStrByVal(ucolumns)//'('//prStrByVal(valShift)//'X,'//trim(ufmt)//'))'

			do i = 1,uM,ucolumns
				if (uM-i.LT.ucolumns) then
					k=uM
				else
					k=i+ucolumns-1
				end if
				fmtCols=k-i+1
				write (iounit,100) (a, a=i,k); call advanceCarriage(2)
				do j = 1,uN
					call advanceCarriage
					if (.NOT.dTranspose) then
						write (iounit,ufmt) j,(array(j,a), a=i,k)
					else
						write (iounit,ufmt) j,(array(a,j), a=i,k)
					endif
				enddo
				write (iounit,101); call advanceCarriage(3)
			enddo
		endif

 99 format (// 2X,A,1X,A<lblLen>/)
199 format (// 2X,A,1X,A<lblLen>,1X,'(transposed)'/)
299 format (// 2X,A,1X,A<lblLen>,1X,'(unrarefied, print tol=',ES7.0,')'/)

100 format (<firstShift>X,<fmtCols>(<indShift>X,i<nmbrColLen>)/)
101	format (//)
105 format (1X,i<nmbrRowLen>,4X,<fmtCols>(<valShift>X,F9.4))
600 format (4X,'Error: ',i3,' occured during array output procedure.')

		return
		end subroutine r8printer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine r4printer(array,iounit,label,fmt,maxwidth,maxcols,sparse,sptol,transpose,iostat)

		implicit none

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4        , intent(in)           :: iounit
		real*4           , intent(in)           :: array(:,:)
		character (len=*), intent(in)           :: label
		real*4           , intent(in), optional :: sptol
		integer*4        , intent(in), optional :: maxwidth,maxcols
		character (len=*), intent(in), optional :: fmt
		logical*1        , intent(in), optional :: transpose,sparse
		integer*4        ,             optional :: iostat

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4                     :: N,M
		real*4                        :: usptol
		integer*4                     :: i,j,k,a,umaxwidth,ucolumns,lblLen,uN,uM
		logical*1                     :: dsparse,diostat,dtranspose
		integer*4                     :: nmbrRowLen,nmbrColLen,fmtCols
		integer*4                     :: firstShift,indShift,valShift

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(array,1); M=UBound(array,2)
		carriageControl=0; arrType=1
		internalIOunit=iounit
		rMid=mid(maxval(abs(array)))

		usptol=defsptol;       if ( present(sptol)    ) usptol=sptol
		dsparse=false;         if ( present(sparse)   ) dsparse=sparse
		dtranspose=false;      if ( present(transpose)) dtranspose=transpose
		umaxwidth=defmaxWidth; if ( present(maxwidth) ) umaxwidth=maxwidth
		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif
		diostat=present(iostat)
		lblLen=len_trim(label)

		if (dsparse) then
			write (iounit,299) prMatrixPrefix,trim(label),usptol; call advanceCarriage(3)
			nmbrRowLen=mid(N); nmbrColLen=mid(M)
			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',1X,i'//prStrByVal(nmbrColLen)//','//prStrByVal(fixedDistance)//'X,'//trim(ufmt)//')'
			do i = 1,N
			do j = 1,M
				if (abs(array(i,j)).GT.usptol) then
					call advanceCarriage
					write (iounit,ufmt) i,j,array(i,j)
				endif
			enddo
			enddo
		else
			if (.NOT.dtranspose) then
				uN=N; uM=M
				nmbrRowLen=mid(N); nmbrColLen=mid(M)
			else
				uN=M; uM=N
				nmbrRowLen=mid(M); nmbrColLen=mid(N)
			endif

			if (present(maxcols)) then
				ucolumns=maxcols
			else
				ucolumns=int( (umaxwidth-nmbrRowLen-1-4)/(vLen+fixedDistance) )
			endif
			firstShift=1+nmbrRowLen+4

			if (vLen-nmbrColLen.LT.0) then
				indShift=fixedDistance
				valShift=nmbrColLen-vLen+fixedDistance
			else
				indShift=vLen-nmbrColLen+fixedDistance
				valShift=fixedDistance
			endif
			
			if (.NOT.dtranspose) then
				write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(4)
			else
				write (iounit,199) prMatrixPrefix,trim(label); call advanceCarriage(4)
			endif

			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',4X,'//prStrByVal(ucolumns)//'('//prStrByVal(valShift)//'X,'//trim(ufmt)//'))'

			do i = 1,uM,ucolumns
				if (uM-i.LT.ucolumns) then
					k=uM
				else
					k=i+ucolumns-1
				end if
				fmtCols=k-i+1
				write (iounit,100) (a, a=i,k); call advanceCarriage(2)
				do j = 1,uN
					call advanceCarriage
					if (.NOT.dTranspose) then
						write (iounit,ufmt) j,(array(j,a), a=i,k)
					else
						write (iounit,ufmt) j,(array(a,j), a=i,k)
					endif

				enddo
				write (iounit,101); call advanceCarriage(3)
			enddo
		endif

 99 format (// 2X,A,1X,A<lblLen>/)
199 format (// 2X,A,1X,A<lblLen>,1X,'(transposed)'/)
299 format (// 2X,A,1X,A<lblLen>,1X,'(unrarefied, print tol=',ES7.0,')'/)

100 format (<firstShift>X,<fmtCols>(<indShift>X,i<nmbrColLen>)/)
101	format (//)
105 format (1X,i<nmbrRowLen>,4X,<fmtCols>(<valShift>X,F9.4))
600 format (4X,'Error: ',i3,' occured during array output procedure.')

		return
		end subroutine r4printer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine i4printer(array,iounit,label,fmt,maxwidth,maxcols,sparse,transpose,iostat)

		implicit none

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4        , intent(in)           :: iounit
		integer*4        , intent(in)           :: array(:,:)
		character (len=*), intent(in)           :: label
		integer*4        , intent(in), optional :: maxwidth,maxcols
		character (len=*), intent(in), optional :: fmt
		logical*1        , intent(in), optional :: transpose,sparse
		integer*4        ,             optional :: iostat

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4                     :: N,M
		integer*4                     :: i,j,k,a,umaxwidth,ucolumns,lblLen,uN,uM
		logical*1                     :: dsparse,diostat,dtranspose
		integer*4                     :: nmbrRowLen,nmbrColLen,fmtCols
		integer*4                     :: firstShift,indShift,valShift

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(array,1); M=UBound(array,2)
		carriageControl=0; arrType=2
		internalIOunit=iounit
		iMid=mid(maxval(abs(array)))

		dsparse=false;         if ( present(sparse)   ) dsparse=sparse
		dtranspose=false;      if ( present(transpose)) dtranspose=transpose
		umaxwidth=defmaxWidth; if ( present(maxwidth) ) umaxwidth=maxwidth
		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif
		diostat=present(iostat)
		lblLen=len_trim(label)

		if (dsparse) then
			write (iounit,299) prMatrixPrefix,trim(label); call advanceCarriage(3)
			nmbrRowLen=mid(N); nmbrColLen=mid(M)
			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',1X,i'//prStrByVal(nmbrColLen)//','//prStrByVal(fixedDistance)//'X,'//trim(ufmt)//')'
			do i = 1,N
			do j = 1,M
				if (abs(array(i,j)).GT.0) then
					call advanceCarriage
					write (iounit,ufmt) i,j,array(i,j)
				endif
			enddo
			enddo
		else
			if (.NOT.dtranspose) then
				uN=N; uM=M
				nmbrRowLen=mid(N); nmbrColLen=mid(M)
			else
				uN=M; uM=N
				nmbrRowLen=mid(M); nmbrColLen=mid(N)
			endif

			if (present(maxcols)) then
				ucolumns=maxcols
			else
				ucolumns=int( (umaxwidth-nmbrRowLen-1-4)/(vLen+fixedDistance) )
			endif
			firstShift=1+nmbrRowLen+4

			if (vLen-nmbrColLen.LT.0) then
				indShift=fixedDistance
				valShift=nmbrColLen-vLen+fixedDistance
			else
				indShift=vLen-nmbrColLen+fixedDistance
				valShift=fixedDistance
			endif
			
			if (.NOT.dtranspose) then
				write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(4)
			else
				write (iounit,199) prMatrixPrefix,trim(label); call advanceCarriage(4)
			endif

			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',4X,'//prStrByVal(ucolumns)//'('//prStrByVal(valShift)//'X,'//trim(ufmt)//'))'

			do i = 1,uM,ucolumns
				if (uM-i.LT.ucolumns) then
					k=uM
				else
					k=i+ucolumns-1
				end if
				fmtCols=k-i+1
				write (iounit,100) (a, a=i,k); call advanceCarriage(2)
				do j = 1,uN
					call advanceCarriage
					if (.NOT.dTranspose) then
						write (iounit,ufmt) j,(array(j,a), a=i,k)
					else
						write (iounit,ufmt) j,(array(a,j), a=i,k)
					endif

				enddo
				write (iounit,101); call advanceCarriage(3)
			enddo
		endif

 99 format (// 2X,A,1X,A<lblLen>/)
199 format (// 2X,A,1X,A<lblLen>,1X,'(transposed)'/)
299 format (// 2X,A,1X,A<lblLen>,1X,'(unrarefied, print tol= 0)'/)

100 format (<firstShift>X,<fmtCols>(<indShift>X,i<nmbrColLen>)/)
101	format (//)
105 format (1X,i<nmbrRowLen>,4X,<fmtCols>(<valShift>X,F9.4))
600 format (4X,'Error: ',i3,' occured during array output procedure.')

		return
		end subroutine i4printer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine i2printer(array,iounit,label,fmt,maxwidth,maxcols,sparse,transpose,iostat)

		implicit none

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4        , intent(in)           :: iounit
		integer*2        , intent(in)           :: array(:,:)
		character (len=*), intent(in)           :: label
		integer*4        , intent(in), optional :: maxwidth,maxcols
		character (len=*), intent(in), optional :: fmt
		logical*1        , intent(in), optional :: transpose,sparse
		integer*4        ,             optional :: iostat

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4                     :: N,M
		integer*4                     :: i,j,k,a,umaxwidth,ucolumns,lblLen,uN,uM
		logical*1                     :: dsparse,diostat,dtranspose
		integer*4                     :: nmbrRowLen,nmbrColLen,fmtCols
		integer*4                     :: firstShift,indShift,valShift

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(array,1); M=UBound(array,2)
		carriageControl=0; arrType=2
		internalIOunit=iounit
		iMid=mid(maxval(abs(array)))

		dsparse=false;         if ( present(sparse)   ) dsparse=sparse
		dtranspose=false;      if ( present(transpose)) dtranspose=transpose
		umaxwidth=defmaxWidth; if ( present(maxwidth) ) umaxwidth=maxwidth
		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif
		diostat=present(iostat)
		lblLen=len_trim(label)

		if (dsparse) then
			write (iounit,299) prMatrixPrefix,trim(label); call advanceCarriage(3)
			nmbrRowLen=mid(N); nmbrColLen=mid(M)
			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',1X,i'//prStrByVal(nmbrColLen)//','//prStrByVal(fixedDistance)//'X,'//trim(ufmt)//')'
			do i = 1,N
			do j = 1,M
				if (abs(array(i,j)).GT.0) then
					call advanceCarriage
					write (iounit,ufmt) i,j,array(i,j)
				endif
			enddo
			enddo
		else
			if (.NOT.dtranspose) then
				uN=N; uM=M
				nmbrRowLen=mid(N); nmbrColLen=mid(M)
			else
				uN=M; uM=N
				nmbrRowLen=mid(M); nmbrColLen=mid(N)
			endif

			if (present(maxcols)) then
				ucolumns=maxcols
			else
				ucolumns=int( (umaxwidth-nmbrRowLen-1-4)/(vLen+fixedDistance) )
			endif
			firstShift=1+nmbrRowLen+4

			if (vLen-nmbrColLen.LT.0) then
				indShift=fixedDistance
				valShift=nmbrColLen-vLen+fixedDistance
			else
				indShift=vLen-nmbrColLen+fixedDistance
				valShift=fixedDistance
			endif
			
			if (.NOT.dtranspose) then
				write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(4)
			else
				write (iounit,199) prMatrixPrefix,trim(label); call advanceCarriage(4)
			endif

			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',4X,'//prStrByVal(ucolumns)//'('//prStrByVal(valShift)//'X,'//trim(ufmt)//'))'

			do i = 1,uM,ucolumns
				if (uM-i.LT.ucolumns) then
					k=uM
				else
					k=i+ucolumns-1
				end if
				fmtCols=k-i+1
				write (iounit,100) (a, a=i,k); call advanceCarriage(2)
				do j = 1,uN
					call advanceCarriage
					if (.NOT.dTranspose) then
						write (iounit,ufmt) j,(array(j,a), a=i,k)
					else
						write (iounit,ufmt) j,(array(a,j), a=i,k)
					endif
				enddo
				write (iounit,101); call advanceCarriage(3)
			enddo
		endif

 99 format (// 2X,A,1X,A<lblLen>/)
199 format (// 2X,A,1X,A<lblLen>,1X,'(transposed)'/)
299 format (// 2X,A,1X,A<lblLen>,1X,'(unrarefied, print tol= 0)'/)

100 format (<firstShift>X,<fmtCols>(<indShift>X,i<nmbrColLen>)/)
101	format (//)
105 format (1X,i<nmbrRowLen>,4X,<fmtCols>(<valShift>X,F9.4))
600 format (4X,'Error: ',i3,' occured during array output procedure.')

		return
		end subroutine i2printer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine i1printer(array,iounit,label,fmt,maxwidth,maxcols,sparse,transpose,iostat)

		implicit none

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4        , intent(in)           :: iounit
		integer*1        , intent(in)           :: array(:,:)
		character (len=*), intent(in)           :: label
		integer*4        , intent(in), optional :: maxwidth,maxcols
		character (len=*), intent(in), optional :: fmt
		logical*1        , intent(in), optional :: transpose,sparse
		integer*4        ,             optional :: iostat

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer*4                     :: N,M
		integer*4                     :: i,j,k,a,umaxwidth,ucolumns,lblLen,uN,uM
		logical*1                     :: dsparse,diostat,dtranspose
		integer*4                     :: nmbrRowLen,nmbrColLen,fmtCols
		integer*4                     :: firstShift,indShift,valShift

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		N=UBound(array,1); M=UBound(array,2)
		carriageControl=0; arrType=2
		internalIOunit=iounit
		iMid=mid(maxval(abs(array)))

		dsparse=false;         if ( present(sparse)   ) dsparse=sparse
		dtranspose=false;      if ( present(transpose)) dtranspose=transpose
		umaxwidth=defmaxWidth; if ( present(maxwidth) ) umaxwidth=maxwidth
		if ( present(fmt) ) then
			void=parsfmt(fmt)
		else
			void=parsfmt()
		endif
		diostat=present(iostat)
		lblLen=len_trim(label)

		if (dsparse) then
			write (iounit,299) prMatrixPrefix,trim(label); call advanceCarriage(3)
			nmbrRowLen=mid(N); nmbrColLen=mid(M)
			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',1X,i'//prStrByVal(nmbrColLen)//','//prStrByVal(fixedDistance)//'X,'//trim(ufmt)//')'
			do i = 1,N
			do j = 1,M
				if (abs(array(i,j)).GT.0) then
					call advanceCarriage
					write (iounit,ufmt) i,j,array(i,j)
				endif
			enddo
			enddo
		else
			if (.NOT.dtranspose) then
				uN=N; uM=M
				nmbrRowLen=mid(N); nmbrColLen=mid(M)
			else
				uN=M; uM=N
				nmbrRowLen=mid(M); nmbrColLen=mid(N)
			endif

			if (present(maxcols)) then
				ucolumns=maxcols
			else
				ucolumns=int( (umaxwidth-nmbrRowLen-1-4)/(vLen+fixedDistance) )
			endif
			firstShift=1+nmbrRowLen+4

			if (vLen-nmbrColLen.LT.0) then
				indShift=fixedDistance
				valShift=nmbrColLen-vLen+fixedDistance
			else
				indShift=vLen-nmbrColLen+fixedDistance
				valShift=fixedDistance
			endif
			
			if (.NOT.dtranspose) then
				write (iounit, 99) prMatrixPrefix,trim(label); call advanceCarriage(4)
			else
				write (iounit,199) prMatrixPrefix,trim(label); call advanceCarriage(4)
			endif

			ufmt='(1X,i'//prStrByVal(nmbrRowLen)//',4X,'//prStrByVal(ucolumns)//'('//prStrByVal(valShift)//'X,'//trim(ufmt)//'))'

			do i = 1,uM,ucolumns
				if (uM-i.LT.ucolumns) then
					k=uM
				else
					k=i+ucolumns-1
				end if
				fmtCols=k-i+1
				write (iounit,100) (a, a=i,k); call advanceCarriage(2)
				do j = 1,uN
					call advanceCarriage
					if (.NOT.dTranspose) then
						write (iounit,ufmt) j,(array(j,a), a=i,k)
					else
						write (iounit,ufmt) j,(array(a,j), a=i,k)
					endif
				enddo
				write (iounit,101); call advanceCarriage(3)
			enddo
		endif

 99 format (// 2X,A,1X,A<lblLen>/)
199 format (// 2X,A,1X,A<lblLen>,1X,'(transposed)'/)
299 format (// 2X,A,1X,A<lblLen>,1X,'(unrarefied, print tol= 0)'/)

100 format (<firstShift>X,<fmtCols>(<indShift>X,i<nmbrColLen>)/)
101	format (//)
105 format (1X,i<nmbrRowLen>,4X,<fmtCols>(<valShift>X,F9.4))
600 format (4X,'Error: ',i3,' occured during array output procedure.')

		return
		end subroutine i1printer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prLongText(str,iounit,adjust,margin,twidth,pwidth,spacer,err)
		implicit none

		character (len=*)          , intent(in)  :: str
		character (len=*), optional, intent(in)  :: adjust,margin
		character (len=1), optional, intent(in)  :: spacer
		integer*4,         optional, intent(in)  :: twidth,pwidth
		integer*4                  , intent(in)  :: iounit
		integer*4,         optional, intent(out) :: err

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=len(str))                 :: ustr,blstr
		character                                :: uadjust*15,umargin*15,uspacer*1
		integer*4                                :: tw,pw

		character, parameter                     :: replChar=char(181)
		character                                :: fmt*3
		logical*1                                :: reverse
		integer*4                                :: sta,sto,i,j,k,l,ln,cnt,numb,scrio

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !


		ustr=str
		uadjust=tpFill(uadjust); uadjust='center'     ; if (present(adjust)) uadjust=trim(adjust)
		umargin=tpFill(umargin); umargin='unjustified'; if (present(margin)) umargin=trim(margin)
		uspacer=tpFill(uspacer); uspacer=' '          ; if (present(spacer)) uspacer=trim(spacer)

		carriageControl=0; internalIOunit=iounit

		tw=60; if (present(twidth)) tw=twidth
		pw=78; if (present(pwidth)) pw=pwidth

		if (len(ustr).LT.tw) tw=len(ustr)

		if (.NOT.tpIsStrInList(uadjust,['left','center','right']) )   then
			if(present(err)) err=-1
			return
		endif
		if (.NOT.tpIsStrInList(umargin,['justified','unjustified']) ) then
			if(present(err)) err=-2
			return
		endif

		fmt='(A)'

		scrio=fcNewID()
		open (scrio,status='scratch')
		do while (index(ustr,newline).GT.0)
			ustr( tpIndex(ustr,newline):tpIndex(ustr,newline)+1 )='  '
		enddo

		ustr=trim(tpReduce(ustr))

		sta=1; sto=sta+tw; cnt=0
		do
			if (sto.GT.len_trim(ustr)) sto=len_trim(ustr)

			if (ustr(sto:sto).NE.' ') then
				if (sto.EQ.len_trim(ustr)) then

					! to be carefull with.
					if ( len(ustr(sta:sto)).GT.tw ) then
						! last change.
						k=tpIndex(ustr(sta:sto),' ',rev=true)
						if (k.EQ.0) then
							if(present(err)) err=-3; close (scrio); void=fcNullID(scrio); return
						endif
						cnt=cnt+1; write (scrio,fmt) ustr(sta:sta+k-1)
						cnt=cnt+1; write (scrio,fmt) ustr(sta+k:sto)
						exit
					endif
					!

					cnt=cnt+1; write (scrio,fmt) ustr(sta:sto)
					exit
				endif
				k=tpIndex(ustr(sta:sto),' ',rev=true)
				if (k.EQ.0) then
					if(present(err)) err=-4; close (scrio); void=fcNullID(scrio); return
				endif
				cnt=cnt+1; write (scrio,fmt) ustr(sta:sta+k-1)
				sta=sta+k
				sto=sta+tw
			else
				cnt=cnt+1; write (scrio,fmt) ustr(sta:sto)
				sta=sto+1
				sto=sto+tw+1
			endif
		enddo
		
		rewind(scrio)

		if (trim(umargin).EQ.'justified') then
			do i = 1,cnt-1
				blstr=tpFill(blstr); read (scrio,fmt) blstr(:tw)

				if (len_trim(blstr).NE.tw) then
					reverse=false; numb=tpCount(trim(blstr),' ')
					numb=numb/2+mod(numb,2); l=0
					do while (len_trim(blstr).NE.tw)
						
						l=l+1; if (l.GT.numb) l=1
						do j = 1,2
							if (len_trim(blstr).EQ.tw) exit

							k=tpIndex(trim(blstr),' ',cnt=l,rev=reverse)
							blstr=trim(tpInsert(blstr,replChar,k))
							reverse=.NOT.reverse
						enddo
					enddo
				endif

				blstr=trim( tpReplace(trim(blstr),replChar,' ') ); ln=pw-tw
				select case (trim(uadjust))
					case('left')  ; write (iounit,fmt) blstr(:tw)//tpFill(ln,uspacer)
					case('center'); write (iounit,fmt) tpFill(ln/2,uspacer)//blstr(:tw)//tpFill(ln/2+mod(ln,2),uspacer)
					case('right') ; write (iounit,fmt) tpFill(ln  ,uspacer)//blstr(:tw)
				end select; call advanceCarriage(1)
			enddo
			read (scrio,fmt) blstr(:tw); ln=pw-tw
			select case (trim(uadjust))
				case('left')  ; write (iounit,fmt) blstr(:tw)//tpFill(ln,uspacer)
				case('center'); write (iounit,fmt) tpFill(ln/2,uspacer)//blstr(:tw)//tpFill(ln/2+mod(ln,2),uspacer)
				case('right') ; write (iounit,fmt) tpFill(ln  ,uspacer)//blstr(:tw)
			end select; call advanceCarriage(1)
		else
			do i = 1,cnt
				blstr=tpFill(blstr); read (scrio,fmt) blstr(:tw); ln=pw-len_trim(blstr)
				select case (trim(uadjust))
					case('left')  ; write (iounit,fmt) trim(blstr)//tpFill(ln,uspacer)
					case('center'); write (iounit,fmt) tpFill(ln/2,uspacer)//trim(blstr)//tpFill(ln/2+mod(ln,2),uspacer)
					case('right') ; write (iounit,fmt) tpFill(ln  ,uspacer)//trim(blstr)
				end select; call advanceCarriage(1)
			enddo
		endif
		close (scrio)
		void=fcNullID(scrio)

		return
		end subroutine prLongText

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prTable(str,iounit,columns,width,indent,adjust,separator,spacer,interColWidth,equalWidth,saveOrder,err)
		implicit none

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=*), intent(in)            :: str
		integer*4        , intent(in)            :: iounit
		integer*4        , intent(in) , optional :: columns,width,indent,interColWidth
		character (len=*), intent(in) , optional :: adjust,separator,spacer
		logical*1        , intent(in) , optional :: saveOrder,equalWidth
		integer*4        , intent(out), optional :: err

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=len(str))                 :: ustr
		character (len=3)                        :: uadjust*6,useparator*3,uspacer*1
		integer*4                                :: ucolumns,uwidth,uindent,uinterColWidth,uerr
		logical*1                                :: usaveOrder,uequalWidth

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		character (len=5000)                     :: pstr
		type(uch), allocatable                   :: hold(:)

		integer*4, allocatable                   :: strLength(:),colLength(:),grid(:,:),lgrid(:,:),ngrid(:,:)

		integer*4                                :: N,wdth,k,l,c,ulines,ln
		logical*1                                :: odd
		character (len=3)                        :: switch

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		ucolumns=10      ; if (present(columns))       ucolumns      =columns
		uwidth  =78      ; if (present(width))         uwidth        =width
		uindent =4       ; if (present(indent))        uindent       =indent
		uinterColWidth=2 ; if (present(interColWidth)) uinterColWidth=interColWidth

		usaveOrder =true ; if (present(saveOrder))     usaveOrder    =saveOrder
		uequalWidth=false; if (present(equalWidth))    uequalWidth   =equalWidth

		uwidth=uwidth-uindent

		ustr      =str               ; ustr=trim(tpReduce(ustr))
		uadjust   =tpFill(uadjust)   ; uadjust='center'; if (present(adjust))    uadjust   =trim(adjust)
		useparator=tpFill(useparator); useparator=','  ; if (present(separator)) useparator=trim(separator)
		uspacer   =' '                                 ; if (present(spacer))    uspacer   =spacer

		if (present(spacer) .AND. (len(spacer).GT.1)) then
			uerr=-1; goto 666
		endif

		if (.NOT.tpIsStrInList(uadjust,['left','center','right']) ) then
			uerr=-2; goto 666
		endif

		if ( present(columns) .AND. present(width) ) then
			uerr=-3; goto 666
		endif

		if (.NOT.tpSplit(ustr,trim(useparator))) then
			uerr=-4; goto 666
		endif
		N=tpSplitLen

		allocate (hold(N),strLength(N)); strLength=0
		do k = 1,N
			hold(k)=uchSet( tpRetSplit(ustr,k) )
			strLength(k)=hold(k)%ln
		enddo

		if (present(columns)) uwidth=huge(uwidth)

		switch='000'
		if (present(columns)) switch(1:1)='1'
		if (usaveOrder)       switch(2:2)='1'
		if (uequalWidth)      switch(3:3)='1'

		select case (switch) ! collected to be more readable (despite code redundancy)
			case ('111') ! ucolumns defined | save order | equal width
				if (N.LT.ucolumns) ucolumns=N
				ulines=int(N/ucolumns)+1
				if (mod(N,ucolumns).EQ.0) ulines=ulines-1
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0
				colLength=maxval(strLength)
				do k = 1,N
					c=mod(k,ucolumns); l=int(k/ucolumns)+1
					if (c.EQ.0) then
						l=l-1; c=ucolumns
					endif

					grid(l,c)=k
				enddo

			case ('110') ! ucolumns defined | save order | NOT equal width
				if (N.LT.ucolumns) ucolumns=N
				ulines=int(N/ucolumns)+1
				if (mod(N,ucolumns).EQ.0) ulines=ulines-1
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0
				do k = 1,N
					c=mod(k,ucolumns)
					if (c.EQ.0) c=ucolumns
					colLength(c)=max(colLength(c),strLength(k))
				enddo
				do k = 1,N
					c=mod(k,ucolumns); l=int(k/ucolumns)+1
					if (c.EQ.0) then
						l=l-1; c=ucolumns
					endif

					grid(l,c)=k
				enddo

			case ('101') ! ucolumns defined | NOT save order | equal width
				if (N.LT.ucolumns) ucolumns=N
				ulines=int(N/ucolumns)+1
				if (mod(N,ucolumns).EQ.0) ulines=ulines-1
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0
				colLength=maxval(strLength)
				do k = 1,N
					c=mod(k,ucolumns); l=int(k/ucolumns)+1
					if (c.EQ.0) then
						l=l-1; c=ucolumns
					endif

					grid(l,c)=k
				enddo

			case ('100') ! ucolumns defined | NOT save order | NOT equal width
				if (N.LT.ucolumns) ucolumns=N
				ulines=int(N/ucolumns)+1
				if (mod(N,ucolumns).EQ.0) ulines=ulines-1
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0
				void=sort(strLength,rev=true)
				void=sort(hold     ,rev=true,attribute='length')
				c=0
				do k = 1,N,ulines
					c=c+1; colLength(c)=strLength(k)
				enddo
				c=1; l=0
				do k = 1,N
					l=l+1
					if (l.GT.ulines) then
						l=1; c=c+1
					endif
					grid(l,c)=k
				enddo

			case ('011') ! NOT ucolumns defined | save order | equal width
				ucolumns=int( uwidth/(maxval(strLength)+uinterColWidth) )

				if (ucolumns.EQ.0) then
					if (present(err)) err=-5
					deallocate (hold,strLength)
					return
				endif
				if (N.LT.ucolumns) ucolumns=N
				ulines=int(N/ucolumns)+1
				if (mod(N,ucolumns).EQ.0) ulines=ulines-1
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0
				colLength=maxval(strLength)
				do k = 1,N
					c=mod(k,ucolumns); l=int(k/ucolumns)+1
					if (c.EQ.0) then
						l=l-1; c=ucolumns
					endif

					grid(l,c)=k
				enddo

			case ('010') ! NOT ucolumns defined | save order | NOT equal width
				ucolumns=0

				allocate (lgrid(0:N,N))
				do ucolumns = N,1,-1
					ulines=int(N/ucolumns)+1
					if (mod(N,ucolumns).EQ.0) ulines=ulines-1

					lgrid=0
					c=0; l=1
					do k = N,1,-1
						c=c+1
						if (c.GT.ucolumns) then
							l=l+1; c=1
						endif
						lgrid(l,c)=strLength(k)
					enddo

					do c = 1,ucolumns
						lgrid(0,c)=maxval(lgrid(1:,c))
					enddo; wdth=sum(lgrid(0,:))+(ucolumns-1)*uinterColWidth

					if (wdth.LE.uwidth) exit
				enddo
				deallocate (lgrid)

				if (ucolumns.EQ.0) then
					uerr=-6; goto 666
				endif
				if (N.LT.ucolumns) ucolumns=N
				ulines=int(N/ucolumns)+1
				if (mod(N,ucolumns).EQ.0) ulines=ulines-1
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0
				do k = 1,N
					c=mod(k,ucolumns)
					if (c.EQ.0) c=ucolumns
					colLength(c)=max(colLength(c),strLength(k))
				enddo
				do k = 1,N
					c=mod(k,ucolumns); l=int(k/ucolumns)+1
					if (c.EQ.0) then
						l=l-1; c=ucolumns
					endif

					grid(l,c)=k
				enddo

			case ('001') ! NOT ucolumns defined | NOT save order | equal width
				ucolumns=int( uwidth/(maxval(strLength)+uinterColWidth) )
				if (ucolumns.EQ.0) then
					uerr=-7; goto 666
				endif
				if (N.LT.ucolumns) ucolumns=N
				ulines=int(N/ucolumns)+1
				if (mod(N,ucolumns).EQ.0) ulines=ulines-1
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0
				colLength=maxval(strLength)
				do k = 1,N
					c=mod(k,ucolumns); l=int(k/ucolumns)+1
					if (c.EQ.0) then
						l=l-1; c=ucolumns
					endif

					grid(l,c)=k
				enddo

			case ('000') ! NOT ucolumns defined | NOT save order | NOT equal width
				void=sort(strLength,rev=true)
				void=sort(hold     ,rev=true,attribute='length')
				ucolumns=0

				allocate (lgrid(0:N,N),ngrid(N,N))
				do ulines = 1,N
					ucolumns=int(N/ulines)+1
					if (mod(N,ulines).EQ.0) ucolumns=ucolumns-1

					lgrid=0; ngrid=0; odd=true
					l=0; c=1
					do k = 1,N
						if (     odd) l=l+1
						if (.NOT.odd) l=l-1

						if (l.GT.ulines) then
							c=c+1; l=ulines; odd=false
						endif

						if (l.LT.1) then
							c=c+1; l=1; odd=true
						endif

						lgrid(l,c)=strLength(k); ngrid(l,c)=k
					enddo

					do c = 1,ucolumns
						lgrid(0,c)=maxval(lgrid(1:,c))
					enddo; wdth=sum(lgrid(0,:))+(ucolumns-1)*uinterColWidth
					
					if (wdth.LE.uwidth) exit
				enddo

				if (ucolumns.EQ.0) then
					uerr=-8; goto 666
				endif
				allocate (colLength(ucolumns),grid(ulines,ucolumns)); colLength=0; grid=0

				grid(1:ulines,1:ucolumns)=ngrid(1:ulines,1:ucolumns)

				do c = 1,ucolumns
					colLength(c)=maxval( lgrid(:,c) )
				enddo
				deallocate (lgrid,ngrid)

		end select

		deallocate(strLength)

		void=collectArray(grid(:,ucolumns))

		carriageControl=0; internalIOunit=iounit
		do l = 1,ulines
			pstr=tpFill(pstr)

			wdth=0
			do c = 1,lenTrimArray( grid(l,:) )
				ln=colLength(c)-hold( grid(l,c) )%ln
				select case (trim(uadjust))
					case ('left')  ; pstr=pstr(1:wdth)//&
					                 uchGet( hold( grid(l,c) ) )//tpFill(uinterColWidth+ln,uspacer)
					case ('center'); pstr=pstr(1:wdth)//tpFill(ln/2,uspacer)//&
					                 uchGet( hold( grid(l,c) ) )//tpFill(uinterColWidth+ln/2+mod(ln,2),uspacer)
					case ('right') ; pstr=pstr(1:wdth)//tpFill(ln  ,uspacer)//&
					                 uchGet( hold( grid(l,c) ) )//tpFill(uinterColWidth,uspacer)
				end select
				wdth=wdth+colLength(c)+uinterColWidth
			enddo
			wdth=wdth-uinterColWidth; write (iounit,'(A)') tpFill(uindent,uspacer)//pstr(:wdth); call advanceCarriage(1)
		enddo

		deallocate(colLength,grid,hold)

		return

666		continue ! Error case.

		if (allocated(strLength)) deallocate(strLength)
		if (allocated(colLength)) deallocate(colLength)
		if (allocated(grid))      deallocate(grid)
		if (allocated(lgrid))     deallocate(lgrid)
		if (allocated(ngrid))     deallocate(ngrid)

		if (present(err)) err=uerr

		return
		end subroutine prTable

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prEchoFile(src,trgt,marker)
		implicit none

		integer*4        , intent(in)           :: src,trgt
		character (len=*), intent(in), optional :: marker

		integer*4, parameter   :: swidth=512,numwidth=4
		character (len=swidth) :: rline
		character (len=64)     :: umarker,fmt

		integer*4              :: k,sta,sto
		logical*1              :: numerate


		umarker=tpFill(umarker); umarker='input file>>'; if (present(marker)) umarker=marker

		numerate='%numeration' .in. umarker
		if (numerate) then
			sta=tpIndex(umarker,'%numeration')-1
			sto=tpIndex(umarker,'%numeration',endp=true)+1
			!fmt=tpFill(fmt)
			!fmt="('"//umarker(:sta)//"'"//",i"//prStrByVal(numwidth)//"."//prStrByVal(numwidth)//",1X,'"//umarker(sto+1:len_trim(umarker))//"',A)"
			!fmt=trim(tpReplace(umarker,'%numeration','i<numwidth>.<numwidth>,1X'))//tpFill(fmt)
		else
			!fmt=tpFill(fmt); fmt='(A)'
			continue
		endif

		!write (*,*) '#'//trim(fmt)//'#'
		!return

		rewind(src)
		k=0
		do
			if (eof(src)) exit
			rline=tpFill(rline)
			k=k+1
			read  (src ,'(A)') rline
			if (numerate) then
				fmt=tpFill(fmt); fmt='(A,i'//prStrByVal(numwidth)//'.'//prStrByVal(numwidth)//',A,A)'
				!write (trgt,trim(fmt)) k,trim(rline)
				write (trgt,fmt) umarker(:sta),k,umarker(sto:len(marker)),trim(rline)
			else
				fmt=tpFill(fmt); fmt='(A,A)'
				write (trgt,fmt) trim(umarker),trim(rline)
			endif
		enddo

		return
		end subroutine prEchoFile

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByLogical_l8(logi) result(ret)
		implicit none
		logical*8, intent(in)     :: logi
		character (len=mid(logi)) :: ret


		if (logi) then; ret='true'
		else;           ret='false'; endif; return
		end function strByLogical_l8
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByLogical_l4(logi) result(ret)
		implicit none
		logical*4, intent(in)     :: logi
		character (len=mid(logi)) :: ret


		if (logi) then; ret='true'
		else;           ret='false'; endif; return
		end function strByLogical_l4
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByLogical_l2(logi) result(ret)
		implicit none
		logical*2, intent(in)     :: logi
		character (len=mid(logi)) :: ret


		if (logi) then; ret='true'
		else;           ret='false'; endif; return
		end function strByLogical_l2
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByLogical_l1(logi) result(ret)
		implicit none
		logical*1, intent(in)     :: logi
		character (len=mid(logi)) :: ret


		if (logi) then; ret='true'
		else;           ret='false'; endif; return
		end function strByLogical_l1
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function nomStrByInt_i8(i) result(ret)
		implicit none
		integer*8, intent(in)  :: i
		character (len=mid(i)) :: ret


		write (ret,'(i<mid(i)>)') i; return
		end function nomStrByInt_i8
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function nomStrByInt_i4(i) result(ret)
		implicit none
		integer*4, intent(in)  :: i
		character (len=mid(i)) :: ret


		write (ret,'(i<mid(i)>)') i; return
		end function nomStrByInt_i4
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function nomStrByInt_i2(i) result(ret)
		implicit none
		integer*2, intent(in)  :: i
		character (len=mid(i)) :: ret


		write (ret,'(i<mid(i)>)') i; return
		end function nomStrByInt_i2
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function nomStrByInt_i1(i) result(ret)
		implicit none
		integer*1, intent(in)  :: i
		character (len=mid(i)) :: ret


		write (ret,'(i<mid(i)>)') i; return
		end function nomStrByInt_i1
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function defStrByInt_i8(i,defln) result(ret)
		implicit none
		integer*8, intent(in) :: i,defln
		character (len=defln) :: ret


		if (defln.LT.mid(i)) then; ret=repeat('*',len(ret)); return; endif
		write (ret,'(i<defln>)') i; return
		end function defStrByInt_i8
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function defStrByInt_i4(i,defln) result(ret)
		implicit none
		integer*4, intent(in) :: i,defln
		character (len=defln) :: ret


		if (defln.LT.mid(i)) then; ret=repeat('*',len(ret)); return; endif
		write (ret,'(i<defln>)') i; return
		end function defStrByInt_i4
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function defStrByInt_i2(i,defln) result(ret)
		implicit none
		integer*2, intent(in) :: i,defln
		character (len=defln) :: ret


		if (defln.LT.mid(i)) then; ret=repeat('*',len(ret)); return; endif
		write (ret,'(i<defln>)') i; return
		end function defStrByInt_i2
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function defStrByInt_i1(i,defln) result(ret)
		implicit none
		integer*1, intent(in) :: i,defln
		character (len=defln) :: ret


		if (defln.LT.mid(i)) then; ret=repeat('*',len(ret)); return; endif

		write (ret,'(i<defln>)') i; return
		end function defStrByInt_i1
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		function strByIntfmt_i8(i,fmt) result(ret)
		implicit none
		integer*8        , intent(in) :: i
		character (len=*), intent(in) :: fmt
		character (len=fmtPrintLen)         :: ret


		iMid=mid(i); arrType=2

		void=parsfmt(fmt); ufmt='('//trim(ufmt)//')'
		write (ret,ufmt) i; return
		end function strByIntfmt_i8
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		function strByIntfmt_i4(i,fmt) result(ret)
		implicit none
		integer*4        , intent(in) :: i
		character (len=*), intent(in) :: fmt
		character (len=fmtPrintLen)         :: ret


		iMid=mid(i); arrType=2

		void=parsfmt(fmt); ufmt='('//trim(ufmt)//')'
		write (ret,ufmt) i; return
		end function strByIntfmt_i4
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		function strByIntfmt_i2(i,fmt) result(ret)
		implicit none
		integer*2        , intent(in) :: i
		character (len=*), intent(in) :: fmt
		character (len=fmtPrintLen)         :: ret


		iMid=mid(i); arrType=2

		void=parsfmt(fmt); ufmt='('//trim(ufmt)//')'
		write (ret,ufmt) i; return
		end function strByIntfmt_i2
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		function strByIntfmt_i1(i,fmt) result(ret)
		implicit none
		integer*1        , intent(in) :: i
		character (len=*), intent(in) :: fmt
		character (len=fmtPrintLen)         :: ret


		iMid=mid(i); arrType=2

		void=parsfmt(fmt); ufmt='('//trim(ufmt)//')'
		write (ret,ufmt) i; return
		end function strByIntfmt_i1
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByRealf_r8(r,indent,accuracy) result(ret)
		implicit none
		real*8   , intent(in)             :: r
		integer*4, intent(in)             :: indent,accuracy
		character (len=indent+accuracy+1) :: ret


		write (ret,'(F<indent+accuracy+1>.<accuracy>)') r; return
		end function strByRealf_r8
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByRealf_r4(r,indent,accuracy) result(ret)
		implicit none
		real*4   , intent(in)             :: r
		integer*4, intent(in)             :: indent,accuracy
		character (len=indent+accuracy+1) :: ret


		write (ret,'(F<indent+accuracy+1>.<accuracy>)') r; return
		end function strByRealf_r4
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByReale_r8(r,indent,accuracy,form) result(ret)
		implicit none
		real*8   , intent(in)             :: r
		integer*4, intent(in)             :: indent,accuracy
		character (len=*), intent(in)     :: form
		character (len=indent+accuracy+7) :: ret

		
		if (form.NE.'exp') then
			ret=''; return
		endif

		write (ret,'(ES<indent+accuracy+7>.<accuracy>)') r; return
		end function strByReale_r8
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function strByReale_r4(r,indent,accuracy,form) result(ret)
		implicit none
		real*4   , intent(in)             :: r
		integer*4, intent(in)             :: indent,accuracy
		character (len=*), intent(in)     :: form
		character (len=indent+accuracy+7) :: ret


		if (form.NE.'exp') then
			ret=''; return
		endif

		write (ret,'(ES<indent+accuracy+7>.<accuracy>)') r; return
		end function strByReale_r4
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		function strByRealfmt_r8(r,fmt) result(ret)
		implicit none
		real*8           , intent(in) :: r
		character (len=*), intent(in) :: fmt
		character (len=fmtPrintLen)   :: ret


		rMid=mid(r); arrType=1

		void=parsfmt(fmt); ufmt='('//trim(ufmt)//')'
		!write (*,*) '$'//trim(ufmt)//'$'
		write (ret,ufmt) r; return
		end function strByRealfmt_r8
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		function strByRealfmt_r4(r,fmt) result(ret)
		implicit none
		real*4           , intent(in) :: r
		character (len=*), intent(in) :: fmt
		character (len=fmtPrintLen)   :: ret


		rMid=mid(r); arrType=1

		void=parsfmt(fmt); ufmt='('//trim(ufmt)//')'
		write (ret,ufmt) r; return
		end function strByRealfmt_r4
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function parsfmt(fmt)
		implicit none

		character (Len=*), optional :: fmt
		character (Len=maxfmtlen)   :: tfmt
		character                   :: keycode*6,prefixPlus*3
		logical*1                   :: fndplus,fndpower,fndexp,fndzero,fnddot,fndund
		integer*4                   :: dPos,ePos,zBeforeDot,zBeforeExp,zBeforeEnd,&
		                               iPlus,fixedPresc,nUnder


		if (present(fmt)) then
			ufmt=repeat(' ',maxfmtlen)
			prefixPlus=repeat(' ',len(prefixPlus))
			tfmt=tpLowerCase(fmt)

			fndund  =Index(tfmt,'_').GT.0
			fndplus =Index(tfmt,'+').GT.0
			fnddot  =Index(tfmt,'.').GT.0
			fndexp  =Index(tfmt,'e').GT.0
			fndpower=Index(tfmt,'^').GT.0
			fndzero =Index(tfmt,'0').GT.0

			iPlus=0
			if (fndplus) then; prefixPlus='sp,'; iPlus=1; endif

			! und|plus|dot|exp|power|zero
			keycode='000000'
			if (fndund)   keycode(1:1)='1'
			if (fndplus)  keycode(2:2)='1'
			if (fnddot)   keycode(3:3)='1'
			if (fndexp)   keycode(4:4)='1'
			if (fndpower) keycode(5:5)='1'
			if (fndzero)  keycode(6:6)='1'

			select case (arrType)
				case (1)
					if (.NOT.fnddot) then
						void=defFmt()
					else
						select case (keycode(4:6))
							case ('000') ! nothing useful in fmt
								if (fndund) then
									dPos=Index(tfmt,'.')
									nUnder=tpCount( tfmt(:dPos),'_' )
									vLen=nUnder+iPlus+1
									ufmt='F'//prStrByVal(vLen)//'.'//prStrByVal(0)
									dotShift=nUnder+iPlus
								else
									void=defFmt()
								endif
							case ('001') ! defined by user fmt (000.000)
								dPos=Index(tfmt,'.')
								nUnder=tpCount( tfmt(:dPos),'_' )
								zBeforeDot=tpCount( tfmt(:dPos),'0' )
								zBeforeExp=tpCount( tfmt(dPos:),'0' )
								vLen=iPlus+zBeforeDot+1+zBeforeExp+nUnder
								ufmt='F'//prStrByVal(vLen)//'.'//prStrByVal(zBeforeExp)
								dotShift=iPlus+zBeforeDot
							case ('010') ! output only integer quantity with a dot
								vLen=iPlus+rMid+1
								ufmt='F'//prStrByVal(vLen)//'.'//prStrByVal(0)
								dotShift=vLen-1
							case ('011') ! defined by user + avoid conversion error
								dPos=Index(tfmt,'.')
								zBeforeExp=tpCount( tfmt(dPos:),'0' )
								vLen=1+rMid+1+zBeforeExp
								ufmt='F'//prStrByVal(vLen)//'.'//prStrByVal(zBeforeExp)
								dotShift=rMid+1
							case ('100') ! unacceptable fmt (because zeros are not present)
								void=defFmt()
							case ('101')
								dPos=Index(tfmt,'.')
								ePos=Index(tfmt,'e')
								zBeforeDot=tpCount( tfmt(:dPos)    ,'0' )
								zBeforeExp=tpCount( tfmt(dPos:ePos),'0' )
								zBeforeEnd=tpCount( tfmt(ePos:)    ,'0' )
								vLen=zBeforeDot+1+zBeforeExp+1+1+zBeforeEnd+1
								fixedPresc=zBeforeExp+zBeforeDot-1
								ufmt=prStrByVal(zBeforeDot)//'PE'//prStrByVal(vLen)//'.'//prStrByVal(fixedPresc)//'E'//prStrByVal(zBeforeEnd)
								dotShift=zBeforeDot+1
							case ('110') ! unacceptable fmt (because zeros are not present)
								void=defFmt()
							case ('111') ! unacceptable combination (ignore power)
								dPos=Index(tfmt,'.')
								ePos=Index(tfmt,'e')
								zBeforeDot=tpCount( tfmt(:dPos)    ,'0' )
								zBeforeExp=tpCount( tfmt(dPos:ePos),'0' )
								zBeforeEnd=tpCount( tfmt(ePos:)    ,'0' )
								vLen=zBeforeDot+1+zBeforeExp+1+1+zBeforeEnd+1
								fixedPresc=zBeforeExp+zBeforeDot-1
								ufmt=prStrByVal(zBeforeDot)//'PE'//prStrByVal(vLen)//'.'//prStrByVal(fixedPresc)//'E'//prStrByVal(zBeforeEnd)
								dotShift=zBeforeDot+1
						end select
						ufmt=trim(prefixPlus)//trim(ufmt)
					endif
				case (2)
					select case (keycode(5:6))
						case ('00')
							if (fndund) then
								nUnder=tpCount( tfmt,'_' )
								vLen=nUnder+iPlus
								ufmt='i'//prStrByVal(vLen)
							else
								void=defFmt()
							endif
						case ('01')
							zBeforeDot=tpCount(tfmt,'0')
							nUnder=tpCount( tfmt,'_' )
							if (zBeforeDot.GT.defiDigits-1) zBeforeDot=defiDigits-1
							vLen=defiDigits+nUnder
							ufmt='i'//prStrByVal(vLen)//'.'//prStrByVal(zBeforeDot)
						case ('10')
							vLen=iMid+1
							ufmt='i'//prStrByVal(vLen)
						case ('11')
							vLen=iMid+1
							ufmt='i'//prStrByVal(iMid+1)//'.'//prStrByVal(vLen)
					end select
					ufmt=trim(prefixPlus)//trim(ufmt)
					dotShift=vLen
			end select
		else
			void=defFmt()
		endif


		parsfmt=0; return
		end function parsfmt

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine advanceCarriage(slen)
		implicit none
		integer*4, optional :: slen


		if (present(slen)) &
		  then; carriageControl=carriageControl+slen
		  else; carriageControl=carriageControl+1   ; endif

		return
		end subroutine advanceCarriage

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function defFmt()
		implicit none

		select case (arrType)
			case (1); ufmt='F'//prStrByVal(defrDigits)//'.'//prStrByVal(defrPresc); vLen=defrDigits; dotShift=defrDigits
			case (2); ufmt='i'//prStrByVal(defiDigits)                            ; vLen=defiDigits
		end select

		defFmt=0; return
		end function defFmt

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prUndoOutput(slen,full)
		implicit none

		integer*4, optional :: slen
		logical*1, optional :: full

		integer*4           :: i

		if (present(full)) then
			if (full) then
				do i = 1,carriageControl
					backspace(internalIOunit)
				enddo
				carriageControl=0
			endif
		else
			if (present(slen)) then
				do i = 1,slen
					backspace(internalIOunit)
				enddo
				carriageControl=carriageControl-slen
			else
				backspace(internalIOunit)
				carriageControl=carriageControl-1
			endif					
		endif

		return
		end subroutine prUndoOutput

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prFinalize
		implicit none


		!continue

		return
		end subroutine prFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module printmod

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!    "0000" number of zeros defines number of digits always to be printed (up to [defiDigits-1] zeros)
!    "+"    will cause always to present the sign of value.
!    "^"    will give guarantee that value will not cause conversion error.
!    "_"    vacancy for additional digits before dot.
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!
!   integer:
!            '+000':
!               1 -> '+001'
!             -20 -> '-020'
!            2000 -> conversion error(63) & '****' as result.
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '+___0':
!               1 -> '   +1'
!             -20 -> '  -20'
!            2000 -> '+2000'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '^':
!               case | maxval(abs(array))=10000
!                   1 -> '     1'
!                 -20 -> '   -20'
!                2000 -> '  2000'
!               10000 -> ' 10000'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '0^': 
!               case | maxval(abs(array))=10000
!                   1 -> ' 00001'
!                 -20 -> '-00001'
!                2000 -> ' 02000'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '+^':
!               case | maxval(abs(array))=10000
!                   1 -> '    +1'
!                 -20 -> '   -20'
!                2000 -> ' +2000'
!               10000 -> '+10000'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '+0^' = '+0000^':
!               case | maxval(abs(array))=10000
!                   1 -> '+00001'
!                 -20 -> '-00020'
!                2000 -> '+02000'
!               10000 -> '+10000'
!               case | maxval(abs(array))=2000
!                   1 -> '+0001'
!                 -20 -> '-0020'
!                2000 -> '+2000'
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   real:
!            '000.0000':
!                            0.1 -> '   0.1000'
!                              1 -> '   1.0000'
!                            100 -> ' 100.0000'
!                           -100 -> '-100.0000'
!                           1000 -> conversion error(63) & '********' as result.
!                          -1000 -> conversion error(63) & '********' as result.
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '_000.0000':
!                            0.1 -> '    0.1000'
!                              1 -> '    1.0000'
!                            100 -> '  100.0000'
!                           -100 -> ' -100.0000'
!                           1000 -> ' 1000.0000'
!                          -1000 -> -> conversion error(63) & '*********' as result.
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!           '+000.0000':
!                            0.1 -> '  +0.1000'
!                              1 -> '  +1.0000'
!                            100 -> '+100.0000'
!                           -100 -> '-100.0000'
!                           1000 -> conversion error(63) & '*********' as result.
!                          -1000 -> conversion error(63) & '*********' as result.
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!           '+_000.0000':
!                            0.1 -> '   +0.1000'
!                              1 -> '   +1.0000'
!                            100 -> ' +100.0000'
!                           -100 -> ' -100.0000'
!                           1000 -> ' 1000.0000'
!                          -1000 -> '+1000.0000'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!           '^.0000' = '^00000.0000' = '^0.0000':
!               case | maxval(abs(array)) = 10000
!                            0.1 -> '     0.1000'
!                              1 -> '     1.0000'
!                            100 -> '   100.0000'
!                           -100 -> '  -100.0000'
!                           1000 -> '  1000.0000'
!                          -1000 -> ' -1000.0000'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!           '^' = '^.':
!               case | maxval(abs(array)) = 10000
!                            0.1 -> '     0.'
!                              1 -> '     1.'
!                            100 -> '   100.'
!                           -100 -> '  -100.'
!                           1000 -> '  1000.'
!                          -1000 -> ' -1000.'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!           '+^0.0000' = '+^00000.0000':
!               case | maxval(abs(array)) = 10000
!                            0.1 -> '    +0.1000'
!                              1 -> '    +1.0000'
!                            100 -> '  +100.0000'
!                           -100 -> '  -100.0000'
!                           1000 -> ' +1000.0000'
!                          -1000 -> ' -1000.0000'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!              '0.0000E00':
!                            0.1 -> ' 1.0000E-01'
!                              1 -> ' 1.0000E+00'
!                            100 -> ' 1.0000E+02'
!                           -100 -> '-1.0000E+02'
!                           1000 -> ' 1.0000E+03'
!                          -1000 -> '-1.0000E+03'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '000.0000E00':
!                            0.1 -> ' 100.0000E-03'
!                              1 -> ' 100.0000E-02'
!                            100 -> ' 100.0000E+00'
!                           -100 -> '-100.0000E+00'
!                           1000 -> ' 100.0000E+01'
!                          -1000 -> '-100.0000E+01'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!            '000.0000E000':
!                            0.1 -> ' 100.0000E-003'
!                              1 -> ' 100.0000E-002'
!                            100 -> ' 100.0000E+000'
!                           -100 -> '-100.0000E+000'
!                           1000 -> ' 100.0000E+001'
!                          -1000 -> '-100.0000E+001'
!          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!           '+000.0000E000':
!                            0.1 -> '+100.0000E-003'
!                              1 -> '+100.0000E-002'
!                            100 -> '+100.0000E+000'
!                           -100 -> '-100.0000E+000'
!                           1000 -> '+100.0000E+001'
!                          -1000 -> '-100.0000E+001'
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
