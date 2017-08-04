	module datablock

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	use glob    , only: uch,uchGet,uchSet,find,collectArray,true,false,void,voidl
	use txtParser
	use fcontrol, only: fcNewID,fcNullID,fcBanID
	use printmod, only: prLongText,prStrByVal,prTable,prEchoFile

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: bdVersion='1.550'
	character (len=*), parameter :: bdDate   ='2017.08.03'
	character (len=*), parameter :: bdAuthor ='Anton B. Zakharov'

	integer*4, parameter :: typeSetSize=12,mxKindNameLen=9,maxListSetLen=100
	integer*4, parameter :: maxSetLen=2500,expectFileWidth=200
	integer*4, parameter :: prepareStringLen=5000
	integer*4, parameter :: storageLen=32,markUpChLen=64

	character (len=*), parameter :: defCommentCh='#',defStartCh='%name{',defEndCh='}',&
	                                defAccordCh='=',defSeparatorCh='|'

	character (len=mxKindNameLen) :: kindLabels(typeSetSize)
	data kindLabels /'real*8','real*4',&
	&                'integer*8','integer*4','integer*2','integer*1',&
	&                'logical*8','logical*4','logical*2','logical*1',&
	&                'character','type(uch)'/

	character (len=2) :: shortKindLabels(typeSetSize)
	data shortKindLabels / 're','re','in','in','in','in','lo','lo','lo','lo','ch','uc' /

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type variable
		type(uch)          :: name
		integer*8          :: address
		integer*4          :: position,kind,vlen,nvals
		logical*1          :: opt,variable,several

		real*8   , pointer :: pntr8
		real*4   , pointer :: pntr4
		integer*8, pointer :: pnti8
		integer*4, pointer :: pnti4
		integer*2, pointer :: pnti2
		integer*1, pointer :: pnti1
		logical*8, pointer :: pntl8
		logical*4, pointer :: pntl4
		logical*2, pointer :: pntl2
		logical*1, pointer :: pntl1
		character, pointer :: pntch
		type(uch), pointer :: pntuch

		character (len=storageLen) :: defstor
		type(uch)                  :: defc,expect,concatCh

		character (len=storageLen) :: varsta,varsto,varste
	end type variable

	type bd
		type(uch)              :: name,bdStr,commentCh,startCh,endCh,accordCh,separatorCh
		integer*4              :: arrayLen,associatedIO
		integer*4, allocatable :: varSet(:)
		logical*1              :: freeBlock
	end type bd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type(variable) :: variableSet(maxSetLen)
	type(bd)       :: bdSet(maxListSetLen)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=storageLen)       :: storeTemplate
	character (len=prepareStringLen) :: holdBlockStr,workBlockStr

	integer*4                        :: varAppend=0,bdAppend=0,iounit=6,errunit=6

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	interface bdShareVariable
		module procedure shareVariable_r8,shareVariable_r4,shareVariable_i8,&
						 shareVariable_i4,shareVariable_i2,shareVariable_i1,&
						 shareVariable_l8,shareVariable_l4,shareVariable_l2,&
						 shareVariable_l1,shareVariable_c ,shareVariable_uc
	end interface bdShareVariable

	interface modify
		module procedure modify_c,modify_r,modify_i,modify_l
	end interface modify

	interface bdCollect
		module procedure bdCollectName,bdCollectAddress
	end interface bdCollect

	interface bdAddVariable
		module procedure bdAddVariableAddress,bdAddVariableName
	end interface bdAddVariable

	interface bdRemoveVariable
		module procedure bdRemoveVariableAddress,bdRemoveVariableName
	end interface bdRemoveVariable

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private

	public :: bdVersion,bdDate,bdAuthor

	public :: bdShareVariable,bdCollect,bdAddVariable,bdRemoveVariable,bdParseFile,&
	          bdPrintBlockData,bdSetIOunit,bdSetERRunit,bdFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdCollectAddress(name,array,startCh,endCh,commentCh,accordCh,separatorCh,freeBlock,associatedIO) result(rcode)
	implicit none

	character (len=*)           :: name
	integer*8, optional         :: array(:)

	character (len=*), optional :: startCh,endCh,commentCh,accordCh,separatorCh
	logical*1        , optional :: freeBlock
	integer*4        , optional :: associatedIO

	character (len=markUpChLen) :: ustartCh,uendCh,ucommentCh,uaccordCh,useparatorCh
	logical*1                   :: ufreeBlock
	integer*4                   :: uassociatedIO

	integer*4                   :: arrayLen,i,j


	ustartCh=tpFill(ustartCh); ustartCh=defStartCh
	if (present(startCh)) ustartCh=trim(startCh)
	
	uendCh=tpFill(uendCh); uendCh=defEndCh
	if (present(endCh)) uendCh=trim(endCh)
	
	ucommentCh=tpFill(ucommentCh); ucommentCh=defCommentCh
	if (present(commentCh)) ucommentCh=trim(commentCh)

	uaccordCh=tpFill(uaccordCh); uaccordCh=defAccordCh
	if (present(accordCh)) uaccordCh=trim(accordCh)

	useparatorCh=tpFill(useparatorCh); useparatorCh=defSeparatorCh
	if (present(separatorCh)) useparatorCh=trim(separatorCh)

	ufreeBlock=false; if (present(freeBlock))    ufreeBlock=freeBlock
	uassociatedIO=0;  if (present(associatedIO)) uassociatedIO=associatedIO

	if (ufreeBlock .AND. .NOT.present(associatedIO) ) then
		write (errunit,*) name//': If present freeBlock option, must present associatedIO value.'
		rcode=-1
		return
	endif

	if (find(bdSet(1:bdAppend)%name,trim(name)).GT.0) then
		write (errunit,*) name//': Block data with such name is already declared.'
		rcode=-2
		return
	endif

	if ( (trim(ucommentCh)).EQ.(trim(uaccordCh)) ) then
		write (errunit,*) name//': Comment and accordance markers must be different.'
		rcode=-3
		return
	endif

	if ( (trim(useparatorCh)).EQ.(trim(uaccordCh)) ) then
		write (errunit,*) name//': Separator and accordance markers must be different.'
		rcode=-4
		return
	endif

	if ( (trim(ucommentCh)).EQ.(trim(useparatorCh)) ) then
		write (errunit,*) name//': Comment and separator markers must be different.'
		rcode=-5
		return
	endif

	if ( (.NOT.present(array)) .AND. (.NOT.freeBlock) ) then
		write (errunit,*) name//': Must contain at least freeBlock or array of addresses.'
		rcode=-6
		return
	endif

	arrayLen=UBound(array,1); bdAppend=bdAppend+1; rcode=bdAppend

	bdSet(bdAppend)%name=uchSet(trim(adjustl(name)))
	bdSet(bdAppend)%startCh=uchSet(trim(tpReplace(ustartCh,'%name',uchGet(bdSet(bdAppend)%name))))
	bdSet(bdAppend)%commentCh=uchSet(trim(ucommentCh))
	bdSet(bdAppend)%accordCh=uchSet(trim(uaccordCh))
	bdSet(bdAppend)%separatorCh=uchSet(trim(useparatorCh))
	bdSet(bdAppend)%endCh=uchSet(trim(tpReplace(uendCh,'%name',uchGet(bdSet(bdAppend)%name))))
	bdSet(bdAppend)%freeBlock=ufreeBlock
	bdSet(bdAppend)%associatedIO=uassociatedIO

	if (ufreeBlock) return

	allocate ( bdSet(bdAppend)%varSet(arrayLen) )

	bdSet(bdAppend)%arrayLen=0

	do i = 1,arrayLen
		j=find(variableSet(1:varAppend)%address,array(i))
		if (j.LT.0) then
			write (errunit,'(A,1X,i,A)') 'Not shared (unknown) variable at adress',array(i),'.'
			stop
		endif
		bdSet(bdAppend)%arrayLen=bdSet(bdAppend)%arrayLen+1
		bdSet(bdAppend)%varSet(bdSet(bdAppend)%arrayLen)=j
	enddo

	return
	end function bdCollectAddress

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdCollectName(name,array,startCh,endCh,commentCh,accordCh,separatorCh,freeBlock,associatedIO) result(rcode)
	implicit none

	character (len=*)           :: name
	character (len=*)           :: array(:)

	character (len=*), optional :: startCh,endCh,commentCh,accordCh,separatorCh
	logical*1        , optional :: freeBlock
	integer*4        , optional :: associatedIO

	character (len=markUpChLen) :: ustartCh,uendCh,ucommentCh,uaccordCh,useparatorCh
	logical*1                   :: ufreeBlock
	integer*4                   :: uassociatedIO

	integer*4                   :: arrayLen,i,j


	ustartCh=tpFill(ustartCh); ustartCh=defStartCh
	if (present(startCh)) ustartCh=trim(startCh)
	
	uendCh=tpFill(uendCh); uendCh=defEndCh
	if (present(endCh)) uendCh=trim(endCh)
	
	ucommentCh=tpFill(ucommentCh); ucommentCh=defCommentCh
	if (present(commentCh)) ucommentCh=trim(commentCh)

	uaccordCh=tpFill(uaccordCh); uaccordCh=defAccordCh
	if (present(accordCh)) uaccordCh=trim(accordCh)

	useparatorCh=tpFill(useparatorCh); useparatorCh=defSeparatorCh
	if (present(separatorCh)) useparatorCh=trim(separatorCh)

	ufreeBlock=false; if (present(freeBlock))    ufreeBlock=freeBlock
	uassociatedIO=0;  if (present(associatedIO)) uassociatedIO=associatedIO

	if (ufreeBlock .AND. .NOT.present(associatedIO) ) then
		write (errunit,*) name//': If present freeBlock option, must present associatedIO value.'
		rcode=-1
		return
	endif

	if (find(bdSet(1:bdAppend)%name,trim(name)).GT.0) then
		write (errunit,*) name//': Block data with such name is already declared.'
		rcode=-2
		return
	endif

	if ( (trim(ucommentCh)).EQ.(trim(uaccordCh)) ) then
		write (errunit,*) name//': Comment and accordance markers must be different.'
		rcode=-3
		return
	endif

	if ( (trim(useparatorCh)).EQ.(trim(uaccordCh)) ) then
		write (errunit,*) name//': Separator and accordance markers must be different.'
		rcode=-4
		return
	endif

	if ( (trim(ucommentCh)).EQ.(trim(useparatorCh)) ) then
		write (errunit,*) name//': Comment and separator markers must be different.'
		rcode=-5
		return
	endif

	arrayLen=UBound(array,1); bdAppend=bdAppend+1; rcode=bdAppend

	bdSet(bdAppend)%name=uchSet(trim(adjustl(name)))
	bdSet(bdAppend)%startCh=uchSet(trim(tpReplace(ustartCh,'%name',uchGet(bdSet(bdAppend)%name))))
	bdSet(bdAppend)%commentCh=uchSet(trim(ucommentCh))
	bdSet(bdAppend)%accordCh=uchSet(trim(uaccordCh))
	bdSet(bdAppend)%separatorCh=uchSet(trim(useparatorCh))
	bdSet(bdAppend)%endCh=uchSet(trim(tpReplace(uendCh,'%name',uchGet(bdSet(bdAppend)%name))))
	bdSet(bdAppend)%freeBlock=ufreeBlock
	bdSet(bdAppend)%associatedIO=uassociatedIO

	allocate ( bdSet(bdAppend)%varSet(arrayLen) )

	bdSet(bdAppend)%arrayLen=0

	do i = 1,arrayLen
		j=find(variableSet(1:varAppend)%name,array(i))
		if (j.LT.0) then
			write (errunit,'(A,1X,i,A)') 'Not shared (unknown) variable with name ',array(i),'.'
			stop
		endif
		bdSet(bdAppend)%arrayLen=bdSet(bdAppend)%arrayLen+1
		bdSet(bdAppend)%varSet(bdSet(bdAppend)%arrayLen)=j
	enddo

	return
	end function bdCollectName

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdAddVariableAddress(bdname,address) result(rcode)
	implicit none

	character (len=*)      :: bdname
	integer*8              :: address
	integer*4              :: j,k,varln
	integer*4, allocatable :: hold(:)


	rcode=0
	j=find(bdSet(1:bdAppend)%name,trim(bdname))
	if (j.LT.0) then
		write (errunit,'(A)') trim(bdname)//': Block data with such name has not been declared.'
		rcode=-1
		return
	endif

	k=find(variableSet(1:varAppend)%address,address)
	if (k.LT.0) then
		write (errunit,'(A,1X,i,A)') 'Not shared (unknown) variable at adress ',address,'.'
		rcode=-2
		return
	endif

	varln=bdSet(j)%arrayLen; allocate (hold(varln))

	hold=bdSet(j)%varSet; deallocate (bdSet(j)%varSet)

	allocate (bdSet(j)%varSet(varln+1))

	bdSet(j)%varSet(1:varln)=hold
	bdSet(j)%varSet(varln+1)=variableSet(k)%address
	bdSet(j)%arrayLen=bdSet(j)%arrayLen+1

	deallocate (hold)

	return
	end function bdAddVariableAddress

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdAddVariableName(bdname,varname) result(rcode)
	implicit none

	character (len=*)      :: bdname,varname
	integer*4              :: j,k,varln
	integer*4, allocatable :: hold(:)


	rcode=0
	j=find(bdSet(1:bdAppend)%name,trim(bdname))
	if (j.LT.0) then
		write (errunit,'(A)') trim(bdname)//': Block data with such name has not been declared.'
		rcode=-1
		return
	endif

	k=find(variableSet(1:varAppend)%name,varname)
	if (k.LT.0) then
		write (errunit,'(A,1X,i,A)') 'Not shared (unknown) variable with name ',varname,'.'
		rcode=-2
		return
	endif

	varln=bdSet(j)%arrayLen; allocate (hold(varln))

	hold=bdSet(j)%varSet; deallocate (bdSet(j)%varSet)

	allocate (bdSet(j)%varSet(varln+1))

	bdSet(j)%varSet(1:varln)=hold
	bdSet(j)%varSet(varln+1)=variableSet(k)%address
	bdSet(j)%arrayLen=bdSet(j)%arrayLen+1

	deallocate (hold)

	return
	end function bdAddVariableName

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdRemoveVariableAddress(bdname,address) result(rcode)
	implicit none

	character (len=*)      :: bdname
	integer*8              :: address
	integer*4              :: j,k,l,varln
	integer*4, allocatable :: hold(:)


	rcode=0
	j=find(bdSet(1:bdAppend)%name,trim(bdname))
	if (j.LT.0) then
		write (errunit,'(A)') trim(bdname)//': Block data with such name has not been declared.'
		rcode=-1
		return
	endif

	k=find(variableSet(1:varAppend)%address,address)
	if (k.LT.0) then
		write (errunit,'(A,1X,i,A)') 'Not shared (unknown) variable at adress ',address,'.'
		rcode=-2
		return
	endif

	varln=bdSet(j)%arrayLen; allocate (hold(varln))

	hold=bdSet(j)%varSet; deallocate (bdSet(j)%varSet)

	l=find(hold,k); if (l.GT.0) hold(l)=0
	varln=collectArray(hold)

	allocate (bdSet(j)%varSet(varln))

	bdSet(j)%varSet(1:varln)=hold(1:varln)
	bdSet(j)%arrayLen=varln

	deallocate (hold)

	return
	end function bdRemoveVariableAddress

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdRemoveVariableName(bdname,varname) result(rcode)
	implicit none

	character (len=*)      :: bdname,varname
	integer*4              :: j,k,l,varln
	integer*4, allocatable :: hold(:)


	rcode=0
	j=find(bdSet(1:bdAppend)%name,trim(bdname))
	if (j.LT.0) then
		write (errunit,'(A)') trim(bdname)//': Block data with such name has not been declared.'
		rcode=-1
		return
	endif

	k=find(variableSet(1:varAppend)%name,varname)
	if (k.LT.0) then
		write (errunit,'(A,1X,i,A)') 'Not shared (unknown) variable with name ',varname,'.'
		rcode=-2
		return
	endif

	varln=bdSet(j)%arrayLen; allocate (hold(varln))

	hold=bdSet(j)%varSet; deallocate (bdSet(j)%varSet)

	l=find(hold,k); if (l.GT.0) hold(l)=0
	varln=collectArray(hold)

	allocate (bdSet(j)%varSet(varln))

	bdSet(j)%varSet(1:varln)=hold(1:varln)
	bdSet(j)%arrayLen=varln

	deallocate (hold)

	return
	end function bdRemoveVariableName

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdParseFile(name,filename,expectWidth) result(rcode)
	implicit none

	character (len=*)   :: name,filename
	integer*4, optional :: expectWidth
	integer*4           :: dexpectWidth

	integer*4           :: bdPos,startFnd,stopFnd,bdStart,bdStop,i,fnd,err,id
	logical*1           :: atTheStart,atTheEnd,ppp
	real*8              :: sta,sto


	rcode=0
	dexpectWidth=expectFileWidth; if (present(expectWidth)) dexpectWidth=expectWidth

	bdPos=find(bdSet(1:bdAppend)%name,trim(name))

	if (bdPos.LT.0) then
		write (errunit,'(A)') trim(name)//': Block data with such name has not been declared.'
		rcode=-1
		return
	endif

	atTheStart=tpIndex( uchGet(bdSet(bdPos)%startCh), uchGet(bdSet(bdPos)%name) ).NE.0
	atTheEnd  =tpIndex( uchGet(bdSet(bdPos)%endCh)  , uchGet(bdSet(bdPos)%name) ).NE.0

	i=fcNewID(); open (i,file=filename,status='old',iostat=err)
	if (err.NE.0) then
		write (errunit,'(A)') 'bd: '//uchGet(bdSet(bdPos)%name)//':'//" File does not exist: "//trim(filename)

		rcode=-2
		return
	else
		close (i); void=fcNullID(i)
	endif

	call tpEnable(dexpectWidth); voidl=tpJoin(filename)
	if (.NOT.voidl) then
		write (errunit,'(A)') 'bd: '//uchGet(bdSet(bdPos)%name)//':'
		write (errunit,'(A)') 'File access denied: '//trim(filename)
		rcode=-3
		return
	endif

	tpCaseSens=false; tpReplTabs=true; tpReduSpac=true;  tpIgnoSpac=false
	tpIgnoTabs=false; tpIgnoComm=true; tpAlloQuot=false

	!write (*,*) 'Start searching'; call cpu_time(sta)
	void=tpSetCommentMark( uchGet( bdSet(bdPos)%commentCh ) )
	startFnd=tpLocate( uchGet(bdSet(bdPos)%startCh),true ); bdStart=tpPointer+1
	stopFnd =tpLocate( uchGet(bdSet(bdPos)%endCh) );        bdStop =tpPointer+1
	!write (*,*) 'End searching'; call cpu_time(sto)
	!write (*,*) 'Spent',sto-sta
	!if (bdSet(bdPos)%freeBlock) stop

	holdBlockStr=tpFill(holdBlockStr)
	if ( (startFnd.EQ.-1).AND.(startFnd.EQ.-1) ) then
		holdBlockStr=uchGet(bdSet(bdPos)%startCh)//uchGet(bdSet(bdPos)%endCh)
	else
		if (bdSet(bdPos)%freeBlock) then
			id=bdSet(bdPos)%associatedIO; void=fcBanID(id)
			open (id,status='scratch') !; write (id,'(A)') uchGet(bdSet(bdPos)%bdStr); rewind(id)

			workBlockStr=tpFill(workBlockStr); workBlockStr=tpGetLineByPointer(bdStart)
			fnd=tpIndex( workBlockStr,uchGet(bdSet(bdPos)%startCh),endp=true )+1
			write (id,'(A)') trim( workBlockStr(fnd:) )

			do i = bdStart+1,bdStop
				workBlockStr=tpFill(workBlockStr); workBlockStr=tpGetLineByPointer(i)

				if (i.EQ.bdStop) then
					fnd=tpIndex(workBlockStr,uchGet(bdSet(bdPos)%endCh),rev=true)
					workBlockStr(fnd:)=''
					if (len_trim(workBlockStr).EQ.0) exit
				endif
				write (id ,'(A)') trim(workBlockStr)
			enddo

			rewind(id)
		else
			do i = bdStart,bdStop

				workBlockStr=tpFill(workBlockStr)
				workBlockStr=tpPrepareString(tpGetLineByPointer(i),ignStr=true)

				if (len_trim(workBlockStr).NE.0) holdBlockStr= trim(holdBlockStr)//' &&&& '//trim(adjustl(workBlockStr))
				if (bdStart.EQ.bdStop) exit
			enddo
			
			holdBlockStr=  adjustl( trim(holdBlockStr) )
			holdBlockStr=tpReduce ( trim(tpReplace(trim(holdBlockStr),'&&&&',' ')),quo=true )

			holdBlockStr=tpReplace( trim(holdBlockStr),uchGet(bdSet(bdPos)%separatorCh),' ',quo=true )
			if (atTheStart) then
				fnd=tpIndex( trim(holdBlockStr), uchGet( bdSet(bdPos)%startCh ),endp=true )
				workBlockStr=tpFill(workBlockStr)
				workBlockStr='{'//trim(adjustl(holdBlockStr(fnd+1:len_trim(holdBlockStr))))

				holdBlockStr=trim(workBlockStr)
			endif

			if (atTheEnd) then
				fnd=tpIndex( trim(holdBlockStr),uchGet( bdSet(bdPos)%endCh ),rev=true )
				workBlockStr=tpFill(workBlockStr)
				workBlockStr=trim(holdBlockStr(1:fnd-1))//'}'
				holdBlockStr=trim(workBlockStr)
			endif

			holdBlockStr=tpReplace( trim(holdBlockStr),' ',uchGet(bdSet(bdPos)%separatorCh),quo=true )
			holdBlockStr=tpReduce( trim(holdBlockStr),uchGet(bdSet(bdPos)%separatorCh),quo=true )

		endif
	endif
	bdSet(bdPos)%bdStr=uchSet( trim(holdBlockStr) )

	rcode=bdSetValues(bdPos)

	call tpDisable()

	return
	end function bdParseFile

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdSetValues(bdPos) result(rcode)
	implicit none

	integer*4              :: bdPos,k,j,l,fnd,pStart,pEnd,id,err,optcount,nvars,seplen,seppos
	integer*4              :: bdsta,bdsto
	logical*1              :: isBlank
	type(uch)              :: ustr
	integer*4, allocatable :: variables(:,:)


	!write (*,'(A)' ) uchGet( bdSet(bdPos)%name  )
	!write (*,'(A)' ) uchGet( bdSet(bdPos)%bdStr )
	!write (*,'(A)' ) '123456789012345678901234567890123456789012345678901234567890123456789012345678'
	!write (*,'(A/)') '         1         2         3         4         5         6         7        '
	!stop

	rcode=0
	if (bdSet(bdPos)%freeBlock) then
		continue
	else
		optcount=0
		do k = 1,bdSet(bdPos)%arrayLen
			j=bdSet(bdPos)%varSet(k)

			if (variableSet(j)%opt) optcount=optcount+1
		enddo; isBlank=optcount.EQ.bdSet(bdPos)%arrayLen

		bdsta=tpIndex( uchGet(bdSet(bdPos)%bdStr) , uchGet(bdSet(bdPos)%startCh) ,endp=true)+1
		bdsto=tpIndex( uchGet(bdSet(bdPos)%bdStr) , uchGet(bdSet(bdPos)%endCh)   ,endp=true)-1

		if ( (tpCount( uchGet(bdSet(bdPos)%bdStr), uchGet(bdSet(bdPos)%accordCh) ).EQ.0) .AND. (.NOT.isBlank)) then
			!write (*,*) uchGet(bdSet(bdPos)%bdStr)
			write (errunit,'(A)') 'bd: '//uchGet(bdSet(bdPos)%name)//':'
			write (errunit,'(A)') 'Format error. i.e. incorrect (accordance,comment) symbols: '//uchGet(bdSet(bdPos)%accordCh)
			rcode=-1
			return
		endif

		nvars=tpCount( uchGet(bdSet(bdPos)%bdStr), uchGet(bdSet(bdPos)%separatorCh) )+1
		allocate (variables(nvars,0:2)); variables=0

		j=2; seplen=len(uchGet(bdSet(bdPos)%separatorCh))
		do k = 1,nvars
			seppos=tpIndex( uchGet(bdSet(bdPos)%bdStr),uchGet(bdSet(bdPos)%separatorCh),cnt=k )
			if (seppos.LE.0) seppos=len(uchGet(bdSet(bdPos)%bdStr))
			variables(k,0)=1
			variables(k,1)=j
			variables(k,2)=seppos-1
			j=seppos+seplen
		enddo
		!do k = 1,nvars
		!	write (*,*) variables(k,1),variables(k,2)
		!enddo

		do k = 1,bdSet(bdPos)%arrayLen
			j=bdSet(bdPos)%varSet(k)

			if (j.LT.-1) then; write (errunit,*) 'Internal error.'; stop; endif

			if (variableSet(j)%opt) void=bdSetDefaultValue(variableSet(j)%address)

			if ( (bdsta.GT.bdsto).AND.(isBlank) ) cycle

			fnd=tpIndex( uchGet(bdSet(bdPos)%bdStr),uchGet( variableSet(j)%name ),endp=true )

			if ( (fnd.EQ.0).AND.variableSet(j)%opt) then
				cycle
			elseif ((fnd.EQ.0).AND. .NOT.variableSet(j)%opt) then
				write (errunit,'(A)') 'bd: '//uchGet(bdSet(bdPos)%name)//':'
				write (errunit,'(A)') 'Non-optional argument '//uchGet( variableSet(j)%name )//' is missing.'
				rcode=-2
				cycle
			endif

			do l = 1,nvars
				!write (*,*) uchGet(variableSet(j)%name),variables(l,1),fnd,variables(l,2)
				if ( (fnd.GE.variables(l,1)).AND.(fnd.LE.variables(l,2)) ) then
					variables(l,0)=0
					exit
				endif
			enddo

			if ( tpQuoted(uchGet(bdSet(bdPos)%bdStr),fnd) ) cycle

			if (shortKindLabels(variableSet(j)%kind).EQ.'lo') then
				pStart=tpIndex( uchGet(bdSet(bdPos)%bdStr), uchGet(bdSet(bdPos)%accordCh)   , start=fnd)
				pEnd  =tpIndex( uchGet(bdSet(bdPos)%bdStr), uchGet(bdSet(bdPos)%separatorCh), start=fnd)
				if ((pStart.EQ.0) .OR. (pStart.GT.pEnd)) then
					void=modify( variableSet(j)%address, true )
					cycle
				endif
			endif

			pStart=tpIndex( uchGet(bdSet(bdPos)%bdStr), uchGet(bdSet(bdPos)%accordCh)   , start=fnd)+1
			pEnd  =tpIndex( uchGet(bdSet(bdPos)%bdStr), uchGet(bdSet(bdPos)%separatorCh), start=fnd)-1

			if (pEnd.GT.0) then
				do while ( tpQuoted(uchGet(bdSet(bdPos)%bdStr),pEnd) )
					pEnd=tpIndex( uchGet(bdSet(bdPos)%bdStr), uchGet(bdSet(bdPos)%separatorCh), start=pEnd+1)
				enddo
			else
				pEnd=len( uchGet(bdSet(bdPos)%bdStr) )-1
			endif
			!write (*,*) uchGet( variableSet(j)%name ),pStart,pEnd

			ustr=uchSet( uchGet(bdSet(bdPos)%bdStr,pStart,pEnd), ustr )
			!write (*,*) uchGet( variableSet(j)%name )//' ===> '//uchGet(bdSet(bdPos)%bdStr,pStart,pEnd)//' ===> '//uchGet(variableSet(j)%expect)

			if ( uchGet(variableSet(j)%expect).NE.'any' ) then
				err=checkOnExpect(j,trim(tpDeQuote(uchGet(ustr))))
				if (err.LT.0) then
					write (errunit,'(A\)') 'bd: '//uchGet(bdSet(bdPos)%name)//':'
					write (errunit,'(A)' ) ' Unexpected value for variable: '//uchGet( variableSet(j)%name )//'='&
					                      //trim(tpDeQuote(uchGet(ustr)))//'. Expect '//uchGet(variableSet(j)%expect)
					rcode=-4
!					return
				endif
			endif

			err=0

			if (rcode.EQ.0) then
				select case ( shortKindLabels(variableSet(j)%kind) )
					case('re'); void=modify( variableSet(j)%address,      tpRealByStr( uchGet(ustr),stat=err )   )
					case('in'); void=modify( variableSet(j)%address,       tpIntByStr( uchGet(ustr),stat=err )   )
					case('lo'); void=modify( variableSet(j)%address,       tpLogByStr( uchGet(ustr),stat=err )   )
					case('ch'); void=modify( variableSet(j)%address, trim(  tpDeQuote( uchGet(ustr),stat=err ) ) )
					case('uc'); void=modify( variableSet(j)%address, trim(  tpDeQuote( uchGet(ustr),stat=err ) ) )
				end select
			endif

			if (err.GE.0) then
				if ( tpIsStrInList(shortKindLabels(variableSet(j)%kind),['ch','uc']) ) then
					if (variableSet(j)%several) then
						!write (*,*) uchGet(variableSet(j)%name)
						variableSet(j)%nvals=tpCount(uchGet(ustr),uchGet(variableSet(j)%concatCh))+1
						!write (*,*) variableSet(j)%nvals
					endif
				endif
			endif

			if (err.LT.0) then
				write (errunit,'(A\)') 'bd: '//uchGet(bdSet(bdPos)%name)//':'
				write (errunit,'(A)') ' Format error. i.e. incorrect (accordance,separator,comment) symbols: '//&
				                     uchGet( variableSet(j)%name )
				rcode=-3
			endif
		enddo

		if ( (bdsta.GT.bdsto).AND.(isBlank) ) then
			deallocate (variables)
			return
		endif

		do k = 1,nvars
			if (variables(k,0).NE.0) then
				if (variables(k,1).GT.variables(k,2)) cycle
				write (errunit,'(A\)') 'bd: '//uchGet(bdSet(bdPos)%name)//':'
				write (errunit,'(A)') ' Unexpected variable: '//uchGet( bdSet(bdPos)%bdStr,variables(k,1),variables(k,2) )
				rcode=-4
			endif
		enddo
		deallocate (variables)
	endif

	return
	end function bdSetValues

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function checkOnExpect(position,str) result(rcode)
	implicit none

	character (len=*) :: str
	real*8, parameter :: tollerance=1D-10
	real*8            :: realbounds(2),realvalue

	integer*4         :: position,k,n,sta,sto,ptype,cnt
	integer*4         :: intbounds(2),intvalue

	character (len=variableSet(position)%concatCh%ln) :: conCh
	character (len=variableSet(position)%expect%ln)   :: expect

	real*8   , allocatable :: realListValues(:)
	integer*4, allocatable :: intListValues(:)


	rcode=-1
	conCh =uchGet(variableSet(position)%concatCh)
	expect=uchGet(variableSet(position)%expect)

	ptype=checkExpectFormat( expect )

	!write (*,*) '#'//expect//'# <=== #'//str//'#'

	if (ptype.EQ.-1) then;           return; endif !error
	if (ptype.EQ. 1) then; rcode= 0; return; endif !any

	sta=tpIndex(expect,'(')+1; sto=tpIndex(expect,')')-1
	select case ( shortKindLabels(variableSet(position)%kind) )

		case('re')
			select case (ptype)
				case (2) !range
					realbounds=tpRealArrayByStr( expect(sta:sto),2 )
					realvalue =tpRealByStr(str)

					if ( (realvalue.GE.realbounds(1)).AND.(realvalue.LE.realbounds(2)) ) rcode=0

				case (3) !list
					n=tpCount(expect,',')+1; allocate ( realListValues(n) )
					realListValues=tpRealArrayByStr(expect(sta:sto),n)

					realvalue=tpRealByStr(str)
					do k = 1,n
						if ( abs(realListValues(k)-realvalue).LT.tollerance ) then; rcode=0; exit; endif
					enddo
					deallocate (realListValues)

			end select

		case('in')
			select case (ptype)
				case (2) !range
					intbounds=tpIntArrayByStr( expect(sta:sto),2 )
					intvalue =tpIntByStr(str)

					!write (*,*) uchGet(variableSet(position)%name),intbounds

					if ( (intvalue.GE.intbounds(1)).AND.(intvalue.LE.intbounds(2)) ) rcode=0

				case (3) !list
					n=tpCount(expect,',')+1; allocate ( intListValues(n) )
					intListValues=tpIntArrayByStr(expect(sta:sto),n)

					if (tpIndex( expect,str ).GT.0) rcode=0
					deallocate ( intListValues )

			end select

		case('ch','uc')
			if (ptype.EQ.3) then !list
				if (variableSet(position)%several) then
					voidl=tpSplit(str,conCh)

					cnt=0
					do k = 1,tpSplitLen
						if (tpIndex(expect,tpRetSplit(str,k)).GT.0) cnt=cnt+1
					enddo; if (cnt.EQ.tpSplitLen) rcode=0
				else
					if (tpIndex( expect,str ).GT.0) rcode=0
				endif
			endif

		case ('lo')
			rcode=0

	end select

	return
	end function checkOnExpect

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function checkExpectFormat(str) result(rcode)
	implicit none

	character (len=*) :: str


	rcode=-1
	if (tpStartsWith(str,'any')                           ) rcode=1
	if (tpStartsWith(str,'range(').AND.tpEndsWith(str,')')) rcode=2
	if (tpStartsWith(str,'list(' ).AND.tpEndsWith(str,')')) rcode=3

	return
	end function checkExpectFormat

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdPrintBlockData(name,trgt,width,namePattern,nameSelector,ignoreFreeBlock) result(rcode)
	implicit none

	character (len=*)           :: name
	integer*4                   :: trgt
	integer*4        , optional :: width
	character (len=*), optional :: namePattern,nameSelector
	logical*1        , optional :: ignoreFreeBlock

	integer*4, parameter        :: stdlen=64,reslen=5000
	real*8   , parameter        :: reUbound=1D+6,reLBound=1D-3

	character (len=stdlen)      :: dnamePattern
	character (len=1)           :: dnameSelector
	integer*4                   :: dwidth
	logical*1                   :: dignoreFreeBlock

	character (len=reslen)      :: resStr
	character (len=256)         :: valueString
	character (len=32)          :: mfmt
	character (len=3)           :: reType

	real*8                      :: uReal
	integer*4                   :: uInt
	logical*4                   :: uLog

	integer*4                   :: bdPos,accuracy,ln,k,j,err,pos


	bdPos=find(bdSet(1:bdAppend)%name,trim(name))

	if (bdPos.LT.0) then
		write (errunit,'(A)') trim(name)//': Block data with such name has not been declared.'
		rcode=-1
		return
	endif

	dwidth=80; if (present(width)) dwidth=width
	dnamePattern=tpFill(dnamePattern); dnamePattern='%name option'
	dignoreFreeBlock=false; if (present(ignoreFreeBlock)) dignoreFreeBlock=ignoreFreeBlock

	if (present(namePattern)) dnamePattern=namePattern

	dnamePattern=tpReplace(dnamePattern,'%name',trim(name))
	!dnamePattern=tpCapitalize(dnamePattern)

	dnameSelector='='; if (present(nameSelector)) dnameSelector=nameSelector

	if ( (.NOT.dignoreFreeBlock) .AND.(bdSet(bdPos)%freeBlock)) then
		write (trgt,'(/A)') tpFill(dwidth,dnameSelector)
		call prEchoFile(bdSet(bdPos)%associatedIO,trgt,'numeration')
		write (trgt,'( A)') tpFill(dwidth,dnameSelector)
		return
	endif

	if (bdSet(bdPos)%freeBlock) return

	ln=len_trim(dnamePattern); if (dwidth.LT.ln) dwidth=ln

	write (trgt,'(/A)') tpAdjustc(tpFill(ln,dnameSelector),dwidth)
	write (trgt,'(A)' ) tpAdjustc(trim(dnamePattern)      ,dwidth)
	write (trgt,'(A)' ) tpAdjustc(tpFill(ln,dnameSelector),dwidth)

	resStr=tpFill(reslen)

	do k = 1,bdSet(bdPos)%arrayLen
		j=bdSet(bdPos)%varSet(k)

		! format presets
		mfmt=tpFill(mfmt)
		select case (variableSet(j)%kind)
			case ( 1); uReal=variableSet(j)%pntr8
			case ( 2); uReal=variableSet(j)%pntr4
			case ( 3); uInt =variableSet(j)%pnti8
			case ( 4); uInt =variableSet(j)%pnti4
			case ( 5); uInt =variableSet(j)%pnti2
			case ( 6); uInt =variableSet(j)%pnti1
			case ( 7); uLog =variableSet(j)%pntl8
			case ( 8); uLog =variableSet(j)%pntl4
			case ( 9); uLog =variableSet(j)%pntl2
			case (10); uLog =variableSet(j)%pntl1
		end select

		select case (shortKindLabels(variableSet(j)%kind))
			case ('re')
				if (variableSet(j)%kind.EQ.1) then
					accuracy=15
				else
					accuracy=7
				endif

				if ((mid(abs(uReal)).GT.6).OR.(abs(uReal).LT.1D-5)) then
					mfmt='0.'//tpFill(accuracy,'0')//'E00'
				else
					mfmt='^.'//tpFill(accuracy-mid(uReal),'0')
				endif

				if (abs(uReal)-abs(int(uReal)).LT.10.**(-accuracy)) then
					mfmt=tpFill(mfmt); mfmt='^.'
				endif
				
			case ('in'); mfmt='^'
		end select

		valueString=tpFill(valueString)

		select case(shortKindLabels(variableSet(j)%kind))
			case ('re'); valueString=trim(adjustl(prStrByVal(uReal,trim(mfmt))))
			case ('in'); valueString=trim(adjustl(prStrByVal(uInt ,trim(mfmt))))
			case ('lo'); valueString=trim(adjustl(prStrByVal(uLog)))
			case ('ch'); valueString=trim(variableSet(j)%pntch(1:variableSet(j)%vlen))
			case ('uc'); valueString=uchGet(variableSet(j)%pntuch)
		end select

		if (shortKindLabels(variableSet(j)%kind).EQ.'re') then
			pos=tpIndex(valueString,'E')
			if (pos.EQ.0) then
				do while (tpEndsWith(trim(valueString),'0'))
					ln=len_trim(valueString)
					valueString(ln:ln)=' '
				enddo
			else
				do while ( tpEndsWith(valueString(:pos-1),'0') )
					ln=len_trim(valueString)
					valueString=valueString(:pos-2)//valueString(pos:ln)//tpFill(len(valueString)-ln+1)
					pos=pos-1
				enddo
			endif
		endif

		resStr=trim(resStr)//trim( uchGet(variableSet(j)%name) )//'='//trim(valueString)
		if (k.NE.bdSet(bdPos)%arrayLen) resStr=trim(resStr)//','
	enddo

	call prTable(trim(resStr),trgt,adjust='left',width=dwidth,indent=2,spacer=' ',interColWidth=2,saveOrder=false,equalWidth=false,err=err)

	return
	end function bdPrintBlockData

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_r8(var,name,opt,def,vary,expect) result(rcode)
	implicit none
	logical*1, optional         :: opt,vary
	character (len=*)           :: name
	logical*1                   :: dopt,dvary
	character (len=*), optional :: expect

	integer*4, parameter  :: vkind=8
	real(vkind), target   :: var
	real(vkind), optional :: def
	real(vkind)           :: ddef


	ddef =real(0,vkind); if (present(def))  ddef =def
	dopt =false;         if (present(opt))  dopt =opt
	dvary=false;         if (present(vary)) dvary=vary

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=dvary
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	variableSet(varAppend)%kind  =1
	variableSet(varAppend)%pntr8 =>var
	variableSet(varAppend)%varsta=transfer(real(0,vkind),storeTemplate)
	variableSet(varAppend)%varsto=transfer(real(0,vkind),storeTemplate)
	variableSet(varAppend)%varste=transfer(real(0,vkind),storeTemplate)

	return
	end function shareVariable_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_r4(var,name,opt,def,vary,expect) result(rcode)
	implicit none
	logical*1, optional         :: opt,vary
	character (len=*)           :: name
	logical*1                   :: dopt,dvary
	character (len=*), optional :: expect

	integer*4, parameter  :: vkind=4
	real(vkind), target   :: var
	real(vkind), optional :: def
	real(vkind)           :: ddef


	ddef =real(0,vkind); if (present(def))  ddef =def
	dopt =false;         if (present(opt))  dopt =opt
	dvary=false;         if (present(vary)) dvary=vary

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=dvary
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	variableSet(varAppend)%kind  =2
	variableSet(varAppend)%pntr4 =>var
	variableSet(varAppend)%varsta=transfer(real(0,vkind),storeTemplate)
	variableSet(varAppend)%varsto=transfer(real(0,vkind),storeTemplate)
	variableSet(varAppend)%varste=transfer(real(0,vkind),storeTemplate)

	return
	end function shareVariable_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_i8(var,name,opt,def,vary,expect) result(rcode)
	implicit none
	logical*1, optional   :: opt,vary
	character (len=*)     :: name
	logical*1             :: dopt,dvary
	character (len=*), optional :: expect

	integer*4, parameter     :: vkind=8
	integer(vkind), target   :: var
	integer(vkind), optional :: def
	integer(vkind)           :: ddef


	ddef =int(0,vkind); if (present(def))  ddef =def
	dopt =false;        if (present(opt))  dopt =opt
	dvary=false;        if (present(vary)) dvary=vary

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=dvary
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	variableSet(varAppend)%kind  =3
	variableSet(varAppend)%pnti8 =>var
	variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

	return
	end function shareVariable_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_i4(var,name,opt,def,vary,expect) result(rcode)
	implicit none
	logical*1, optional   :: opt,vary
	character (len=*)     :: name
	logical*1             :: dopt,dvary
	character (len=*), optional :: expect

	integer*4, parameter     :: vkind=4
	integer(vkind), target   :: var
	integer(vkind), optional :: def
	integer(vkind)           :: ddef


	ddef =int(0,vkind); if (present(def))  ddef =def
	dopt =false;        if (present(opt))  dopt =opt
	dvary=false;        if (present(vary)) dvary=vary

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=dvary
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	variableSet(varAppend)%kind  =4
	variableSet(varAppend)%pnti4 =>var
	variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

	return
	end function shareVariable_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_i2(var,name,opt,def,vary,expect) result(rcode)
	implicit none
	logical*1, optional   :: opt,vary
	character (len=*)     :: name
	logical*1             :: dopt,dvary
	character (len=*), optional :: expect

	integer*4, parameter     :: vkind=2
	integer(vkind), target   :: var
	integer(vkind), optional :: def
	integer(vkind)           :: ddef


	ddef =int(0,vkind); if (present(def))  ddef =def
	dopt =false;        if (present(opt))  dopt =opt
	dvary=false;        if (present(vary)) dvary=vary

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=dvary
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	variableSet(varAppend)%kind  =5
	variableSet(varAppend)%pnti2 =>var
	variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

	return
	end function shareVariable_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_i1(var,name,opt,def,vary,expect) result(rcode)
	implicit none
	logical*1, optional   :: opt,vary
	character (len=*)     :: name
	logical*1             :: dopt,dvary
	character (len=*), optional :: expect

	integer*4, parameter     :: vkind=1
	integer(vkind), target   :: var
	integer(vkind), optional :: def
	integer(vkind)           :: ddef


	ddef =int(0,vkind); if (present(def))  ddef =def
	dopt =false;        if (present(opt))  dopt =opt
	dvary=false;        if (present(vary)) dvary=vary

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=dvary
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	variableSet(varAppend)%kind  =6
	variableSet(varAppend)%pnti1 =>var
	variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
	variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

	return
	end function shareVariable_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_l8(var,name,opt,def) result(rcode)
	implicit none
	logical*1, optional   :: opt
	character (len=*)     :: name
	logical*1             :: dopt

	integer*4, parameter     :: vkind=8
	logical(vkind), target   :: var
	logical(vkind), optional :: def
	logical(vkind)           :: ddef


	ddef =logical(false,vkind); if (present(def))  ddef =def
	dopt =false;                if (present(opt))  dopt =opt

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=false
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%several =false
	variableSet(varAppend)%expect  =uchSet('any')
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	variableSet(varAppend)%kind =7
	variableSet(varAppend)%pntl8=>var

	return
	end function shareVariable_l8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_l4(var,name,opt,def) result(rcode)
	implicit none
	logical*1, optional   :: opt
	character (len=*)     :: name
	logical*1             :: dopt

	integer*4, parameter     :: vkind=4
	logical(vkind), target   :: var
	logical(vkind), optional :: def
	logical(vkind)           :: ddef


	ddef =logical(false,vkind); if (present(def))  ddef =def
	dopt =false;                if (present(opt))  dopt =opt

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=false
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%several =false
	variableSet(varAppend)%expect  =uchSet('any')
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	variableSet(varAppend)%kind =8
	variableSet(varAppend)%pntl4=>var

	return
	end function shareVariable_l4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_l2(var,name,opt,def) result(rcode)
	implicit none
	logical*1, optional   :: opt
	character (len=*)     :: name
	logical*1             :: dopt

	integer*4, parameter     :: vkind=2
	logical(vkind), target   :: var
	logical(vkind), optional :: def
	logical(vkind)           :: ddef


	ddef =logical(false,vkind); if (present(def))  ddef =def
	dopt =false;                if (present(opt))  dopt =opt

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=false
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%several =false
	variableSet(varAppend)%expect  =uchSet('any')
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0

	variableSet(varAppend)%kind =9
	variableSet(varAppend)%pntl2=>var

	return
	end function shareVariable_l2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_l1(var,name,opt,def) result(rcode)
	implicit none
	logical*1, optional   :: opt
	character (len=*)     :: name
	logical*1             :: dopt

	integer*4, parameter     :: vkind=1
	logical(vkind), target   :: var
	logical(vkind), optional :: def
	logical(vkind)           :: ddef


	ddef =logical(false,vkind); if (present(def))  ddef =def
	dopt =false;                if (present(opt))  dopt =opt

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%vlen    =vkind
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=false
	variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
	variableSet(varAppend)%several =false
	variableSet(varAppend)%expect  =uchSet('any')
	variableSet(varAppend)%concatCh=uchSet('')
	variableSet(varAppend)%nvals   =0


	variableSet(varAppend)%kind =10
	variableSet(varAppend)%pntl1=>var

	return
	end function shareVariable_l1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_c(var,name,opt,def,expect,several,concat) result(rcode)
	implicit none
	character (len=*), target   :: var
	character (len=*), optional :: def
	logical*1, optional         :: opt,several
	character (len=*), optional :: expect,concat

	character (len=*)           :: name
	logical*1                   :: dopt,dseveral


	dopt =false;    if (present(opt)) dopt =opt
	dseveral=false; if (present(several)) dseveral=several

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%kind    =11
	variableSet(varAppend)%vlen    =len(var)
	variableSet(varAppend)%pntch   =>var
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=false
	variableSet(varAppend)%several =dseveral
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	if (present(concat)) then
		variableSet(varAppend)%concatCh=uchSet(concat)
	else
		variableSet(varAppend)%concatCh=uchSet('+')
	endif

	if (present(def)) variableSet(varAppend)%defc=uchSet(trim(def),variableSet(varAppend)%defc)

	return
	end function shareVariable_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*8 function shareVariable_uc(var,name,opt,def,expect,several,concat) result(rcode)
	implicit none
	type(uch), target           :: var
	character (len=*), optional :: def
	logical*1, optional         :: opt,several
	character (len=*), optional :: expect,concat

	character (len=*)           :: name
	logical*1                   :: dopt,dseveral


	dopt =false;    if (present(opt)) dopt =opt
	dseveral=false; if (present(several)) dseveral=several

	rcode=loc(var)
	if (find(variableSet(1:varAppend)%address,int(loc(var),8)).GT.0) then
		write (errunit,*) 'Variable '//name//' at address',loc(var),' is already shared.'
		rcode=-1; return
	endif

	if ( dopt .AND. .NOT.present(def) ) then
		write (errunit,*) name//': If variable is optional, call must contain default value.'
		rcode=-2; return
	endif

	if ( .NOT.dopt .AND. present(def) ) then
		write (errunit,*) name//': Call must contain default value only in case variable is optional.'
		rcode=-3; return
	endif

	if (present(expect)) then
		if (checkExpectFormat(expect).LT.0) then
			write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
			rcode=-4; return
		endif
	endif

	varAppend=varAppend+1

	variableSet(varAppend)%name    =uchSet( trim(adjustl(name)),variableSet(varAppend)%name )
	variableSet(varAppend)%address =loc(var)
	variableSet(varAppend)%position=varAppend
	variableSet(varAppend)%kind    =12
	variableSet(varAppend)%vlen    =0 !len(var)
	variableSet(varAppend)%pntuch  =>var
	variableSet(varAppend)%opt     =dopt
	variableSet(varAppend)%variable=false
	variableSet(varAppend)%several =dseveral
	variableSet(varAppend)%nvals   =0

	if (present(expect)) then
		variableSet(varAppend)%expect=uchSet(expect)
	else
		variableSet(varAppend)%expect=uchSet('any')
	endif

	if (present(concat)) then
		variableSet(varAppend)%concatCh=uchSet(concat)
	else
		variableSet(varAppend)%concatCh=uchSet('+')
	endif

	if (present(def)) variableSet(varAppend)%defc=uchSet(trim(def),variableSet(varAppend)%defc)

	return
	end function shareVariable_uc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function bdSetDefaultValue(address) result(ret)
	implicit none

	integer*8 :: address
	integer*4 :: i


	i=find(variableSet(:)%address,address)
	if (i.GT.0) then
		if (variableSet(i)%opt) then
			select case(variableSet(i)%kind)
				case ( 1); variableSet(i)%pntr8=transfer(variableSet(i)%defstor,variableSet(i)%pntr8)
				case ( 2); variableSet(i)%pntr4=transfer(variableSet(i)%defstor,variableSet(i)%pntr4)
				case ( 3); variableSet(i)%pnti8=transfer(variableSet(i)%defstor,variableSet(i)%pnti8)
				case ( 4); variableSet(i)%pnti4=transfer(variableSet(i)%defstor,variableSet(i)%pnti4)
				case ( 5); variableSet(i)%pnti2=transfer(variableSet(i)%defstor,variableSet(i)%pnti2)
				case ( 6); variableSet(i)%pnti1=transfer(variableSet(i)%defstor,variableSet(i)%pnti1)
				case ( 7); variableSet(i)%pntl8=transfer(variableSet(i)%defstor,variableSet(i)%pntl8)
				case ( 8); variableSet(i)%pntl4=transfer(variableSet(i)%defstor,variableSet(i)%pntl4)
				case ( 9); variableSet(i)%pntl2=transfer(variableSet(i)%defstor,variableSet(i)%pntl2)
				case (10); variableSet(i)%pntl1=transfer(variableSet(i)%defstor,variableSet(i)%pntl1)

				! Please, let me 'shoot myself in the foot'.
				! Will raise exception 'array bounds exceeded', but works fine without boundary check.
				! Source file should be compiled without "-CB" (for ifort) and "-fbounds-check" for GNU fortran
				case (11)
					variableSet(i)%pntch(1:variableSet(i)%vlen)=tpFill(variableSet(i)%vlen)
					variableSet(i)%pntch(1:variableSet(i)%defc%ln)=uchGet(variableSet(i)%defc)

				case (12)
					variableSet(i)%pntuch=uchSet( uchGet(variableSet(i)%defc) )
					
			end select
		endif
	endif

	ret=0; return
	end function bdSetDefaultValue

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function modify_c(address,val) result(ret)
	implicit none

	character (len=*), intent(in) :: val
	integer*8                     :: address
	integer*4                     :: i


	i=find(variableSet(:)%address,address)
	if (i.GT.0) then

		if (variableSet(i)%kind.EQ.12) then
			variableSet(i)%pntuch=uchSet( val )
			ret=0
			return
		endif

		if (variableSet(i)%kind.NE.11) then
			ret=-1
			return
		endif

		! Please, let me 'shoot myself in the foot'.
		! Will raise exception 'array bounds exceeded', but works fine without boundary check.
		! Source file should be compiled without "-CB" (for ifort) and "-fbounds-check" (for GNU fortran) keys.
		variableSet(i)%pntch(1:variableSet(i)%vlen)=tpFill(variableSet(i)%vlen)
		variableSet(i)%pntch(1:variableSet(i)%vlen)=val
	else
		ret=-2
		return
	endif

	ret=0; return
	end function modify_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function modify_r(address,val) result(ret)
	implicit none

	real*8, intent(in) :: val
	integer*8          :: address
	integer*4          :: i


	i=find(variableSet(:)%address,address)
	if (i.GT.0) then
		if ((variableSet(i)%kind.NE.1).AND.(variableSet(i)%kind.NE.2)) then
			ret=-1
			return
		endif
		select case (variableSet(i)%vlen)
			case (8); variableSet(i)%pntr8=real(val,kind=8)
			case (4); variableSet(i)%pntr4=real(val,kind=4)
		end select
	else; ret=-2; return
	endif

	ret=0; return
	end function modify_r

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function modify_i(address,val) result(ret)
	implicit none

	integer*8, intent(in) :: val
	integer*8             :: address
	integer*4             :: i


	i=find(variableSet(:)%address,address)
	if (i.GT.0) then
		if ((variableSet(i)%kind.NE.3).AND.(variableSet(i)%kind.NE.4).AND.&
		    (variableSet(i)%kind.NE.5).AND.(variableSet(i)%kind.NE.6)) then
			ret=-1
			return
		endif
		select case (variableSet(i)%vlen)
			case (8); variableSet(i)%pnti8=int(val,kind=8)
			case (4); variableSet(i)%pnti4=int(val,kind=4)
			case (2); variableSet(i)%pnti2=int(val,kind=2)
			case (1); variableSet(i)%pnti1=int(val,kind=1)
		end select
	else; ret=-2; return
	endif

	ret=0; return
	end function modify_i

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer*4 function modify_l(address,val) result(ret)
	implicit none

	logical*1, intent(in) :: val
	integer*8             :: address
	integer*4             :: i


	i=find(variableSet(:)%address,address)
	if (i.GT.0) then
		if ((variableSet(i)%kind.NE.6).AND.(variableSet(i)%kind.NE.7).AND.&
			(variableSet(i)%kind.NE.8).AND.(variableSet(i)%kind.NE.9)) then
			ret=-1
			return
		endif
		select case (variableSet(i)%vlen)
			case (8); variableSet(i)%pntl8=logical(val,kind=8)
			case (4); variableSet(i)%pntl4=logical(val,kind=4)
			case (2); variableSet(i)%pntl2=logical(val,kind=2)
			case (1); variableSet(i)%pntl1=val
		end select
	else; ret=-2; return	
	endif

	ret=0; return
	end function modify_l

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine bdSetIOunit(iunt)
	implicit none

	integer*4, intent(in) :: iunt
	integer*4             :: unt


	unt=iunt
	if (unt.EQ.5) unt=6 !stdin  (forbidden)
	if (unt.EQ.0) unt=6 !stdout (for compatibility)
	iounit=unt

	return
	end subroutine bdSetIOunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine bdSetERRunit(iunt)
	implicit none

	integer*4, intent(in) :: iunt
	integer*4             :: unt


	unt=iunt
	if (unt.EQ.5) unt=6 !stdin  (forbidden)
	if (unt.EQ.0) unt=6 !stdout (for compatibility)
	errunit=unt

	return
	end subroutine bdSetERRunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine bdFinalize
	implicit none


	!continue

	return
	end subroutine bdFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module datablock