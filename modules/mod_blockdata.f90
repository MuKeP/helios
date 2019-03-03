    module datablock

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,      only: assignment (=)
    use glob,      only: true,false,void,voidl,NaN,mid,uch,find,collectArray
    use fcontrol,  only: fcNewID,fcNullID,fcBanID
    use printmod,  only: prLongText,prStrByVal,prTable,prEchoFile

    ! almost everything :(
    use txtParser, only: operator(.in.)
    use txtParser, only: tpIgnoSpac,tpIgnoTabs,tpCaseSens,tpReplTabs,tpReduSpac,tpIgnoComm,tpAlloQuot
    use txtParser, only: tpSetCommentMark,tpPrepareString
    use txtParser, only: tpJoin,tpEnable,tpDisable,tpLocate,tpPointer,tpGetLineByPointer
    use txtParser, only: tpSplit,tpSplitLen,tpSplitHold,tpRetSplit
    use txtParser, only: tpAdjustl,tpAdjustc,tpFill,tpUpperCase
    use txtParser, only: tpIndex,tpCount,tpReplace,tpReduce,tpDeQuote
    use txtParser, only: tpQuoted,tpStartsWith,tpEndsWith
    use txtParser, only: tpRealByStr,tpIntByStr,tpLogByStr,tpRealArrayByStr,tpIntArrayByStr
    use txtParser, only: tpRealNumber,tpLetters,tpSetAccordance

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: bdVersion='1.711'
    character (len=*), parameter :: bdDate   ='2019.01.03'
    character (len=*), parameter :: bdAuthor ='Anton B. Zakharov'

    integer*4, parameter :: typeSetSize=12,mxKindNameLen=9,maxListSetLen=100
    integer*4, parameter :: maxSetLen=2500,expectFileWidth=200
    integer*4, parameter :: prepareStringLen=5000
    integer*4, parameter :: storageLen=32,markUpChLen=64

    character (len=*), parameter  :: defCommentCh='#',defStartCh='%name{',defEndCh='}',&
                                     defAccordCh='=',defSeparatorCh='|'

    character (len=mxKindNameLen) :: kindLabels(0:typeSetSize)
    data kindLabels /'real*16',  'real*8',   'real*4',&
    &                'integer*8','integer*4','integer*2','integer*1',&
    &                'logical*8','logical*4','logical*2','logical*1',&
    &                'character','type(uch)'/

    character (len=2) :: shortKindLabels(0:typeSetSize)
    data shortKindLabels / 're','re','re','in','in','in','in','lo','lo','lo','lo','ch','uc' /

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type expectHolder
        type(uch)         :: kind ! [range,list,any]
        type(uch)         :: value ! [ge,le,in,ins]
        type(uch)         :: listStore ! store all values
        type(uch)         :: listType ! [all,none]
        character (len=1) :: even,positive ![+,-,n]
        real*8            :: ub,lb,step

        contains

        procedure :: init => expectInitiate
    end type expectHolder

    type variable
        type(uch)          :: name
        integer*8          :: address
        integer*4          :: position,kind,vlen,nvals
        logical*1          :: opt,variable,several

        character (len=1)  :: potentiate

        real*16  , pointer :: pntr16
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

        type(expectHolder) :: exph

        character (len=storageLen) :: defstor
        type(uch)                  :: defc,expect,concatCh,description

        character (len=storageLen) :: varsta,varsto,varste
    end type variable

    type bd
        type(uch)              :: name,bdStr,commentCh,startCh,endCh,accordCh,separatorCh
        type(uch)              :: description
        integer*4              :: arrayLen,associatedIO
        integer*4, allocatable :: varSet(:)
        logical*1              :: freeBlock
    end type bd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type(variable) :: variableSet(maxSetLen)
    type(bd)       :: bdSet(maxListSetLen)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=storageLen)       :: storeTemplate
    character (len=prepareStringLen) :: holdBlockStr,workBlockStr

    integer*4                        :: varAppend=0,bdAppend=0,iounit=6,errunit=6

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    interface bdShareVariable
        module procedure shareVariable_r16,shareVariable_r8,shareVariable_r4,&
                         shareVariable_i8 ,shareVariable_i4,shareVariable_i2,&
                         shareVariable_i1 ,shareVariable_l8,shareVariable_l4,&
                         shareVariable_l2 ,shareVariable_l1,shareVariable_c ,&
                         shareVariable_uc
    end interface bdShareVariable

    interface modify
        module procedure modify_c,modify_r,modify_i,modify_l
    end interface modify

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    public :: bdVersion,bdDate,bdAuthor

    public :: bdShareVariable,bdCollect,bdAddVariable,bdRemoveVariable,bdParseFile,&
              bdPrintBlockData,bdSetIOunit,bdSetERRunit,bdFinalize

    public :: bdAddDescription,bdVariableAddDescription,bdPrintHelp,bdCheckExternalParameter

    !public :: checkExpectFormat,checkOnExpect,variableSet,varAppend

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdPrintHelp(bdname,unit,width) result(rcode)
    implicit none

    character (len=*), intent(in) :: bdname
    integer*4        , intent(in) :: unit,width
    integer*4                     :: i,j,k,l,m


    rcode=0
    k=find(bdSet(1:bdAppend)%name,trim(bdname))
    if (k.LT.0) then
        write (errunit,'(A)') trim(bdname)//': Block data with such name has not been declared.'
        rcode=-1
        return
    endif

    write (unit,'(A)') tpFill(width,'()')
    write (unit,'(A/A/)') tpUpperCase('Block: '),bdSet(k)%name%get()
    write (unit,'(A)') tpUpperCase('Block description:')
    write (unit,'(A)') bdSet(k)%description%get()
    write (unit,*)
    if (bdSet(k)%freeBlock) then
        write (unit,'(A)') tpFill(width,'()')
        return
    endif

    write (unit,'(A/)') tpUpperCase(tpFill(5,'=')//' Block contains variables '//tpFill(5,'='))
    do i = 1,bdSet(k)%arrayLen
        j=bdSet(k)%varSet(i)

        write (unit,'(A,1X,A,1X,A,1X,A/)') tpFill(5,'*'),tpUpperCase('Variable'),prStrByVal(i),tpFill(5,'*')
        write (unit,'(A/A/)') tpUpperCase('Name:'),variableSet(j)%name%get()
        write (unit,'(A/A/)') tpUpperCase('Variable description:'),variableSet(j)%description%get()

        select case(variableSet(j)%kind)
            case (0)    ; write (unit,100) 'float128'
            case (1)    ; write (unit,100) 'float64'
            case (2)    ; write (unit,100) 'float32'
            case (3)    ; write (unit,100) 'int64'
            case (4)    ; write (unit,100) 'int32'
            case (5)    ; write (unit,100) 'int16'
            case (6)    ; write (unit,100) 'int8'
            case ( 7:10); write (unit,100) 'boolean'
            case (11:12); write (unit,100) 'string'
        end select

100     format ('TYPE:'/A/)

        if (variableSet(j)%several) then
            write (unit,'(A/A/)') tpUpperCase('Multivalued:'),'true (cancatination symbol: '//variableSet(j)%concatCh%get()//')'
        endif

        write (unit,'(A/A/)') tpUpperCase('Optional:'),prStrByVal(variableSet(j)%opt)

        if (variableSet(j)%potentiate .in. '+-') then
            write (unit,'(A/A/)') tpUpperCase('Exponent:'),'true (sign: 10^'//variableSet(j)%potentiate//'var)'
        endif

        if (variableSet(j)%opt) then
            write (unit,'(A)') tpUpperCase('Default value:')
            select case(variableSet(j)%kind)
                case ( 0); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pntr16),3,7,'exp')
                case ( 1); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pntr8 ),3,7,'exp')
                case ( 2); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pntr4 ),3,7,'exp')
                case ( 3); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pnti8))
                case ( 4); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pnti4))
                case ( 5); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pnti2))
                case ( 6); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pnti1))
                case ( 7); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pntl8))
                case ( 8); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pntl4))
                case ( 9); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pntl2))
                case (10); write (unit,101) prStrByVal(transfer(variableSet(j)%defstor,variableSet(j)%pntl1))
                case (11); write (unit,101) prStrByVal(variableSet(j)%defc)
                case (12); write (unit,101) prStrByVal(variableSet(j)%defc)
            end select
101         format (A)
            write (unit,*)
        endif

        if (variableSet(j)%expect%get().NE.'any') then
            write (unit,'(A/A/)') tpUpperCase('Expected values preset:'),variableSet(j)%expect%get()
        endif
    enddo
    write (unit,'(A)') tpFill(width,'()')

    return
    end function bdPrintHelp

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdAddDescription(bdname,description) result(rcode)
    implicit none

    character (len=*) :: bdname,description
    integer*4         :: k


    rcode=0
    k=find(bdSet(1:bdAppend)%name,trim(bdname))
    if (k.LT.0) then
        write (errunit,'(A)') trim(bdname)//': Block data with such name has not been declared.'
        rcode=-1
        return
    endif

    bdSet(k)%description=description

    return
    end function bdAddDescription

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdVariableAddDescription(id,description) result(rcode)
    implicit none

    integer*8         :: id
    character (len=*) :: description
    integer*4         :: k


    k=find(variableSet(1:varAppend)%address,id)
    if (k.LT.0) then
        write (errunit,'(A,1X,i,A)') 'Not shared (unknown) variable at address ',id,'.'
        rcode=-1
        return
    endif

    variableSet(k)%description=description

    rcode=0; return
    end function bdVariableAddDescription

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine bdAttachString(bdname,str,stat)
    implicit none

    character (len=*)  , intent(in)  :: bdname,str
    integer*4, optional, intent(out) :: stat
    integer*4                        :: j


    j=find(bdSet(1:bdAppend)%name,trim(bdname))
    if (j.LT.0) then
        write (errunit,'(A)') bdname//': Block data with such name has not been declared.'
        if (present(stat)) stat=-1
        return
    endif

    bdSet(j)%bdStr=str
    if (present(stat)) stat=0

    return
    end subroutine bdAttachString

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdCollect(name,array,startCh,endCh,commentCh,accordCh,separatorCh,freeBlock,associatedIO) result(rcode)
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

    bdSet(bdAppend)%name=tpAdjustl(name)

    bdSet(bdAppend)%startCh=tpReplace(trim(ustartCh),'%name',bdSet(bdAppend)%name%get())
    bdSet(bdAppend)%commentCh=trim(ucommentCh)
    bdSet(bdAppend)%accordCh=trim(uaccordCh)
    bdSet(bdAppend)%separatorCh=trim(useparatorCh)
    bdSet(bdAppend)%endCh=tpReplace(trim(uendCh),'%name',bdSet(bdAppend)%name%get())
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
    end function bdCollect

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdAddVariable(bdname,address) result(rcode)
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
    end function bdAddVariable

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdRemoveVariable(bdname,address) result(rcode)
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
    end function bdRemoveVariable

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdParseFile(name,filename,expectWidth) result(rcode)
    implicit none

    character (len=*)   :: name,filename
    integer*4, optional :: expectWidth
    integer*4           :: dexpectWidth

    integer*4           :: bdPos,startFnd,stopFnd,bdStart,bdStop,i,fnd,err,id
    logical*1           :: atTheStart,atTheEnd


    rcode=0
    dexpectWidth=expectFileWidth; if (present(expectWidth)) dexpectWidth=expectWidth

    bdPos=find(bdSet(1:bdAppend)%name,trim(name))

    if (bdPos.LT.0) then
        write (errunit,'(A)') trim(name)//': Block data with such name has not been declared.'
        rcode=-1
        return
    endif

    atTheStart=tpIndex( bdSet(bdPos)%startCh%get(), bdSet(bdPos)%name%get() ).NE.0
    atTheEnd  =tpIndex( bdSet(bdPos)%endCh%get()  , bdSet(bdPos)%name%get() ).NE.0

    i=fcNewID(); open (i,file=filename,status='old',iostat=err)
    if (err.NE.0) then
        write (errunit,'(A)') 'bd: '//bdSet(bdPos)%name%get()//':'//" File does not exist: "//trim(filename)

        rcode=-2; void=fcNullID(i)
        return
    else
        close (i); void=fcNullID(i)
    endif

    call tpEnable(dexpectWidth); voidl=tpJoin(filename)
    if (.NOT.voidl) then
        write (errunit,'(A)') 'bd: '//bdSet(bdPos)%name%get()//':'
        write (errunit,'(A)') 'File access denied: '//trim(filename)
        rcode=-3
        return
    endif

    tpCaseSens=false; tpReplTabs=true; tpReduSpac=true;  tpIgnoSpac=false
    tpIgnoTabs=false; tpIgnoComm=true; tpAlloQuot=false

    !write (*,*) 'Start searching'; call cpu_time(sta)

    !write (*,*) 'Looking for #'//bdSet(bdPos)%startCh%get()//'# and #'//bdSet(bdPos)%endCh%get()//'#'

    void=tpSetCommentMark( bdSet(bdPos)%commentCh%get() )
    startFnd=tpLocate( bdSet(bdPos)%startCh%get(),true ); bdStart=tpPointer+1
    stopFnd =tpLocate( bdSet(bdPos)%endCh%get() );        bdStop =tpPointer+1

    !write (*,*) 'End searching'; call cpu_time(sto)
    !write (*,*) 'Spent',sto-sta
    !if (bdSet(bdPos)%freeBlock) stop

    !write (*,*) 'Found ',startFnd,stopFnd
    !write (*,*) 'Found ',bdStart,bdStop

    holdBlockStr=tpFill(holdBlockStr)
    if ( (startFnd.EQ.-1).AND.(startFnd.EQ.-1) ) then
        holdBlockStr=bdSet(bdPos)%startCh%get()//bdSet(bdPos)%endCh%get()
    else
        if (bdSet(bdPos)%freeBlock) then
            id=bdSet(bdPos)%associatedIO; void=fcBanID(id)
            open (id,status='scratch') !; write (id,'(A)') bdSet(bdPos)%bdStr%get(); rewind(id)

            workBlockStr=tpFill(workBlockStr); workBlockStr=tpGetLineByPointer(bdStart)
            fnd=tpIndex( workBlockStr,bdSet(bdPos)%startCh%get(),endp=true )+1
            write (id,'(A)') trim( workBlockStr(fnd:) )

            do i = bdStart+1,bdStop
                workBlockStr=tpFill(workBlockStr); workBlockStr=tpGetLineByPointer(i)

                if (i.EQ.bdStop) then
                    fnd=tpIndex(workBlockStr,bdSet(bdPos)%endCh%get(),rev=true)
                    workBlockStr(fnd:)=''
                    if (len_trim(workBlockStr).EQ.0) exit
                endif
                write (id ,'(A)') trim(workBlockStr)
            enddo

            rewind(id)
        else
            do i = bdStart,bdStop

                workBlockStr=tpFill(workBlockStr)
                !write (*,*) 'Line by line ',trim(tpGetLineByPointer(i))
                workBlockStr=tpPrepareString(tpGetLineByPointer(i),ignStr=true)

                if (len_trim(workBlockStr).NE.0) holdBlockStr= trim(holdBlockStr)//' &&&& '//trim(adjustl(workBlockStr))
                if (bdStart.EQ.bdStop) exit
            enddo

            holdBlockStr=  adjustl( trim(holdBlockStr) )
            holdBlockStr=tpReduce( trim(tpReplace(trim(holdBlockStr),'&&&&',' ')),quo=true )

            holdBlockStr=tpReplace( trim(holdBlockStr),bdSet(bdPos)%separatorCh%get(),' ',quo=true )
            if (atTheStart) then
                fnd=tpIndex( trim(holdBlockStr),bdSet(bdPos)%startCh%get(),endp=true )
                workBlockStr=tpFill(workBlockStr)
                workBlockStr='{'//trim(adjustl(holdBlockStr(fnd+1:len_trim(holdBlockStr))))

                holdBlockStr=trim(workBlockStr)
            endif

            if (atTheEnd) then
                fnd=tpIndex( trim(holdBlockStr),bdSet(bdPos)%endCh%get(),rev=true )
                workBlockStr=tpFill(workBlockStr)
                workBlockStr=trim(holdBlockStr(1:fnd-1))//'}'
                holdBlockStr=trim(workBlockStr)
            endif

            holdBlockStr=tpReplace( trim(holdBlockStr),' ',bdSet(bdPos)%separatorCh%get(),quo=true )
            holdBlockStr=tpReduce( trim(holdBlockStr),bdSet(bdPos)%separatorCh%get(),quo=true )

        endif
    endif

    !write (*,*) trim(holdBlockStr)
    bdSet(bdPos)%bdStr=trim(holdBlockStr)

    ! if (trim(name).EQ.'rdm') then
    !     write (*,*) bdSet(bdPos)%bdStr%get()
    !     stop
    ! endif

    rcode=bdSetValues(bdPos)

    call tpDisable()

    return
    end function bdParseFile

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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

    real*8                      :: uReal
    integer*4                   :: uInt
    logical*4                   :: uLog

    integer*4                   :: bdPos,accuracy,ln,k,j,err,pos


    rcode=0
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

    dnamePattern=tpReplace(trim(dnamePattern),'%name',trim(name))
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
            case ( 0); uReal=variableSet(j)%pntr16
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
                select case (variableSet(j)%kind)
                    case (0); accuracy=31
                    case (1); accuracy=15
                    case (2); accuracy=7
                end select

                if ((abs(uReal).GT.1D+5).OR.(uReal.LT.1D-5)) then
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
            case ('uc'); valueString=variableSet(j)%pntuch%get()
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

        resStr=trim(resStr)//trim( variableSet(j)%name%get() )//'='//trim(valueString)
        if (k.NE.bdSet(bdPos)%arrayLen) resStr=trim(resStr)//'&&&'
    enddo

    !write (*,'(A)') '#'//trim(resStr)//'#'
    call prTable(trim(resStr),trgt,adjust='left',width=dwidth,indent=2,separator='&&&',spacer=' ',interColWidth=2,saveOrder=false,equalWidth=false)

    return
    end function bdPrintBlockData

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdCheckExternalParameter(bdentry, paramentry, value, autoset) result(rcode)
    implicit none

    type(uch), intent(in)           :: bdentry, paramentry, value
    logical*1, intent(in), optional :: autoset
    logical*1                       :: dautoset,found
    integer*4                       :: i,j,k,l,position


    dautoset=false; if (present(autoset)) dautoset=autoset

    rcode=0
    j=find(bdSet(1:bdAppend)%name,bdentry%get())
    if (j.LT.0) then
        write (errunit,'(A)') 'External parameter check.'
        write (errunit,'(A)') trim(bdentry%get())//': Block data with such name has not been declared.'
        rcode=-1
        return
    endif

    do i = 1,bdSet(j)%arrayLen
        if (variableSet(bdSet(j)%varset(i))%name%get().EQ.paramentry%get()) then
            found=true
            position=bdSet(j)%varset(i)
            exit
        endif
    enddo

    if (.NOT.found) then
        write (errunit,'(A)') 'External parameter check.'
        write (errunit,'(A,1X,A)') 'Not found variable in bd: '//bdSet(j)%name%get(),&
                                   ' with name '//paramentry%get()//'.'
        rcode=-2
        return
    endif

    ! write(*,*) 'BD',bdSet(j)%name%get()
    ! write(*,*) 'PA',variableSet(position)%name%get()
    ! write(*,*) 'VA',value%get()

    rcode=checkOnExpect(position,value%get())

    if (dautoset) then
        rcode=bdSetValue(position,value%get())
    endif

    return
    end function bdCheckExternalParameter

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdSetValues(bdPos) result(rcode)
    implicit none

    integer*4              :: bdPos,k,j,l,fnd,pStart,pEnd,err,optcount,nvars,seplen,seppos
    integer*4              :: bdsta,bdsto,quoshift
    logical*1              :: isBlank
    type(uch)              :: ustr
    integer*4, allocatable :: variables(:,:)


    ! if (bdSet(bdPos)%name%get().EQ.'iteration') then
    ! call system('CLS')
    ! write (*,*)
    ! write (*,'(A)' ) bdSet(bdPos)%name%get()
    ! write (*,'(A)' ) bdSet(bdPos)%bdStr%get()
    ! write (*,'(A)' ) '123456789012345678901234567890123456789012345678901234567890123456789012345678'
    ! write (*,'(A/)') '         1         2         3         4         5         6         7        '
    ! endif

    rcode=0
    if (bdSet(bdPos)%freeBlock) then
        continue
    else
        optcount=0
        do k = 1,bdSet(bdPos)%arrayLen
            j=bdSet(bdPos)%varSet(k)

            if (variableSet(j)%opt) optcount=optcount+1
        enddo
        isBlank=optcount.EQ.bdSet(bdPos)%arrayLen

        bdsta=tpIndex( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%startCh%get() ,endp=true)+1
        bdsto=tpIndex( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%endCh%get()   ,endp=true)-1

        if ( (tpCount( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%accordCh%get() ).EQ.0) .AND. (.NOT.isBlank)) then
            write (errunit,'(A\)') 'bd '//bdSet(bdPos)%name%get()//':'
            write (errunit,'(A)' ) ' Format error. e.g. incorrect accordance/comment symbols: '//bdSet(bdPos)%accordCh%get()
            rcode=-1
            return
        endif

        nvars=tpCount( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%separatorCh%get() )+1
        allocate (variables(nvars,0:2)); variables=0

        j=2; seplen=bdSet(bdPos)%separatorCh%ln; quoshift=0
        do k = 1,nvars
            seppos=tpIndex( bdSet(bdPos)%bdStr%get(),bdSet(bdPos)%separatorCh%get(),cnt=k+quoshift )

            ! ====> one of modifications. need time to be tested <====
            if ( tpQuoted(bdSet(bdPos)%bdStr%get(), seppos) ) then
                quoshift=quoshift+1
                cycle
            endif

            if (seppos.LE.0) seppos=bdSet(bdPos)%bdStr%ln
            variables(k,0)=1
            variables(k,1)=j
            variables(k,2)=seppos-1
            j=seppos+seplen
        enddo
        !do k = 1,nvars
        !    write (*,*) variables(k,1),variables(k,2)
        !enddo

        do k = 1,bdSet(bdPos)%arrayLen
            j=bdSet(bdPos)%varSet(k)

            if (j.LT.0) then; write (errunit,*) 'Internal error.'; stop; endif

            ! write (*,*) 'Variable: ', variableSet(j)%name%get()

            if (variableSet(j)%opt) void=bdSetDefaultValue(variableSet(j)%address)

            if ( (bdsta.GT.bdsto).AND.(isBlank) ) cycle

            fnd=tpIndex( bdSet(bdPos)%bdStr%get(),variableSet(j)%name%get(),endp=true )

            if ( (fnd.EQ.0).AND.variableSet(j)%opt) then
                cycle
            elseif ((fnd.EQ.0).AND. .NOT.variableSet(j)%opt) then
                write (errunit,'(A\)') 'bd '//bdSet(bdPos)%name%get()//':'
                write (errunit,'(A )') ' Non-optional variable '//variableSet(j)%name%get()//' missing.'
                rcode=-2
                cycle
            endif

            do l = 1,nvars
                !write (*,*) variableSet(j)%name%get(),variables(l,1),fnd,variables(l,2)
                if ( (fnd.GE.variables(l,1)).AND.(fnd.LE.variables(l,2)) ) then
                    variables(l,0)=0
                    exit
                endif
            enddo

            if ( tpQuoted(bdSet(bdPos)%bdStr%get(),fnd) ) cycle

            if (shortKindLabels(variableSet(j)%kind).EQ.'lo') then

                ! ====> one of modifications. need time to be tested <====
                pEnd=tpIndex( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%separatorCh%get(), start=fnd)
                if (pEnd.EQ.0) pEnd=bdSet(bdPos)%bdStr%ln-1

                pStart=tpIndex( bdSet(bdPos)%bdStr%get(ustop=pEnd),bdSet(bdPos)%accordCh%get(), start=fnd)

                if (pStart.EQ.0) then
                    void=modify( variableSet(j)%address, true )
                    !write (*,*) variableSet(j)%name%get(), ' modified'
                    cycle
                endif
            endif

            pStart=tpIndex( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%accordCh%get()   , start=fnd, endp=true)+1
            pEnd  =tpIndex( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%separatorCh%get(), start=fnd)-1

            ! ====> one of modifications. need time to be tested <====
            if (tpQuoted(bdSet(bdPos)%bdStr%get(),pEnd)) then
                do while ( tpQuoted(bdSet(bdPos)%bdStr%get(),pEnd) )
                    pEnd=tpIndex( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%separatorCh%get(), start=pEnd+1)
                enddo
                pEnd=pEnd-1
            endif

            ! if (bdSet(bdPos)%name%get().EQ.'rdm') then
            !     write (*,*) 'Accord ',bdSet(bdPos)%accordCh%get()
            !     write (*,*) 'Separator ',bdSet(bdPos)%separatorCh%get()
            !     write (*,*) 'Bounds ',pStart,pEnd,bdSet(bdPos)%bdStr%ln
            !     write (*,*) 'Statement ',bdSet(bdPos)%bdStr%get(pStart,pEnd)
            ! endif

            if (pEnd.GT.0) then
                do while ( tpQuoted(bdSet(bdPos)%bdStr%get(),pEnd) )
                    pEnd=tpIndex( bdSet(bdPos)%bdStr%get(), bdSet(bdPos)%separatorCh%get(), start=pEnd+1)
                enddo
            else
                pEnd=bdSet(bdPos)%bdStr%ln-1
            endif
            ! if (bdSet(bdPos)%name%get().EQ.'rdm') then
            !     write (*,*) variableSet(j)%name%get(),pStart,pEnd
            ! endif

            if ((pStart.EQ.1) .OR. (pstart.GT.pEnd)) then
                write (errunit,'(A\)') 'bd '//bdSet(bdPos)%name%get()//':'
                write (errunit,'(A )') ' Non-boolean variable '//variableSet(j)%name%get()//' requires value.'
                rcode=-3
                cycle
            end if

            ustr=bdSet(bdPos)%bdStr%get(pStart,pEnd)
            ! if (bdSet(bdPos)%name%get().EQ.'rdm') then
            !     write (*,*) '#'//ustr%get()//'#'
            !     write (*,*) variableSet(j)%name%get()//' ===> '//bdSet(bdPos)%bdStr%get(pStart,pEnd)//' ===> '//variableSet(j)%expect%get()
            ! endif

            if ( variableSet(j)%expect%get().NE.'any' ) then
                err=checkOnExpect(j,trim(tpDeQuote(ustr%get())))
                !write (*,*) '#'//trim(tpDeQuote(ustr%get()))//'#',err
                if (err.LT.0) then
                    write (errunit,'(A\)') 'bd '//bdSet(bdPos)%name%get()//':'
                    write (errunit,'(A )') ' Unexpected value for variable: '//variableSet(j)%name%get()//'='&
                                          //trim(tpDeQuote(ustr%get()))//'. Expect '//variableSet(j)%expect%get()
                    rcode=-4
                    cycle
                endif
            endif

            err=0
            if (rcode.EQ.0) then
                select case ( shortKindLabels(variableSet(j)%kind) )
                    case('re'); void=modify( variableSet(j)%address,      tpRealByStr( ustr%get(),stat=err )   )
                    case('in'); void=modify( variableSet(j)%address,       tpIntByStr( ustr%get(),stat=err )   )
                    case('lo'); void=modify( variableSet(j)%address,       tpLogByStr( ustr%get(),stat=err )   )
                    case('ch'); void=modify( variableSet(j)%address, trim(  tpDeQuote( ustr%get(),stat=err ) ) )
                    case('uc'); void=modify( variableSet(j)%address, trim(  tpDeQuote( ustr%get(),stat=err ) ) )
                end select
            endif

            if (err.GE.0) then
                if ( shortKindLabels(variableSet(j)%kind) .in. ['ch','uc'] ) then
                    if (variableSet(j)%several) then
                        !write (*,*) variableSet(j)%name%get()
                        variableSet(j)%nvals=tpCount(ustr%get(),variableSet(j)%concatCh%get())+1
                        !write (*,*) variableSet(j)%nvals
                    endif
                endif
            endif

            if (err.LT.0) then
                write (errunit,'(A\)') 'bd '//bdSet(bdPos)%name%get()//':'
                write (errunit,'(A )') ' Format error. i.e. incorrect accordance/separator/comment symbols: '//&
                                     variableSet(j)%name%get()
                rcode=-5
                cycle
            endif
        enddo

        if ( (bdsta.GT.bdsto).AND.(isBlank) ) then
            deallocate (variables)
            return
        endif

        do k = 1,nvars
            if (variables(k,0).NE.0) then
                if (variables(k,1).GT.variables(k,2)) cycle
                write (errunit,'(A\)') 'bd '//bdSet(bdPos)%name%get()//':'
                write (errunit,'(A )') ' Unexpected variable '//bdSet(bdPos)%bdStr%get(variables(k,1),variables(k,2))
                rcode=-6
                cycle
            endif
        enddo
        deallocate (variables)
    endif

    ! if (bdSet(bdPos)%name%get().EQ.'iteration') stop

    return
    end function bdSetValues

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdSetValue(variable, value) result(rcode)
    implicit none

    integer*4        :: variable,err
    character(len=*) :: value


    rcode=0
    select case ( shortKindLabels(variableSet(variable)%kind) )
        case('re'); void=modify( variableSet(variable)%address,      tpRealByStr( value,stat=err )   )
        case('in'); void=modify( variableSet(variable)%address,       tpIntByStr( value,stat=err )   )
        case('lo'); void=modify( variableSet(variable)%address,       tpLogByStr( value,stat=err )   )
        case('ch'); void=modify( variableSet(variable)%address, trim(  tpDeQuote( value,stat=err ) ) )
        case('uc'); void=modify( variableSet(variable)%address, trim(  tpDeQuote( value,stat=err ) ) )
    end select

    return
    end function bdSetValue

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function checkOnExpect(position,str) result(rcode)
    implicit none

    character (len=*), intent(in) :: str
    integer*4        , intent(in) :: position

    character (len=variableSet(position)%concatCh%ln) :: conCh
    character (len=variableSet(position)%expect%ln)   :: expect

    integer*4              :: k,n

    real*8   , allocatable :: realListValues(:)
    integer*4, allocatable :: intListValues(:)
    real*8                 :: realvalue,rvar,rub,rlb,rstep
    integer*4              :: intvalue,cnt


    rcode=-1
    conCh =variableSet(position)%concatCh%get()
    expect=variableSet(position)%expect%get()

    !ptype=checkExpectFormat( expect,position )

    !write (*,*) '#'//expect//'# <=== #'//str//'#'


    if (variableSet(position)%exph%kind%get().EQ.'any') then; rcode= 0; return; endif !any

!    type expectHolder
!        type(uch)         :: kind !range,list,any
!        type(uch)         :: value !ge,le,in,ins
!        type(uch)         :: listStore !(x,y,z)
!        type(uch)         :: listType !all,none
!        character (len=1) :: even,positive
!        real*8            :: ub,lb,step
!   end type expectHolder


!    sta=tpIndex(expect,'(')+1; sto=tpIndex(expect,')')-1
    select case ( shortKindLabels(variableSet(position)%kind) )

        case ('re')
            realvalue=tpRealByStr(str)
            rlb  =variableSet(position)%exph%lb
            rub  =variableSet(position)%exph%ub
            rstep=variableSet(position)%exph%step
            select case (variableSet(position)%exph%kind%get())
                case ('range')
                    select case (variableSet(position)%exph%value%get())
                        case ('any','norange') !ignoring even/odd
                            select case (variableSet(position)%exph%positive)
                                case('n'); rcode=0; return
                                case('+'); if (realvalue.GE.0) then; rcode=0; return; endif
                                case('-'); if (realvalue.LE.0) then; rcode=0; return; endif
                            end select

                        case ('ge'); if (realvalue.GE.rlb) then; rcode=0; return; endif

                        case ('le'); if (realvalue.LE.rub) then; rcode=0; return; endif

                        case ('in')
                            if ((realvalue.GE.rlb).AND.(realvalue.LE.rub)) then
                                rcode=0; return
                            endif

                        case ('ins')
                            do rvar = rlb,rub,rstep
                                if ( abs(rvar-realvalue).LT.epsilon(rvar)*1000 ) then
                                    rcode=0; return
                                endif
                            enddo
                    end select

                case ('list')
                    n=tpCount(variableSet(position)%exph%listStore%get(),',')+1
                    allocate ( realListValues(n) )
                    realListValues=tpRealArrayByStr(variableSet(position)%exph%listStore%get(),n)

                    select case (variableSet(position)%exph%listType%get())
                        case('any')
                             do k = 1,n
                                if (abs(realListValues(k)-realvalue).LT.epsilon(realvalue)*10) then
                                    rcode=0; deallocate (realListValues); return
                                endif
                            enddo
                            rcode=-1; deallocate (realListValues); return

                        case ('none')
                             do k = 1,n
                                if (abs(realListValues(k)-realvalue).LT.epsilon(realvalue)*10) then
                                    rcode=-1; deallocate (realListValues); return
                                endif
                            enddo
                            rcode=0; deallocate (realListValues); return

                    end select

            end select

        case ('in')
            intvalue=tpIntByStr(str)
            rlb  =variableSet(position)%exph%lb
            rub  =variableSet(position)%exph%ub
            rstep=variableSet(position)%exph%step
            select case (variableSet(position)%exph%kind%get())
                case ('range')
                    !write (*,*) variableSet(position)%exph%value%get()
                    select case (variableSet(position)%exph%value%get())
                        case ('any')
                            !write (*,*) variableSet(position)%exph%positive,variableSet(position)%exph%even
                            select case (variableSet(position)%exph%positive)
                                case('n'); continue
                                case('+'); if (intvalue.LT.0) then; rcode=-1; return; endif
                                case('-'); if (intvalue.GE.0) then; rcode=-1; return; endif
                            end select

                            select case (variableSet(position)%exph%even)
                                case('n'); rcode=0; return
                                case('+'); if (mod(intvalue,2).EQ.0) then; rcode=0; return; endif
                                case('-'); if (mod(intvalue,2).EQ.1) then; rcode=0; return; endif
                            end select

                        case ('norange')
                            select case (variableSet(position)%exph%positive)
                                case('n'); continue
                                case('+'); if (intvalue.LT.0) then; rcode=-1; return; endif
                                case('-'); if (intvalue.GE.0) then; rcode=-1; return; endif
                            end select

                            select case (variableSet(position)%exph%even)
                                case('n'); rcode=0; return
                                case('+'); if (mod(intvalue,2).EQ.0) then; rcode=0; return; endif
                                case('-'); if (mod(intvalue,2).EQ.1) then; rcode=0; return; endif
                            end select

                        case ('ge')
                            if (intvalue.GE.int(rlb)) then
                                select case (variableSet(position)%exph%even)
                                    case('n'); rcode=0; return
                                    case('+'); if (mod(intvalue,2).EQ.0) then; rcode=0; return; endif
                                    case('-'); if (mod(intvalue,2).EQ.1) then; rcode=0; return; endif
                                end select
                            endif

                        case ('le')
                            if (intvalue.LE.int(rub)) then
                                select case (variableSet(position)%exph%even)
                                    case('n'); rcode=0; return
                                    case('+'); if (mod(intvalue,2).EQ.0) then; rcode=0; return; endif
                                    case('-'); if (mod(intvalue,2).EQ.1) then; rcode=0; return; endif
                                end select
                            endif

                        case ('in')
                            if ((intvalue.GE.int(rlb)).AND.(intvalue.LE.int(rub))) then
                                select case (variableSet(position)%exph%even)
                                    case('n'); rcode=0; return
                                    case('+'); if (mod(intvalue,2).EQ.0) then; rcode=0; return; endif
                                    case('-'); if (mod(intvalue,2).EQ.1) then; rcode=0; return; endif
                                end select
                            endif

                        case ('ins')
                            do k = int(rlb),int(rub),int(rstep)
                                if (k.EQ.intvalue) then
                                    select case (variableSet(position)%exph%even)
                                        case('n'); rcode=0; return
                                        case('+'); if (mod(intvalue,2).EQ.0) then; rcode=0; return; endif
                                        case('-'); if (mod(intvalue,2).EQ.1) then; rcode=0; return; endif
                                    end select
                                endif
                            enddo

                    end select

                case ('list')
                    n=tpCount(variableSet(position)%exph%listStore%get(),',')+1
                    allocate ( intListValues(n) )
                    intListValues=tpIntArrayByStr(variableSet(position)%exph%listStore%get(),n)

                    select case (variableSet(position)%exph%listType%get())
                        case('any')
                            do k = 1,n
                                if (intvalue.EQ.intListValues(k)) then
                                   rcode=0; deallocate (intListValues); return
                                endif
                            enddo
                            rcode=-1; deallocate (intListValues); return

                        case('none')
                            do k = 1,n
                                if (intvalue.EQ.intListValues(k)) then
                                   rcode=-1; deallocate (intListValues); return
                                endif
                            enddo
                            rcode=0; deallocate (intListValues); return

                    end select

            end select

        case ('uc','ch')
            !write (*,*) 'Kind '//variableSet(position)%exph%kind%get()
            select case (variableSet(position)%exph%kind%get())
                case ('range'); rcode=-1; return

                case ('list')
                    !write (*,*) 'List type '//variableSet(position)%exph%listType%get()
                    select case (variableSet(position)%exph%listType%get())
                        case('any')
                            if (variableSet(position)%several) then
                                voidl=tpSplit(str,conCh); cnt=0
                                do k = 1,tpSplitLen
                                    if (tpRetSplit(str,k) .in. variableSet(position)%exph%listStore%get()) then
                                        cnt=cnt+1
                                    endif
                                enddo
                                if (cnt.EQ.tpSplitLen) rcode=0
                            else
                                !write (*,*) 'String #'//str//'#'
                                if (str .in. variableSet(position)%exph%listStore%get()) then
                                    rcode=0; return
                                endif
                            endif

                        case('none')
                            if (variableSet(position)%several) then
                                voidl=tpSplit(str,conCh)
                                do k = 1,tpSplitLen
                                    if (tpRetSplit(str,k) .in. variableSet(position)%exph%listStore%get()) then
                                        return
                                    endif
                                enddo
                            else
                                if ( .NOT.(str .in. variableSet(position)%exph%listStore%get()) ) then
                                    rcode=0; return
                                endif
                            endif
                    end select

            end select

        case ('lo')
            rcode=0

    end select

    return
    end function checkOnExpect

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function checkExpectFormat(str,position) result(rcode)
    implicit none

    character (len=*), intent(in) :: str
    integer*4        , intent(in) :: position
    character (len=len(str))      :: astr,bstr,cstr,dstr
    character (len=3)             :: switch

    character (len=1)             :: posit,even

    integer*4                     :: k,err,sta,sto,cnt,cnt2
    logical*1                     :: fndEven,fndOdd,fndNegative,fndPositive
    logical*1                     :: fndAny,fndNone
    real*8                        :: slice(3)


!    write (*,*) '#'//str//'#'
!    write (*,*) position

    if ( (.NOT. ('any' .in. str)) .AND. (.NOT. ('range' .in. str)) .AND. (.NOT. ('list' .in. str)) ) then
        rcode=-1; return
    endif

    astr=str; bstr=tpFill(bstr); cstr=tpFill(cstr); dstr=tpFill(dstr)

    rcode=-1
    if (tpStartsWith(str,'any')                           ) rcode=1
    if (tpStartsWith(str,'range(').AND.tpEndsWith(str,')')) rcode=2
    if (tpStartsWith(str,'list(' ).AND.tpEndsWith(str,')')) rcode=3

    if (rcode.EQ.1) then
        call variableSet(position)%exph%init()
        variableSet(position)%exph%kind='any'
    endif

    if (rcode.EQ.2) then
        switch='111'; slice=NaN; posit='n'; even='n'
        sta=tpIndex(str,'range(',endp=true)+1
        sto=tpIndex(str,')'     ,rev =true)-1

        fndEven=false; fndOdd=false; fndPositive=false; fndNegative=false

        if (sto.LT.sta) then
            rcode=-2; return
        endif

        if ('|' .in. str) sto=tpIndex(str,'|')-1

        dstr=str(sta:sto); cnt=tpCount(str,':')
        select case (cnt)
            case (1:2)

                voidl=tpSplit(trim(dstr),':')
                do k = 1,tpSplitLen
                    if (.NOT.tpSetAccordance( tpSplitHold(k)%get(), tpRealNumber() )) then
                        rcode=-3; return
                    endif
                enddo

                if (cnt.EQ.2) slice(1:3)=tpRealArrayByStr(trim(dstr),3,':')
                if (cnt.EQ.1) slice(1:2)=tpRealArrayByStr(trim(dstr),2,':')

                do k = 1,3
                    if (isNan(slice(k))) switch(k:k)='0'
                enddo

                select case (switch)
                    case ('000','010','100'); continue
                    case ('001','101','011'); rcode=-4; return
                    case ('110')
                        if (slice(1).GT.slice(2)) then
                            rcode=-5; return
                        endif
                    case ('111')
                        if ( (slice(3).GT.0) .AND. (slice(1).GT.slice(2)) ) then
                            rcode=-6; return
                        endif
                        if ( (slice(3).LT.0) .AND. (slice(2).GT.slice(1)) ) then
                            rcode=-7; return
                        endif
                end select

            case (0)
                switch='000'
                if (sto.GE.sta) then
                    rcode=-8; return
                endif

            case default
                rcode=-9; return

        end select

        if ('|' .in. astr) then
            voidl=tpSplit(astr,'|')
            bstr=tpSplitHold(2)%get()
            bstr(len_trim(bstr):len_trim(bstr))=' '

            fndEven    =tpIndex(bstr,'even'    ).GT.0
            fndOdd     =tpIndex(bstr,'odd'     ).GT.0
            fndPositive=tpIndex(bstr,'positive').GT.0
            fndNegative=tpIndex(bstr,'negative').GT.0

            if (fndEven.AND.fndOdd) then
                rcode=-11; return
            endif

            if (fndPositive.AND.fndNegative) then
                rcode=-12; return
            endif

            if (fndPositive) posit='+'
            if (fndNegative) posit='-'
            if (fndEven)     even='+'
            if (fndOdd)      even='-'

            if (.NOT.tpSetAccordance( bstr, tpLetters()//'; ' ) ) then
                rcode=-13; return
            endif

            if (';' .in. bstr) then
                voidl=tpSplit(bstr,';')
                do k = 1,tpSplitLen
                    cstr=tpRetSplit(bstr,k)
                    if ( .NOT. (cstr .in. ['even','odd','positive','negative']) ) then
                        rcode=-14; return
                    endif
                enddo
            endif

            select case (switch)
                case ('000'); continue
                case ('010','100','110','111')
                    if (fndPositive.OR.fndNegative) then
                        rcode=-15; return
                    endif
            end select
        endif
        !write (*,*) switch
        !write (*,*) 'slice',slice
        !write (*,'(1X,A,4(1X,L)\)') 'EOPN',fndEven,fndOdd,fndPositive,fndNegative

        call variableSet(position)%exph%init()

        variableSet(position)%exph%kind='range'
        variableSet(position)%exph%lb  =slice(1)
        variableSet(position)%exph%ub  =slice(2)
        variableSet(position)%exph%step=slice(3)

        variableSet(position)%exph%positive=posit
        variableSet(position)%exph%even    =even

        select case (switch)
            case ('000')
                if ( (posit.EQ.'n') .AND. (even.EQ.'n') ) then
                    variableSet(position)%exph%kind='any'
                endif
                variableSet(position)%exph%value='norange'

            case ('010'); variableSet(position)%exph%value='le'
            case ('100'); variableSet(position)%exph%value='ge'
            case ('110'); variableSet(position)%exph%value='in'
            case ('111'); variableSet(position)%exph%value='ins'
        end select
    endif

    if (rcode.EQ.3) then
        sta=tpIndex(str,'list(',endp=true)+1
        sto=tpIndex(str,')'    ,rev =true)-1

        fndAny=false; fndNone=false

        if (sta.GT.sto) then
            rcode=-16; return
        endif

        bstr=str(sta:sto)

        cnt=tpCount(bstr,':')
        select case (cnt)
            case (1)
                voidl=tpSplit(trim(bstr),':')

                if (tpSplitHold(1)%ln.EQ.0) then
                    fndAny=true
                else
                    if ( .NOT. (tpSplitHold(1)%get() .in. ['any','none']) ) then
                        rcode=-17; return
                    endif
                    fndAny = 'any'  .in. tpSplitHold(1)%get()
                    fndNone= 'none' .in. tpSplitHold(1)%get()
                endif

                if (tpCount( tpSplitHold(2)%get() , ',' ).EQ.0) then
                    rcode=-18; return
                endif

                if (fndAny .AND. fndNone) then
                    rcode=-19; return
                endif

                voidl=tpSplit( tpSplitHold(2)%get() , ',' ); cnt2=0
                do k = 1,tpSplitLen
                    cnt2=cnt2+tpSplitHold(k)%ln
                enddo
                if (cnt2.EQ.0) then
                    rcode=-20; return
                endif

            case (0)
                if (tpCount( trim(bstr), ',' ).EQ.0) then
                    rcode=-21; return
                endif

                fndAny=true

                voidl=tpSplit(trim(bstr),','); cnt2=0
                do k = 1,tpSplitLen
                    cnt2=cnt2+tpSplitHold(k)%ln
                enddo
                if (cnt2.EQ.0) then
                    rcode=-22; return
                endif

            case default
                rcode=-23; return

        end select

        call variableSet(position)%exph%init()

        variableSet(position)%exph%kind='list'
        if (fndAny)  variableSet(position)%exph%listType='any'
        if (fndNone) variableSet(position)%exph%listType='none'

        select case (cnt)
            case(1)
                sta=tpIndex(bstr,':')+1
                variableSet(position)%exph%listStore='{'//bstr(sta:len_trim(bstr))//'}'

            case(0)
                variableSet(position)%exph%listStore='{'//trim(bstr)//'}'

        end select
    endif

    return
    end function checkExpectFormat

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_r16(var,name,opt,def,vary,expect,potentiate,description) result(rcode)
    implicit none

    integer*4, parameter        :: vkind=16
    real(vkind), target         :: var
    character (len=*)           :: name
    logical*1        , optional :: opt,vary
    real(vkind)      , optional :: def
    character (len=*), optional :: expect,description
    character (len=1), optional :: potentiate

    logical*1                   :: dopt,dvary
    integer*4                   :: err
    real(vkind)                 :: ddef


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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name      =tpAdjustl(name)
    variableSet(varAppend)%address   =loc(var)
    variableSet(varAppend)%position  =varAppend
    variableSet(varAppend)%vlen      =vkind
    variableSet(varAppend)%opt       =dopt
    variableSet(varAppend)%variable  =dvary
    !variableSet(varAppend)%defstor   =transfer(ddef,storeTemplate)
    variableSet(varAppend)%concatCh  =''
    variableSet(varAppend)%nvals     =0
    variableSet(varAppend)%potentiate='n'

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(potentiate)) then
        variableSet(varAppend)%potentiate=potentiate
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    select case (variableSet(varAppend)%potentiate)
        case('+'); variableSet(varAppend)%defstor   =transfer(10.**(+ddef),storeTemplate)
        case('-'); variableSet(varAppend)%defstor   =transfer(10.**(-ddef),storeTemplate)
        case('n'); variableSet(varAppend)%defstor   =transfer(ddef,storeTemplate)
    end select

    variableSet(varAppend)%kind  =0
    variableSet(varAppend)%pntr16=>var
    variableSet(varAppend)%varsta=transfer(real(0,vkind),storeTemplate)
    variableSet(varAppend)%varsto=transfer(real(0,vkind),storeTemplate)
    variableSet(varAppend)%varste=transfer(real(0,vkind),storeTemplate)

    return
    end function shareVariable_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_r8(var,name,opt,def,vary,expect,potentiate,description) result(rcode)
    implicit none

    integer*4, parameter        :: vkind=8
    real(vkind), target         :: var
    character (len=*)           :: name
    logical*1        , optional :: opt,vary
    real(vkind)      , optional :: def
    character (len=*), optional :: expect,description
    character (len=1), optional :: potentiate

    logical*1                   :: dopt,dvary
    integer*4                   :: err
    real(vkind)                 :: ddef


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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name      =tpAdjustl(name)
    variableSet(varAppend)%address   =loc(var)
    variableSet(varAppend)%position  =varAppend
    variableSet(varAppend)%vlen      =vkind
    variableSet(varAppend)%opt       =dopt
    variableSet(varAppend)%variable  =dvary
    !variableSet(varAppend)%defstor   =transfer(ddef,storeTemplate)
    variableSet(varAppend)%concatCh  =''
    variableSet(varAppend)%nvals     =0
    variableSet(varAppend)%potentiate='n'

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(potentiate)) then
        variableSet(varAppend)%potentiate=potentiate
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    select case (variableSet(varAppend)%potentiate)
        case('+'); variableSet(varAppend)%defstor   =transfer(10.**(+ddef),storeTemplate)
        case('-'); variableSet(varAppend)%defstor   =transfer(10.**(-ddef),storeTemplate)
        case('n'); variableSet(varAppend)%defstor   =transfer(ddef,storeTemplate)
    end select

    variableSet(varAppend)%kind  =1
    variableSet(varAppend)%pntr8 =>var
    variableSet(varAppend)%varsta=transfer(real(0,vkind),storeTemplate)
    variableSet(varAppend)%varsto=transfer(real(0,vkind),storeTemplate)
    variableSet(varAppend)%varste=transfer(real(0,vkind),storeTemplate)

    return
    end function shareVariable_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_r4(var,name,opt,def,vary,expect,potentiate,description) result(rcode)
    implicit none

    integer*4, parameter        :: vkind=4
    real(vkind), target         :: var
    character (len=*)           :: name
    logical*1        , optional :: opt,vary
    real(vkind)      , optional :: def
    character (len=*), optional :: expect,description
    character (len=1), optional :: potentiate

    logical*1                   :: dopt,dvary
    integer*4                   :: err
    real(vkind)                 :: ddef


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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name      =tpAdjustl(name)
    variableSet(varAppend)%address   =loc(var)
    variableSet(varAppend)%position  =varAppend
    variableSet(varAppend)%vlen      =vkind
    variableSet(varAppend)%opt       =dopt
    variableSet(varAppend)%variable  =dvary
    !variableSet(varAppend)%defstor   =transfer(ddef,storeTemplate)
    variableSet(varAppend)%concatCh  =''
    variableSet(varAppend)%nvals     =0
    variableSet(varAppend)%potentiate='n'

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(potentiate)) then
        variableSet(varAppend)%potentiate=potentiate
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    select case (variableSet(varAppend)%potentiate)
        case('+'); variableSet(varAppend)%defstor   =transfer(10.**(+ddef),storeTemplate)
        case('-'); variableSet(varAppend)%defstor   =transfer(10.**(-ddef),storeTemplate)
        case('n'); variableSet(varAppend)%defstor   =transfer(ddef,storeTemplate)
    end select

    variableSet(varAppend)%kind  =2
    variableSet(varAppend)%pntr4 =>var
    variableSet(varAppend)%varsta=transfer(real(0,vkind),storeTemplate)
    variableSet(varAppend)%varsto=transfer(real(0,vkind),storeTemplate)
    variableSet(varAppend)%varste=transfer(real(0,vkind),storeTemplate)

    return
    end function shareVariable_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_i8(var,name,opt,def,vary,expect,description) result(rcode)
    implicit none
    logical*1, optional         :: opt,vary
    character (len=*)           :: name
    logical*1                   :: dopt,dvary
    character (len=*), optional :: expect,description

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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=dvary
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind  =3
    variableSet(varAppend)%pnti8 =>var
    variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

    return
    end function shareVariable_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_i4(var,name,opt,def,vary,expect,description) result(rcode)
    implicit none
    logical*1, optional         :: opt,vary
    character (len=*)           :: name
    logical*1                   :: dopt,dvary
    character (len=*), optional :: expect,description

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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=dvary
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind  =4
    variableSet(varAppend)%pnti4 =>var
    variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

    return
    end function shareVariable_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_i2(var,name,opt,def,vary,expect,description) result(rcode)
    implicit none
    logical*1, optional         :: opt,vary
    character (len=*)           :: name
    logical*1                   :: dopt,dvary
    character (len=*), optional :: expect,description

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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=dvary
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind  =5
    variableSet(varAppend)%pnti2 =>var
    variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

    return
    end function shareVariable_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_i1(var,name,opt,def,vary,expect,description) result(rcode)
    implicit none
    logical*1, optional         :: opt,vary
    character (len=*)           :: name
    logical*1                   :: dopt,dvary
    character (len=*), optional :: expect,description

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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=dvary
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind  =6
    variableSet(varAppend)%pnti1 =>var
    variableSet(varAppend)%varsta=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varsto=transfer(int(0,vkind),storeTemplate)
    variableSet(varAppend)%varste=transfer(int(0,vkind),storeTemplate)

    return
    end function shareVariable_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_l8(var,name,opt,def,description) result(rcode)
    implicit none
    logical*1, optional         :: opt
    character (len=*)           :: name
    logical*1                   :: dopt

    character (len=*), optional :: description

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

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=false
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%several =false
    variableSet(varAppend)%expect  ='any'
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind =7
    variableSet(varAppend)%pntl8=>var

    return
    end function shareVariable_l8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_l4(var,name,opt,def,description) result(rcode)
    implicit none
    logical*1, optional         :: opt
    character (len=*)           :: name
    logical*1                   :: dopt

    character (len=*), optional :: description

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

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=false
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%several =false
    variableSet(varAppend)%expect  ='any'
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind =8
    variableSet(varAppend)%pntl4=>var

    return
    end function shareVariable_l4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_l2(var,name,opt,def,description) result(rcode)
    implicit none
    logical*1, optional         :: opt
    character (len=*)           :: name
    logical*1                   :: dopt

    character (len=*), optional :: description

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

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=false
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%several =false
    variableSet(varAppend)%expect  ='any'
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind =9
    variableSet(varAppend)%pntl2=>var

    return
    end function shareVariable_l2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_l1(var,name,opt,def,description) result(rcode)
    implicit none
    logical*1, optional         :: opt
    character (len=*)           :: name
    logical*1                   :: dopt

    character (len=*), optional :: description

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

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%vlen    =vkind
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=false
    variableSet(varAppend)%defstor =transfer(ddef,storeTemplate)
    variableSet(varAppend)%several =false
    variableSet(varAppend)%expect  ='any'
    variableSet(varAppend)%concatCh=''
    variableSet(varAppend)%nvals   =0

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    variableSet(varAppend)%kind =10
    variableSet(varAppend)%pntl1=>var

    return
    end function shareVariable_l1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_c(var,name,opt,def,expect,several,concat,description) result(rcode)
    implicit none
    character (len=*), target   :: var
    character (len=*), optional :: def
    logical*1, optional         :: opt,several
    character (len=*), optional :: expect,concat,description

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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name    =tpAdjustl(name)
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
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(concat)) then
        variableSet(varAppend)%concatCh=concat
    else
        variableSet(varAppend)%concatCh='+'
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    if (present(def)) variableSet(varAppend)%defc=trim(def)

    return
    end function shareVariable_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*8 function shareVariable_uc(var,name,opt,def,expect,several,concat,description) result(rcode)
    implicit none
    type(uch), target           :: var
    character (len=*), optional :: def
    logical*1, optional         :: opt,several
    character (len=*), optional :: expect,concat,description

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
        if (checkExpectFormat(expect,varAppend+1).LT.0) then
            write (errunit,*) name//': Incorrect format for expectation values: '//expect//'.'
            rcode=-4; return
        endif
    endif

    varAppend=varAppend+1

    variableSet(varAppend)%name    =tpAdjustl(name)
    variableSet(varAppend)%address =loc(var)
    variableSet(varAppend)%position=varAppend
    variableSet(varAppend)%kind    =12
    variableSet(varAppend)%vlen    =var%ln
    variableSet(varAppend)%pntuch  =>var
    variableSet(varAppend)%opt     =dopt
    variableSet(varAppend)%variable=false
    variableSet(varAppend)%several =dseveral
    variableSet(varAppend)%nvals   =0

    if (present(expect)) then
        variableSet(varAppend)%expect=expect
    else
        variableSet(varAppend)%expect='any'
    endif

    if (present(concat)) then
        variableSet(varAppend)%concatCh=concat
    else
        variableSet(varAppend)%concatCh='+'
    endif

    if (present(description)) then
        variableSet(varAppend)%description=description
    else
        variableSet(varAppend)%description=''
    endif

    if (present(def)) variableSet(varAppend)%defc=trim(def)

    return
    end function shareVariable_uc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function bdSetDefaultValue(address) result(ret)
    implicit none

    integer*8 :: address
    integer*4 :: i


    i=find(variableSet(:)%address,address)
    if (i.GT.0) then
        if (variableSet(i)%opt) then
            select case(variableSet(i)%kind)
                case ( 0); variableSet(i)%pntr16=transfer(variableSet(i)%defstor,variableSet(i)%pntr16)
                case ( 1); variableSet(i)%pntr8 =transfer(variableSet(i)%defstor,variableSet(i)%pntr8 )
                case ( 2); variableSet(i)%pntr4 =transfer(variableSet(i)%defstor,variableSet(i)%pntr4 )
                case ( 3); variableSet(i)%pnti8 =transfer(variableSet(i)%defstor,variableSet(i)%pnti8 )
                case ( 4); variableSet(i)%pnti4 =transfer(variableSet(i)%defstor,variableSet(i)%pnti4 )
                case ( 5); variableSet(i)%pnti2 =transfer(variableSet(i)%defstor,variableSet(i)%pnti2 )
                case ( 6); variableSet(i)%pnti1 =transfer(variableSet(i)%defstor,variableSet(i)%pnti1 )
                case ( 7); variableSet(i)%pntl8 =transfer(variableSet(i)%defstor,variableSet(i)%pntl8 )
                case ( 8); variableSet(i)%pntl4 =transfer(variableSet(i)%defstor,variableSet(i)%pntl4 )
                case ( 9); variableSet(i)%pntl2 =transfer(variableSet(i)%defstor,variableSet(i)%pntl2 )
                case (10); variableSet(i)%pntl1 =transfer(variableSet(i)%defstor,variableSet(i)%pntl1 )

                ! Please, let me 'shoot myself in the foot'.
                ! Will raise 'array bounds exceeded', but works fine without boundary check.
                ! Source file should be compiled without "-CB" for ifort and "-fbounds-check" for gfortran
                case (11)
                    variableSet(i)%pntch(1:variableSet(i)%vlen)=tpFill(variableSet(i)%vlen)
                    variableSet(i)%pntch(1:variableSet(i)%defc%ln)=variableSet(i)%defc%get()

                case (12)
                    variableSet(i)%pntuch=variableSet(i)%defc

            end select
        endif
    endif

    ret=0; return
    end function bdSetDefaultValue

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine expectInitiate(this)
    implicit none

    class(expectHolder), intent(out) :: this


    this%kind='None'; this%value='None';
    this%listStore='None'; this%listType='None'
    this%even='n'; this%positive='n'
    this%ub=NaN; this%lb=NaN; this%step=NaN

    return
    end subroutine expectInitiate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function modify_c(address,val) result(ret)
    implicit none

    character (len=*), intent(in) :: val
    integer*8                     :: address
    integer*4                     :: i


    i=find(variableSet(:)%address,address)
    if (i.GT.0) then

        if (variableSet(i)%kind.EQ.12) then
            variableSet(i)%pntuch=val
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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function modify_r(address,val) result(ret)
    implicit none

    real*8, intent(in) :: val
    real*8             :: uval
    integer*8          :: address
    integer*4          :: i


    i=find(variableSet(:)%address,address)
    if (i.GT.0) then
        if ((variableSet(i)%kind.NE.0).AND.(variableSet(i)%kind.NE.1).AND.(variableSet(i)%kind.NE.2)) then
            ret=-1; return
        endif

        !write (*,*) variableSet(i)%name%get(),variableSet(i)%potentiate

        select case (variableSet(i)%potentiate)
            case ('+'); uval=10.**(+val)
            case ('-'); uval=10.**(-val)
            case ('n'); uval=val
        end select

        select case (variableSet(i)%vlen)
            case (16); variableSet(i)%pntr16=real(uval,kind=16)
            case ( 8); variableSet(i)%pntr8 =real(uval,kind=8 )
            case ( 4); variableSet(i)%pntr4 =real(uval,kind=4 )
        end select
    else
        ret=-2; return
    endif

    ret=0; return
    end function modify_r

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function modify_l(address,val) result(ret)
    implicit none

    logical*1, intent(in) :: val
    integer*8             :: address
    integer*4             :: i


    i=find(variableSet(:)%address,address)
    if (i.GT.0) then
        if ((variableSet(i)%kind.NE.7).AND.(variableSet(i)%kind.NE.8).AND.&
            (variableSet(i)%kind.NE.9).AND.(variableSet(i)%kind.NE.10)) then
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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine bdFinalize
    implicit none


    !continue

    return
    end subroutine bdFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module datablock
