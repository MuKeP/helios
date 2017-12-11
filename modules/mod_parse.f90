    module argsParser

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,      only: assignment(=)
    use glob,      only: mid,uch,find,true,false,osSeparator,ivarVector
    use txtParser, only: tpRealByStr,tpIntByStr,tpDeQuote,tpFill,tpLowerCase
    use txtParser, only: operator(.in.)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: apVersion='4.300'
    character (len=*), parameter :: apDate   ='2017.12.10'
    character (len=*), parameter :: apAuthor ='Anton B. Zakharov'

    integer*4, parameter :: maxArgsListLen=128  ! maximum number of defined arguments
    integer*4, parameter :: maxGroupListLen=16  ! maximum number of groups
    integer*4, parameter :: maxSynonymListLen=5 ! maximum number of the synonyms for arguments
    integer*4, parameter :: argLen=512          ! length of the "getarg" return string
    integer*4, parameter :: appNameLen=512      ! length of the application name
    integer*4, parameter :: maxRulesListLen=512 ! length of the rules (xor,and) list

    character (len=*), parameter :: splitSymbol='|',concatSynonymsSymbol=', '

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type variable
        integer*8          :: address
        integer*4          :: position,kind,vlen

                                     !kind
        real*16  , pointer :: pntr16 ! 0
        real*8   , pointer :: pntr8  ! 1
        real*4   , pointer :: pntr4  ! 2
        integer*8, pointer :: pnti8  ! 3
        integer*4, pointer :: pnti4  ! 4
        integer*2, pointer :: pnti2  ! 5
        integer*1, pointer :: pnti1  ! 6
        logical*8, pointer :: pntl8  ! 7
        logical*4, pointer :: pntl4  ! 8
        logical*2, pointer :: pntl2  ! 9
        logical*1, pointer :: pntl1  !10
        character, pointer :: pntch  !11
        type(uch), pointer :: pntuch !12
    end type variable

    type cArgument
        type(uch) :: str(maxSynonymListLen)

        integer*4 :: group,var,synonymsListLen,synonym,position
        logical*1 :: csensitive,default,required,expect,found,value
        type(uch) :: description,expectValue
    end type

    type rule
        type(ivarVector) :: argsList
        type(uch)        :: type
    end type rule

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type(variable)  :: variableSet(maxArgsListLen)
    type(cArgument) :: arg(0:maxArgsListLen)
    type(uch)       :: groupSet(-1:maxGroupListLen)
    type(rule)      :: ruleSet(maxRulesListLen)

    type(uch), allocatable :: holdArguments(:),holdUnexpectArguments(:)
    logical*1, allocatable :: argumentUsed(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type(uch) :: appName,preamble,comment

    integer*4 :: void,iounit=6,errunit=6,argsListStart=1,lastError=0
    integer*4 :: argsListLen=0,groupListLen=0,varListLen=0
    logical*1 :: alreadyParsed=false,appNameSet=false
    logical*1 :: preambleSet=false,commentSet=false
    integer*4 :: argsCount=0,unexpectArgsCount=0,outputArgumentLen=1
    integer*4 :: ruleListLen=0

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    interface apShare
        module procedure shareVariable_r16,shareVariable_r8,shareVariable_r4,&
                         shareVariable_i8 ,shareVariable_i4,shareVariable_i2,&
                         shareVariable_i1 ,shareVariable_l8,shareVariable_l4,&
                         shareVariable_l2 ,shareVariable_l1,shareVariable_c ,&
                         shareVariable_uc
    end interface apShare

    interface modify
        module procedure modify_c,modify_r,modify_i,modify_l
    end interface modify

    interface getArgumentValue
        module procedure getArgumentValue_i, getArgumentValue_c
    end interface getArgumentValue

    interface apArgumentFound
        module procedure apArgumentFound_c,apArgumentFound_i
    end interface apArgumentFound

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    public :: addArg,addGroup,parseArgs,setApplicationName,setPreamble,     &
              setComment,showHelp,getDefaultArgument,getArgumentValue,      &
              ap_errDescription,setGroupMembership,apSetIOunit,apSetERRunit,&
              ap_popLastError,apFinalize,apShare,apArgumentFound,           &
              apGetCommandLine,appendRule

    public :: apVersion,apDate,apAuthor

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function addArg(argument,description,var,expect,default,casesensitive,required,group,expectvalue) result(rcode)
    implicit none

    character (len=*)                :: argument,description

    character (len=*), optional      :: expectvalue
    logical*1, optional              :: expect,default,casesensitive,required
    integer*4, optional              :: group,var

    character (len=len(expectvalue)) :: dexpectvalue
    logical*1                        :: dexpect,ddefault,dcasesensitive,drequired
    integer*4                        :: dgroup,dvar

    logical*1                        :: dvalue

    type(cArgument)                  :: recarg
    integer*4                        :: err,ind,k


    ! parse function arguments
    if (argsListLen.EQ.maxArgsListLen)   then; rcode=-1; lastError=rcode; return; endif

    dexpect=false;        if (present(expect))        dexpect=expect
    ddefault=false;       if (present(default))       ddefault=default
    dcasesensitive=false; if (present(casesensitive)) dcasesensitive=casesensitive
    drequired=false;      if (present(required))      drequired=required
    dgroup=0;             if (present(group))         dgroup=group
    dvar=0;               if (present(var))           dvar=var
    dexpectvalue='';      if (present(expectvalue))   dexpectvalue=trim(expectvalue)

    if (dexpect.AND.ddefault.AND.len_trim(dexpectvalue).EQ.0) then; rcode=-4; lastError=rcode; return; endif
    err=splitArgument(argument,recarg)

    if (err.EQ.-1) then; rcode=-5; lastError=rcode; return; endif
    if (err.EQ.-2) then; rcode=-6; lastError=rcode; return; endif

    if (.NOT.dcasesensitive) then
        do k = 1,recarg%synonymsListLen
            recarg%str(k)=tpLowerCase( recarg%str(k)%get() )
        enddo
    endif

    err=checkExistence(recarg)

    if (err.EQ.-1)                                            then; rcode= -7; lastError=rcode; return; endif
    if (ddefault.AND.drequired)                               then; rcode= -8; lastError=rcode; return; endif
    if ( (dgroup.LT.0) .OR. (dgroup.GT.groupListLen) )        then; rcode= -9; lastError=rcode; return; endif
    if ( (.NOT.dexpect)  .AND. (len_trim(dexpectvalue).GT.0)) then; rcode=-10; lastError=rcode; return; endif
    if ( (.NOT.ddefault) .AND. (len_trim(dexpectvalue).GT.0)) then; rcode=-11; lastError=rcode; return; endif

    dvalue=(.NOT.dexpect).AND.ddefault

    argsListLen=argsListLen+1; ind=argsListLen

    arg(ind)%var=0; if (present(var)) arg(ind)%var=var

    select case (arg(ind)%var)
        case(-1);     rcode=-12; lastError=rcode; return
        case(-2);     rcode=-13; lastError=rcode; return
        case( 0);     continue
        case default; continue
    end select

    arg(ind)%csensitive=dcasesensitive
    arg(ind)%default=ddefault
    arg(ind)%required=drequired
    arg(ind)%expect=dexpect
    arg(ind)%value=dvalue
    arg(ind)%expectValue=dexpectvalue
    arg(ind)%group=dgroup
    arg(ind)%str=recarg%str
    arg(ind)%synonymsListLen=recarg%synonymsListLen
    arg(ind)%description=trim(description)
    arg(ind)%found=false
    arg(ind)%synonym=0
    arg(ind)%position=0

    rcode=ind; return
    end function addArg

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function addGroup(group) result(rcode)
    implicit none

    character (len=*) :: group
    integer*4         :: i


    if (groupListLen.EQ.maxGroupListLen) then; rcode=-101; lastError=rcode; return; endif

    do i = 1,groupListLen
        if (trim(group).EQ.trim( groupSet(i)%get() )) then
            rcode=i; return
        endif
    enddo

    groupListLen=groupListLen+1
    groupSet(groupListLen)=trim(adjustl(group))

    rcode=groupListLen; return
    end function addGroup

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function parseArgs(onlyDeclared) result(rcode)
    implicit none

    logical*1, optional    :: onlyDeclared
    logical*1              :: donlyDeclared
    character (len=argLen) :: tempStr

    integer*4              :: fnd,i,j


    donlyDeclared=false; if (present(onlyDeclared)) donlyDeclared=onlyDeclared

    if (alreadyParsed) then; rcode=-201; return; lastError=rcode; endif

    void=addHelpArgument()
    if (.NOT.appNameSet) void=getApplicationName()

    argsCount=nargs()-1

    outputArgumentLen=0
    do i = 0,argsListLen
        outputArgumentLen=max(outputArgumentLen,len_trim(collectSynonyms(arg(i))))
    enddo

    if (groupListLen.EQ.0) then
        do i = 1,argsListLen; arg(i)%group=0; enddo
    endif

    if (argsCount.EQ.0) then
        if (checkDependencies(true).NE.0) then; rcode=-202; lastError=rcode; return; endif
    endif

    alreadyParsed=true
    allocate (holdArguments(1:argsCount))
    allocate ( argumentUsed(1:argsCount)); argumentUsed=false

    do i = 1,argsCount
        ! store all arguments
        call getarg(i,tempStr); holdArguments(i)=trim(adjustl(tempStr))

        ! if found -h/--help, show help
        if (argumentConforms(holdArguments(i)%get(),arg(0)).GT.0) then; void=showHelp(); stop; endif
    enddo

    ! get matching arguments
    do i = 1,argsCount
        do j = 1,argsListLen
            fnd=argumentConforms(holdArguments(i)%get(),arg(j))
            if (fnd.GT.0) then
                arg(j)%found=true; arg(j)%position=i; arg(j)%synonym=fnd; arg(j)%value=true
                argumentUsed(i)=true
            endif
        enddo
    enddo

    ! check on "expect" arguments. stop if "expect" argument is used as "unexpect"
    do j = 1,argsListLen
        if (arg(j)%expect.AND.arg(j)%found) then
            i=arg(j)%position+1
            if ( (i.GT.argsCount) .OR. (argumentUsed(i)) ) then
                write (iounit,99) arg(j)%str(arg(j)%synonym)%get()
                void=showHelp(); stop
            endif
            argumentUsed(i)=true

            arg(j)%expectValue=holdArguments(i)
        endif
    enddo

    ! append all non-ruled required arguments to the rule list
    void=appendSingles()

    ! check dependencies
    void=checkDependencies(true)

    ! collect unexpected arguments
    unexpectArgsCount=0
    do i = 1,argsCount
        if (.NOT.argumentUsed(i)) unexpectArgsCount=unexpectArgsCount+1
    enddo

    allocate (holdUnexpectArguments(unexpectArgsCount))

    j=0 ! count all unexpected arguments
    do i = 1,argsCount
        if (.NOT.argumentUsed(i)) then
            !write (*,*) holdArguments(i)%get()
            j=j+1; holdUnexpectArguments(j)=holdArguments(i)
        endif
    enddo

    ! show all unexpected arguments
    if (unexpectArgsCount.GT.0) then
        do j = 1,unexpectArgsCount
            write (iounit,100) holdUnexpectArguments(j)%get()
        enddo
        if (donlyDeclared) then; void=showHelp(); stop; endif
    endif

    ! free array
    deallocate (argumentUsed)

    ! set values
    void=setValues()

 99 format (2X,'Expected value for argument',1X,A,1X,'missing.')
100 format (2X,'Unexpected argument found:',1X,A)

    rcode=0; return
    end function parseArgs

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function apArgumentFound_c(str,err) result(ret)
    implicit none
    character (len=*)   :: str
    integer*4, optional :: err
    integer*4           :: position


    ret=false
    if (.NOT.alreadyParsed) then; lastError=-304; if (present(err)) err=lastError; return; endif
    position=getArgumentIndex(trim(str))

    if ( (position.LT.0) ) then
        lastError=-307
        if (present(err)) err=lastError ! Argument not declared.
        ret=false
        return
    endif

    ret=arg(position)%found

    return
    end function apArgumentFound_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function apArgumentFound_i(position,err) result(ret)
    implicit none
    integer*4           :: position
    integer*4, optional :: err


    ret=false
    if (.NOT.alreadyParsed) then; lastError=-304; if (present(err)) err=lastError; return; endif
    if (present(err)) err=0

    if ( (position.LE.0).OR.(position.GT.argsListLen) ) then
        lastError=-306
        if (present(err)) err=lastError ! Incorrect position.
        return
    endif

    ret=arg(position)%found

    return
    end function apArgumentFound_i

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function setApplicationName(str) result(ret)
    implicit none
    character (len=*) :: str


    appNameSet=true; appName=trim(str)

    ret=0; return
    end function setApplicationName

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function setPreamble(str) result(ret)
    implicit none
    character (len=*) :: str


    preambleSet=true; preamble=trim(str)

    ret=0; return
    end function setPreamble

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function setComment(str) result(ret)
    implicit none
    character (len=*) :: str


    commentSet=true; comment=trim(str)

    ret=0; return
    end function setComment

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function showHelp() result(ret)
    implicit none
    integer*4 :: k


    write (iounit,100) appName%get()
    if (preamble%ln.NE.0) write (iounit,101) preamble%get()

    do k = 1,groupListLen
        void=outputGroupHelp(k)
    enddo

    void=outputGroupHelp( 0)
    void=outputGroupHelp(-1)

    if (comment%ln.NE.0) write (iounit,101) comment%get()

100 format (/'syntax:',1X,A)
101 format (A/)

    ret=0; return
    end function showHelp

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function getDefaultArgument(position,err) result(ret)
    implicit none
    character (len=argLen) :: ret
    integer*4, optional    :: position,err


    ret=repeat(' ',argLen)
    if (.NOT.alreadyParsed) then; lastError=-304; if (present(err)) err=lastError; return; endif

    if (present(err)) err=0

    if (unexpectArgsCount.GT.0) then
        ret=holdUnexpectArguments(1)%get()
        return
    endif

    if (present(position)) then
        if ( (position.LE.0).OR.(position.GT.unexpectArgsCount) ) then
            lastError=-305; if (present(err)) err=lastError ! Incorrect position.
            return
        endif

        ret=holdUnexpectArguments(position)%get()
    endif

    return
    end function getDefaultArgument

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function getArgumentValue_i(position,exist,err) result(ret)
    implicit none

    integer*4              :: position
    logical*1, optional    :: exist
    integer*4, optional    :: err
    character (len=argLen) :: ret


    ret=repeat(' ',argLen)
    if (.NOT.alreadyParsed) then; lastError=-304; if (present(err)) err=lastError; return; endif

    if (present(err)) err=0

    if ( (position.LE.0).OR.(position.GT.argsListLen) ) then
        lastError=-306
        if (present(err)) err=lastError ! Incorrect position.
        return
    endif

    ret=arg(position)%expectValue%get()
    if (present(exist)) exist=arg(position)%found

    return
    end function getArgumentValue_i

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function getArgumentValue_c(str,exist,err) result(ret)
    implicit none

    character (len=*)      :: str
    logical*1, optional    :: exist
    integer*4, optional    :: err
    character (len=argLen) :: ret

    integer*4              :: position


    ret=repeat(' ',argLen)

    if (.NOT.alreadyParsed) then; lastError=-304; if (present(err)) err=lastError; return; endif
    if (present(err)) err=0

    position=getArgumentIndex(trim(str))

    if ( (position.LT.0) ) then
        lastError=-307
        if (present(err)) err=lastError ! Argument not declared.
        if (present(exist)) exist=false
        return
    endif

    ret=arg(position)%expectValue%get()
    if (present(exist)) exist=arg(position)%found

    return
    end function getArgumentValue_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function setGroupMembership(str,group) result(ret)
    implicit none

    character (len=*) :: str
    integer*4         :: group,position


    ret=0
    if ( (group.LT.0) .OR. (group.GT.groupListLen) ) then; ret=-9; lastError=ret; return; endif

    position=getArgumentIndex(trim(str))

    if ( (position.LT.0) ) then
        ret=-307; lastError=ret; return
    endif

    arg(position)%group=group

    return
    end function setGroupMembership

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine apFinalize
    implicit none

    integer*4 :: i


    alreadyParsed=false
    if ( allocated(holdArguments) ) then
        !do i = 1,UBound(holdArguments,1)
            call holdArguments%del()
        !enddo
        deallocate(holdArguments)
    endif

    if ( allocated(holdUnexpectArguments) ) then
        !do i = 1,UBound(holdUnexpectArguments,1)
            call holdUnexpectArguments%del()
        !enddo
        deallocate(holdUnexpectArguments)
    endif

    iounit=6; lastError=0
    argsListStart=1; outputArgumentLen=1
    argsListLen=0; argsCount=0; groupListLen=0; unexpectArgsCount=0
    alreadyParsed=false; appNameSet=false; preambleSet=false; commentSet=false

    return
    end subroutine apFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine apSetIOunit(iunt)
    implicit none

    integer*4, intent(in) :: iunt
    integer*4             :: unt


    unt=iunt
    if (unt.EQ.5) unt=6 !stdin  (forbidden)
    if (unt.EQ.0) unt=6 !stdout (for compatibility)
    iounit=unt

    return
    end subroutine apSetIOunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine apSetERRunit(iunt)
    implicit none

    integer*4, intent(in) :: iunt
    integer*4             :: unt


    unt=iunt
    if (unt.EQ.5) unt=6 !stdin  (forbidden)
    if (unt.EQ.0) unt=6 !stdout (for compatibility)
    errunit=unt

    return
    end subroutine apSetERRunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function ap_errDescription(err) result(ret)
    implicit none

    integer*4 :: err


    select case(err)
        case(   0);   write (errunit,100) err,'No errors occured.'
        case(  -1);   write (errunit,100) err,'The size of argument list exceeded.'
        case(  -2);   write (errunit,100) err,'The length of "description" is too big.'
        case(  -3);   write (errunit,100) err,'The length of "expectation value" is too big.'
        case(  -4);   write (errunit,100) err,'If contain "expect" and "default", must contain "expectation value".' !" hello sublime
        case(  -5);   write (errunit,100) err,'Too many synonyms for one argument.'
        case(  -6);   write (errunit,100) err,'The length for one of the synonyms is too big.'
        case(  -7);   write (errunit,100) err,'Argument has already been declared.'
        case(  -8);   write (errunit,100) err,'Can not contain "default" and "required" at once.' !" hello sublime
        case(  -9);   write (errunit,100) err,'Group ID does not exist.'
        case( -10);   write (errunit,100) err,'If contain "expectation value", must contain "expect".'
        case( -11);   write (errunit,100) err,'If contain "expectation value", must contain "default".' !" hello sublime
        case( -12);   write (errunit,100) err,'Unable to link variable with the argument. Already linked.'
        case( -13);   write (errunit,100) err,'Unable to link variable with the argument. List size exceeded.'
        case(-101);   write (errunit,100) err,'The size of the group list exceeded.'
        case(-102);   write (errunit,100) err,'The length of the "groupname" is too big.'
        case(-201);   write (errunit,100) err,'Arguments have already been parsed.'
        case(-202);   write (errunit,100) err,'One or more required arguments missing.'
        case(-301);   write (errunit,100) err,'The length of the "application name" is too big.'
        case(-302);   write (errunit,100) err,'The length of "preamble" is too big.'
        case(-303);   write (errunit,100) err,'The length of "comment" is too big.'
        case(-304);   write (errunit,100) err,'Arguments have not been parsed yet.'
        case(-305);   write (errunit,100) err,'Incorrect ID in the unexpected arguments list.'
        case(-306);   write (errunit,100) err,'Incorrect ID for recieved argument.'
        case(-307);   write (errunit,100) err,'Argument has not been declared.'
        case(-308);   write (errunit,100) err,'Argument has not been declared.'
        case(-309);   write (errunit,100) err,'Incorrect parameter format (Type mismatch).'
        case default; write (errunit,100) err,'Uknown error ID. Internal error.'
    end select

100 format (1X,'Error',1X,i3,1X,A)

    ret=0; return
    end function ap_errDescription

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function ap_popLastError() result(ret)
    implicit none


    void=ap_errDescription(lastError)

    ret=lastError; return
    end function ap_popLastError

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function apGetCommandLine() result(ret)
    implicit none

    type(uch)           :: ret
    character (len=512) :: tempStr
    integer*4           :: i


    ret=appName; i=1
    do
        tempStr=repeat(' ',len(tempStr))
        call getarg(i,tempStr)
        ret=ret%get()//' '//trim(tempStr)
        if (len_trim(tempStr).EQ.0) exit
        i=i+1
    enddo

    return
    end function apGetCommandLine

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function appendRule(arr,ruletype) result(rcode)
    implicit none

    character (len=*), intent(in) :: arr(:),ruletype
    integer*4, allocatable        :: iarr(:)
    integer*4                     :: ub,k,j


    rcode=0; ub=Ubound(arr,1)
    if (ruleListLen+1.GT.maxRulesListLen) then
        rcode=-1; return
    endif

    if ( .NOT. (ruletype .in. ['xor','and','single']) ) then
        rcode=-2; return
    endif

    do k = 1,ub
        !write (*,*) 'Calling from appendRule (preliminary)'
        j=getArgumentIndex(arr(k))

        if (j.EQ.-1) then
            rcode=-3; return
        endif
    enddo

    allocate (iarr(ub)); iarr=0

    do k = 1,ub
        !write (*,*) 'Calling from appendRule (final)'
        iarr(k)=getArgumentIndex(arr(k))
    enddo

    ruleListLen=ruleListLen+1
    ruleSet(ruleListLen)%argsList=iarr
    ruleset(ruleListLen)%type=ruletype

    deallocate (iarr)

    return
    end function appendRule

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function appendSingles() result(rcode)
    implicit none

    integer*4 :: k,i,j
    logical*1 :: notUsed


    rcode=0
    do k = 1,argsListLen
        if (arg(k)%required) then
            notUsed=true
            do i = 1,ruleListLen
            do j = 1,ruleSet(i)%argsList%ln
                if (k.EQ.ruleSet(i)%argsList%v(j)) then
                    notUsed=false; exit
                endif
            enddo
            enddo
            if (notUsed) rcode=appendRule([arg(k)%str(1)%get()],'single')
        endif
    enddo

    return
    end function appendSingles

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine clearRules
    implicit none

    integer*4 :: k


    do k = 1,ruleListLen
        call ruleSet(k)%argsList%del()
        ruleSet(k)%type=''
    enddo
    ruleListLen=0

    return
    end subroutine clearRules

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function setValues() result(rcode)
    implicit none

    integer*4 :: err,k,j
    type(uch) :: ustr


    rcode=0
    do k = 1,argsListLen
        j=arg(k)%var; if (j.EQ.0) cycle

        err=0

        if (.NOT.arg(k)%expect) then
            select case (variableSet(j)%vlen)
                case (8); variableSet(j)%pntl8=logical(arg(k)%value,kind=8)
                case (4); variableSet(j)%pntl4=logical(arg(k)%value,kind=4)
                case (2); variableSet(j)%pntl2=logical(arg(k)%value,kind=2)
                case (1); variableSet(j)%pntl1=logical(arg(k)%value,kind=1)
            end select
        else
            ustr=arg(k)%expectValue
            select case ( variableSet(j)%kind )
                case(0:2)  ; void=modify( variableSet(j)%address,tpRealByStr( ustr%get(),stat=err ) )
                case(3:6)  ; void=modify( variableSet(j)%address, tpIntByStr( ustr%get(),stat=err ) )
                case(11:12); void=modify( variableSet(j)%address,  tpDeQuote( ustr%get(),stat=err ) )
                !write (*,*) arg(k)%str(1)%get(),err,ustr%get()
            end select
        endif
    enddo

    if (err.NE.0) then; rcode=-309; lastError=rcode; endif

    return
    end function setValues

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function outputGroupHelp(gnumber) result(ret)
    implicit none

    logical*1 :: belong
    integer*4 :: i,gnumber,ln


    belong=false
    do i = 0,argsListLen
        if (arg(i)%group.NE.gnumber) cycle

        if (gnumber.GE.0) then
            if (.NOT.belong) then
                write (iounit,101) groupSet(gnumber)%get()
            endif
        endif
        belong=true

        ln=len_trim(collectSynonyms(arg(i)))
        if (arg(i)%required) then
            write (iounit,102) collectSynonyms(arg(i)),arg(i)%description%get(),' (required)'
        else
            write (iounit,102) collectSynonyms(arg(i)),arg(i)%description%get(),''
        endif
    enddo
    if (belong) write (iounit,103)

101 format (A,':')
102 format (7X,A<ln>,<outputArgumentLen+3-ln>X,A<arg(i)%description%ln>,A)
103 format ( )

    ret=0; return
    end function outputGroupHelp

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function addHelpArgument() result (rcode)
    implicit none

    type(cArgument)         :: recarg
    integer*4               :: err


    recarg%csensitive=true
    recarg%default=false
    recarg%required=false
    recarg%expect=false
    recarg%str(1)='-h'
    recarg%str(2)='--help'
    recarg%description='show help message'
    recarg%group=-1
    recarg%synonymsListLen=2
    recarg%synonym=0
    recarg%position=0
    recarg%expectValue=''
    recarg%found=false
    recarg%value=false
    groupSet(-1)=''
    groupSet(0)='Others'

    err=checkExistence(recarg)

    if (err.EQ.0) then
        arg(0)=recarg
        argsListStart=0
    endif

    rcode=0; return
    end function addHelpArgument

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function splitArgument(str,arg) result(ret)
    implicit none

    character (len=*) :: str
    type(cArgument)   :: arg

    type(uch)         :: split(maxSynonymListLen)
    integer*4         :: i,sta,sto,splCount


    splCount=0; sta=1; sto=len_trim(str)
    do
        if (splCount.GE.maxSynonymListLen) then; ret=-1; return; endif

        i=index(str(sta:sto),splitSymbol)
        if (i.EQ.0) then
            if (sto-sta+1.GT.argLen) then; ret=-2; return; endif
            splCount=splCount+1; split(splCount)=str(sta:sto)
            exit
        endif

        if (i-1.GT.argLen) then; ret=-2; return; endif
        splCount=splCount+1; split(splCount)=str(sta:sta+i-2)
        sta=sta+i
    enddo

    do i = 1,splCount
        arg%str(i)=split(i)
    enddo
    arg%synonymsListLen=splCount

    ret=splCount; return
    end function splitArgument

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function checkExistence(recarg) result(ret)
    implicit none

    type(cArgument) :: recarg
    integer*4       :: i,k,l


    do i = 1,argsListLen
    do k = 1,recarg%synonymsListLen
    do l = 1,arg(i)%synonymsListLen
        if ( (len_trim( recarg%str(k)%get() ).GT.0) .AND.&
             (len_trim( arg(i)%str(l)%get() ).GT.0) .AND.&
             (trim( recarg%str(k)%get() ).EQ.trim( arg(i)%str(l)%get() )) ) then

            ret=-1; return

        endif
    enddo
    enddo
    enddo

    ret=0; return
    end function checkExistence

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function checkDependencies(stopOnMissing) result(rcode)
    implicit none

    logical*1, optional    :: stopOnMissing
    logical*1              :: dstopOnMissing,missing,finally
    integer*4              :: mx,ln,i,j,k
    logical*1, allocatable :: larr(:)


    rcode=0
    dstopOnMissing=false; if (present(stopOnMissing)) dstopOnMissing=stopOnMissing

    mx=1
    do i = 1,ruleListLen
        mx=max(mx,ruleSet(i)%argsList%ln)
    enddo

    allocate (larr(mx))

    ! xor,and rules.
    missing=false
    do i = 1,ruleListLen
        ln=ruleSet(i)%argsList%ln; larr=false

        do j = 1,ln
            k=ruleSet(i)%argsList%v(j)
            larr(j)=arg(k)%found
        enddo

        finally=larr(1)
        select case(ruleSet(i)%type%get())
            case('xor')
                do j = 2,ln
                    finally=finally .XOR. larr(j)
                enddo

            case('and')
                do j = 2,ln
                    finally=finally .AND. larr(j)
                enddo

            case('single')
                continue

        end select

        if (.NOT.finally) then
            missing=true
            write (iounit,100) ruleSet(i)%type%get()
            do j = 1,ln
                k=ruleSet(i)%argsList%v(j)
                write (iounit,101) trim(collectSynonyms(arg(k)))
            enddo
        endif

    enddo

    deallocate (larr)

    if (missing) then
        void=showHelp(); rcode=-1
        if (dstopOnMissing) stop
    endif

100 format ('Rule type: ',A)
101 format (4X,'Expected argument [',A,'] missing.')

    return
    end function checkDependencies

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function getApplicationName() result(ret)
    implicit none

    character (len=appNameLen) :: tempStr
    integer*4                  :: sta,sto


    tempStr=repeat(' ',appNameLen); call getarg(0,tempStr)
    sta=index(tempStr,osSeparator,true)+1; sto=len_trim(tempStr)

    appName=tempStr(sta:sto)

    ret=0; return
    end function getApplicationName

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function collectSynonyms(recarg) result(string)
    implicit none

    character (len=argLen*maxSynonymListLen) :: string
    type(cArgument)                          :: recarg
    integer*4                                :: k


    string=repeat(' ',len(string)); string=recarg%str(1)%get()
    do k = 2,recarg%synonymsListLen
        string=trim(string)//concatSynonymsSymbol//recarg%str(k)%get()
    enddo

    return
    end function collectSynonyms

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function argumentConforms(str,recarg) result(ret)
    implicit none

    character (len=*)        :: str
    type(cArgument)          :: recarg
    character (len=len(str)) :: tempStr
    integer*4                :: k


    tempStr=repeat(' ',len(str)); tempStr=trim(str)
    if (.NOT.recarg%csensitive) tempStr=tpLowerCase(tempStr)
    do k = 1,recarg%synonymsListLen
        !write (*,*) trim(tempStr),' =====> ',recarg%str(k)%get()
        if ( trim(tempStr).EQ.recarg%str(k)%get() ) then
            ret=k; return
        endif
    enddo

    ret=0; return
    end function argumentConforms

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function getArgumentIndex(str) result(ret)
    implicit none

    character (len=*)        :: str
    character (len=len(str)) :: tempStr
    integer*4                :: k,i


    ! write (*,*) 'All arguments',argsListLen
    ! do i = 1,argsListLen
    !     do k = 1,arg(i)%synonymsListLen
    !         write (*,*) '#'//arg(i)%str(k)%get()//'#'
    !     enddo
    ! enddo
    ! write (*,*) 'Proceeding...'

    do i = 1,argsListLen
        tempStr=repeat(' ',len(str)); tempStr=trim(str)
        if (.NOT.arg(i)%csensitive) tempStr=tpLowerCase(tempStr)
        do k = 1,arg(i)%synonymsListLen
            !write (*,*) '#'//trim(tempStr)//'# #'//arg(i)%str(k)%get()//'#'
            if ( trim(tempStr).EQ.arg(i)%str(k)%get() ) then
                ret=i; return
            endif
        enddo
    enddo

    ret=-1; return
    end function getArgumentIndex

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_r16(var) result(rcode)
    implicit none

    integer*4       :: vkind=0,vlen=16
    real*16, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntr16  =>var

    return
    end function shareVariable_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_r8(var) result(rcode)
    implicit none

    integer*4      :: vkind=1,vlen=8
    real*8, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntr8   =>var

    return
    end function shareVariable_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_r4(var) result(rcode)
    implicit none

    integer*4      :: vkind=2,vlen=4
    real*4, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntr4   =>var

    return
    end function shareVariable_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_i8(var) result(rcode)
    implicit none

    integer*4         :: vkind=3,vlen=8
    integer*8, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pnti8   =>var

    return
    end function shareVariable_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_i4(var) result(rcode)
    implicit none

    integer*4         :: vkind=4,vlen=4
    integer*4, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pnti4   =>var

    return
    end function shareVariable_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_i2(var) result(rcode)
    implicit none

    integer*4         :: vkind=5,vlen=2
    integer*2, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pnti2   =>var

    return
    end function shareVariable_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_i1(var) result(rcode)
    implicit none

    integer*4         :: vkind=6,vlen=1
    integer*1, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pnti1   =>var

    return
    end function shareVariable_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_l8(var) result(rcode)
    implicit none

    integer*4         :: vkind=7,vlen=8
    logical*8, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntl8   =>var

    return
    end function shareVariable_l8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_l4(var) result(rcode)
    implicit none

    integer*4         :: vkind=8,vlen=4
    logical*4, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntl4   =>var

    return
    end function shareVariable_l4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_l2(var) result(rcode)
    implicit none

    integer*4         :: vkind=9,vlen=2
    logical*2, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntl2   =>var

    return
    end function shareVariable_l2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_l1(var) result(rcode)
    implicit none

    integer*4         :: vkind=9,vlen=1
    logical*1, target :: var


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =vlen
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntl1   =>var

    return
    end function shareVariable_l1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_c(var) result(rcode)
    implicit none

    character (len=*), target :: var
    integer*4                 :: vkind=11


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =len(var)
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntch   =>var

    return
    end function shareVariable_c

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function shareVariable_uc(var) result(rcode)
    implicit none

    type(uch), target :: var
    integer*4         :: vkind=12


    if (find(variableSet(1:varListLen)%address,int(loc(var),8)).GT.0) then
        write (errunit,*) 'Variable at address',loc(var),' is already shared.'
        rcode=-1; return
    endif

    if (varListLen+1.GT.maxArgsListLen) then
        write (errunit,*) 'Max number of shared variables exceeded.'
        rcode=-2; return
    endif

    varListLen=varListLen+1; rcode=varListLen
    variableSet(varListLen)%address =loc(var)
    variableSet(varListLen)%position=varListLen
    variableSet(varListLen)%vlen    =0
    variableSet(varListLen)%kind    =vkind
    variableSet(varListLen)%pntuch  =>var

    return
    end function shareVariable_uc

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
    integer*8          :: address
    integer*4          :: i


    i=find(variableSet(:)%address,address)
    if (i.GT.0) then
        if ((variableSet(i)%kind.NE.1).AND.(variableSet(i)%kind.NE.1).AND.(variableSet(i)%kind.NE.2)) then
            ret=-1
            return
        endif
        select case (variableSet(i)%vlen)
            case (16); variableSet(i)%pntr16=real(val,kind=16)
            case ( 8); variableSet(i)%pntr8 =real(val,kind=8)
            case ( 4); variableSet(i)%pntr4 =real(val,kind=4)
        end select
    else; ret=-2; return
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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module argsParser
