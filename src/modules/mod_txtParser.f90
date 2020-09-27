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

    module txtParser

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHANGELOG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!
!   *** 2019.07.27 (v3.700):
!     > added tpIsAnyInList(flist, slist).
!       If any in the first list present in the second list, returns true.
!     > added tpGetSplit([str | uch], delim) + overrides .in. operator.
!       Returns array of stings splited from str.
!       !!!IMPORTANT, returns with trailing spaces,
!       needs trim() for every array element: drawback of Fortran LHS procedure.
!
!   *** 2019.06.01 (v3.621):
!     > changed behaviour of tpSplit. Now if separator is absent - return the whole string.
!
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob    , only: assignment(=)
    use glob    , only: r8kind,void,voidl,true,false,NaN,uch
    use fcontrol, only: fcNewID,fcBanID,fcUnBanID,fcNullID

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: tpVersion='3.700'
    character (len=*), parameter :: tpDate   ='2019.07.27'
    character (len=*), parameter :: tpAuthor ='Anton B. Zakharov'

    integer*4, parameter         :: maxStrLen=1024,maxCommentDefLen=5

    integer*4, parameter         :: maxReplaceStrLen=4*maxStrLen

    character (len=1), parameter :: spaceChar=char(32),tabChar=char(9),nullchar=char(0)

    character (len=*), parameter :: sigset='+-',&
                                    expset='edqEDQ',&
                                    sepset='.,',&
                                    odnset='0123456789edqEDQ+-.,',&
                                    numset='0123456789abcdefqABCDEFQ+-.,',&
                                    hexset='0123456789abcdefABCDEF',&
                                    decset='0123456789',&
                                    octset='01234547',&
                                    binset='01',&
                                    abcset='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPRSTUVWXYZ',&
                                    bancomset=':,./_{}+-'//"'"//'"'//spaceChar//tabChar//nullchar,&
                                    quoset='"'//"'"

#   if (__OS==1)
        character (len=*), parameter :: newline=char(13)//char(10)
#   else
        character (len=*), parameter :: newline=char(10)
#   endif


!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type(uch), allocatable :: tpSplitHold(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=maxStrLen)        :: tpholdCurrentString
    character (len=maxCommentDefLen) :: commentMarkup='!    '

    integer*8 :: tpSplitAdress
    integer*4 :: iount,iostat,uStrLen,tpPointer
    integer*4 :: tpSplitLen
    logical*1 :: tpEnabled=false,fAttached

    logical*1 :: tpCaseSens,tpReplTabs,tpIgnoSpac,tpReduSpac,&
                 tpIgnoTabs,tpAlloQuot,tpIgnoComm

    !tpCaseSens if true, does not lowerCase strings during search.  Default: true.
    !tpIgnoTabs if true, remove tabs.                               Default: false.
    !tpReplTabs if true, replace tabs by spaces.                    Default: false.
    !tpIgnoSpac if true, remove spaces.                             Default: false.
    !tpReduSpac if true, reduce spaces.                             Default: false.
    !tpAlloQuot if true, allow search in strings.                   Default: true.
    !tpIgnoComm if true, remove chars after comment mark.           Default: false.

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    interface tpAdjustc
        module procedure tpDefCentrString,tpNomCentrString
    end interface tpAdjustc

    interface tpFill
        module procedure tpFillbc_blank,tpFillbl_blank,tpFillbc_any,tpFillbl_any
    end interface tpFill

    interface tpLocate
        module procedure tpLocateOne,tpLocateOneOf
    end interface tpLocate

    interface tpIsStrInList
        module procedure tpIsStrInList_ch,tpIsStrInList_uc
    end interface tpIsStrInList

    interface tpIsIn
        module procedure tpIsIn_ch_ch,tpIsIn_ch_uc,tpIsIn_uc_ch,tpIsIn_uc_uc
    end interface tpIsIn

    interface tpSplit
        module procedure tpSplit_uc,tpSplit_ch
    end interface tpSplit

    interface tpRetSplit
        module procedure tpRetSplit_ch,tpRetSplit_uc
    end interface tpRetSplit

    interface operator (.in.)
        module procedure tpIsIn_ch_ch,tpIsIn_ch_uc,tpIsIn_uc_ch,tpIsIn_uc_uc,&
                         tpIsStrInList_ch,tpIsStrInList_uc,tpIsAnyInList
    end interface

    interface tpGetSplit
        module procedure tpGetSplit_ch, tpGetSplit_uc
    end interface tpGetSplit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    public

    private   :: readstr,newline
    private   :: tpFillbc_blank,tpFillbl_blank,tpFillbc_any,tpFillbl_any,&
                 tpNomCentrString,tpDefCentrString,defIntBySystem,       &
                 tpIsStrInList_ch,tpIsStrInList_uc,tpLocateOne,          &
                 tpLocateOneOf,tpIsIn_ch_ch,tpIsIn_ch_uc,tpIsIn_uc_ch,   &
                 tpIsIn_uc_uc,tpSplit_uc,tpSplit_ch,tpRetSplit_ch,       &
                 tpRetSplit_uc,tpIsAnyInList,tpGetSplit_ch,tpGetSplit_uc

    private   :: fcNewID,fcBanID,fcUnBanID,fcNullID
    private   :: uch,void,voidl,true,false

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~ General Functions ~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine tpEnable(strlen)
    implicit none
    integer*4, optional :: strlen


    uStrLen=maxStrLen; fAttached=false
    iount=fcNewID() !; void=fcBanID(iount,true)

    if (present(strlen)) then
        if (strlen.LE.maxStrLen) then
            uStrLen=strlen
        else
            !write (*,100) maxStrLen
            !100 format (1X,'Max possible length is',1X,i4)
        endif
    endif

    void=setDefaultMask(); tpEnabled=true

    return
    end subroutine tpEnable

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine tpDisable
    implicit none


    if (fAttached) voidl=tpRelease()
    uStrLen=maxStrLen; tpEnabled=false

    return
    end subroutine tpDisable

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function setDefaultMask() result(ret)
    implicit none


    tpCaseSens=true;  tpReplTabs=false; tpReduSpac=false
    tpIgnoSpac=false; tpIgnoTabs=false; tpIgnoComm=false
    tpAlloQuot=true
    commentMarkup=tpFill(commentMarkup); commentMarkup(1:1)='!'

    ret=0; return
    end function setDefaultMask

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function tpSetCommentMark(str) result(ret)
    implicit none

    character (len=*)        :: str
    character (len=len(str)) :: ustr
    integer*4                :: i


    ret=0; ustr=adjustl(trim(str))
    if (len_trim(ustr).GT.len(commentMarkup))                  then; ret=-1; return; endif
    if (len(ustr).EQ.0)                                        then; ret=-2; return; endif
    if (tpOneof(ustr(1:1),banComset//tpDigits()//tpLetters())) then; ret=-3; return; endif

    commentMarkup=tpFill(commentMarkup)
    do i = 1,len_trim(ustr); commentMarkup(i:i)=ustr(i:i); enddo

    return
    end function tpSetCommentMark

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~ Atomic File Functions ~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpJoin(path) result(ret)
    implicit none
    character (len=*) :: path


    ret=tpEnabled; if (.NOT.tpEnabled) return

    open(iount,file=trim(path), action='read', status='old', iostat=iostat)
    if (iostat.NE.0) ret=false

    fAttached=ret; tpPointer=0

    return
    end function tpJoin

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpRelease() result(ret)
    implicit none


    fAttached=false
    !void=fcUnBanID(iount)
    void=fcNullID(iount)
    !fAttached=false; write (*,*) 'uban',fcUnBanID(iount); write (*,*) 'null',fcNullID(iount)

    ret=true; return
    end function tpRelease

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function readStr(again) result(ret)
    implicit none

    logical*1, optional :: again


    if (present(again)) then
        if (again) voidl=tpBackSpace()
    endif

    ret=true; tpholdCurrentString=tpFill(len(tpholdCurrentString))
    read (iount,100,iostat=iostat) tpholdCurrentString; tpPointer=tpPointer+1
100 format (A<uStrLen>)

    if (iostat.NE.0) ret=false

    return
    end function readStr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpRewind() result(ret)
    implicit none


    rewind(iount); tpPointer=0

    ret=true; return
    end function tpRewind

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpBackSpace() result(ret)
    implicit none
    integer*4 :: err


    backspace(iount,iostat=err)
    if (err.EQ.0) then
        tpPointer=tpPointer-1
        ret=true
    else
        ret=false
    endif

    return
    end function tpBackSpace

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function tpSetFilePointer(pntr) result(ret)
    implicit none
    integer*4 :: pntr,i


    ret=-1; if (.NOT.tpEnabled) return
    ret=-1; if (.NOT.fAttached) return

    ret=0
    if (tpPointer.GT.pntr-1) then
        do i = 1,tpPointer-pntr+1
            if (.NOT.tpBackSpace()) exit
        enddo
    else
        do i = tpPointer,pntr-2
            if (eof(iount)) then; ret=-1; return; endif
            voidl=readStr()
        enddo
    endif

    return
    end function tpSetFilePointer

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpGetLineByPointer(pntr) result(ret)
    implicit none

    character (len=uStrLen) :: ret
    integer*4               :: pntr


    ret=tpFill(ret)
    if (.NOT.tpEnabled) return
    if (.NOT.fAttached) return

    void=tpSetFilePointer(pntr)
    if (void.EQ.0) then
        if (readStr()) ret=tpholdCurrentString
    endif

    return
    end function tpGetLineByPointer

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpReadAdvance() result(ret)
    implicit none
    character (len=uStrLen) :: ret


    voidl=readstr(); ret=tpFill(len(ret)); ret=tpHoldCurrentString
    return
    end function tpReadAdvance

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~ Complex File Functions ~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function tpLocateOne(substr,doRewind,atTheStart) result(ret)
    implicit none

    character (len=*)  , intent(in) :: substr
    logical*1, optional, intent(in) :: doRewind,atTheStart
    character (len=len(substr))     :: usubstr

    logical*1 :: fnd,dts
    integer*4 :: pos


    ret=0
    if (.NOT.tpEnabled) ret=-1; if (.NOT.tpEnabled) return
    if (.NOT.fAttached) ret=-2; if (.NOT.fAttached) return

    usubstr=tpFill(usubstr); usubstr=tpPrepareString(substr)

    dts=false; if (present(atTheStart)) dts=atTheStart
    if (present(doRewind)) then; if (doRewind) voidl=tpRewind(); endif

    do
        if (    eof(iount)) then; ret=-3; return; endif
        if (.NOT.readstr()) then; ret=-4; return; endif

        tpholdCurrentString=tpPrepareString(tpholdCurrentString)

        !write (*,*) '$'//trim(tpholdCurrentString)//'$'
        pos=index(tpholdCurrentString,usubstr); fnd=pos.GT.0

        !write (*,*) 'Locate: '//trim(usubstr)//' in '//trim(tpHoldCurrentString)

        if (fnd) then
            if (.NOT.tpAlloQuot) then
                if (tpQuoted(tpholdCurrentString,pos)) cycle
            endif
            if (dts.AND.(pos.NE.1)) cycle
            !write (*,*) 'Found.', tpPointer

            ! 4.07.17 change
            voidl=tpBackSpace()
            exit
        endif

    enddo

    return
    end function tpLocateOne

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function tpLocateOneOf(arr,doRewind,atTheStart) result(ret)
    implicit none

    character (len=*)  , intent(in) :: arr(:)
    logical*1, optional, intent(in) :: doRewind,atTheStart

    character (len=len(arr))        :: uarr(1:UBound(arr,1))

    logical*1 :: fnd,dts
    integer*4 :: pos,i,alen


    ret=tpEnabled; if (.NOT.tpEnabled) return
    ret=fAttached; if (.NOT.fAttached) return

    dts=false; if (present(atTheStart)) dts=atTheStart

    if (present(doRewind)) then
        if (doRewind) voidl=tpRewind()
    endif

    alen=UBound(arr,1)

    do i = 1,alen; uarr(i)=tpFill(uarr(i));                enddo
    do i = 1,alen; uarr(i)=trim(arr(i));                   enddo
    do i = 1,alen; uarr(i)=trim(tpPrepareString(uarr(i))); enddo

    do
        if (    eof(iount)) then; ret=-1; return; endif
        if (.NOT.readstr()) then; ret=-1; return; endif

        tpholdCurrentString=tpPrepareString(tpholdCurrentString)

        do i = 1,alen
            pos=index(trim(tpholdCurrentString),trim(uarr(i))); fnd=pos.GT.0

            if (fnd) then
                if (.NOT.tpAlloQuot) then
                    if (tpQuoted(tpholdCurrentString,pos)) cycle
                endif
                if (dts.AND.(pos.NE.1)) cycle

                ! 4.07.17 change
                voidl=tpBackSpace()
                ret=i; return
            endif
        enddo
    enddo

    return
    end function tpLocateOneOf

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer*4 function tpStringOccurance(substr,doRewind) result(ret)
    implicit none

    character (len=*)   :: substr
    logical*1, optional :: doRewind
    logical*1           :: fnd


    ret=0
    if (.NOT.tpEnabled) return
    if (.NOT.fAttached) return

    if (present(doRewind)) then
        if (doRewind) voidl=tpRewind()
    endif

    do
        fnd=tpLocate(substr); fnd=readstr()
        if (fnd) ret=ret+1

        if (.NOT.fnd) exit
    enddo

    return
    end function tpStringOccurance

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpReplaceInFile(file,src,trg) result(ret)
    implicit none

    character (len=*), intent(in) :: file,trg,src
    logical*1                     :: fnd
    integer*4                     :: isrc,itrg
    character (len=maxStrLen)     :: tstr


    isrc=fcNewID(); open (isrc,file=trim(file))
    itrg=fcNewID(); open (itrg,status='scratch')
    do
        if (eof(isrc)) exit

        tstr=tpFill(len(tstr))
        read  (isrc,100) tstr
        fnd=tpIndex(trim(adjustl(tpLowerCase(tstr))),trim(tpLowerCase(src)))
        if (fnd) tstr=tpReplace(trim(tpLowerCase(tstr)),trim(tpLowerCase(src)),trim(trg))
        write (itrg,100) trim(tstr)
    enddo
    close (isrc)

    open (isrc,file=trim(file),status='replace')
    rewind(itrg)

    do
        if (eof(itrg)) exit
        tstr=tpFill(len(tstr))
        read  (itrg,100) tstr
        write (isrc,100) trim(tstr)
    enddo
100 format (A)

    close (isrc); close (itrg)

    void=fcNullID(isrc); void=fcNullID(itrg)

    ret=true; return
    end function tpReplaceInFile

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~ Atomic String Functions ~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpLowerCase(str,quo) result(ret)
    implicit none
    character (len=*), intent(in)   :: str
    character (len=len(str))        :: ret
    logical*1, optional, intent(in) :: quo
    logical*1                       :: dquo
    integer*4                       :: i,ascii


    dquo=false; if (present(quo)) dquo=quo

    ret=tpFill(ret)
    do i = 1,len_trim(str)
        ascii=iachar(str(i:i))
        if (dquo) then
            if ((ascii.GE.65).AND.(ascii.LE.90) .AND.(.NOT.(dquo.AND.tpQuoted(str,i))) ) then
                ret(i:i)=char(ascii+32)
            else
                ret(i:i)=char(ascii)
            endif
        else
            if ( (ascii.GE.65).AND.(ascii.LE.90) ) then
                ret(i:i)=char(ascii+32)
            else
                ret(i:i)=char(ascii)
            endif
        endif
    enddo

    return
    end function tpLowerCase

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpUpperCase(str,quo) result(ret)
    implicit none
    character (len=*), intent(in) :: str
    character (len=len(str))      :: ret
    logical*1, optional, intent(in) :: quo
    logical*1                       :: dquo
    integer*4                     :: i,ascii


    dquo=false; if (present(quo)) dquo=quo
    do i = 1,len(str)
        ascii=iachar(str(i:i))
        if (dquo) then
            if ((ascii.GE.97).AND.(ascii.LE.122) .AND.(.NOT.(dquo.AND.tpQuoted(str,i))) ) then
                ret(i:i)=char(ascii-32)
            else
                ret(i:i)=char(ascii)
            endif
        else
            if ((ascii.GE.97).AND.(ascii.LE.122) ) then
                ret(i:i)=char(ascii-32)
            else
                ret(i:i)=char(ascii)
            endif
        endif
    enddo

    return
    end function tpUpperCase

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpSwitchCase(str) result(ret)
    implicit none
    character (len=*), intent(in) :: str
    character (len=len(str))      :: ret
    integer*4                     :: i,ascii


    do i = 1,len(str)
        ascii=ichar(str(i:i))
        if ((ascii.GE.65).AND.(ascii.LE.90)) then
            ret(i:i)=char(ascii+32)
        elseif ((ascii.GE.97).AND.(ascii.LE.122)) then
            ret(i:i)=char(ascii-32)
        else
            ret(i:i)=char(ascii)
        endif
    enddo

    return
    end function tpSwitchCase

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpCapitalize(str) result(ret)
    implicit none
    character (len=*), intent(in) :: str
    character (len=len(str))      :: ret


    ret=str; ret(1:1)=tpUpperCase(ret(1:1))
    return
    end function tpCapitalize

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpAdjustl(str,spacer) result(ret)
    implicit none
    character (len=*), intent(in)           :: str
    character (len=:), allocatable          :: ret
    character (len=1), intent(in), optional :: spacer
    integer*4                               :: ln


    ln=len_trim(adjustl(str))
    if (present(spacer)) then
        ret=trim(adjustl(str))//tpFill(len(str)-ln,spacer)
    else
        ret=trim(adjustl(str))
    endif

    return
    end function tpAdjustl

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpAdjustr(str,spacer) result(ret)
    implicit none
    character (len=*), intent(in)           :: str
    character (len=len(str))                :: ret
    character (len=1), intent(in), optional :: spacer
    integer*4                               :: ln


    ln=len_trim(adjustl(str)); ret=adjustr(str)

    if (present(spacer)) then
        ret(1:len(str)-ln)=tpFill(len(str)-ln,spacer)
    endif

    return
    end function tpAdjustr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpFillbl_blank(ln) result(ret)
    implicit none
    integer*4, intent(in) :: ln
    character (len=ln)    :: ret


    ret=repeat(' ',ln); return
    end function tpFillbl_blank

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpFillbc_blank(str) result(ret)
    implicit none
    character (len=*), intent(in) :: str
    character (len=len(str))      :: ret


    ret=repeat(' ',len(str)); return
    end function tpFillbc_blank

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpFillbl_any(ln,spacer) result(ret)
    implicit none
    integer*4        , intent(in) :: ln
    character (len=*), intent(in) :: spacer
    character (len=ln)            :: ret

    integer*4                     :: splen


    splen=len(spacer)
    ret=repeat( spacer,int(ln/splen) )//spacer(1:mod(ln,splen) )

    return
    end function tpFillbl_any

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpFillbc_any(str,spacer) result(ret)
    implicit none
    character (len=*), intent(in) :: str
    character (len=*), intent(in) :: spacer
    character (len=len(str))      :: ret

    integer*4                     :: ln,splen


    splen=len(spacer); ln=len(str)
    ret=repeat( spacer,int(ln/splen) )//spacer(1:mod(ln,splen) )

    return
    end function tpFillbc_any

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpNomCentrString(str,spacer) result(ret) !tpAdjustc
    implicit none
    character (len=*)          , intent(in) :: str
    character (len=1), optional, intent(in) :: spacer

    character (len=len(str)) :: ret,tmp
    integer*4                :: ln,rst,lft,rght
    character (len=1)        :: defSpacer


    defSpacer=' '; if (present(spacer)) defSpacer=spacer
    if (defSpacer.EQ.' ') then
        tmp=tpFill(len(tmp)); tmp=trim(adjustl(str)); ln=len_trim(tmp)
    else
        tmp=str; ln=len(tmp)
    endif

    rst=len(str)-ln
    lft=rst/2; rght=rst/2

    lft=lft+ibits(rst,0,1)

    if (defSpacer.EQ.' ') then
        ret=tpFill(lft,defSpacer)//trim(tmp)//tpFill(rght,defSpacer)
    else
        ret=tpFill(lft,defSpacer)//tmp//tpFill(rght,defSpacer)
    endif

    return
    end function tpNomCentrString

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpDefCentrString(str,defLen,spacer) result(ret) !tpAdjustc
    implicit none
    integer*4, intent(in)                   :: defLen
    character (len=*)          , intent(in) :: str
    character (len=1), optional, intent(in) :: spacer

    character (len=len(str)) :: tmp
    character (len=defLen)   :: ret
    integer*4                :: ln,rst,lft,rght
    character (len=1)        :: defSpacer


    defSpacer=' '; if (present(spacer)) defSpacer=spacer
    if (defSpacer.EQ.' ') then
        tmp=tpFill(len(tmp)); tmp=trim(adjustl(str)); ln=len_trim(tmp)
    else
        tmp=str; ln=len(tmp)
    endif

    if (ln.GT.defLen) then
        ret=repeat('*',defLen)
        return
    endif

    rst=defLen-ln; lft=rst/2; rght=rst/2

    lft=lft+ibits(rst,0,1)

    if (defSpacer.EQ.' ') then
        ret=tpFill(lft,defSpacer)//trim(tmp)//tpFill(rght,defSpacer)
    else
        ret=tpFill(lft,defSpacer)//tmp//tpFill(rght,defSpacer)
    endif

    return
    end function tpDefCentrString

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer*4 function tpCount(str,sub,overlap) result(ret)
    implicit none
    character (len=*)  , intent(in) :: str,sub
    logical*1, optional, intent(in) :: overlap
    logical*1                       :: doverlap
    integer*4                       :: i,k


    ret=0; i=1
    doverlap=false; if (present(overlap)) doverlap=overlap
    if (doverlap) then
        do
            if (i.GT.len(str)) exit
            k=index(str(i:),sub); if (k.EQ.0) exit
            ret=ret+1
            i=i+k
        enddo
    else
        do
            if (i.GT.len(str)) exit
            k=index(str(i:),sub); if (k.EQ.0) exit
            ret=ret+1
            i=i+k+len(sub)-1
        enddo
    endif

    return
    end function tpCount

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpReverse(str) result(ret)
    implicit none
    character (len=*), intent(in) :: str
    character (len=len(str))      :: ret
    integer*4                     :: ln,i


    ln=len(str); ret=tpFill(ret)
    do i = 1,ln
        ret(ln-i+1:ln-i+1)=str(i:i)
    enddo

    return
    end function tpReverse

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~ Complex String Functions ~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpPrepareString(str,ignStr) result(ret)
    implicit none

    character (len=*)        :: str
    character (len=len(str)) :: ret
    logical*1, optional      :: ignStr
    logical*1                :: dignStr


    dignStr=false; if (present(ignStr)) dignStr=ignStr
    ret=tpFill(len(ret)); ret=str
    if (.NOT.tpCaseSens) ret=tpLowerCase(ret                ,quo=dignStr)
    if (     tpReplTabs) ret=tpReplace(ret,tabChar,spaceChar,quo=dignStr)
    if (     tpReduSpac) ret=trim( tpReduce(ret,spaceChar   ,quo=dignStr) )
    if (     tpIgnoSpac) ret=trim( tpRemove(ret,spaceChar   ,quo=dignStr) )
    if (     tpIgnoTabs) ret=trim( tpRemove(ret,tabChar     ,quo=dignStr) )
    if (     tpIgnoComm) ret=tpRemoveComment(ret            ,quo=dignStr)

    return
    end function tpPrepareString

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! be cyka very careful. if something goes wrong, it is due to this shit.
    pure integer*4 function tpIndex(str,sub,cnt,start,rev,endp,casens) result(ret)
    implicit none
    character (len=*)  , intent(in) :: str,sub
    integer*4, optional, intent(in) :: cnt,start
    logical*1, optional, intent(in) :: rev,endp,casens
    logical*1                       :: drev,dendp,dcasens
    integer*4                       :: pos,fnsh,nrml,ucnt,fnd,endpcor
    character (len=len(str))        :: ustr
    character (len=len(sub))        :: usub


    drev=false;    if (present(rev))    drev=rev
    dendp=false;   if (present(endp))   dendp=endp
    ucnt=1;        if (present(cnt))    ucnt=cnt
    dcasens=false; if (present(casens)) dcasens=casens

    if (.NOT.dcasens) then
        ustr=tpLowerCase(str); usub=tpLowerCase(sub)
    else
        ustr=str; usub=sub
    endif

    if (present(start)) then
        if (start.LE.0) then;         ret=0; return; endif
        if (start.GE.len(ustr)) then; ret=0; return; endif
    endif

    pos=1; fnsh=len(ustr)
    if (drev) then
        if (present(start)) fnsh=len(ustr)-start
    else
        if (present(start)) pos=start
    endif

    endpcor=0; if (dendp) endpcor=len_trim(sub)-1

    if (drev) then
        ustr=tpReverse(ustr)
        usub=tpReverse(usub)
    endif

    ! write (*,*) 'Looking for '//trim(usub)//' in '//trim(ustr)

    nrml=0
    do
        if (pos+len(usub)-1.GT.fnsh) then; ret=0; exit; endif
        fnd=index(ustr(pos:fnsh),usub)

        if (fnd.GT.0) then; nrml=nrml+1
        else;               ret=0; exit
        endif

        if (nrml.EQ.ucnt) then
            if (drev) then; ret=len(ustr)-(pos+fnd-1)-(len(usub)-2)+endpcor
            else;           ret=           pos+fnd-1               +endpcor
            endif
            return
        endif

        pos=pos+fnd
    enddo

    return
    end function tpIndex

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpStartsWith(str,sub) result(ret)
    implicit none

    character (len=*), intent(in) :: str,sub


    ret=tpIndex(str,sub).EQ.1

    return
    end function tpStartsWith

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpEndsWith(str,sub) result(ret)
    implicit none

    character (len=*), intent(in) :: str,sub


    ret=tpIndex(str,sub,rev=true,endp=true).EQ.len(str)

    return
    end function tpEndsWith

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpGetSplit_ch(str,delim) result(ret)
    implicit none

    character (len=*)               , intent(in)  :: str,delim
    character (len=:),   allocatable              :: ret(:)
    type(uch), allocatable                        :: arr(:)
    integer*4                                     :: k,ln,sta,sto


    ln=tpCount(str,delim,overlap=false)+1
    if (ln.EQ.1) then
        allocate(ret, source=[str])
    else
        allocate (arr(ln))
        sta=1
        do k = 1,ln-1
            sto=tpIndex(str,delim,cnt=k)
            arr(k)=str(sta:sto-1)

            sta=sto+len(delim)
        enddo
        arr(ln)=str(sta:)

        allocate(ret, source=[(arr(k)%get(), k=1,ln)])
    endif

    return
    end function tpGetSplit_ch

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpGetSplit_uc(str,delim) result(ret)
    implicit none

    type(uch)                       , intent(in)  :: str
    character(len=*)                , intent(in)  :: delim
    character (len=:),   allocatable              :: ret(:)
    type(uch), allocatable                        :: arr(:)
    integer*4                                     :: k,ln,sta,sto


    ln=tpCount(str%get(),delim,overlap=false)+1
    if (ln.EQ.1) then
        allocate(ret, source=[str%get()])
    else
        allocate (arr(ln))
        sta=1
        do k = 1,ln-1
            sto=tpIndex(str%get(),delim,cnt=k)
            arr(k)=str%get(sta,sto-1)

            sta=sto+len(delim)
        enddo
        arr(ln)=str%get(sta)

        allocate(ret, source=[(arr(k)%get(), k=1,ln)])
    endif

    return
    end function tpGetSplit_uc

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpSplit_ch(str,delim,arr) result(ret)
    implicit none

    character (len=*)               , intent(in)  :: str,delim
    type(uch), allocatable, optional, intent(out) :: arr(:)
    integer*4                                     :: k,ln,sta,sto


    if (allocated(tpSplitHold)) then
        call tpSplitHold%del(); deallocate (tpSplitHold)
        tpSplitLen=0; tpSplitAdress=0
    endif

    ln=tpCount(str,delim,overlap=false)+1
    tpSplitLen=ln; tpSplitAdress=int(loc(str),8)
    allocate (tpSplitHold(ln))

    ! write (*,*)
    ! write (*,'(A," looking for ",A)') trim(str)//'#','#'//delim//'#'
    ! write (*,'(A)') '         1         2         3         4         5'
    ! write (*,'(A)') '12345678901234567890123456789012345678901234567890'

    if (tpSplitLen.EQ.1) then
        ret=true; tpSplitHold(1)=str
    else
        sta=1
        do k = 1,ln-1
            sto=tpIndex(str,delim,cnt=k)
            tpSplitHold(k)=str(sta:sto-1)

            sta=sto+len(delim)
        enddo
        tpSplitHold(ln)=str(sta:)
    endif

    if (present(arr)) then
        if (allocated(arr)) then
            call arr%del(); deallocate (arr)
        endif
        allocate (arr(ln)); arr=tpSplitHold
    endif

    ret=true; return
    end function tpSplit_ch

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical*1 function tpSplit_uc(str,delim,arr) result(ret)
    implicit none

    type(uch)                       , intent(in)  :: str
    character(len=*)                , intent(in)  :: delim
    type(uch), allocatable, optional, intent(out) :: arr(:)
    integer*4                                     :: k,ln,sta,sto


    if (allocated(tpSplitHold)) then
        call tpSplitHold%del(); deallocate (tpSplitHold)
        tpSplitLen=0; tpSplitAdress=0
    endif

    ln=tpCount(str%get(),delim,overlap=false)+1
    tpSplitLen=ln; tpSplitAdress=int(loc(str),8)
    allocate (tpSplitHold(ln))

    ! write (*,*)
    ! write (*,'(A," looking for ",A)') trim(str)//'#','#'//delim//'#'
    ! write (*,'(A)') '         1         2         3         4         5'
    ! write (*,'(A)') '12345678901234567890123456789012345678901234567890'

    if (tpSplitLen.EQ.1) then
        ret=true; tpSplitHold(1)=str%get()
    else
        sta=1
        do k = 1,ln-1
            sto=tpIndex(str%get(),delim,cnt=k)
            tpSplitHold(k)=str%get(sta,sto-1)

            sta=sto+len(delim)
        enddo
        tpSplitHold(ln)=str%get(sta)
    endif

    if (present(arr)) then
        if (allocated(arr)) then
            call arr%del(); deallocate (arr)
        endif
        allocate (arr(ln)); arr=tpSplitHold
    endif

    ret=true; return
    end function tpSplit_uc

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpRetSplit_ch(str,nret,fcall,stat) result(ret)
    implicit none
    integer*4                                        :: nret
    character (len=*)                                :: str
    integer*4, optional                              :: stat
    logical*1, optional                              :: fcall

    character (len=tpSplitHold(splitCheck(nret))%ln) :: ret
    integer*4                                        :: ln


    ret=tpFill(ret)
    if (present(fcall)) then
        if (.NOT.fcall) then
            if (present(stat)) stat=-1
            return
        endif
    endif

    if (int(loc(str),8).NE.tpSplitAdress) then
        if (present(stat)) stat=-1
        return
    endif

    ln=len(ret)
    if (ln.LE.0) then
        if (present(stat)) stat=-1
        return
    endif

    ret=tpSplitHold(splitCheck(nret))%get()

    if (present(stat)) stat=ln; return
    end function tpRetSplit_ch

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpRetSplit_uc(str,nret,fcall,stat) result(ret)
    implicit none
    integer*4                                        :: nret
    type(uch)                                        :: str
    integer*4, optional                              :: stat
    logical*1, optional                              :: fcall

    character (len=tpSplitHold(splitCheck(nret))%ln) :: ret
    integer*4                                        :: ln


    ret=tpFill(ret)
    if (present(fcall)) then
        if (.NOT.fcall) then
            if (present(stat)) stat=-1
            return
        endif
    endif

    if (int(loc(str),8).NE.tpSplitAdress) then
        if (present(stat)) stat=-1
        return
    endif

    ln=len(ret)
    if (ln.LE.0) then
        if (present(stat)) stat=-1
        return
    endif

    ret=tpSplitHold(splitCheck(nret))%get()

    if (present(stat)) stat=ln; return
    end function tpRetSplit_uc

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer*4 function splitCheck(pos) result(ret)
    implicit none
    integer*4, intent(in) :: pos


    ret=pos
    if (pos.GT.tpSplitLen) ret=tpSplitLen
    if (tpSplitLen.EQ.0)   ret=0

    return
    end function splitCheck

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpTranslateEscapes(str,stat) result(ret)
    implicit none

    character (len=*), intent(in)    :: str
    character (len=:), allocatable   :: ret
    integer*4, intent(out), optional :: stat

    character (len=len(str))         :: tmp
    character (len=1)                :: chr
    character (len=*), parameter     :: backslash=char(92)
    character (len=*), parameter     :: keys     ='btnvfr'//backslash
    integer*4        , parameter     :: values(7)= (/ 8,9,10,11,12,13,92 /)
    integer*4                        :: k,l,pos


    if (present(stat)) stat=0

    tmp=tpFill(tmp); k=0; l=0
    do
        k=k+1; if (k.GT.len(str)) exit
        if (str(k:k).EQ.backslash) then
            k=k+1
            pos=Index(keys,str(k:k))

            ! account fucking difference of CRLF & LF for OSs
#           if (__OS==1)
                if (pos.EQ.3) then
                    l=l+1; tmp(l:l)=char(13)
                endif
#           endif

            if (pos.GT.0) then
                chr=char( values(pos) )
            else
                ! char is not escape symbol, but backslashed
                if (present(stat)) stat=-1
                ret=''; return
            endif
        else
            chr=str(k:k)
        endif
        l=l+1; tmp(l:l)=chr
    enddo

    allocate (character (len=l) :: ret); ret=tmp(:l)

    return
    end function tpTranslateEscapes

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpReplace(str,fsub,ssub,rev,cnt,quo) result(ret)
    implicit none
    character (len=*)  , intent(in) :: str,fsub,ssub
    logical*1, optional, intent(in) :: rev,quo
    integer*4, optional, intent(in) :: cnt

    character (len=len(str))        :: ustr
    character (len=len(fsub))       :: ufsub
    character (len=len(ssub))       :: ussub
    logical*1                       :: drev,dquo
    integer*4                       :: ucnt

    character (len=:), allocatable  :: tmp,ret

    integer*4                       :: rcnt,i,tlen,fnd,fprev,sprev,diff,fdiff,sdiff


    if (len(fsub).EQ.0) then
        allocate (character (len=len(str)) :: ret); ret=str
        return
    endif

    drev=false; if(present(rev)) drev=rev
    dquo=false; if(present(quo)) dquo=quo

    ucnt=len(str)+1; if (present(cnt)) ucnt=cnt
    rcnt=tpCount(str,fsub)
    diff=len(ssub)-len(fsub)

    if (.NOT.present(cnt).OR.(ucnt.GT.rcnt)) drev=false

    ucnt=min(ucnt,rcnt)

    if (ucnt.LE.0) then
        allocate (character (len=len(str)) :: ret); ret=str
        return
    endif

    if (drev) then
        ustr =tpReverse(str)
        ufsub=tpReverse(fsub)
        ussub=tpReverse(ssub)
    else
        ustr=str
        ufsub=fsub
        ussub=ssub
    endif

    if (diff.LT.0) then
        tlen=len(ustr)
    else
        tlen=len(ustr)+diff*ucnt
    endif
    allocate (character (len=tlen) :: tmp)

    sprev=1; fprev=1; i=0; rcnt=0
    if (dquo) then
        do ! fortran cannot into short-circuit evaluation
            i=i+1; fnd=tpIndex(ustr,ufsub,cnt=i)-1

            if (tpQuoted(ustr,fnd+1)) cycle

            if ((fnd.EQ.-1).OR.(rcnt.EQ.ucnt)) exit

            rcnt=rcnt+1

            fdiff=len(ustr(fprev:fnd)//ufsub); sdiff=len(ustr(fprev:fnd)//ussub)
            tmp(sprev:)=ustr(fprev:fnd)//ussub
            fprev=fprev+fdiff; sprev=sprev+sdiff
        enddo
    else
        do
            i=i+1; fnd=tpIndex(ustr,ufsub,cnt=i)-1
            if ((fnd.EQ.-1).OR.(i-1.EQ.ucnt)) exit

            rcnt=rcnt+1

            fdiff=len(ustr(fprev:fnd)//ufsub); sdiff=len(ustr(fprev:fnd)//ussub)
            tmp(sprev:)=ustr(fprev:fnd)//ussub
            fprev=fprev+fdiff; sprev=sprev+sdiff
        enddo
    endif
    i=i-1; tmp(sprev:)=ustr(fprev:)

    allocate (character (len=len(ustr)+diff*rcnt) :: ret)
    ret=tmp(1:len(ustr)+diff*rcnt)

    if (drev) ret=tpReverse(ret)

    return
    end function tpReplace

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpRemove(str,sub,rev,cnt,quo) result(ret)
    implicit none
    character (len=*),           intent(in) :: str
    character (len=*), optional, intent(in) :: sub
    logical*1,         optional, intent(in) :: rev,quo
    integer*4,         optional, intent(in) :: cnt

    character (len=:), allocatable          :: ret,usub
    character (len=len(str))                :: ustr
    logical*1                               :: drev,dquo
    integer*4                               :: ucnt

    character (len=len(str))                :: tmp
    integer*4                               :: k,ln,red,fnd,rcnt


    drev=false; if(present(rev)) drev=rev
    dquo=false; if(present(quo)) dquo=quo

    if (present(sub)) then
        ln=len(sub); allocate (character (len=ln) :: usub); usub=sub
    else
        ln=1       ; allocate (character (len=ln) :: usub); usub=' '
    endif
    ustr=str

    ucnt=len(ustr)+1; if (present(cnt)) ucnt=cnt
    tmp=tpFill(tmp); tmp=tpReplace(str,usub,'',drev,ucnt,dquo)

    rcnt=max(0,min(ucnt,tpCount(str,usub)))

    red=0
    if (dquo) then
        do k = 1,rcnt
            fnd=tpIndex(str,usub)
            if (tpQuoted(str,fnd)) cycle
            red=red+1
        enddo
    else
        red=rcnt
    endif

    allocate (character (len=len(str)-red*ln) :: ret); ret=tmp(1:len(str)-red*ln)

    return
    end function tpRemove

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpInsert(str,sub,position) result(ret)
    implicit none

    character (len=*), intent(in)     :: str,sub
    integer*4        , intent(in)     :: position
    character (len=len(str)+len(sub)) :: ret


    ret=repeat(' ',len(ret))
    if ((position.LE.0).OR.(position.GE.len(str))) return

    ret=str(:position-1)//sub//str(position:)
    return
    end function tpInsert

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpReduce(str,sub,quo) result(ret)
    implicit none
    character (len=*),           intent(in) :: str
    character (len=*), optional, intent(in) :: sub
    logical*1        , optional, intent(in) :: quo
    logical*1                               :: dquo
    character (len=:), allocatable          :: ret,usub

    character (len=len(str))                :: tmp
    integer*4                               :: i,ln,fpos,spos,rp,rlen,red


    if (present(sub)) then
        ln=len(sub); allocate (character (len=ln) :: usub); usub=sub
    else
        ln=1       ; allocate (character (len=ln) :: usub); usub=' '
    endif
    dquo=false; if (present(quo)) dquo=quo

    if (len(usub).EQ.0) then
        allocate (character (len=len(str)) :: ret); ret=str
        return
    endif

    tmp=tpFill(tmp); fpos=1; spos=1; rp=0; red=0
    if (dquo) then
        do
            if (fpos+ln-1.GT.len(str)) then
                do i = 0,len(str)-fpos
                    tmp(spos+i:spos+i)=str(fpos+i:fpos+i)
                enddo
                exit
            endif

            if ((str(fpos:fpos+ln-1).EQ.usub).AND.(.NOT.tpQuoted(str,fpos))) then
                rp=rp+1
                if (rp.EQ.1) then
                    tmp(spos:spos+ln-1)=usub
                    fpos=fpos+ln; spos=spos+ln
                else
                    fpos=fpos+ln; red=red+1
                endif
            else
                tmp(spos:spos)=str(fpos:fpos); rp=0
                fpos=fpos+1; spos=spos+1
            endif
        enddo
    else
        do
            if (fpos+ln-1.GT.len(str)) then
                do i = 0,len(str)-fpos
                    tmp(spos+i:spos+i)=str(fpos+i:fpos+i)
                enddo
                exit
            endif

            if ((str(fpos:fpos+ln-1).EQ.usub)) then
                rp=rp+1
                if (rp.EQ.1) then
                    tmp(spos:spos+ln-1)=usub
                    fpos=fpos+ln; spos=spos+ln
                else
                    fpos=fpos+ln; red=red+1
                endif
            else
                tmp(spos:spos)=str(fpos:fpos); rp=0
                fpos=fpos+1; spos=spos+1
            endif
        enddo
    endif

    allocate (character (len=len(str)-red*ln) :: ret); ret=tmp(1:len(ret))

    return
    end function tpReduce

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer*4 function tpFindKet(str,pos) result(ret)
    implicit none
    character (len=*), intent(in) :: str
    character (len=*), parameter  :: sibra='(',siket=')'
    character (len=*), parameter  :: sqbra='[',sqket=']'
    character (len=*), parameter  :: fibra='{',fiket='}'
    character (len=1)             :: bra,ket
    integer*4,         intent(in) :: pos
    integer*4                     :: i,incl


    bra=str(pos:pos)
    select case (bra)
        case (sibra); ket=siket
        case (sqbra); ket=sqket
        case (fibra); ket=fiket
        case default; ret=-1; return
    end select

    ret=0; incl=1
    do i = pos+1,len_trim(str)

        if(    str(i:i).EQ.bra) then
            incl=incl+1
        elseif(str(i:i).EQ.ket) then
            incl=incl-1
        endif

        if (incl.EQ.0) then
            ret=i
            exit
        endif

    enddo

    return
    end function tpFindKet

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpRemoveComment(str,quo) result(ret)
    implicit none
    character (len=*)              :: str
    character (len=len(str))       :: ustr
    character (len=:), allocatable :: ret
    integer*4                      :: cPos,sta
    logical*1, optional            :: quo
    logical*1                      :: dquo


    dquo=false; if (present(quo)) dquo=quo
    ustr=str; sta=0
    do
        cPos=   Index(ustr(sta+1:),trim(commentMarkup)) ! is there a comment?
        if (cPos.EQ.0) exit ! there is no comment in the string.

        if (dquo) then
            if (tpQuoted(ustr,cPos+sta)) then ! is it quoted?
                sta=sta+cPos ! it is a part of the string, it is not a comment.
                cycle
            endif
        endif

        ustr(sta+cPos:)=repeat(' ',len(ustr)-cPos-sta) ! clear the comment.
        exit
    enddo

    if (cPos.EQ.0) then
        allocate (character (len=len(ustr) ) :: ret); ret=ustr
    else
        allocate (character (len=sta+cPos-1) :: ret); ret=ustr(1:len(ret))
    endif

    return
    end function tpRemoveComment

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpQuoted(str,whereIsIt) result(ret) ! recognize whether symbol is quoted (is string)
    implicit none
    character (len=*), intent(in) :: str
    integer*4, intent(in)         :: whereIsIt
    integer*1, allocatable        :: cstatus(:)
    integer*4                     :: i,k,qslen
    integer*4, allocatable        :: num(:)
    logical*1, allocatable        :: opened(:),bcopened(:),cquo(:)
    logical*1                     :: nopened


    if (whereIsIt.GT.len_trim(str)) then; ret=false; return; endif
    if (whereIsIt.LT.1            ) then; ret=false; return; endif

    allocate (cstatus(len_trim(str)))

    qslen=len(quoset)

    allocate (num(qslen),opened(qslen),bcopened(qslen),cquo(qslen))

    do k = 1,qslen
        num(k)=tpCount(trim(str),quoset(k:k)) !how many quos of each type
    enddo
    opened=false; i=0 !at zero positions all quos are closed

    do
        i=i+1; if (i.GT.len_trim(str)) exit !exit on the end of the line

        nopened=true
        do k = 1,qslen
            nopened= nopened .AND. (.NOT.opened(k))
        enddo !none quos are opened

        do k = 1,qslen
            cquo(k)=str(i:i).EQ.quoset(k:k) !is it quo?
            if (cquo(k)) num(k)=num(k)-1    !reduce number of quos left
        enddo

        if ( nopened ) then
            do k = 1,qslen
                opened(k)=cquo(k) !is current char a quo?
            enddo
            cstatus(i)=0   !cause none are opened, this char is not quoted
            cycle
        endif

        do k = 1,qslen
            if ( opened(k).AND.cquo(k)) then  !we are ready to close quos, that are opened
                opened(k)=false; cstatus(i)=0 !cause all are closed, this char is not quoted
                exit
            endif

            if ( opened(k).AND.( .NOT.cquo(k) )) then
                if (num(k).EQ.0) then
                    cstatus(i)=0; opened(k)=false !if it is the last quo of this type (in case of odd num(k))
                else
                    cstatus(i)=k ! still opened and we expect to see it closed, chars are quoted by "k" type quos.
                endif
                exit
            endif
        enddo
    enddo

    ret=cstatus(whereIsIt).NE.0 !is quoted by none of quotes.

100 format (i3,3A,1X,L,1X,L,i2,1X,i2,1X,L,1X,L,1X,i1)

    deallocate (cstatus,num,opened,bcopened,cquo)

    return
    end function tpQuoted

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpIsIn_ch_ch(sub,str) result(ret)
    implicit none

    character (len=*), intent(in) :: str,sub


    ret=false
    if (tpIndex(str,sub).GT.0) ret=true

    return
    end function tpIsIn_ch_ch

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpIsIn_ch_uc(sub,str) result(ret)
    implicit none

    character (len=*), intent(in) :: str
    type(uch)        , intent(in) :: sub


    ret=false
    if (tpIndex(str,sub%get()).GT.0) ret=true

    return
    end function tpIsIn_ch_uc

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpIsIn_uc_ch(sub,str) result(ret)
    implicit none

    type(uch)        , intent(in) :: str
    character (len=*), intent(in) :: sub


    ret=false
    if (tpIndex(str%get(),sub).GT.0) ret=true

    return
    end function tpIsIn_uc_ch

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpIsIn_uc_uc(sub,str) result(ret)
    implicit none

    type(uch)        , intent(in) :: str,sub


    ret=false
    if (tpIndex(str%get(),sub%get()).GT.0) ret=true

    return
    end function tpIsIn_uc_uc

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpIsStrInList_ch(str,list) result(ret)
    implicit none

    character (len=*), intent(in) :: str,list(:)
    integer*4                     :: i


    ret=false
    do i = 1,UBound(list,1)
        if (trim(str).EQ.trim(list(i))) then
            ret=true
            return
        endif
    enddo

    return
    end function tpIsStrInList_ch

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpIsAnyInList(slist,tlist) result(ret)
    implicit none

    character (len=*), intent(in) :: slist(:),tlist(:)
    integer*4                     :: i,j


    ret=false
    do i = 1,UBound(slist,1)
        do j = 1,UBound(tlist,1)
            if (trim(slist(i)).EQ.trim(tlist(j))) then
                ret=true
                return
            endif
        enddo
    enddo

    return
    end function tpIsAnyInList

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpIsStrInList_uc(str,list) result(ret)
    implicit none

    type(uch), intent(in) :: str,list(:)
    integer*4             :: i


    ret=false
    do i = 1,UBound(list,1)
        if ( trim(str%get()) .EQ. trim(list(i)%get()) ) then
            ret=true
            return
        endif
    enddo

    return
    end function tpIsStrInList_uc

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpOneof(chr,set) result(ret)
    implicit none
    character (len=*), intent(in) :: set
    character (len=1), intent(in) :: chr


    ret=index(set,chr).GT.0; return
    end function tpOneof

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer*4 function tpWhichOne(chr,set) result(ret)
    implicit none
    character (len=*), intent(in) :: set
    character (len=1), intent(in) :: chr


    ret=index(set,chr); return
    end function tpWhichOne

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure logical*1 function tpSetAccordance(str,set) result(ret)
    implicit none
    character (len=*), intent(in) :: set
    character (len=*), intent(in) :: str
    integer*4                     :: i


    ret=true
    do i = 1,len(str)
        if (.NOT.tpOneof(str(i:i),set)) then
            ret=false
            return
        endif
    enddo

    return
    end function tpSetAccordance

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~ Preset Info String Functions ~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpDigits() result(ret)
    implicit none
    character (len=len(decset)) :: ret
    ret=decset; return
    end function tpDigits

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpRealNumber() result(ret)
    implicit none
    character (len=len(odnset)) :: ret
    ret=odnset; return
    end function tpRealNumber

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpLetters() result(ret)
    implicit none
    character (len=len(abcset)) :: ret
    ret=abcset; return
    end function tpLetters

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure function tpNewLine() result(ret)
    implicit none
    character (len=len(newline)) :: ret
    ret=newline; return
    end function tpNewLine

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine tpShowRuler(ln,iounit)
    implicit none

    integer*4, intent(in)       :: ln,iounit
    integer*4                   :: k,l


    do k = int(log10(float(ln))),0,-1
        do l = 10**k,ln,10**k
            write (iounit,'(A,i1\)') repeat(' ',10**k-1),mod(l/10**k,10)
        enddo
        write (iounit,*)
    enddo

    return
    end subroutine tpShowRuler

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpDeQuote(str,stat) result(ret)
    implicit none

    character (len=*)        :: str
    integer*4, optional      :: stat
    character (len=len(str)) :: ret
    integer*4                :: k,ln
    logical*1                :: fsym,lsym


    if (present(stat)) stat=0
    if (len(str).LT.1) then
        if (present(stat)) stat=-1
        return
    endif

    ln=len(str); ret=tpFill(ret)
    do k = 1,2
        fsym=str( 1: 1).EQ.quoset(k:k)
        lsym=str(ln:ln).EQ.quoset(k:k)

        if (fsym.AND.lsym) then
            if (present(stat)) stat=1
            ret=str(2:ln-1)
            return
        endif

        if ( (fsym.AND. .NOT.lsym) .OR. (.NOT.fsym .AND.lsym) ) then
            if (present(stat)) stat=-1
            ret=str
            return
        endif
    enddo

    ret=str

    return
    end function tpDeQuote

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpRealByStr(str,stat) result(ret)
    implicit none
    integer*4, parameter          :: vkind=8
    real (kind=vkind)             :: ret
    character (len=*), intent(in) :: str
    integer*4, optional           :: stat
    character (len=len(str))      :: ustr
    integer*4                     :: istat,tlen


    ret=NaN; if (present(stat)) stat=0
    istat=0; ustr=tpFill(ustr); ustr=tpLowerCase( trim( adjustl(str) ) ); tlen=len_trim(ustr)

    ! symbol set is not suitable for datatype
    if (.NOT.tpSetAccordance(trim(ustr),odnset)) goto 666

    ! in case of ',' separator
    if (tpIndex(trim(ustr),',').NE.0) ustr=tpReplace(ustr,',','.')

    read (ustr(1:tlen),*,err=666) ret

    return
666 istat=-1; if (present(stat)) stat=istat; return

    end function tpRealByStr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    complex*8 function tpComplexByStr(str,stat) result(ret)
    implicit none

    character (len=*)   :: str
    integer*4, optional :: stat
    integer*4           :: istat

    istat=0

    ret=(0._r8kind,0._r8kind); read (str,*,err=666) ret
    return
666 istat=-1; if (present(stat)) stat=istat; return

    end function tpComplexByStr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpIntByStr(str,stat) result(ret)
    implicit none
    integer*4, parameter          :: vkind=8
    integer (kind=vkind)          :: ret
    character (len=*), intent(in) :: str
    integer*4, optional           :: stat
    character (len=len(str))      :: ustr,ustr2
    integer*4                     :: istat,tlen,chck,pos,ppos,mpos,sig
    logical*4                     :: lfsemi,lfund,lf0x,lfh,isdec


    ret=int(0,vkind); if (present(stat)) stat=0; isdec=false
    istat=0; ustr2=tpFill(ustr2); ustr=tpFill(ustr); ustr2=tpLowerCase( trim( adjustl(str) ) ); tlen=len_trim(ustr2)
    ppos=tpIndex(trim(ustr2),'+'); mpos=tpIndex(trim(ustr2),'-')

    chck=0
    lfsemi=tpIndex(ustr2,':').GT.0;                     if (lfsemi) chck=chck+1
    lfund =tpIndex(ustr2,'_').GT.0;                     if (lfund ) chck=chck+1
    lf0x  =tpIndex(ustr2,'0x').EQ.1;                    if (lf0x  ) chck=chck+1
    lfh   =tpIndex(ustr2(1:tlen),'h',rev=true).EQ.tlen; if (lfh   ) chck=chck+1
    if (chck.GT.1) goto 666

    if (ppos+mpos.GT.0) then
        pos=tpIndex(ustr2,':')
        if (pos.GT.0) then
            if ((ppos.NE.0).AND.(mpos.EQ.0)) then
                if (ppos.EQ.pos+1) then
                    sig=1 ; ustr=ustr2(1:pos)//ustr2(ppos+1:tlen)
                else
                    goto 666
                endif
            elseif ((ppos.EQ.0).AND.(mpos.NE.0)) then
                if (mpos.EQ.pos+1) then
                    sig=-1; ustr=ustr2(1:pos)//ustr2(mpos+1:tlen)
                else
                    goto 666
                endif
            else
                goto 666
            endif
        endif

        if (lf0x) then
            if (ppos+mpos.EQ.2) then
                if (ppos.EQ.2) then
                    sig=1 ; ustr=ustr2(1:2)//ustr2(4:tlen)
                else
                    sig=-1; ustr=ustr2(1:2)//ustr2(4:tlen)
                endif
            else
                goto 666
            endif
        endif

        if (ppos+mpos.EQ.1) then
            if (ppos.EQ.1) then
                sig=1 ; ustr=ustr2(2:tlen)
            else
                sig=-1; ustr=ustr2(2:tlen)
            endif
        endif
    else
        sig=1; ustr=adjustl(trim(ustr2))
    endif
    tlen=len_trim(ustr)

    if (chck.EQ.1) then
        if (lfsemi) then
            pos=tpIndex(ustr,':')

            select case (ustr(1:pos-1))

                case ('bin','2','b')
                    if (.NOT.tpSetAccordance(ustr(pos+1:tlen),binset)) goto 666
                    ret=int(defIntBySystem(ustr(pos+1:tlen),'bin'),vkind)

                case ('oct','8','o')
                    if (.NOT.tpSetAccordance(ustr(pos+1:tlen),octset)) goto 666
                    ret=int(defIntBySystem(ustr(pos+1:tlen),'oct'),vkind)

                case ('dec','10','d')
                    if (.NOT.tpSetAccordance(ustr(pos+1:tlen),decset)) goto 666
                    ret=int(defIntBySystem(ustr(pos+1:tlen),'dec'),vkind)

                case ('hex','16','h')
                    if (.NOT.tpSetAccordance(ustr(pos+1:tlen),hexset)) goto 666
                    ret=int(defIntBySystem(ustr(pos+1:tlen),'hex'),vkind)

                case default; goto 666

            end select

        endif !lfsemi

        if (lfund) then
            pos=tpIndex(ustr,'_')

            select case (ustr(pos+1:tlen))

                case ('2')
                    if (.NOT.tpSetAccordance(ustr(1:pos-1),binset)) goto 666
                    ret=int(defIntBySystem(ustr(1:pos-1),'bin'),vkind)

                case ('8')
                    if (.NOT.tpSetAccordance(ustr(1:pos-1),octset)) goto 666
                    ret=int(defIntBySystem(ustr(1:pos-1),'oct'),vkind)

                case ('10'); isdec=true
                    if (.NOT.tpSetAccordance(ustr(1:pos-1),decset)) goto 666
                    ret=int(defIntBySystem(ustr(1:pos-1),'dec'),vkind)

                case ('16')
                    if (.NOT.tpSetAccordance(ustr(1:pos-1),hexset)) goto 666
                    ret=int(defIntBySystem(ustr(1:pos-1),'hex'),vkind)

                case default; goto 666

            end select
        end if

        if (lf0x) then
            if (.NOT.tpSetAccordance(ustr(3:tlen),binset)) goto 666
            ret=int(defIntBySystem(ustr(3:tlen),'bin'),vkind)
        endif

        if (lfh) then
            if (.NOT.tpSetAccordance(ustr(1:tlen-1),hexset)) goto 666
            ret=int(defIntBySystem(ustr(1:tlen-1),'hex'),vkind)
        endif
    else
        if (.NOT.tpSetAccordance(ustr(1:tlen),decset)) goto 666
        ret=int(defIntBySystem(ustr(1:tlen),'dec'),vkind)
        isdec=true
    endif

    if (isdec) ret=ret*int(sig,vkind)

    return
666 istat=-1; if (present(stat)) stat=istat; return

    end function tpIntByStr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpLogByStr(str,stat) result(ret)
    implicit none
    integer*4, parameter          :: vkind=1
    logical (kind=vkind)          :: ret
    character (len=*), intent(in) :: str
    integer*4, optional           :: stat
    character (len=len(str))      :: ustr
    integer*4                     :: istat,tlen,i
    character (len=*), parameter  :: trueset ='on|yes|true|.true.|.t.|t',&
                                         falseset='off|no|false|.false.|.f.|f'


    ret=logical(false,vkind); if (present(stat)) stat=0
    istat=0; ustr=tpFill(ustr); ustr=tpLowerCase( trim( adjustl(str) ) ); tlen=len_trim(ustr)


    voidl=tpSplit(falseset,'|')
    do i = 1,tpSplitLen
        if (ustr(1:tlen).EQ.tpRetSplit(falseset,i)) return
    enddo

    voidl=tpSplit(trueset,'|')
    do i = 1,tpSplitLen
        if (ustr(1:tlen).EQ.tpRetSplit(trueset,i)) then
            ret=logical(true,vkind)
            return
        endif
    enddo

666 istat=-1; if (present(stat)) stat=istat; return
    end function tpLogByStr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure integer*4 function defIntBySystem(str,stype) result(ret)
    implicit none

    character (len=*), intent(in) :: str,stype
    integer*4                     :: i,j,id,lstr
    integer*4, allocatable        :: arr(:)


    ret=0; lstr=len_trim(str); allocate (arr(1:lstr))
    select case (stype)
        case ('bin'); id=2
            do i = 1,lstr
            do j = 1,len(binset)
                if (str(i:i).EQ.binset(j:j)) then
                    arr(lstr-i+1)=j-1; exit
                endif
            enddo
            enddo

        case ('oct'); id=8
            do i = 1,lstr
            do j = 1,len(octset)
                if (str(i:i).EQ.octset(j:j)) then
                    arr(lstr-i+1)=j-1; exit
                endif
            enddo
            enddo

        case ('dec'); id=10
            do i = 1,lstr
            do j = 1,len(decset)
                if (str(i:i).EQ.decset(j:j)) then
                    arr(lstr-i+1)=j-1; exit
                endif
            enddo
            enddo

        case ('hex'); id=16
            do i = 1,lstr
            do j = 1,len(hexset)
                if (str(i:i).EQ.hexset(j:j)) then
                    arr(lstr-i+1)=j-1; exit
                endif
            enddo
            enddo

        case default
            deallocate (arr)
            return

    end select

    do i = 1,lstr; ret=ret+arr(i)*id**(i-1); enddo
    deallocate (arr)

    return
    end function defIntBySystem

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpRealArrayByStr(str,n,symbol,stat) result(ret)
    implicit none

    character (len=*)           :: str
    character (len=1), optional :: symbol
    character (len=1)           :: dsymbol
    integer*4, optional         :: stat
    integer*4                   :: n,istat
    real*8                      :: ret(n)
    integer*4                   :: i,sta,sto,cnt,ln


    dsymbol=','; if (present(symbol)) dsymbol=symbol

    ln=len_trim(str); sta=1; sto=ln
    if ( (str(1:1).EQ.'{') .AND. (str(ln:ln).EQ.'}') ) then
        sta=2; sto=ln-1
    endif

    if ( (str(1:1).EQ.'(') .AND. (str(ln:ln).EQ.')') ) then
        sta=2; sto=ln-1
    endif

    if ( (str(1:1).EQ.'[') .AND. (str(ln:ln).EQ.']') ) then
        sta=2; sto=ln-1
    endif

    if ( (str(1:1).EQ.'<') .AND. (str(ln:ln).EQ.'>') ) then
        sta=2; sto=ln-1
    endif

    ret=NaN; cnt=0
    do
        if (cnt.GE.n) exit

        i=tpIndex(str(sta:sto),dsymbol)

        if (sta.GT.sto) exit
        if (i.EQ.0) then
            cnt=cnt+1
            ret(cnt)=tpRealByStr(str(sta:sto),istat)
            exit
        endif

        cnt=cnt+1
        if (sta.GT.sta+i-2) then
            sta=sta+1
            cycle
        endif
        ret(cnt)=tpRealByStr(str(sta:sta+i-2),istat)
        sta=sta+i
    enddo

    if (present(stat)) stat=istat; return
    end function tpRealArrayByStr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function tpIntArrayByStr(str,n,symbol,stat) result(ret)
    implicit none

    character (len=*)           :: str
    character (len=1), optional :: symbol
    character (len=1)           :: dsymbol
    integer*4, optional         :: stat
    integer*4                   :: n,istat
    integer*4                   :: ret(n)
    integer*4                   :: i,sta,sto,cnt,ln


    dsymbol=','; if (present(symbol)) dsymbol=symbol

    ln=len_trim(str); sta=1; sto=ln
    if ( (str(1:1).EQ.'{') .AND. (str(ln:ln).EQ.'}') ) then
        sta=2; sto=ln-1
    endif

    if ( (str(1:1).EQ.'(') .AND. (str(ln:ln).EQ.')') ) then
        sta=2; sto=ln-1
    endif

    if ( (str(1:1).EQ.'[') .AND. (str(ln:ln).EQ.']') ) then
        sta=2; sto=ln-1
    endif

    if ( (str(1:1).EQ.'<') .AND. (str(ln:ln).EQ.'>') ) then
        sta=2; sto=ln-1
    endif

    ret=0; cnt=0
    do
        if (cnt.GE.n) exit

        i=index(str(sta:sto),dsymbol)
        if (i.EQ.0) then
            cnt=cnt+1
            ret(cnt)=tpIntByStr(str(sta:sto),istat)
            exit
        endif

        cnt=cnt+1
        ret(cnt)=tpIntByStr(str(sta:sta+i-2),istat)
        sta=sta+i
    enddo

    if (present(stat)) stat=istat; return
    end function tpIntArrayByStr

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine tpFinalize
    implicit none


    !continue

    return
    end subroutine tpFinalize

    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module txtParser
