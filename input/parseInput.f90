    subroutine parseInput

    use glob,       only: assignment(=)
    use glob,       only: iglu,rglu,true,void,false,i8kind,glShareMemory
    use glob,       only: uch,mid
    use hdb,        only: bdNames,eu,ou,init,ouWidth,controlFileName
    use hdb,        only: mol,generalbd,geometrybd,systembd,densitybd,scfbd,methodNames
    use hdb,        only: iterationbd,lrbd
    use hdb,        only: scfg,lrg,ierror,outOfMemory,changeMemoryState,cparamstring
    use datablock,  only: bdParseFile,bdPrintBlockData,bdCheckExternalParameter
    use datablock,  only: bdReplaceRegisteredPatterns
    use argsParser, only: parseArgs,apArgumentFound
    use printmod,   only: prEchoFile,prJoin
    use txtParser,  only: tpAdjustc,tpIndex,operator(.in.)
    use txtParser,  only: tpSplit,tpSplitLen,tpSplitHold,tpRetSplit

    implicit none

    integer(kind=iglu)     :: err,k,spl,rwidth=150
    type(uch), allocatable :: cpar(:),holder(:),holder2(:)
    type(uch)              :: ebd, eparam, evalue


    void=parseArgs(true)

    if (apArgumentFound('-i')) call inputHelp
    if (apArgumentFound('-t')) call printTemplate

    if ( tpIndex( generalbd%infile%get(),'.inp').GT.0 ) then
        generalbd%outfile    =generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.inp',rev=true)-1)//'.out'
        generalbd%harvestfile=generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.inp',rev=true)-1)//'.hrv'
        controlFileName      =generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.inp',rev=true)-1)//'.ctrl'
    elseif (tpIndex( generalbd%infile%get(),'.',rev=true).GT.0) then
        generalbd%outfile    =generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.'   ,rev=true)-1)//'.out'
        generalbd%harvestfile=generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.'   ,rev=true)-1)//'.hrv'
        controlFileName      =generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.'   ,rev=true)-1)//'.ctrl'
    else
        generalbd%outfile    =generalbd%infile%get()//'.out'
        generalbd%harvestfile=generalbd%infile%get()//'.hrv'
        controlFileName      =generalbd%infile%get()//'.ctrl'
    endif

    do k = 1,UBound(bdNames,1)
        err=bdParseFile(bdNames(k),generalbd%infile%get(),rwidth)
        if (err.NE.0) then
            ierror=bdNames(k)//': Error occured while parsing input file.'
            call primaryInformation('error')
        endif
    enddo

    if (apArgumentFound('-p')) then
        void=tpSplit(cparamstring,'|',cpar)
        spl=tpSplitLen
        do k = 1,spl
            void=tpSplit(cpar(k),':',holder)
            if (tpSplitLen.NE.2) then
                ierror='Incorrect syntax of --change-parameter.'
                call primaryInformation('error')
            endif
            ebd = holder(1)
            void=tpSplit(holder(2),'=',holder2)
            if (tpSplitLen.NE.2) then
                ierror='Expected to get parameter=value for bd '//ebd%get()
                call primaryInformation('error')
            endif
            eparam=holder2(1)
            evalue=holder2(2)
            void=bdCheckExternalParameter(ebd,eparam,evalue,true)
        enddo
    endif

    open (ou  ,file=generalbd%outfile%get())
    open (init,file=generalbd%infile%get())

    call primaryInformation('init')
    write (ou,'(A/)') tpAdjustc('Input file',ouWidth,'=')
    call prEchoFile(init,ou,'input %numeration >>'); close(init)
    write (ou,'(/A)') tpAdjustc('Computation settings',ouWidth,'=')
    do k = 1,UBound(bdNames,1)
        if ((bdNames(k).EQ.'general').AND.(generalbd%methods%ln+10.GT.ouWidth)) then
            void=bdPrintBlockData(bdNames(k),ou,width=generalbd%methods%ln+20,ignoreFreeBlock=true)
        else
            void=bdPrintBlockData(bdNames(k),ou,width=ouWidth,ignoreFreeBlock=true)
        endif
    enddo
    write (ou,'(/A/)') tpAdjustc('Molecule information',ouWidth,'=')

    void=glShareMemory(int(systembd%memory*1024*1024,kind=i8kind),outOfMemory,changeMemoryState)

    void=bdReplaceRegisteredPatterns()

    ! preset default values
    geometrybd%searchPlanar(2)=false
    geometrybd%searchLinear(2)=false
    systembd%throughEnable(2)=false
    iterationbd%chkStopIteration(2)=false
    lrbd%storeSolution(2)=false

    return
    end subroutine parseInput