    subroutine parseInput

    use glob,       only: assignment(=)
    use glob,       only: iglu,rglu,true,void,false,i8kind,glShareMemory
    use glob,       only: uch,mid
    use hdb,        only: bdNames,eu,ou,init,ouWidth,controlFileName
    use hdb,        only: mol,generalbd,geometrybd,systembd,densitybd,scfbd,methodNames
    use hdb,        only: iterationbd,lrbd,throughbd
    use hdb,        only: scfg,lrg,ierror,outOfMemory,changeMemoryState,cparamstring
    use datablock,  only: bdParseFile,bdPrintBlockData,bdCheckExternalParameter
    use datablock,  only: bdReplaceRegisteredPatterns
    use argsParser, only: parseArgs,apArgumentFound
    use printmod,   only: prEchoFile,prJoin
    use txtParser,  only: tpAdjustc,tpIndex,operator(.in.)
    use txtParser,  only: tpSplit,tpSplitLen,tpSplitHold,tpRetSplit,tpReplace

    implicit none

    integer(kind=iglu)     :: err,k,spl,rwidth=150
    type(uch), allocatable :: cpar(:),holder(:),holder2(:)
    type(uch)              :: ebd, eparam, evalue, inp_part


    void=parseArgs(true)

    if (apArgumentFound('-i')) call inputHelp
    if (apArgumentFound('-t')) call printTemplate

    do k = 1,UBound(bdNames,1)
        err=bdParseFile(bdNames(k),generalbd%infile%get(),rwidth)
        if (err.NE.0) then
            ierror=bdNames(k)//': Error occured while parsing input file.'
            call primaryInformation('error')
        endif
    enddo

    if ( tpIndex( generalbd%infile%get(),'.inp').GT.0 ) then
        inp_part=generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.inp',rev=true)-1)
    else
        inp_part=generalbd%infile%get()
    endif

    controlFileName=inp_part%get()//'.ctrl'
    if (tpIndex(generalbd%outfile%get(), '%input%').GT.0) then
        generalbd%outfile    =tpReplace(generalbd%outfile%get(),'%input%',inp_part%get(),cnt=1)//'.out'
        generalbd%harvestfile=tpReplace(generalbd%outfile%get(),'%input%',inp_part%get(),cnt=1)//'.hrv'
    else
        generalbd%outfile    =generalbd%outfile%get()//'.out'
        generalbd%harvestfile=generalbd%outfile%get()//'.hrv'
    endif

    if (apArgumentFound('-p')) then
        void=tpSplit(cparamstring,'|',cpar)
        spl=tpSplitLen
        do k = 1,spl
            void=tpSplit(cpar(k),':',holder)
            if (tpSplitLen.NE.2) then
                ierror='Incorrect syntax of [-p|--change-parameter].'
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
    void=bdReplaceRegisteredPatterns()

    open (init,file=generalbd%infile%get())
    open (ou  ,file=generalbd%outfile%get())

    call primaryInformation('init')
    write (ou,'(A/)') tpAdjustc(' Input file ',ouWidth,'=')
    call prEchoFile(init,ou,'input %numeration >>'); close(init)
    write (ou,'(/A)') tpAdjustc(' Computation settings ',ouWidth,'=')
    do k = 1,UBound(bdNames,1)
        if ((bdNames(k).EQ.'general').AND.(generalbd%methods%ln+10.GT.ouWidth)) then
            void=bdPrintBlockData(bdNames(k),ou,width=generalbd%methods%ln+20,ignoreFreeBlock=true) ! fix?
        else
            void=bdPrintBlockData(bdNames(k),ou,width=ouWidth,ignoreFreeBlock=true)
        endif
    enddo
    write (ou,'(/A/)') tpAdjustc(' Molecule information ',ouWidth,'=')

    void=glShareMemory(int(systembd%memory*1024*1024,kind=i8kind),outOfMemory,changeMemoryState)

    ! preset default values
    geometrybd%searchPlanar(2)=false
    geometrybd%searchLinear(2)=false
    throughbd%enabled(2)=false
    iterationbd%chkStopIteration(2)=false
    lrbd%storeSolution(2)=false

    return
    end subroutine parseInput