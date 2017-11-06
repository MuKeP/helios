    subroutine parseInput

    use glob,       only: assignment(=)
    use glob,       only: iglu,rglu,true,void,false
    use glob,       only: uch
    use hdb,        only: bdNames,eu,ou,init,ouWidth
    use hdb,        only: mol,generalbd,geometrybd,densitybd,fcibd,scfbd,pipekbd,ccbd,lrbd,iterationbd
    use hdb,        only: scfg,lrg
    use datablock,  only: bdParseFile,bdPrintBlockData
    use argsParser, only: parseArgs,apArgumentFound
    use printmod,   only: prEchoFile
    use txtParser,  only: tpAdjustc,tpIndex

    implicit none

    integer(kind=iglu) :: err,k,rwidth=150


    void=parseArgs(true)

    if (apArgumentFound('-i')) then
        !call ShowInputHelp
        stop
    endif

    if ( tpIndex( generalbd%infile%get(),'.inp').GT.0 ) then
        generalbd%outfile=generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.inp',rev=true)-1)//'.out'
    elseif (tpIndex( generalbd%infile%get(),'.',rev=true).GT.0) then
        generalbd%outfile=generalbd%infile%get(1,tpIndex( generalbd%infile%get(),'.'   ,rev=true)-1)//'.out'
    else
        generalbd%outfile=generalbd%infile%get()//'.out'
    endif


    do k = 1,UBound(bdNames,1)
        err=bdParseFile(bdNames(k),generalbd%infile%get(),rwidth)
        if (err.NE.0) then
            write (eu,*) bdNames(k),': Error occured while parsing input file.'; stop
        endif
    enddo

    generalbd%bondsAlternated=abs(generalbd%alternation).GT.0
    geometrybd%searchLinear(2)=false
    geometrybd%searchPlanar(2)=false

    !geometrybd%symmetryTolerance=real(10,rglu)**-geometrybd%symmetryTolerance
    !geometrybd%searchTolerance=  real(10,rglu)**-geometrybd%searchTolerance
    !densitybd%prntAccuracy=      real(10,rglu)**-densitybd%prntAccuracy
    !fcibd%accuracy=              real(10,rglu)**-fcibd%accuracy
    !fcibd%zeroThreshold=         real(10,rglu)**-fcibd%zeroThreshold
    !scfbd%accuracy=              real(10,rglu)**-scfbd%accuracy
    !pipekbd%accuracy=            real(10,rglu)**-pipekbd%accuracy
    !ccbd%accuracy=               real(10,rglu)**-ccbd%accuracy
    !lrbd%accuracy=               real(10,rglu)**-lrbd%accuracy

    !write (*,*) fcibd%accuracy,fcibd%zeroThreshold
    !stop

    !call infiniteVoidLoop

    open (ou  ,file=generalbd%outfile%get())
    open (init,file=generalbd%infile%get())


    call primaryInformation('init')
    write (ou,'(A/)') tpAdjustc('Input file',ouWidth,'=')
    call prEchoFile(init,ou,'input %numeration >>')
    write (ou,'(/A)') tpAdjustc('Computation settings',ouWidth,'=')
    do k = 1,UBound(bdNames,1)
        void=bdPrintBlockData(bdNames(k),ou,width=ouWidth,ignoreFreeBlock=true)
    enddo
    write (ou,'(/A/)') tpAdjustc('Molecule information',ouWidth,'=')

    return
    end subroutine parseInput