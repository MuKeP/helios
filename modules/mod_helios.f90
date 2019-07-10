#   if defined(__GFORTRAN__)
#       define __unix 1
#   endif

#    if defined(__unix)
#       define __OS 2
#    else
#       define __OS 1
#   endif

    module hdb

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,       only: assignment (=)
    use glob,       only: rglu,iglu,lglu,r16kind,i4kind,i8kind,void,true,false,gluCompare
    use glob,       only: SIGABRT,SIGINT,SIGTERM,system
    use glob,       only: uch,mid,find,glMemoryLeft,glFromBytes
    use glob,       only: setThreadsNumber,timeControl,definePi,getpid,glSetIOunit
    use glob,       only: glSharedMemory
    use fcontrol,   only: fcNewID,fcSetIOunit,fcNullID
    use printmod,   only: prMatrix,prStrByVal

    use argsParser, only: apSetIOunit,apSetERRunit
    use txtParser,  only: tpFill
    use datablock,  only: bdSetIOunit,bdSetERRunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! ~~~~~ Module information ~~~~~ !

    character (len=*), parameter :: heVersion   ='1.610a'
    character (len=*), parameter :: heDate      ='11-Dec-2018'
    character (len=*), parameter :: heAuthor    ='Anton B. Zakharov'
    character (len=*), parameter :: heCompilDate='11-Dec-2018'
    ! Memory: 45 bytes

    ! ~~~~~ Fundamental constants ~~~~~ !

    !doi:10.5281/zenodo.22826  2014
    real(kind=rglu), parameter :: PlankConstant     =6.626070040e-34_rglu  ! J*s
    real(kind=rglu), parameter :: BohrRadius        =0.52917721067_rglu    ! Angstrom
    real(kind=rglu), parameter :: HartreeEnergy     =27.21138602_rglu      ! eV
    real(kind=rglu), parameter :: ElementaryCharge  =1.6021766208e-19_rglu ! C
    real(kind=rglu), parameter :: VacuumSpeedOfLight=299792458._rglu       ! m/s
    ! Memory: 5*rglu bytes

    ! ~~~~~ Conversion factors ~~~~~ !

    real(kind=rglu), parameter :: dipoleToDeby=4.8031_rglu
    real(kind=rglu), parameter :: gammaToesu  =5.036238e-4_rglu
    ! Memory: 2*rglu bytes

    ! ~~~~~ cue constants ~~~~~ !

    ! compiler does not allow using "sqrt" in constant declaration statements.
    real(kind=rglu) :: cueConstant1,cueConstant2,cueConstant4
    ! Memory: 3*rglu bytes

    ! ~~~~~ Service ~~~~~ !

    ! preset of the line length for output file
    integer(kind=iglu)  , parameter :: setfLineLen=79
    ! signals
    integer(kind=i4kind), parameter :: SIGHUP=1,SIGUSR1=10,SIGCONT=18,SIGSTOP=19
    ! Memory: iglu+16

    ! ~~~~~ Block data names ~~~~~ !
    integer(kind=iglu)  , parameter :: bdListLen=19
    character (len=*), dimension(bdListLen), parameter :: bdNames=['general','system','geometry',&
                                                            'polariz','rdm','scf','hypercharges',&
                                                            'states','coulson','coupled-cluster',&
                                                            'linear-response','molecule','local',&
                                                            'lrguess','iteration','scfguess',    &
                                                            'fci','cue','field']

    ! ~~~~~ Methods names ~~~~~ !
    integer(kind=iglu)  , parameter :: MethodListLen=15
    character (len=*), dimension(MethodListLen), parameter :: methodNames=['huckel','rhf','cue-ccs',&
                                                              'r-ccd','u-ccd','mp2','mp3','r-ccsd' ,&
                                                              'cue-ccsd','u-ccsd','r-ccsd(t)','fci',&
                                                              'cue-ccsdt','u-ccsdt','r-ccsdt']

    ! Memory: ~19*10+15*10 bytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA BLOCKS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type bdgeneral
        type(uch)          :: methods,task,infile,outfile,harvestfile,coulombType
        real(kind=rglu)    :: alternation
        logical(kind=lglu) :: bondsAlternated
    end type bdgeneral

    type bdgeometry
        logical(kind=lglu) :: symmetryAccount,searchPlanar(2),searchLinear(2)
        real(kind=rglu)    :: symmetryTolerance,searchTolerance,randomDisplacement(2)
    end type bdgeometry

    type bdsystem
        type(uch)            :: muttDestination,memoryUnits,throughHeader,throughFile,throughPrefix
        real(kind=rglu)      :: memory,memoryThreshold
        integer(kind=iglu)   :: nNodes,verboselvl
        integer(kind=i8kind) :: imemory
        logical(kind=lglu)   :: allowRestart,allowMutt,muttSendTared,ignoreSIGHUP,harvest
        logical(kind=lglu)   :: memoryReport,throughEnable(2)
    end type bdsystem

    type bditeration
        logical(kind=lglu) :: chkStagnation,chkDivergence,chkStopIteration(2),afterPause
        real(kind=rglu)    :: feelDivergence,feelStagnation,thresholdStagnation
        real(kind=rglu)    :: printFrequency,printNotRearly

        ! 0 - iteration, 1 - converged, -1 - failed, -2 - diverged, -3 - stagnated, -4 - iterrupted
        integer(kind=iglu) :: lastProcedureStatus

        logical(kind=lglu) :: doRestart

        logical(kind=lglu) :: energy
    end type bditeration

    type bdstates
        integer(kind=iglu) :: nStates,spin
    end type bdstates

    type bdpolariz
        type(uch)          :: scales
        integer(kind=iglu) :: nPoints,maxPower
        real(kind=rglu)    :: derivStep
    end type bdpolariz

    type bddensity
        type(uch)          :: dtype,scharges,sorders,gcharges
        integer(kind=iglu) :: nPoints,prntAccuracy
        real(kind=rglu)    :: derivStep
        logical(kind=lglu) :: NOAnalize,forceNumerical
    end type bddensity

    type bdcoulson
        type(uch)          :: ctype,selected,sdiagonals,soffdiagonals
        integer(kind=iglu) :: nPoints,derivPower,prntAccuracy
        real(kind=rglu)    :: derivStep
    end type bdcoulson

    type bdhypercharges
        type(uch)          :: scales,scharges,gcharges,dtype
        integer(kind=iglu) :: naPoints,nfPoints,derivPower
        real(kind=rglu)    :: derivaStep,derivfStep
    end type bdhypercharges

    type bdcue
        integer(kind=iglu) :: radius(0:3)
        logical(kind=lglu) :: sparse,showBasis,local(0:3)
    end type bdcue

    type bdfci
        integer(kind=iglu) :: nSteps,maxiters
        real(kind=rglu)    :: accuracy,zeroThreshold
    end type bdfci

    type bdscf
        type(uch)          :: guess,exctype
        integer(kind=iglu) :: maxiters
        real(kind=rglu)    :: accuracy,iterStep,iterStepVariation,iterStepChange
        logical(kind=lglu) :: keep,achieveSolution

        logical(kind=lglu) :: exceeded
    end type bdscf

    type bdlocal
        type(uch)          :: method
        integer(kind=iglu) :: maxiters
        real(kind=rglu)    :: accuracy
        logical(kind=lglu) :: enabled
    end type bdlocal

    type bdcc
        type(uch)          :: projType,diisStorage
        integer(kind=iglu) :: maxiters,diisSteps,wfSwitches
        real(kind=rglu)    :: accuracy,iterStep(3),printThreshold
        logical(kind=lglu) :: dcue,forceSpin,storeIntegrals,diisEnabled
    end type bdcc

    type bdlr
        type(uch)          :: guess,diisStorage
        integer(kind=iglu) :: maxiters,diisSteps
        real(kind=rglu)    :: accuracy,iterStep(2),guessThreshold
        logical(kind=lglu) :: orthogonalize,diisEnabled
    end type bdlr

    type bdfield
        real(kind=rglu)    :: strength(3)
    end type bdfield

    type cparam
        type(uch)          :: bd,param,value
    end type

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type atom
        integer(kind=iglu) :: number,atype
        real   (kind=rglu) :: coords(3),alpha,gamma,nels
        character  (len=2) :: symbol

        !integer(kind=iglu) :: mem=3*iglu+6*rglu+2
    end type atom

    type bond
        integer(kind=iglu) :: atoms(2)
        real(kind=rglu)    :: coords(3),beta,length
        character(len=6)   :: kind

        !integer(kind=iglu) :: mem=3*iglu+5*rglu+6
    end type bond

    type cueorb
        integer(kind=iglu) :: atoms(2),nels,ova !ova - occ-vac accordance
        real(kind=rglu)    :: coords(3),coef(2)

        !integer(kind=iglu) :: mem=6*iglu+3*rglu
    end type

    type molecule
        type(uch)          :: name
        integer(kind=iglu) :: nAtoms,nBonds,nEls,nOrbs,naEls,nbEls,nEthylenes
        integer(kind=iglu) :: uniqueAtoms,uniqueBonds,chBonds
        integer(kind=iglu) :: nCUELayers

        real   (kind=rglu) :: clearance(3)
        real   (kind=rglu) :: cueLevel(0:3)

        type(atom)        , allocatable :: atm(:)
        type(bond)        , allocatable :: bnd(:)
        type(cueorb)      , allocatable :: orb(:)
        real(kind=rglu)   , allocatable :: atmdist(:,:),cuedist(:,:)
        real(kind=rglu)   , allocatable :: g(:,:),core(:,:),huckelCore(:,:)
        real(kind=rglu)   , allocatable :: connect(:,:),coreImage(:,:),perturbation(:,:)
        real(kind=rglu)   , allocatable :: cueLayers(:)
        integer(kind=iglu), allocatable :: inverse(:)
    end type molecule ! Memory: uch%ln+na*atom+nb*bond+no*cueorb+rglu*(4*na*na+3)+(9+na)*iglu

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type (molecule)        :: mol
    type (bdgeneral)       :: generalbd
    type (bdsystem)        :: systembd
    type (bdgeometry)      :: geometrybd
    type (bdstates)        :: statesbd
    type (bdpolariz)       :: polarizbd
    type (bddensity)       :: densitybd
    type (bdcoulson)       :: coulsonbd
    type (bdhypercharges)  :: hyperchargesbd
    type (bdcue)           :: cuebd
    type (bdfci)           :: fcibd
    type (bdscf)           :: scfbd
    type (bdlocal)         :: localbd
    type (bdcc)            :: ccbd
    type (bdlr)            :: lrbd
    type (bditeration)     :: iterationbd
    type (bdfield)         :: fieldbd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: GlEt(:,:,:),MMEt(:,:,:),MEt(:,:),Et(:)
    real   (kind=rglu), allocatable :: gEnergyHolder(:,:,:,:,:)
    integer(kind=iglu), allocatable :: pointAccordance(:,:,:),pointSet(:,:),atomEqu(:,:),bondEqu(:,:)

    type(cparam)      , allocatable :: cparHolder(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! ~~~~~ Derivatives ~~~~~ !
    integer(kind=iglu)   :: pointToPut,pointToCalc,perturbationID=0
    logical(kind=lglu)   :: noPerturbation=true
    ! Memory: 3*iglu

    ! ~~~~~ File IDs ~~~~~ !
    ! su=screen; in=input; ou=output; rf=restart file; eu=error unit
    ! debug=debug file; lsmf=for lsm check; unul=/dev/null or NUL
    integer(kind=iglu)   :: su,init,in,ou,eu,unul,lrg,scfg
    ! Memory: 8*iglu

    ! ~~~~~ Service ~~~~~ !

    integer(kind=i4kind) :: SIGNSET
    integer(kind=iglu)   :: appPID,ouWidth,ouIndent,nMethods
    type(uch)            :: showInputHelp
    real   (kind=rglu)   :: timeSpent(2,3)
    integer(kind=i8kind) :: store_memory
    logical(kind=lglu)   :: afterPause
    ! Memory: 4+4*iglu+lglu+6*rglu+8

    type(uch)            :: controlFileName,holdMethods(30)
    ! Memory: based on the length

    type(uch)            :: cparamstring
    type(uch)            :: iheader,ierror
    logical(kind=lglu)   :: ipFailed

    integer(kind=iglu)   :: CATOM

    ! ~~~~~ Single Session Output ~~~~~ !

    integer(kind=iglu)   :: single_io=-1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    public
    private :: iglu,rglu,lglu,i4kind,void
    private :: uch
    private :: fcSetIOunit,apSetIOunit,bdSetIOunit,bdSetERRunit,apSetERRunit
    private :: prMatrix,tpFill,fcNewID
    private :: setThreadsNumber,timeControl,definePi,getpid
    private :: single_io

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setParams
    implicit none

    integer(kind=i4kind) :: fid,err


    void=setThreadsNumber(systembd%nNodes)
    iterationbd%afterPause=false
    ipFailed=false

    ! only for POSIX systems
#   if(__OS==2)
        fid=fcNewID()
        open(fid,file=controlFileName%get()); write(fid,100) trim(adjustl(prStrByVal(appPID))); close(fid)
        void=fcNullID(fid)
        err=system('chmod +x '//controlFileName%get())

100 format('MPID=',A/&
           'case $1 in'/&
           4X,'stop)'/4X         ,'kill -2 $MPID'/4X,';;'/&
           4X,'kill)'/4X         ,'kill -9 $MPID'/4X,';;'/&
           4X,'stopiteration)'/4X,'kill -10 $MPID'/4X,';;'/&
           4X,'continue)'/4X     ,'kill -18 $MPID'/4X,';;'/&
           4X,'pause)'/4X        ,'kill -19 $MPID'/4X,';;'/&
           'esac')
#   endif

    return
    end subroutine setParams

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine trapSignals
    implicit none

    integer(kind=i4kind) :: err


    ! only for POSIX systems
#   if(__OS==2)
        call PXFSTRUCTCREATE('sigset',SIGNSET,err)
        call PXFSIGADDSET(SIGNSET,SIGCONT,err)
        call PXFSIGADDSET(SIGNSET,SIGSTOP,err)
        call PXFSIGADDSET(SIGNSET,SIGHUP ,err)
        call PXFSIGADDSET(SIGNSET,SIGUSR1,err)
#   endif

    return
    end subroutine trapSignals

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeHelios
    implicit none

    integer(kind=iglu) :: err


    ! only for POSIX systems
#   if(__OS==2)
        err=system('rm '//controlFileName%get())
#   endif

    end subroutine finalizeHelios

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=i4kind) function onTrap(SIGNUM) result(rcode)
    implicit none
    integer(kind=i4kind) :: SIGNUM
    integer(kind=i4kind) :: err
    character            :: ddate*9,ttime*8


    rcode=0

    ! only for POSIX systems
#   if(__OS==2)
        call time(ttime); call date(ddate); write (ou,99) ddate,ttime
        select case (SIGNUM)
            case (SIGINT) ; write (ou,100) 'SIGINT  received. Stopping.'; call finalizeHelios; stop
            case (SIGABRT); write (ou,100) 'SIGABRT received. Stopping.'; call finalizeHelios; stop
            case (SIGTERM); write (ou,100) 'SIGTERM received. Stopping.'; call finalizeHelios; stop

            case (SIGSTOP)
                write (ou,100) 'SIGSTOP received. Suspending.'
                call PXFPAUSE(err)
                rcode=1
                return

            case (SIGHUP)
                write (ou,100) 'SIGHUP  received. Ignoring.'
                rcode=1
                return

            case (SIGCONT)
                write (ou,100) 'SIGCONT received. Waking up.'
                iterationbd%afterPause=true
                rcode=1
                return

            case (SIGUSR1)
                write (ou,101) SIGUSR1; rcode=1
                if (iterationbd%chkStopIteration(1)) iterationbd%chkStopIteration(2)=true
                return

        end select
     99 format (/2X,A9,1X,A8,1X\)
    100 format (2X,A/)
    101 format ( 2X,'Defined StopIteration signal received (',i2,').'/&
                22X,'Iteration procedure will be interrupted.'/)
#   endif

    return
    end function onTrap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine onLoad
    implicit none


    ! time started.
    timeSpent(1,1)=timeControl(timeSpent(2,1))

    ! some actions.
    call random_seed
    void=definePi()

    appPID=getpid()

    ! some constants.
    ouWidth=setfLineLen; ouIndent=4
    cueConstant1=1._rglu/sqrt(2._rglu)
    cueConstant2=cueConstant1**2
    cueConstant4=cueConstant1**4

    su=6 ! screen/console unit
    eu=6 ! error unit

    ! setting i/o units for modules.
    ! glob
    call glSetIOunit(su)

    ! args parser
    call apSetIOunit(su); call apSetERRunit(eu)

    ! file control
    call fcSetIOunit(su)

    ! data block
    call bdSetIOunit(su); call bdSetERRunit(eu)

    ! prepare i/o units.
    init =fcNewID()
    in   =fcNewID()
    ou   =fcNewID()
    !rf   =fcNewID() !restart file
    !debug=fcNewID()
    !lsmf =fcNewID()
    lrg  =fcNewID()
    scfg =fcNewID()

    unul =fcNewID()
#   if (__OS==2)
        open(unul,file='/dev/null')
#   else
        open(unul,file='NUL')
#   endif

    call definebd
    call defineArguments

    call resetIterationState

    return
    end subroutine onLoad

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine outOfMemory(label, req, left)
    implicit none

    character (len=*)   , intent(in) :: label
    integer(kind=i8kind), intent(in) :: req,left

    character (len=:), allocatable   :: units,creq,cleft


    units=systembd%memoryUnits%get()
    creq=trim(prStrByVal(glFromBytes(req, units), '^.00'))
    cleft=trim(prStrByVal(glFromBytes(left, units), '^.00'))

    ierror=label//': not enough memory. Required '//creq//' '//units//', left '//cleft//' '//units//'.'
    call primaryInformation('error')

    return
    end subroutine outOfMemory

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine changeMemoryState(label, size)
    implicit none

    character (len=*)   , intent(in) :: label
    integer(kind=i8kind), intent(in) :: size
    real(kind=r16kind)               :: rsize,rleft
    character (len=:), allocatable   :: units,csize,cleft


    if (systembd%memoryReport) then
        units=systembd%memoryUnits%get()

        rsize=glFromBytes(size, units)
        rleft=glFromBytes(glSharedMemory, units)

        csize=prStrByVal(rsize, '^.00')
        cleft=prStrByVal(rleft, '^.00')

        if (rleft.LT.0) return

        if (abs(rsize).GT.systembd%memoryThreshold) then
            if (rsize.GT.0) then
                write (ou,'(/A)') ' ***** '//label//'. Allocate: '//csize//' '//units//&
                                                    ', Left: '//cleft//' '//units//' ***** '
            elseif(rsize.LT.0) then
                write (ou,'(/A)') ' ***** '//label//'. Free: '//csize//' '//units//&
                                                    ', Left: '//cleft//' '//units//' ***** '
            endif
        endif

    endif

    return
    end subroutine changeMemoryState

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine perturbate(output)
    implicit none

    logical(kind=lglu), optional :: output
    integer(kind=iglu)           :: i,j


    noPerturbation=maxval(abs(mol%perturbation)).LE.gluCompare

    perturbationID=perturbationID+1
    do i = 1,mol%nAtoms
        do j = 1,mol%nAtoms
            mol%core(i,j)      =mol%coreImage(i,j)+mol%perturbation(i,j)
            mol%huckelCore(i,j)=mol%connect(i,j)  +mol%perturbation(i,j)
        enddo
    enddo
    ! call prMatrix(mol%huckelCore,100,'Huckel core.','^.0000',maxwidth=200)

    if (present(output)) then
        if (output) call prMatrix(mol%perturbation,ou,'Perturbation ('//&
                                  prStrByVal(perturbationID)//')','^.0000',maxwidth=ouWidth)
    endif

    return
    end subroutine perturbate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setIterationHeader(str)
    implicit none

    character(len=*), intent(in), optional :: str


    if (present(str)) then
        iheader=str
    else
        iheader=''
    endif

    return
    end subroutine setIterationHeader

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function getMethodNumber(method) result(ret)
    implicit none

    character(len=*) :: method


    ret=find(methodNames,method)
    if (ret.GT.0) return

    stop 'Internal error (hdb::getMethodNumber): Unknown method.'
    end function getMethodNumber

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine singleSession(str)
    implicit none

    character (len=*), intent(in) :: str
    integer(kind=iglu)            :: err


    ! not enabled
    if (.NOT.systembd%throughEnable(1)) return

    ! need to open file
    if (.NOT.systembd%throughEnable(2)) then
        single_io=fcNewID()

        ! check wether file exists
        open(single_io, file=systembd%throughFile%get(), status='old', iostat=err)

        ! file does not exist
        if (err.NE.0) then
            systembd%throughEnable(2)=true
            open (single_io, file=systembd%throughFile%get())
            write(single_io, '(A)') systembd%throughHeader%get()
        else ! file exists
            close(single_io)
            systembd%throughEnable(2)=true
            open (single_io, file=systembd%throughFile%get(), access='append')
            write(single_io, *)
        endif
    endif

    write(single_io, '(A\)') str

    return
    end subroutine singleSession

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine resetIterationState(state)
    implicit none

    logical(kind=lglu), optional :: state


    if (present(state)) then
        iterationbd%energy=state
    else
        iterationbd%energy=true
    endif

    return
    end subroutine resetIterationState

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module hdb
