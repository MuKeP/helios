    !DEC$if defined(__unix)
        !DEC$define OS=2
    !DEC$else
        !DEC$define OS=1
    !DEC$endif

    module hdb

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob      , only: rglu,iglu,lglu,i4kind,void,true,false
    use glob      , only: SIGABRT,SIGINT,SIGTERM,system
    use glob      , only: uch,uch_set
    use glob      , only: setThreadsNumber,timeControl,definePi,getpid
    use fcontrol  , only: fcNewID,fcSetIOunit,fcNullID
    use printmod  , only: prMatrix

    use argsParser, only: apSetIOunit,apSetERRunit
    use txtParser , only: tpFill
    use datablock , only: bdSetIOunit,bdSetERRunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! ~~~~~ Module information ~~~~~ !

    character (len=*), parameter :: heVersion   ='1.500a'
    character (len=*), parameter :: heDate      ='18-Aug-2017'
    character (len=*), parameter :: heAuthor    ='Anton B. Zakharov'
    character (len=*), parameter :: heCompilDate='18-Aug-2017'
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
                                                            'linear-responce','molecule','local',&
                                                            'fci','diis','iteration','scfguess', &
                                                            'lrguess','cue']

    ! ~~~~~ Methods names ~~~~~ !
    integer(kind=iglu)  , parameter :: MethodListLen=15
    character (len=*), dimension(MethodListLen), parameter :: methodNames=['huckel','hf','cue-ccs' ,&
                                                              'r-ccd','u-ccd','mp2','mp3','r-ccsd' ,&
                                                              'cue-ccsd','u-ccsd','r-ccsd(t)','fci',&
                                                              'cue-ccsdt','u-ccsdt','r-ccsdt']

    ! Memory: ~19*10+15*10 bytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DATA BLOCKS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type bdgeneral
        type(uch)          :: methods,task,infile,outfile,coulombType
        real(kind=rglu)    :: alternation
        logical(kind=lglu) :: bondsAlternated
    end type bdgeneral

    type bdgeometry
        logical(kind=lglu) :: symmetryAccount,searchPlanar(2),searchLinear(2)
        real(kind=rglu)    :: symmetryTolerance,searchTolerance,randomDisplacement(2)
    end type bdgeometry

    type bdsystem
        type(uch)          :: muttDestination
        real(kind=rglu)    :: memory
        integer(kind=iglu) :: nNodes,verboselvl
        logical(kind=lglu) :: allowRestart,allowMutt,muttSendTared,ignoreSIGHUP
    end type bdsystem

    type bditeration
        logical(kind=lglu) :: chkStagnation,chkDivergence,chkStopIteration(2)
        real(kind=rglu)    :: feelDivergence,feelStagnation,thresholdStagnation
        real(kind=rglu)    :: printFrequency,printNotRearly
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
        type(uch)          :: dtype
        integer(kind=iglu) :: nPoints,prntAccuracy
        real(kind=rglu)    :: derivStep
        logical(kind=lglu) :: NOAnalize
    end type bddensity

    type bdcoulson
        type(uch)          :: ctype
        integer(kind=iglu) :: nPoints,derivPower
        real(kind=rglu)    :: derivStep
    end type bdcoulson

    type bdhypercharges
        type(uch)          :: scales
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
        real(kind=rglu)    :: accuracy,iterStep
    end type bdscf

    type bdpipek
        integer(kind=iglu) :: maxiters
        real(kind=rglu)    :: accuracy
        logical(kind=lglu) :: enabled
    end type bdpipek

    type bdcc
        type(uch)          :: projType
        integer(kind=iglu) :: maxiters
        real(kind=rglu)    :: accuracy,iterStep(3)
        logical(kind=lglu) :: dcue,forceSpin,storeIntegrals
    end type bdcc

    type bdlr
        type(uch)          :: guess
        integer(kind=iglu) :: maxiters
        real(kind=rglu)    :: accuracy,iterStep,guessThreshold
        logical(kind=lglu) :: orthogonalize
    end type bdlr

    type bddiis
        type(uch)          :: storage
        integer(kind=iglu) :: steps
        logical(kind=lglu) :: enabled
    end type bddiis

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

        type(atom)     , allocatable :: atm(:)
        type(bond)     , allocatable :: bnd(:)
        type(cueorb)   , allocatable :: orb(:)
        real(kind=rglu), allocatable :: atmdist(:,:),cuedist(:,:)
        real(kind=rglu), allocatable :: g(:,:),core(:,:),huckelCore(:,:)
        real(kind=rglu), allocatable :: connect(:,:),coreImage(:,:),perturbation(:,:)
        real(kind=rglu), allocatable :: cueLayers(:)
    end type molecule ! Memory: uch%ln+na*atom+nb*bond+no*cueorb+rglu*(4*na*na+3)+9*iglu

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
    type (bdpipek)         :: pipekbd
    type (bdcc)            :: ccbd
    type (bdlr)            :: lrbd
    type (bddiis)          :: diisbd
    type (bditeration)     :: iterationbd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: GlEt(:,:,:),MMEt(:,:,:),MEt(:,:),Et(:)
    real   (kind=rglu), allocatable :: gEnergyHolder(:,:,:,:,:)
    integer(kind=iglu), allocatable :: pointAccordance(:,:,:),pointSet(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! ~~~~~ Derivatives ~~~~~ !
    integer(kind=iglu)   :: pointToPut,pointToCalc,perturbationID=0
    ! Memory: 2*iglu

    ! ~~~~~ File IDs ~~~~~ !
    ! su=screen; in=input; ou=output; rf=restart file; eu=error unit
    ! debug=debug file; lsmf=for lsm check; unul=/dev/null or NUL
    integer(kind=iglu)   :: su,init,in,ou,eu,unul,lrg,scfg
    ! Memory: 6*iglu

    ! ~~~~~ Service ~~~~~ !

    integer(kind=i4kind) :: SIGNSET
    integer(kind=iglu)   :: appPID,ouWidth,ouIndent,nMethods
    logical(kind=lglu)   :: showInputHelp
    real   (kind=rglu)   :: timeSpent(2,3)
    ! Memory: 4+4*iglu+lglu+6*rglu

    type(uch)            :: pauseFileName,continueFileName,holdMethods(30)
    ! Memory: based on the length

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    public
    private :: iglu,rglu,lglu,i4kind,void
    private :: uch
    private :: fcSetIOunit,apSetIOunit,bdSetIOunit,bdSetERRunit,apSetERRunit
    private :: prMatrix,tpFill,fcNewID
    private :: setThreadsNumber,timeControl,definePi,getpid

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setParams
    implicit none


    void=setThreadsNumber(systembd%nNodes)

    return
    end subroutine setParams

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine trapSignals
    implicit none

    integer(kind=i4kind) :: fid,err


    !MS$if(OS.EQ.2)
        fid=fcNewID()
        open  (fid,file=pauseFileName%get())   ; write (fid,"('kill -19 ',i5)") appPID; close (fid)
        open  (fid,file=continueFileName%get()); write (fid,"('kill -18 ',i5)") appPID; close (fid)
        void=fcNullID(fid)
        err=system('chmod +x '//pauseFileName%get()); err=system('chmod +x '//continueFileName%get())

        call PXFSTRUCTCREATE("sigset",SIGNSET,err)
        call PXFSIGADDSET(SIGNSET,SIGCONT,err)
        call PXFSIGADDSET(SIGNSET,SIGSTOP,err)
        call PXFSIGADDSET(SIGNSET,SIGHUP ,err)
        call PXFSIGADDSET(SIGNSET,SIGUSR1,err)
    !MS$endif

    return
    end subroutine trapSignals

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeHelios
    implicit none

    integer(kind=i4kind) :: fid


    !MS$if(OS.EQ.2)
        fid=fcNewID()
        open (fid,file=pauseFileName%get())   ; close (fid,status='delete')
        open (fid,file=continueFileName%get()); close (fid,status='delete')
        void=fcNullID(fid)
    !MS$endif

    return
    end subroutine finalizeHelios

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=i4kind) function onTrap(SIGNUM) result(rcode)
    implicit none
    integer(kind=i4kind) :: SIGNUM
    integer(kind=i4kind) :: err
    character            :: ddate*9,ttime*8


    rcode=0
    !MS$if(OS.EQ.2)
        call time(ttime); call date(ddate); write (ou,99) ddate,ttime
        select case (SIGNUM)
            case (SIGHUP) ; write (ou,100) 'SIGHUP  received. Ignoring.  '; rcode=1; return
            case (SIGINT) ; write (ou,100) 'SIGINT  received. Stopping.  '; call finalizeHelios; stop
            case (SIGABRT); write (ou,100) 'SIGABRT received. Stopping.  '; call finalizeHelios; stop
            case (SIGTERM); write (ou,100) 'SIGTERM received. Stopping.  '; call finalizeHelios; stop
            case (SIGCONT); write (ou,100) 'SIGCONT received. Waking up. '; rcode=1; return
            case (SIGSTOP); write (ou,100) 'SIGSTOP received. Suspending.'; call PXFPAUSE(err)
            case (SIGUSR1)
                write (ou,101) SIGUSR1; rcode=1
                if (iterationbd%chkStopIteration(1)) iterationbd%chkStopIteration(2)=true
                return
        end select
     99 format (/2X,A9,1X,A8,1X\)
    100 format (2X,A/)
    101 format (2X,'Defined StopIteration signal received (',i2,'). Iteration procedure will be interrupted.'/)
    !MS$endif

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

    !define names.
    pauseFileName   =uch_set('pauseproc.sh')
    continueFileName=uch_set('continueproc.sh')

    su=6 ! screen/console unit
    eu=6 ! error unit

    ! setting i/o units for modules.
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
    !rf   =fcNewID()
    !debug=fcNewID()
    !lsmf =fcNewID()
    lrg  =fcNewID()
    scfg =fcNewID()

    unul =fcNewID()
    !MS$if(OS.EQ.2)
        open(unul,file='/dev/null')
    !MS$else
        open(unul,file='NUL')
    !MS$endif

    call definebd
    call defineArguments

    return
    end subroutine onLoad

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine perturbate(output)
    implicit none

    logical(kind=lglu), optional :: output
    integer(kind=iglu)           :: i,j


    perturbationID=perturbationID+1
    do i = 1,mol%nAtoms
        do j = 1,mol%nAtoms
            mol%core(i,j)      =mol%coreImage(i,j)+mol%perturbation(i,j)
            mol%huckelCore(i,j)=mol%connect(i,j)  +mol%perturbation(i,j)
        enddo
    enddo

    if (present(output)) then
        if (output) call prMatrix(mol%perturbation,ou,'Perturbation matrix','^.0000',maxwidth=ouWidth)
    endif

    return
    end subroutine perturbate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module hdb
