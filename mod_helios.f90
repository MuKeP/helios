	!DEC$if defined(__unix)
		!DEC$define OS=2
	!DEC$else
		!DEC$define OS=1
	!DEC$endif

	module hdb

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	use glob
	use argsParser, only: apSetIOunit,apSetERRunit
	use fcontrol
	use printmod
	!use approx
	use txtParser
	use datablock, only: bdSetIOunit,bdSetERRunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	character (len=*), parameter :: heVersion   ='1.500'
	character (len=*), parameter :: heDate      ='18-Aug-2017'
	character (len=*), parameter :: heAuthor    ='Anton B. Zakharov'
	character (len=*), parameter :: heCompilDate='18-Aug-2017'

	! Memory: 32 bytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! glu = real for global use  (storage and non-accuracy-demanding procedures).
	! spu = real for special use (accuracy-demanding procedures).

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! ~~~~~ Common constants ~~~~~ !

	real(kind=rglu), parameter :: diisTol=real(1d-10,rglu) !control value for the sum of DIIS coefficients
	real(kind=rglu), parameter :: cordTol=real(1d-5 ,rglu) !tolerance for coordinate
	!real(kind=rglu), parameter :: symmTol=real(1d-10,rglu) !tolerance for the symmetry search

	! Memory: 4*rglu+rspu bytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! ~~~~~ Fundamental constants ~~~~~ !

	!doi:10.5281/zenodo.22826   2014
	real(kind=rglu), parameter :: PlankConstant     =real(6.626070040d-34 ,rglu) ! J*s
	real(kind=rglu), parameter :: BohrRadius        =real(0.52917721067d0 ,rglu) ! Angstrom
	real(kind=rglu), parameter :: HartreeEnergy     =real(27.21138602d0   ,rglu) ! eV
	real(kind=rglu), parameter :: ElementaryCharge  =real(1.6021766208d-19,rglu) ! C
	real(kind=rglu), parameter :: VacuumSpeedOfLight=real(299792458.d0    ,rglu) ! m/s

	! ~~~~~ Conversion factor ~~~~~ !

	real(kind=rglu), parameter :: dipoleToDeby=real(4.8031d0,rglu)
	real(kind=rglu), parameter :: gammaToesu  =real(5.036238d-4,rglu)

	! ~~~~~ CUE Constants ~~~~~ !

	real(kind=rglu) :: cueConstant1,cueConstant2,cueConstant4

	! Memory: 8*rglu bytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! ~~~~~ Service ~~~~~ !
	character (len=*) , parameter :: inputType='inp'

	! preset of the line length for s=screen, f=default file
	integer(kind=iglu), parameter :: setsLineLen=78,setfLineLen=79
	integer(kind=iglu), parameter :: fnLen=128 !length of the file name.

	! Memory: 3+3*iglu

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	type atom
		integer(kind=iglu) :: number,atype
		real(kind=rglu)    :: coords(3),alpha,gamma,nels
		character (len=2)  :: symbol

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
		real(kind=rglu), allocatable :: g(:,:),connect(:,:),core(:,:),coreImage(:,:),perturbation(:,:)
		real(kind=rglu), allocatable :: cueLayers(:)
	end type molecule ! Memory: uch%ln+na*atom+nb*bond+no*cueorb+rglu*(4*na*na+3)+9*iglu

	type(molecule) :: mol

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	character (len=*), dimension(16), parameter :: bdNames=['general','system','geometry','diis',&
	                                                        'polariz','rdm','scf','hypercharges',&
	                                                        'states','coulson','coupled-cluster',&
	                                                        'linear-responce','molecule','local',&
	                                                        'fci','cue']

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	type bdgeneral
		type(uch)                   :: methods,task,fname,outfile,coulombType
		real(kind=rglu)             :: alternation
		logical(kind=lglu)          :: bondsAlternated
	end type bdgeneral

	type bdgeometry
		logical(kind=lglu)          :: symmetryAccount,searchPlanar(2),searchLinear(2)
		real(kind=rglu)             :: symmetryTolerance,searchTolerance,randomDisplacement(2)
	end type bdgeometry
	
	type bdsystem
		type(uch)                   :: muttDestination
		real(kind=rglu)             :: memory
		integer(kind=iglu)          :: nNodes,verboselvl
		logical(kind=lglu)          :: allowRestart,allowMutt,muttSendTared,ignoreSIGHUP
	end type bdsystem

	type bdstates
		integer(kind=iglu)          :: nStates,spin
	end type bdstates

	type bdpolariz
		type(uch)                   :: scales
		integer(kind=iglu)          :: nPoints,maxPower
		real(kind=rglu)             :: derivStep
	end type bdpolariz

	type bddensity
		type(uch)                   :: dtype
		integer(kind=iglu)          :: nPoints,prntAccuracy
		real(kind=rglu)             :: derivStep
		logical(kind=lglu)          :: NOAnalize
	end type bddensity

	type bdcoulson
		type(uch)                   :: ctype
		integer(kind=iglu)          :: nPoints,derivPower
		real(kind=rglu)             :: derivStep
	end type bdcoulson

	type bdhypercharges
		type(uch)                   :: scales
		integer(kind=iglu)          :: naPoints,nfPoints,derivPower
		real(kind=rglu)             :: derivaStep,derivfStep
	end type bdhypercharges

	type bdcue
		integer(kind=iglu)          :: radius(0:3)
		logical(kind=lglu)          :: sparse,showBasis,local(0:3)
	end type bdcue

	type bdfci
		integer(kind=iglu)          :: nSteps,maxiters
		real(kind=rglu)             :: accuracy,zeroThreshold
	end type bdfci

	type bdscf
		type(uch)                   :: guess
		integer(kind=iglu)          :: maxiters
		real(kind=rglu)             :: accuracy,iterStep
	end type bdscf

	type bdpipek
		integer(kind=iglu)          :: maxiters
		real(kind=rglu)             :: accuracy
		logical(kind=lglu)          :: enabled
	end type bdpipek

	type bdcc
		type(uch)                   :: projType
		integer(kind=iglu)          :: maxiters
		real(kind=rglu)             :: accuracy,iterStep(3)
		logical(kind=lglu)          :: dcue,forceSpin
	end type bdcc

	type bdlr
		type(uch)                   :: guess
		integer(kind=iglu)          :: maxiters
		real(kind=rglu)             :: accuracy,iterStep
		logical(kind=lglu)          :: orthogonalize
	end type bdlr

	type bddiis
		type(uch)                   :: storage
		integer(kind=iglu)          :: steps
		logical(kind=lglu)          :: enabled
	end type bddiis

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	type (bdgeneral)      :: generalbd
	type (bdsystem)       :: systembd
	type (bdgeometry)     :: geometrybd
	type (bdstates)       :: statesbd
	type (bdpolariz)      :: polarizbd
	type (bddensity)      :: densitybd
	type (bdcoulson)      :: coulsonbd
	type (bdhypercharges) :: hyperchargesbd
	type (bdcue)          :: cuebd
	type (bdfci)          :: fcibd
	type (bdscf)          :: scfbd
	type (bdpipek)        :: pipekbd
	type (bdcc)           :: ccbd
	type (bdlr)           :: lrbd
	type (bddiis)         :: diisbd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! ~~~~~ Signals ~~~~~ !
	integer(kind=i4kind) :: SIGNSET
	integer(kind=i4kind), parameter :: SIGHUP=1,SIGCONT=18,SIGSTOP=19

	! ~~~~~ File IDs ~~~~~ !
	! su=screen; in=input; ou=output; rf=restart file; eu=error unit
	! debug=debug file; lsmf=for lsm check; unul=/dev/null or NUL
	integer(kind=iglu)    :: su,init,in,ou,eu,rf,debug,lsmf,unul
	character (len=fnLen) :: inFile,outFile,restartFile,debugFile

	! ~~~~~ Service ~~~~~ !

	integer(kind=iglu)   :: appPID,ouWidth,ouIndent
	logical(kind=lglu)   :: showInputHelp

	! ~~~~~ Time control ~~~~~ !

	real(kind=rglu) :: timeSpent(2,3)

	! Memory: 4+7*iglu

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	real   (kind=rglu), allocatable :: GlEt(:,:,:),Et(:),MEt(:,:),MMEt(:,:,:)
	integer(kind=iglu), allocatable :: pointAccordance(:,:,:),pointSet(:,:)
	integer(kind=iglu)              :: pointToPut,pointToCalc

	character (len=64) :: pauseFileName,continueFileName
	type(uch)          :: holdMethods(30)
	integer(kind=iglu) :: nMethods

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine setParams
		implicit none


		void=setThreadsNumber(systembd%nNodes)

		return
		end subroutine setParams

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine trapSignals
		implicit none

		integer(kind=i4kind) :: fid,err


		!MS$if(OS.EQ.2)
			fid=fcNewID()
			open  (fid,file=trim(pauseFileName))   ; write (fid,"('kill -19 ',i5)") appPID; close (fid)
			open  (fid,file=trim(continueFileName)); write (fid,"('kill -18 ',i5)") appPID; close (fid)
			void=fcNullID(fid)
			err=system('chmod +x '//trim(pauseFileName)); err=system('chmod +x '//trim(continueFileName))

			call PXFSTRUCTCREATE("sigset",SIGNSET,err)
			call PXFSIGADDSET(SIGNSET,SIGCONT,err)
			call PXFSIGADDSET(SIGNSET,SIGSTOP,err)
			call PXFSIGADDSET(SIGNSET,SIGHUP ,err)
		!MS$endif

		return
		end subroutine trapSignals

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine finalizeHelios
		implicit none

		integer(kind=i4kind) :: fid


		!MS$if(OS.EQ.2)
			fid=fcNewID()
			open (fid,file=trim(pauseFileName))   ; close (fid,status='delete')
			open (fid,file=trim(continueFileName)); close (fid,status='delete')
			void=fcNullID(fid)
		!MS$endif

		return
		end subroutine finalizeHelios

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		integer(kind=i4kind) function onTrap(sig_num) result(rcode)
		implicit none
		integer(kind=i4kind) :: sig_num,err
		
		
		!MS$if(OS.EQ.2)
			character :: ddate*9,ttime*8
			
			call time(ttime); call date(ddate); write (ou,99) ddate,ttime
			select case (sig_num)
				case (SIGHUP) ; write (ou,100) 'SIGHUP  accured. Ignoring.  '; return
				case (SIGINT) ; write (ou,100) 'SIGINT  accured. Stopping.  '; call finalizeHelios; stop
				case (SIGABRT); write (ou,100) 'SIGABRT accured. Stopping.  '; call finalizeHelios; stop
				case (SIGTERM); write (ou,100) 'SIGTERM accured. Stopping.  '; call finalizeHelios; stop
				case (SIGCONT); write (ou,100) 'SIGCONT accured. Waking up. '
				case (SIGSTOP); write (ou,100) 'SIGSTOP accured. Suspending.'; call PXFPAUSE(err)
				case default  ; write (ou,101) 'Signal ',sig_num,'recived. Ignoring.'; return
			end select
		 99 format (/2X,A9,1X,A8,1X\)
		100	format (2X,A/)
		101 format (2X,A,i2,A/)
		!MS$endif

		rcode=0; return
		end function onTrap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

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
		cueConstant1=real(1,rglu)/sqrt( real(2,rglu) )
		cueConstant2=cueConstant1**2
		cueConstant4=cueConstant1**4

		!define names.
		pauseFileName   =tpFill(pauseFileName)   ; pauseFileName   ='pauseproc.sh' 
		continueFileName=tpFill(continueFileName); continueFileName='continueproc.sh'

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
		rf   =fcNewID()
		debug=fcNewID()
		lsmf =fcNewID()

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine perturbate(output)
		implicit none

		logical(kind=lglu), optional :: output
		integer(kind=iglu)           :: i,j


		do i = 1,mol%nAtoms
			do j = 1,mol%nAtoms
				mol%core(i,j)=mol%coreImage(i,j)+mol%perturbation(i,j)
			enddo
		enddo

		if (present(output)) then
			if (output) call prMatrix(mol%perturbation,ou,'Perturbation matrix','^.0000',maxwidth=ouWidth)
		endif

		return
		end subroutine perturbate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine Null
		implicit none


		continue

		return
		end subroutine Null

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module hdb
