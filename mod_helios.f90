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

	character (len=*), parameter :: dbVersion='1.000'
	character (len=*), parameter :: dbDate   ='2017.03.28'
	character (len=*), parameter :: dbAuthor ='Anton B. Zakharov'

	! Memory: 32 bytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! glu = real for global use  (storage and non-accuracy-demanding procedures).
	! spu = real for special use (accuracy-demanding procedures).

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! ~~~~~ Common constants ~~~~~ !

	!tolerance in comparison of two real variables.
	real(kind=rglu), parameter :: gluCompare=10**floor(log10(epsilon(gluUnity))+real(2,rglu))
	real(kind=rspu), parameter :: spuCompare=10**floor(log10(epsilon(spuUnity))+real(2,rspu))

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

	! ~~~~~ CUE Constants ~~~~~ !

	! Unfortunately, compiler does not allow to use "sqrt" in this statement.
	real(kind=rglu) :: cueConstant1,cueConstant2,cueConstant4 !=real(1.,rglu)/sqrt( real(2.,rglu) )

	! Memory: 8*rglu bytes

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! ~~~~~ Service ~~~~~ !
	character (len=*) , parameter :: inputType='inp'

	! preset of the line length for s=screen, f=default file
	integer(kind=iglu), parameter :: setsLineLen=78,setfLineLen=150
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
		real   (kind=rglu) :: cueLevel(3)

		type(atom)     , allocatable :: atm(:)
		type(bond)     , allocatable :: bnd(:)
		type(cueorb)   , allocatable :: orb(:)
		real(kind=rglu), allocatable :: atmdist(:,:),cuedist(:,:)
		real(kind=rglu), allocatable :: g(:,:),connect(:,:),core(:,:)
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
		integer(kind=iglu)          :: radius(3)
		logical(kind=lglu)          :: sparse,showBasis
	end type bdcue

	type bdfci
		integer(kind=iglu)          :: nSteps,maxiters
		real(kind=rglu)             :: accuracy,zeroThreshold
	end type bdfci

	type bdscf
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
		logical(kind=lglu)          :: dcue
	end type bdcc

	type bdlr
		type(uch)                   :: guess
		integer(kind=iglu)          :: maxiters
		real(kind=rglu)             :: accuracy,iterStep
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
	integer(kind=i4kind) :: SIGNSET,SIGCONT,SIGSTOP,SIGHUP

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

	character (len=64) :: pauseFileName,continueFileName

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

		integer(kind=i4kind) :: fid,pid,err


		!MS$if(OS.EQ.2)
			pid=getpid(); fid=fcNewID()
			open  (fid,file=trim(pauseFileName))   ; write (fid,"('kill -19 ',i5)") pid; close (fid)
			open  (fid,file=trim(continueFileName)); write (fid,"('kill -18 ',i5)") pid; close (fid)
			void=fcNullID(fid)
			err=system('chmod +x '//trim(pauseFileName)); err=system('chmod +x '//trim(continueFileName))

			SIGHUP =1
			SIGCONT=18
			SIGSTOP=19
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
			
			call time(ttime); call date(ddate); write (ou,99) ddate,ttime; write (su,99) ddate,ttime
			select case (sig_num)

				case ( 1)
					write (su,100) 'SIGHUP  accured. Ignoring.  '
					write (ou,100) 'SIGHUP  accured. Ignoring.  '; return
				case ( 2)
					write (su,100) 'SIGINT  accured. Stopping.  '
					write (ou,100) 'SIGINT  accured. Stopping.  '; call finalizeHelios; stop
				case ( 6)
					write (su,100) 'SIGABRT accured. Stopping.  '
					write (ou,100) 'SIGABRT accured. Stopping.  '; call finalizeHelios; stop
				case (15)
					write (su,100) 'SIGTERM accured. Stopping.  '
					write (ou,100) 'SIGTERM accured. Stopping.  '; call finalizeHelios; stop
				case (18)
					write (su,100) 'SIGCONT accured. Waking up. '
					write (ou,100) 'SIGCONT accured. Waking up. '
				case (19)
					write (su,100) 'SIGSTOP accured. Suspending.'
					write (ou,100) 'SIGSTOP accured. Suspending.'; call PXFPAUSE(err)

			end select
		 99 format (/2X,A9,1X,A8,1X\)
		100	format (2X,A/)
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

		subroutine Null
		implicit none


		continue

		return
		end subroutine Null

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module hdb
