	!MS$if defined(__unix) ! 1=windows; 2=unix
		!MS$DEFINE OS=2
	!MS$else
		!MS$DEFINE OS=1
	!MS$endif

	program HELIOS

	use hdb
	use coupledCluster, only: setCCParameters,initCC,iterationCC,stepCC,energyCC

	implicit none

	external :: parseInput

	integer*4 :: k
	real*8    :: sta


	call onLoad
	!MS$if(OS.EQ.2)
		call trapSignals
		void=signal(SIGHUP , ontrap, -1)
		void=signal(SIGABRT, ontrap, -1)
		void=signal(SIGINT , ontrap, -1)
		void=signal(SIGTERM, ontrap, -1)
		void=signal(SIGCONT, ontrap, -1)
		void=signal(SIGSTOP, ontrap, -1)
	!MS$endif

	call parseInput
	call setParams
	!call primaryInformation('init')
	call readMoleculeInformation

	call setCCParameters('spare-cue-ccsd')
	!call setCCParameters('cue-ccsdt')
	call initCC

	do k = 1,250
		sta=timeControl()
		write (*,'(i3\)') k
		call iterationCC
		call stepCC
		write (*,'(1X,F8.4)') timeControl()-sta
	enddo

!	call energyCC



	select case( uchGet(generalbd%task) )
		case ('polarizability'); call Null
		case ('density')       ; call Null
		case ('coulson')       ; call Null
		case ('hypercharges')  ; call Null
		case ('energy')        ; call Null
		case ('wf-analize')    ; call Null
	end select

	stop
	end program HELIOS
