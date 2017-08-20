	!MS$if defined(__unix) ! 1=windows; 2=unix
		!MS$DEFINE OS=2
	!MS$else
		!MS$DEFINE OS=1
	!MS$endif

	program HELIOS

	use glob    , only: uchGet,void,signal,true,false
	use hdb     , only: onLoad,trapSignals,generalbd,setParams,onTrap,Null
	use hdb     , only: sighup,sigabrt,sigint,sigterm,sigcont,sigstop
	use hdb     , only: ccbd
	use property, only: getPolarizability
	use coupledCluster
	use lrccsdModule, only: setLRParameters,initLR

	implicit none


	call onLoad; call trapSignals
	void=signal(SIGHUP , ontrap, -1)
	void=signal(SIGABRT, ontrap, -1)
	void=signal(SIGINT , ontrap, -1)
	void=signal(SIGTERM, ontrap, -1)
	void=signal(SIGCONT, ontrap, -1)
	void=signal(SIGSTOP, ontrap, -1)

	call parseInput
	call setParams
	call readMoleculeInformation

	call setCCParameters('spin-cue-ccsd')
	call initCC
	call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false)

	call setLRParameters('spin-cue-ccsd')
	call initLR

	stop


	select case( uchGet(generalbd%task) )
		case ('polarizability'); call getPolarizability
		case ('density')       ; call Null
		case ('coulson')       ; call Null
		case ('hypercharges')  ; call Null
		case ('energy')        ; call Null
		case ('wf-analize')    ; call Null
	end select

	call primaryInformation('end')

	stop
	end program HELIOS
