	program HELIOS

	use glob    , only: uchGet,void,signal,false,rglu,nullSub,pi
	use hdb     , only: onLoad,trapSignals,generalbd,setParams,onTrap,ccbd,ou
	use hdb     , only: sighup,sigabrt,sigint,sigterm,sigcont,sigstop,sigusr1
	use property, only: getPolarizability,getEnergy,getRDM

	use coupledCluster
	use lrccsdmodule

	implicit none

	real(kind=rglu) :: Ax,Eref


	call onLoad; call trapSignals
	void=signal(SIGHUP , ontrap, -1)
	void=signal(SIGABRT, ontrap, -1)
	void=signal(SIGINT , ontrap, -1)
	void=signal(SIGTERM, ontrap, -1)
	void=signal(SIGCONT, ontrap, -1)
	void=signal(SIGSTOP, ontrap, -1)
	void=signal(SIGUSR1, ontrap, -1)

	call parseInput
	call setParams
	call readMoleculeInformation

!	ou=6

!	call setCCParameters('spin-cue-ccsd')
!	call initCC
!	call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false)
!
!	call setLRParameters('spin-cue-ccsd')
!	call initLR
!
!	stop

!	Ax=getEnergy('huckel')   ; write (*,*) 'huckel   ',Ax
!	Ax=getEnergy('cue-ccs')  ; write (*,*) 'cue-ccs  ',Ax
!	Ax=getEnergy('hf')       ; write (*,*) 'hf       ',Ax
!	Ax=getEnergy('mp2')      ; write (*,*) 'mp2      ',Ax
!	Ax=getEnergy('mp3')      ; write (*,*) 'mp3      ',Ax
!	Ax=getEnergy('r-ccd')    ; write (*,*) 'r-ccd    ',Ax
!	Ax=getEnergy('u-ccd')    ; write (*,*) 'u-ccd    ',Ax
!	Ax=getEnergy('cue-ccsd') ; write (*,*) 'cue-ccsd ',Ax
!	Ax=getEnergy('r-ccsd')   ; write (*,*) 'r-ccsd   ',Ax
!	Ax=getEnergy('u-ccsd')   ; write (*,*) 'u-ccsd   ',Ax
!	Ax=getEnergy('r-ccsd(t)'); write (*,*) 'r-ccsd(t)',Ax
	Ax=getEnergy('cue-ccsdt'); write (*,*) 'cue-ccsdt',Ax
	Ax=getEnergy('u-ccsdt')  ; write (*,*) 'u-ccsdt  ',Ax
	Ax=getEnergy('r-ccsdt')  ; write (*,*) 'r-ccsdt  ',Ax
	Ax=getEnergy('fci')      ; write (*,*) 'fci      ',Ax !-getEnergy('fci')
	stop

	!Ax=getEnergy('huckel',1)   ; write (*,*) 'huckel   ',Ax
	Ax=getEnergy('hf',1)       ; write (*,*) '1 hf       ',Ax
	Ax=getEnergy('hf',2)       ; write (*,*) '2 hf       ',Ax

	Ax=getEnergy('cue-ccsd',1) ; write (*,*) '1 cue-ccsd ',Ax
	Ax=getEnergy('cue-ccsd',2) ; write (*,*) '2 cue-ccsd ',Ax

	Ax=getEnergy('r-ccsd',1)   ; write (*,*) '1 r-ccsd   ',Ax
	Ax=getEnergy('r-ccsd',2)   ; write (*,*) '2 r-ccsd   ',Ax

	Ax=getEnergy('u-ccsd',1)   ; write (*,*) '1 u-ccsd   ',Ax
	Ax=getEnergy('u-ccsd',2)   ; write (*,*) '2 u-ccsd   ',Ax

	Ax=getEnergy('fci',1)      ; write (*,*) '1 fci      ',Ax
	Ax=getEnergy('fci',2)      ; write (*,*) '2 fci      ',Ax

	stop
	select case( uchGet(generalbd%task) )
		case ('polarizability'); call getPolarizability
		case ('density')       ; call getRDM
		case ('coulson')       ; call nullSub
		case ('hypercharges')  ; call nullSub
		case ('energy')        ; call nullSub
		case ('wf-analize')    ; call nullSub
	end select

	call primaryInformation('end')

	stop
	end program HELIOS
