	module property

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	use glob          , only: rglu,iglu,lglu,true,false,uch,uchGet,uchSet,void
	use txtParser     , only: tpSplit,tpSplitLen,tpRetSplit,tpFill,tpAdjustc,operator(.in.)
	use printmod      , only: prStrByVal
	use hdb           , only: mol,scfbd,ccbd,statesbd,generalbd,polarizbd,ou,ouWidth
	use hdb           , only: pointToCalc,pointSet,GlEt,perturbate,perturbationID
	use huckel        , only: getHuckelResult
	use scf           , only: setSCFParameters,initSCF,iterationSCF,energySCF
	use scf           , only: getSCFResult,finalizeSCF
	use mbpt          , only: setMBPTParameters,initMBPT,energyMBPT,finalizeMBPT
	use coupledCluster, only: setCCParameters,initCC,iterationCC,energyCC,finalizeCC
	use fci           , only: setFCIParameters,initFCI,energyFCI,finalizeFCI,fciHoldStateEnergy
	use excitedStates , only: getExcitationEnergy
	use lrccsdModule  , only: setLRParameters,initLR,energyLR,finalizeLR,lrHoldStateEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	logical(kind=lglu)              :: firstRun=true,excitedReady
	type(uch)                       :: lastMethod
	integer(kind=iglu)              :: lastPerturbationID=-1
	real   (kind=rglu), allocatable :: lastEnergyHolder(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private
	public :: getEnergy,getPolarizability,getRDM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu) function getEnergy(method,state) result(ret)
	implicit none

	character  (len=*), intent(in)            :: method
	integer(kind=iglu), intent(in), optional  :: state
	integer(kind=iglu)                        :: dstate

	integer(kind=iglu)                        :: N
!	real   (kind=rglu)                        :: allExcitedStates(1:statesbd%nStates)

	real   (kind=rglu)                        :: excEnergy
	real   (kind=rglu)                        :: holdEnergy(5)
	real   (kind=rglu), allocatable           :: V(:,:),E(:)


	!write (*,*) '*** Method before '//uchGet(lastMethod)
	!write (*,*) '*** New method '//method
	!write (*,*) '*** First run? ',firstRun

	if (.NOT.allocated(lastEnergyHolder)) then
		allocate (lastEnergyHolder(0:statesbd%nStates))
		lastEnergyHolder=0
	endif

	dstate=0; if (present(state)) dstate=state

	if (dstate.GT.statesbd%nStates) then
		!error
		stop
	endif

	if ((lastPerturbationID.EQ.perturbationID).AND.(method.EQ.uchGet(lastMethod))) then
		if (dstate.EQ.0) then
			ret=lastEnergyHolder(0)
			return
		elseif ((dstate.GT.0).AND.(excitedReady)) then
			ret=lastEnergyHolder(dstate)
			return
		else
			continue
		endif
	endif

	N=mol%nAtoms
	allocate (V(N,N),E(N)); V=0; E=0
	do
		if (firstRun) then
			lastMethod=uchSet(method)
			lastEnergyHolder=0

			select case (method)

				case ('huckel') !done
					call getHuckelResult(holdEnergy,vectors=V,energies=E)

					lastEnergyHolder(0)=holdEnergy(1)

					if (dstate.GT.0) then
						ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
						excitedReady=true
					else
						ret=holdEnergy(1)
						excitedReady=false
					endif
					lastPerturbationID=perturbationID

				case ('hf') !done
					call setSCFParameters
					call initSCF
					call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,false)
					call energySCF(holdEnergy)
					call getSCFResult(vectors=V,energies=E)

					lastEnergyHolder(0)=holdEnergy(1)

					if (dstate.GT.0) then
						ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
						excitedReady=true
					else
						ret=holdEnergy(1)
						excitedReady=false
					endif
					lastPerturbationID=perturbationID

				case ('mp2','mp3') !done
					call setMBPTParameters(method)
					call initMBPT
					call energyMBPT(holdEnergy)

					if (dstate.GT.0) then
						ret=0; lastEnergyHolder(1:statesbd%nStates)=0
						excitedReady=true
					else
						ret=holdEnergy(1)
						excitedReady=false
					endif
					lastPerturbationID=perturbationID

				case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd',&
					  'r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt') !done
					call setCCParameters(method)
					call initCC
					call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false)
					call energyCC(holdEnergy)

					lastEnergyHolder(0)=holdEnergy(1)

					if ((method .in. ['cue-ccsd','r-ccsd','u-ccsd']).AND.(dstate.GT.0)) then
						call setLRParameters(method)
						call initLR

						ret=lrHoldStateEnergy(dstate)
						lastEnergyHolder(1:statesbd%nStates)=lrHoldStateEnergy(1:statesbd%nStates)
						excitedReady=true
					else
						ret=holdEnergy(1)
						excitedReady=false
					endif
					lastPerturbationID=perturbationID

				case ('fci') !done
					call setFCIParameters
					call initFCI
					call energyFCI(holdEnergy)

					lastEnergyHolder(0)=holdEnergy(1)

					excitedReady=true
					if (dstate.GT.0) then
						lastEnergyHolder(1:)=fciHoldStateEnergy(1:)-fciHoldStateEnergy(0)
						ret=lastEnergyHolder(dstate)
					else
						ret=lastEnergyHolder(0)
					endif
					lastPerturbationID=perturbationID

				case default
					ret=0
					lastEnergyHolder=0
					excitedReady=false
					lastPerturbationID=perturbationID

			end select
			firstRun=false
			exit
		else
			if (method.EQ.uchGet(lastMethod)) then
				select case (method)

					case ('huckel') !done
						call getHuckelResult(holdEnergy,vectors=V,energies=E)

						lastEnergyHolder(0)=holdEnergy(1)

						if (dstate.GT.0) then
							ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
							excitedReady=true
						else
							ret=holdEnergy(1)
							excitedReady=false
						endif
						lastPerturbationID=perturbationID

					case ('hf') !done
						!call setSCFParameters
						call initSCF
						call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,false)
						call energySCF(holdEnergy)
						call getSCFResult(vectors=V,energies=E)

						lastEnergyHolder(0)=holdEnergy(1)

						if (dstate.GT.0) then
							ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
							excitedReady=true
						else
							ret=holdEnergy(1)
							excitedReady=false
						endif
						lastPerturbationID=perturbationID

					case ('mp2','mp3') !done
						!call setMBPTParameters(method)
						call initMBPT
						call energyMBPT(holdEnergy)

						if (dstate.GT.0) then
							ret=0; lastEnergyHolder(1:statesbd%nStates)=0
							excitedReady=true
						else
							ret=holdEnergy(1)
							excitedReady=false
						endif
						lastPerturbationID=perturbationID

					case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd',&
						  'r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt') !done
						!call setCCParameters(method)
						call initCC
						call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false)
						call energyCC(holdEnergy)

						lastEnergyHolder(0)=holdEnergy(1)

						if ((method .in. ['cue-ccsd','r-ccsd','u-ccsd']).AND.(dstate.GT.0)) then
							call setLRParameters(method)
							call initLR

							ret=lrHoldStateEnergy(dstate)
							lastEnergyHolder(1:statesbd%nStates)=lrHoldStateEnergy(1:statesbd%nStates)
							excitedReady=true
						else
							ret=holdEnergy(1)
							excitedReady=false
						endif
						lastPerturbationID=perturbationID

					case ('fci') !done
						!call setFCIParameters
						call initFCI
						call energyFCI(holdEnergy)

						lastEnergyHolder(0)=holdEnergy(1)

						excitedReady=true
						if (dstate.GT.0) then
							lastEnergyHolder(1:)=fciHoldStateEnergy(1:)-fciHoldStateEnergy(0)
							ret=lastEnergyHolder(dstate)
						else
							ret=lastEnergyHolder(0)
						endif
						lastPerturbationID=perturbationID

					case default
						ret=0
						lastEnergyHolder=0
						excitedReady=false
						lastPerturbationID=perturbationID

				end select
				exit
			else
				select case (uchGet(lastMethod))

					case ('huckel')
						continue

					case ('hf')
						call finalizeSCF

					case ('mp2','mp3')
						call finalizeMBPT

					case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd',&
						  'r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt')
						call finalizeCC
						if (excitedReady) call finalizeLR

					case ('fci')
						call finalizeFCI

					case default
						continue

				end select
				firstRun=true
				cycle
			endif
		endif
	enddo
	deallocate (V,E)

	return
	end function getEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine getPolarizability
	implicit none

	type(uch), allocatable  :: methodSet(:)
	type(uch)               :: sendString

	integer(kind=iglu)      :: pp,pX,pY,pZ,mu,i,state,ln
	real   (kind=rglu)      :: field,rez
	character (len=ouWidth) :: chLine


	field=polarizbd%derivStep
	void=tpSplit(generalbd%methods,'+')

	allocate (methodSet(tpSplitLen))

	do i = 1,tpSplitLen
		methodSet(i)=uchSet( tpRetSplit(generalbd%methods,i) )
	enddo

	do state = 0,statesbd%nStates
		do i = 1,UBound(methodSet,1)
			do pp = 1,pointToCalc
				pX=pointSet(1,pp); pY=pointSet(2,pp); pZ=pointSet(3,pp)

				ln=len(uchGet(methodSet(i)))+19
				chLine=tpFill(chLine); write (chLine,100) uchGet(methodSet(i)),pp,pointToCalc,pX,pY,pZ

				write (ou,*)
				write (ou,99) tpAdjustc(tpFill(ln+10,'*'),ouWidth)
				write (ou,99) tpAdjustc(trim(chLine),ouWidth)
				write (ou,99) tpAdjustc(tpFill(ln+10,'*'),ouWidth)

				mol%perturbation=0
				do mu = 1,mol%nAtoms
					mol%perturbation(mu,mu)=+pX*field*mol%atm(mu)%coords(1)&
											+pY*field*mol%atm(mu)%coords(2)&
											+pZ*field*mol%atm(mu)%coords(3)
				enddo; call perturbate

				GlEt(pX,pY,pZ)=getEnergy(uchGet(methodSet(i)))
			enddo
			rez=getEnergy('finalize')
			
			if (state.EQ.0) then
				sendString=uchSet( uchGet(methodSet(i)) )
			else
				sendString=uchSet( uchGet(methodSet(i))//' State #'//prStrByVal(state) )
			endif

			call putPoint
			call showPolarizability( uchGet(sendString) ,ou,ouWidth)
		enddo
	enddo

 99 format (A)
100 format (1X,A,1X,i3,'/',i3,3(1X,i2),1X)

	return
	end subroutine getPolarizability

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine getRDM

	implicit none




	return
	end subroutine getRDM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module property