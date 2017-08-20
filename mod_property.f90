	module property

	use glob, only: rglu,iglu,lglu,true,false,uch,uchGet,uchSet,void
	use hdb , only: mol,scfbd,ccbd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	type(uch)          :: cmethod
	logical(kind=lglu) :: firstRun=true
	real   (kind=rglu) :: gsEnergy(5)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	private

	public :: getGroundStateEnergy,getPolarizability !,gsEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	real(kind=rglu) function getGroundStateEnergy(method) result(ret)

	use scf
	use mbpt
	use coupledCluster
	use huckel
	use fci

	implicit none

	character (len=*) :: method


	if (method.EQ.'finalize') then
		select case (uchGet(cmethod))

			case ('huckel')
				continue

			case ('hf')
				call finalizeSCF

			case ('mp2','mp3')
				call finalizeMBPT

			case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd','r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt')
				call finalizeCC

			case ('fci')
				call finalizeFCI

		end select
		firstRun=true; cmethod=uchSet(''); ret=0
		return
	endif

	do
		if (firstRun) then
			cmethod=uchSet(method)

			select case (method)

				case ('huckel')
					call energyHuckel(gsEnergy)

				case ('hf')
					call setSCFParameters
					call initSCF
					call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,false)
					call energySCF(gsEnergy)
					!call finalizeSCF

				case ('mp2','mp3')
					call setMBPTParameters(method)
					call initMBPT
					call energyMBPT(gsEnergy)
					!call finalizeMBPT

				case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd','r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt')
					call setCCParameters(method)
					call initCC
					call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false)
					call energyCC(gsEnergy)
					!call finalizeCC

				case ('fci')
					call setFCIParameters
					call initFCI
					call energyFCI(gsEnergy)
					!call finalizeFCI

			end select
			firstRun=false
			exit
		else
			if (method.EQ.uchGet(cmethod)) then
				select case (method)

					case ('huckel')
						call energyHuckel(gsEnergy)

					case ('hf')
						!call setSCFParameters
						call initSCF
						call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,false)
						call energySCF(gsEnergy)
						!call finalizeSCF

					case ('mp2','mp3')
						!call setMBPTParameters(method)
						call initMBPT
						call energyMBPT(gsEnergy)
						!call finalizeMBPT

					case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd','r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt')
						!call setCCParameters(method)
						call initCC
						call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false)
						call energyCC(gsEnergy)
						!call finalizeCC

					case ('fci')
						!call setFCIParameters
						call initFCI
						call energyFCI(gsEnergy)
						!call finalizeFCI

				end select
				exit
			else
				select case (method)

					case ('huckel')
						continue

					case ('hf')
						call finalizeSCF

					case ('mp2','mp3')
						call finalizeMBPT

					case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd','r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt')
						call finalizeCC

					case ('fci')
						call finalizeFCI

				end select
				firstRun=true
				cycle
			endif
		endif
	enddo
	ret=gsEnergy(1)

	return
	end function getGroundStateEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine getPolarizability

	use hdb      , only: generalbd,polarizbd,ou,ouWidth
	use hdb      , only: pointToCalc,pointSet,GlEt,perturbate
	use txtParser, only: tpSplit,tpSplitLen,tpRetSplit,tpFill,tpAdjustc

	implicit none

	type(uch), allocatable :: methodSet(:)

	integer(kind=iglu)      :: pp,pX,pY,pZ,mu,i,ln
	real   (kind=rglu)      :: field,rez
	character (len=ouWidth) :: chLine


	field=polarizbd%derivStep
	void=tpSplit(generalbd%methods,'+')

	allocate (methodSet(tpSplitLen))

	do i = 1,tpSplitLen
		methodSet(i)=uchSet( tpRetSplit(generalbd%methods,i) )
	enddo

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

			GlEt(pX,pY,pZ)=getGroundStateEnergy(uchGet(methodSet(i)))
		enddo
		rez=getGroundStateEnergy('finalize')
		
		call putPoint
		call showPolarizability(uchGet(methodSet(i)),ou,ouWidth)
	enddo

 99 format (A)
100 format (1X,A,1X,i3,'/',i3,3(1X,i2),1X)

	return
	end subroutine getPolarizability

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine getRDM

	implicit none




	return
	end subroutine getRDM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module property