	subroutine iterator(callsub,energysub,maxiter,eps,dnull)

	use glob     , only: mid,rglu,iglu,lglu,true,false,void,uch,uchSet,uchGet,uchDel,timeControl,gluCompare
	use printmod , only: prStrByVal
	use txtParser, only: tpSplit,tpSplitLen,tpRetSplit,operator(.in.),tpFill,tpAdjustc,tpIndex,tpNewLine
	use hdb      , only: ou

	implicit none

	interface
		subroutine callsub(iteration,epsilon,accuracy)
		use glob, only: iglu,rglu
		integer(kind=iglu), intent(in)  :: iteration
		real   (kind=rglu), intent(in)  :: epsilon
		real   (kind=rglu), intent(out) :: accuracy(5)
		end subroutine callsub
	end interface

	interface
		subroutine energysub(energy)
		use glob, only: rglu
		real   (kind=rglu), intent(out) :: energy(5)
		end subroutine energysub
	end interface

	real   (kind=rglu), parameter  :: feelDivergence=1000,frequency=-5,notRearly=20

	integer(kind=iglu), intent(in) :: maxiter
	real   (kind=rglu), intent(in) :: eps
	logical(kind=lglu)             :: dnull

	character (len=1)              :: success
	integer(kind=iglu)             :: iteration,successfulIterations,printRate,headLen
	real   (kind=rglu)             :: accuracy(5),energy(5)
	real   (kind=rglu)             :: currentAccuracy,previousAccuracy,successPercent
	real   (kind=rglu)             :: improvementPercent,meanImprovement,saveEstimatedIters
	real   (kind=rglu)             :: estimatedIters
	real   (kind=rglu)             :: timePerIter,estimatedTime,sTime(2),cTime(2),fTime(2)


	!ou=6

	sTime(1)=timeControl(sTime(2))
	iteration=0; successfulIterations=1
	currentAccuracy=0; previousAccuracy=10
	meanImprovement=0; printRate=1

	if (.NOT.dnull) call printInfo('header')
	do
		iteration=iteration+1
		if (iteration.GT.maxiter) then
			iteration=iteration-1
			call printInfo('failed')
			exit
		endif

		call callSub(iteration,eps,accuracy)
		call energysub(energy)

		cTime(1)=timeControl(cTime(2))

		timePerIter=(cTime(1)-sTime(1))/iteration

		previousAccuracy=currentAccuracy
		currentAccuracy=maxval(accuracy)

		success='u'
		if (currentAccuracy.LT.previousAccuracy) then
			successfulIterations=successfulIterations+1
			success=' '
		endif

		if (iteration.EQ.1) then
			if (timePerIter.GT.frequency) then
				printRate=1
			else
				if (timePerIter.LT.10**6*gluCompare) then
					printRate=notRearly
				else
					printRate=int(frequency/timePerIter)
				endif
			endif
			if (printRate.GT.notRearly) printRate=notRearly

			improvementPercent=0
			successPercent=100
			estimatedIters=maxiter
			estimatedTime=0
			success=' '
		else
			improvementPercent=100.*(previousAccuracy-currentAccuracy)/previousAccuracy
			successPercent=100.*(successfulIterations)/iteration
			if (improvementPercent.GT.1D-5) meanImprovement=meanImprovement+improvementPercent

			if (success.EQ.'u') then
				improvementPercent=abs(improvementPercent)
				estimatedIters=saveEstimatedIters
			else
				!estimatedIters=int( (log(eps)-log(currentAccuracy) )/log(1.-improvementPercent/100.) )
				estimatedIters=int( (log(eps)-log(currentAccuracy) )/log(1.-(meanImprovement/iteration)/100.) )

				if (estimatedIters.LT.0) estimatedIters=0

				saveEstimatedIters=estimatedIters
			endif

			if (estimatedIters.GT.maxiter) estimatedIters=maxiter-iteration
			estimatedTime=estimatedIters*timePerIter
		endif

		if (.NOT.dnull) then
			if (mod(iteration-1,printRate).EQ.0) then
				write (*,*) iteration,currentAccuracy,energy(1)
				call printInfo('iteration')
			endif
		endif

		if (currentAccuracy.GT.feelDivergence) then
			call printInfo('diverged')
			exit
		endif

		if (currentAccuracy.LT.eps) then
			call printInfo('converged')
			exit
		endif
	enddo

	return

	contains

		subroutine printInfo(string)
		implicit none

		character (len=*) :: string
		type(uch)         :: printLine


		select case (string)
			case ('header')
				printLine=uchSet(&
				'Iteration'         //tpFill(1)//'|'//tpFill(3)//'Accuracy'   //tpFill(3)//'|'//tpFill(1)//&
				'Success'           //tpFill(1)//'|'//tpFill(1)//'Improvement'//tpFill(1)//'|'//tpFill(1)//&
				'Time per iteration'//tpFill(1)//'|'//tpFill(7)//'Energy'     //tpFill(8)//'|'//tpFill(1)//&
				'Iters left'        //tpFill(1)//'|'//tpFill(1)//'Time left'  //tpFill(1)&
				)
				headLen=len( uchGet(printLine) )

				write (ou,'(/A)') tpFill(headLen,'=')
				write (ou,'(A)' ) uchGet(printLine)
				write (ou,'(A)' ) tpFill(headLen,'=')

				void=uchDel(printLine)

			case ('iteration')
				write (ou,'(A)') &

				tpFill(1)//prStrByVal(iteration,5)//&
				tpFill(4)//'|'//prStrByVal(currentAccuracy,0,6,'exp')//&
				tpFill(1)//'|'//tpFill(1)//prStrByVal(successPercent,3,3)//&
				tpFill(1)//'|'//tpFill(1)//success//' '//prStrByVal(improvementPercent,3,2)//&
				tpFill(4)//'|'//tpFill(3)//prStrByVal(timePerIter,0,4,'exp')//&
				tpFill(6)//'|'//tpFill(1)//prStrByVal(energy(1),0,12,'exp')//&
				tpFill(1)//'|'//tpFill(2)//prStrByVal(int(estimatedIters),5)//&
				tpFill(5)//'|'//prStrByVal(estimatedTime,0,4,'exp')

			case ('converged')
				fTime(1)=timeControl(fTime(2))

				write (ou,'(A)') tpFill(headLen,'=')

				write (ou,'(A)') &

				tpAdjustc(&
				'Procedure converged after '//prStrByVal(iteration)//'('//prStrByVal(successfulIterations)//') iterations'//&
				' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'. Resulting energy '//prStrByVal(energy(1),0,15,'exp')//' eV.'&
				,headLen)//tpNewLine()//&

				tpAdjustc(&
				'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent. Average iteration efficiency'//&
				prStrByVal(meanImprovement/iteration,3,2)//'%. Average iteration time '//trim(adjustl(prStrByVal(timePerIter,8,5)))//' seconds.'&
				,headLen)

				write (ou,'(A)') tpFill(headLen,'=')

			case ('failed')
				fTime(1)=timeControl(fTime(2))

				write (ou,'(A)') tpFill(headLen,'=')

				write (ou,'(A)') &

				tpAdjustc(&
				'Required accuracy was not reached after '//prStrByVal(iteration)//' iterations'//&
				' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'. Energy '//prStrByVal(energy(1),0,13,'exp')//' eV.'&
				,headLen)//tpNewLine()//&

				tpAdjustc(&
				'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent. Average iteration efficiency'//&
				prStrByVal(meanImprovement/iteration,3,2)//'%. Average iteration time '//trim(adjustl(prStrByVal(timePerIter,8,5)))//' seconds.'&
				,headLen)

				write (ou,'(A)') tpFill(headLen,'=')

			case ('diverged')
				fTime(1)=timeControl(fTime(2))

				write (ou,'(A)') tpFill(headLen,'=')

				write (ou,'(A)') &

				tpAdjustc(&
				'Procedure diverged and stopped after '//prStrByVal(iteration)//' iterations'//&
				' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'. Energy '//prStrByVal(energy(1),0,13,'exp')//' eV.'&
				,headLen)//tpNewLine()//&

				tpAdjustc(&
				'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent. Average iteration efficiency'//&
				prStrByVal(meanImprovement/iteration,3,2)//'%. Average iteration time '//trim(adjustl(prStrByVal(timePerIter,8,5)))//' seconds.'&
				,headLen)

				write (ou,'(A)') tpFill(headLen,'=')

		end select

		return
		end subroutine printInfo

	end subroutine iterator
