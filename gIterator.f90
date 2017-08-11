	subroutine iterator(callsub,energysub,maxiter,eps,dnull)

	use hdb, only: rglu,iglu,lglu,true,false,ou,uch,uchSet,uchGet,timeControl
	use hdb, only: tpSplit,tpSplitLen,tpRetSplit,operator(.in.),mid,prStrByVal
	use hdb, only: tpFill,tpAdjustc,tpIndex

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

	real   (kind=rglu), parameter  :: feelDivergence=1000

	integer(kind=iglu), intent(in) :: maxiter
	real   (kind=rglu), intent(in) :: eps
	logical(kind=lglu)             :: dnull

	character (len=1)              :: success
	real   (kind=rglu)             :: accuracy(5),frequency=-10
	integer(kind=iglu)             :: iteration,successfulIterations,printRate
	real   (kind=rglu)             :: currentAccuracy,previousAccuracy
	real   (kind=rglu)             :: improvementPercent,meanImprovement
	real   (kind=rglu)             :: estimatedIters,estimatedTime,energy(5),saveEstimatedIters
	real   (kind=rglu)             :: timePerIter,sTime(2),cTime(2),successPercent
	type(uch)                      :: printLine,header


	!ou=6

	header=uchSet( 'Iteration'         //tpFill(1)//'|'//tpFill(3)//'Accuracy'   //tpFill(3)//'|'//tpFill(1)//&
	               'Success'           //tpFill(1)//'|'//tpFill(1)//'Improvement'//tpFill(1)//'|'//tpFill(1)//&
				   'Time per iteration'//tpFill(1)//'|'//tpFill(7)//'Energy'     //tpFill(8)//'|'//tpFill(1)//&
				   'Iters left'        //tpFill(1)//'|'//tpFill(1)//'Time left'  //tpFill(1) )

	sTime(1)=timeControl(sTime(2))
	iteration=0; successfulIterations=1
	currentAccuracy=0; previousAccuracy=10
	meanImprovement=0

	if (.NOT.dnull) then
		write (ou,'(/A)') tpFill (len( uchGet(header) ),'=' )
		write (ou,'( A)') uchGet(header)
		write (ou,'( A)') tpFill (len( uchGet(header) ),'=' )
	endif
	do
		iteration=iteration+1
		if (iteration.GT.maxiter) then
			write (*,*) 'Iterations failed.'
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
				printRate=int(frequency/timePerIter)
			endif
			printRate=1
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
				!estimatedIters=int( (log(eps) - log(currentAccuracy))/log(1.-improvementPercent/100.) )
				estimatedIters=int( (log(eps) - log(currentAccuracy))/log(1.-(meanImprovement/iteration)/100.) )

				saveEstimatedIters=estimatedIters
			endif

			if (estimatedIters.GT.maxiter) estimatedIters=maxiter-iteration
			estimatedTime=estimatedIters*timePerIter
		endif

		write (*,*) iteration, currentAccuracy,improvementPercent

		if (.NOT.dnull) then
			if (mod(iteration-1,printRate).EQ.0) then
				call output
			endif
		endif

		if (currentAccuracy.GT.feelDivergence) then
			write (*,*) 'Divergance occured.'
			exit
		endif

		if (currentAccuracy.LT.eps) then
			write (*,*) 'Iterations done.'
			exit
		endif
	enddo

	return

	contains

		subroutine output
		implicit none

		printLine=uchSet( tpFill(1)//prStrByVal(iteration,5) )
		printLine=uchSet( uchGet(printLine)//tpFill(4)//'|'//prStrByVal(currentAccuracy,0,6,'exp') )
		printLine=uchSet( uchGet(printLine)//tpFill(1)//'|'//tpFill(1)//prStrByVal(successPercent,3,3) )
		printLine=uchSet( uchGet(printLine)//tpFill(1)//'|'//tpFill(1)//success//' '//prStrByVal(improvementPercent,3,2) )
		printLine=uchSet( uchGet(printLine)//tpFill(4)//'|'//tpFill(3)//prStrByVal(timePerIter,0,4,'exp') )
		printLine=uchSet( uchGet(printLine)//tpFill(6)//'|'//tpFill(1)//prStrByVal(energy(1),0,12,'exp') )
		printLine=uchSet( uchGet(printLine)//tpFill(1)//'|'//tpFill(2)//prStrByVal(int(estimatedIters),5) )
		printLine=uchSet( uchGet(printLine)//tpFill(5)//'|'//prStrByVal(estimatedTime,0,4,'exp') )

		write (ou,'(A)') uchGet(printLine)

		return
		end subroutine output

	end subroutine iterator
