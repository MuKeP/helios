	subroutine iterator(callsub,energysub,maxiter,eps,dnull)

	use glob     , only: mid,rglu,iglu,lglu,true,false,void,uch,uchSet,uchGet,uchDel,timeControl,gluCompare
	use printmod , only: prStrByVal
	use txtParser, only: tpSplit,tpSplitLen,tpRetSplit,operator(.in.),tpFill,tpAdjustc,tpIndex,tpNewLine
	use hdb      , only: ou,ouWidth

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

!	ou=6

	sTime(1)=timeControl(sTime(2))
	iteration=0; successfulIterations=1
	currentAccuracy=0; previousAccuracy=10
	meanImprovement=0; printRate=1

	if (.NOT.dnull) call printInfo('header',ou)
	do
		iteration=iteration+1
		if (iteration.GT.maxiter) then
			iteration=iteration-1
			if (.NOT.dnull) call printInfo('failed',ou)
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
				call printInfo('iteration',ou)
			endif
		endif

		if (currentAccuracy.GT.feelDivergence) then
			if (.NOT.dnull) call printInfo('diverged',ou)
			exit
		endif

		if (currentAccuracy.LT.eps) then
			if (.NOT.dnull) call printInfo('converged',ou)
			exit
		endif
	enddo

	return

	contains

		subroutine printInfo(string,iou)
		implicit none

		character (len=*)  :: string
		integer(kind=iglu) :: i,iou


		headLen=79
		select case (string)
			case ('header')
				write (iou,'(/A)') tpFill(headLen,'=')
				write (iou, '(A)') tpAdjustc('1 - Iteration, 2 - Accuracy, 3 - Step efficiency, 4 - Time per iteration',headLen)
				write (iou, '(A)') tpAdjustc('5 - Energy, 6 - Estimated iterations left, 7 - Estimated time left',headLen)
				write (iou, '(A)') tpFill(headLen,'=')
				write (iou, '(2X,i1,2X,"|",6X,i1,6X,"|",4X,i1,3X,"|",6X,i1,5X,"|",9X,i1,9X,"|",2X,i1,2X,"|",5X,i1,5X)' ) (i, i = 1,7)
				write (iou, '(A)') tpFill(headLen,'=')

			case ('iteration')
				write (iou,'(A,"|",A,"|",A,1X,A,"|",A,1X,"|",A,"|",A,"|",A)') &
				                tpAdjustc(trim(adjustl(prStrByVal(iteration,5))),5),&
								prStrByVal(currentAccuracy,0,6,'exp'),&
				                success,prStrByVal(improvementPercent,3,2),&
								prStrByVal(timePerIter,0,4,'exp'),&
				                prStrByVal(energy(1),0,12,'exp'),&
								tpAdjustc(trim(adjustl(prStrByVal(int(estimatedIters),5))),5),&
				                prStrByVal(estimatedTime,0,4,'exp')

			case ('converged')
				fTime(1)=timeControl(fTime(2))

				write (iou,'(A)') tpFill(headLen,'=')
				write (iou,'(A)') tpAdjustc('Procedure converged after '//prStrByVal(iteration)//'('//prStrByVal(successfulIterations)//&
				                            ') iterations'//' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
				write (iou,'(A)') tpAdjustc('Resulting energy '//prStrByVal(energy(1),0,15,'exp')//' eV. '//&
				                            'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
				write (iou,'(A)') tpAdjustc('Avg iteration efficiency'//prStrByVal(meanImprovement/iteration,3,2)//&
				                            '%. Avg iteration time '//trim(adjustl(prStrByVal(timePerIter,8,5)))//' seconds.',headLen)
				write (iou,'(A)') tpFill(headLen,'=')

			case ('failed')
				fTime(1)=timeControl(fTime(2))

				write (iou,'(A)') tpFill(headLen,'=')
				write (iou,'(A)') tpAdjustc('Accuracy was not reached after '//prStrByVal(iteration)//' iterations'//&
				                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
				write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
				                            'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
				write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(meanImprovement/iteration,3,2)))//'%. '//&
				                            'Avg iteration time '//trim(adjustl(prStrByVal(timePerIter,8,5)))//' seconds.',headLen)
				write (iou,'(A)') tpFill(headLen,'=')

			case ('diverged')
				fTime(1)=timeControl(fTime(2))

				write (iou,'(A)') tpFill(headLen,'=')
				write (iou,'(A)') tpAdjustc('Procedure diverged and stopped after '//prStrByVal(iteration)//' iterations'//&
				                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
				write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
				                            'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
				write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(meanImprovement/iteration,3,2)))//'%. '//&
				                            'Avg iteration time '//trim(adjustl(prStrByVal(timePerIter,8,5)))//' seconds.',headlen)
				write (iou,'(A)') tpFill(headLen,'=')

		end select

		return
		end subroutine printInfo

	end subroutine iterator
