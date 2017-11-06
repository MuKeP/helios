    subroutine iterator(callsub,energysub,maxiter,eps,dnull)

    use glob     , only: rglu,iglu,lglu,true,false,timeControl,gluCompare
    use printmod , only: prStrByVal
    use txtParser, only: tpFill,tpAdjustc
    use hdb      , only: iterationbd,ou

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

    integer(kind=iglu), intent(in) :: maxiter
    real   (kind=rglu), intent(in) :: eps
    logical(kind=lglu)             :: dnull

    character (len=1)              :: success
    integer(kind=iglu)             :: iteration,successfulIterations,printRate,headLen
    real   (kind=rglu)             :: accuracy(5),energy(5)
    real   (kind=rglu)             :: currentAccuracy,previousAccuracy,successPercent
    real   (kind=rglu)             :: improvementPercent,meanImprovement,saveEstimatedIters
    real   (kind=rglu)             :: estimatedIters,timePerIter
    real   (kind=rglu)             :: avgTimePerIter,estimatedTime,sTime(2),cTime(2),fTime(2)

    !ou=6

    sTime(1)=timeControl(sTime(2))
    iteration=0; successfulIterations=1
    currentAccuracy=0; previousAccuracy=10
    meanImprovement=0; printRate=1
    avgTimePerIter=0
    iterationbd%chkStopIteration(2)=false

    if (.NOT.dnull) call printInfo('header',ou)
    do
        iteration=iteration+1

        call callSub(iteration,eps,accuracy)
        call energysub(energy)

        cTime(1)=timeControl(cTime(2))

        timePerIter=(cTime(1)-sTime(1))-avgTimePerIter*(iteration-1)
        if (timePerIter.LT.gluCompare) timePerIter=0
        avgTimePerIter=(cTime(1)-sTime(1))/iteration

        previousAccuracy=currentAccuracy; currentAccuracy=maxval(accuracy)

        success='u'
        if (currentAccuracy.LT.previousAccuracy) then
            successfulIterations=successfulIterations+1
            success=' '
        endif

        if (iteration.EQ.1) then
            if (avgTimePerIter.GT.iterationbd%printFrequency) then
                printRate=1
            else
                if (avgTimePerIter.LT.10**6*gluCompare) then
                    printRate=iterationbd%printNotRearly
                else
                    printRate=int(iterationbd%printFrequency/avgTimePerIter)
                endif
            endif
            if (printRate.GT.iterationbd%printNotRearly) printRate=iterationbd%printNotRearly

            improvementPercent=0
            successPercent=100
            estimatedIters=maxiter
            estimatedTime=0
            success=' '
        else
            improvementPercent=100.*(previousAccuracy-currentAccuracy)/previousAccuracy
            successPercent=100.*(successfulIterations)/iteration
            if (improvementPercent.GT.1E-5_rglu) meanImprovement=meanImprovement+improvementPercent

            if (success.EQ.'u') then
                improvementPercent=abs(improvementPercent)

                if (improvementPercent.GT.999.99_rglu) improvementPercent=999.99_rglu
                estimatedIters=saveEstimatedIters
            else
                !estimatedIters=int( (log(eps)-log(currentAccuracy) )/log(1.-improvementPercent/100.) )
                estimatedIters=int( (log(eps)-log(currentAccuracy) )/log(1.-(meanImprovement/iteration)/100.) )

                if (estimatedIters.LT.0) estimatedIters=0

                saveEstimatedIters=estimatedIters
            endif

            if (estimatedIters.GT.maxiter) estimatedIters=maxiter-iteration
            estimatedTime=estimatedIters*avgTimePerIter
        endif

        !write (*,*) iteration,maxval(accuracy)

        if (iterationbd%chkStopIteration(2)) then
            iterationbd%chkStopIteration(2)=false
            if (.NOT.dnull) call printInfo('interrupted',ou)
            exit
        endif

        if (.NOT.dnull) then
            if (mod(iteration-1,printRate).EQ.0) then
                call printInfo('iteration',ou)
            endif
        endif

        if (currentAccuracy.LT.eps) then
            if (.NOT.dnull) call printInfo('converged',ou)
            exit
        endif

        if (iterationbd%chkDivergence) then
            if (currentAccuracy.GT.iterationbd%feelDivergence) then
                if (.NOT.dnull) call printInfo('diverged',ou)
                exit
            endif
        endif

        if (iteration+1.GT.maxiter) then
            if (.NOT.dnull) call printInfo('failed',ou)
            exit
        endif

        if (iterationbd%chkStagnation) then
            if (iteration.GT.iterationbd%thresholdStagnation) then
                if (successPercent.LE.iterationbd%feelStagnation) then
                    if (.NOT.dnull) call printInfo('stagnated',ou)
                    exit
                endif
            endif
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
                                !prStrByVal(avgTimePerIter,0,4,'exp'),&
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
                                            '%. Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headLen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('failed')
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Accuracy was not reached after '//prStrByVal(iteration)//' iterations'//&
                                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
                write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                            'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(meanImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headLen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('diverged')
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Procedure diverged and stopped after '//prStrByVal(iteration)//' iterations'//&
                                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
                write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                            'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(meanImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headlen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('stagnated')
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Procedure stagnated and stopped after '//prStrByVal(iteration)//' iterations'//&
                                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
                write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                            'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(meanImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headlen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('interrupted')
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Procedure was interrupted due to StopIteration after '//prStrByVal(iteration)//' iterations with',headLen)
                write (iou,'(A)') tpAdjustc('accuracy '//prStrByVal(currentAccuracy,0,2,'exp')//'. '//&
                                            'Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                            'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(meanImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headlen)
                write (iou,'(A)') tpFill(headLen,'=')

        end select

        return
        end subroutine printInfo

    end subroutine iterator
