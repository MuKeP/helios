    subroutine iterator(callsub,energysub,maxiter,eps,dnull,callback,docallback,converged)

    use glob,      only: rglu,iglu,lglu,true,false,timeControl,gluCompare
    use printmod,  only: prStrByVal
    use txtParser, only: tpFill,tpAdjustc
    use hdb,       only: iterationbd,ou,iheader,ipFailed,resetIterationState

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

    interface
        subroutine callback
        end subroutine callback
    end interface

    integer(kind=iglu), intent(in)  :: maxiter
    real   (kind=rglu), intent(in)  :: eps
    logical(kind=lglu), intent(in)  :: dnull,docallback
    logical(kind=lglu), intent(out) :: converged

    character (len=1)              :: success
    integer(kind=iglu)             :: iteration,lastPrinted,successfulIterations,printRate,headLen
    real   (kind=rglu)             :: accuracy(5),energy(5)
    real   (kind=rglu)             :: currentAccuracy,previousAccuracy,successPercent
    real   (kind=rglu)             :: improvementPercent,cummImprovement,saveEstimatedIters
    real   (kind=rglu)             :: estimatedIters,timePerIter
    real   (kind=rglu)             :: avgTimePerIter,estimatedTime,sTime(2),cTime(2),fTime(2)


    !ou=6

    iterationbd%doRestart=true
    iterationbd%lastProcedureStatus=0

    do

        sTime(1)=timeControl(sTime(2))
        iteration=0; successfulIterations=1
        currentAccuracy=0; previousAccuracy=10
        cummImprovement=0; printRate=1
        avgTimePerIter=0
        iterationbd%chkStopIteration(2)=false
        lastPrinted=1

        converged=false

        call performIterationProcedure

        if ((iterationbd%lastProcedureStatus.LE.-1).AND.&
            (iterationbd%lastProcedureStatus.GE.-3)) then
            if (docallback) then
                call callback; if (iterationbd%doRestart) cycle
            endif
        endif

        exit
    enddo

    if (iterationbd%lastProcedureStatus.LE.0) then
        ipFailed=true
    endif

    call resetIterationState

    return

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine performIterationProcedure

        implicit none


        call printInfo('header',ou)
        do
            iteration=iteration+1

            call callSub(iteration,eps,accuracy)
            call energysub(energy)

            cTime(1)=timeControl(cTime(2))

            if (iterationbd%afterPause) then
                timePerIter=avgTimePerIter
                iterationbd%afterPause=false
            else
                timePerIter=(cTime(1)-sTime(1))-avgTimePerIter*(iteration-1)
            endif

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
                        printRate=iterationbd%printNotRarely
                    else
                        printRate=int(iterationbd%printFrequency/avgTimePerIter)
                    endif
                endif
                if (printRate.GT.iterationbd%printNotRarely) printRate=iterationbd%printNotRarely

                improvementPercent=0
                successPercent=100
                estimatedIters=maxiter
                saveEstimatedIters=maxiter
                estimatedTime=0
                success=' '
            else
                improvementPercent=100.*(previousAccuracy-currentAccuracy)/previousAccuracy
                successPercent=100.*(successfulIterations)/iteration
                if (improvementPercent.GT.1E-5_rglu) cummImprovement=cummImprovement+improvementPercent

                if (success.EQ.'u') then
                    ! try
                    !improvementPercent=100.*(previousAccuracy-currentAccuracy)/previousAccuracy
                    improvementPercent=100.*(previousAccuracy-currentAccuracy)/currentAccuracy

                    improvementPercent=abs(improvementPercent)

                    if (improvementPercent.GT.999.99_rglu) improvementPercent=999.99_rglu
                    estimatedIters=saveEstimatedIters
                else
                    !estimatedIters=int( (log(eps)-log(currentAccuracy) )/log(1.-improvementPercent/100.) )
                    estimatedIters=int( (log(eps)-log(currentAccuracy) )/log(1.-(cummImprovement/iteration)/100.) )

                    if (estimatedIters.LT.0) estimatedIters=0

                    saveEstimatedIters=estimatedIters
                endif

                if (estimatedIters.GT.maxiter) estimatedIters=maxiter-iteration
                estimatedTime=estimatedIters*avgTimePerIter
            endif

            if (iterationbd%chkStopIteration(2)) then
                iterationbd%chkStopIteration(2)=false
                if (lastPrinted.NE.iteration) call printInfo('iteration',ou)
                call printInfo('interrupted',ou)
                exit
            endif

            if (mod(iteration-1,printRate).EQ.0) then
                lastPrinted=iteration
                call printInfo('iteration',ou)
            endif

            if (currentAccuracy.LT.eps) then
                if (lastPrinted.NE.iteration) call printInfo('iteration',ou)
                call printInfo('converged',ou)
                converged=true
                exit
            endif

            if (iterationbd%chkDivergence) then
                if (currentAccuracy.GT.iterationbd%feelDivergence) then
                    if (lastPrinted.NE.iteration) call printInfo('iteration',ou)
                    call printInfo('diverged',ou)
                    exit
                endif
            endif

            if (iteration+1.GT.maxiter) then
                            if (lastPrinted.NE.iteration) call printInfo('iteration',ou)
                call printInfo('failed',ou)
                exit
            endif

            if (iterationbd%chkStagnation) then
                if (iteration.GT.iterationbd%thresholdStagnation) then
                    if (successPercent.LE.iterationbd%feelStagnation) then
                        if (lastPrinted.NE.iteration) call printInfo('iteration',ou)
                        call printInfo('stagnated',ou)
                        exit
                    endif
                endif
            endif
        enddo

        return
        end subroutine performIterationProcedure

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine printInfo(string,iou)
        implicit none

        character (len=*)  :: string
        integer(kind=iglu) :: i,iou
        character (len=1)  :: separator


        headLen=79
        select case (string)
            case ('header')
                iterationbd%lastProcedureStatus=0

                if (dnull) return

                if (iheader%ln.GT.0) then
                    write (iou,'(/A)') tpFill(headLen,'=')
                    write (iou,'( A)') tpAdjustc(iheader%get(),headlen,'*')
                    write (iou,'( A)') tpFill(headLen,'=')
                else
                    write (iou,'(/A)') tpFill(headLen,'=')
                endif
                write (iou, '(A)') tpAdjustc('1 - Iteration, 2 - Accuracy, 3 - Step efficiency, 4 - Time per iteration',headLen)
                if (iterationbd%energy) then
                    write (iou, '(A)') tpAdjustc('5 - Energy, 6 - Estimated iterations left, 7 - Estimated time left',headLen)
                else
                    write (iou, '(A)') tpAdjustc('5 - Fitness, 6 - Estimated iterations left, 7 - Estimated time left',headLen)
                endif
                write (iou, '(A)') tpFill(headLen,'=')
                write (iou, '(2X,i1,2X,"|",6X,i1,6X,"|",4X,i1,3X,"|",6X,i1,5X,"|",9X,i1,9X,"|",2X,i1,2X,"|",5X,i1,5X)' ) (i, i = 1,7)
                write (iou, '(A)') tpFill(headLen,'=')

            case ('iteration')
                iterationbd%lastProcedureStatus=0

                if (dnull) return
                separator=" "
                write (iou,'(A,A1,A,A1,A,1X,A,A1,A,1X,A1,A,A1,A,A1,A)') &
                                tpAdjustc(trim(adjustl(prStrByVal(iteration,5))),5),separator,&
                                prStrByVal(currentAccuracy,0,6,'exp'),separator,&
                                success,prStrByVal(improvementPercent,3,2),separator,&
                                prStrByVal(timePerIter,0,4,'exp'),separator,&
                                !prStrByVal(avgTimePerIter,0,4,'exp'),&
                                prStrByVal(energy(1),0,12,'exp'),separator,&
                                tpAdjustc(trim(adjustl(prStrByVal(min(int(estimatedIters),maxiter-iteration+1),5))),5),separator,&
                                prStrByVal(estimatedTime,0,4,'exp')

            case ('converged')
                iterationbd%lastProcedureStatus=1

                if (dnull) return
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Procedure converged after '//prStrByVal(iteration)//'('//prStrByVal(successfulIterations)//&
                                            ') iterations'//' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)

                if (iterationbd%energy) then
                    write (iou,'(A)') tpAdjustc('Resulting energy '//prStrByVal(energy(1),0,15,'exp')//' eV. '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                else
                    write (iou,'(A)') tpAdjustc('Resulting fitness '//prStrByVal(energy(1),0,15,'exp')//' '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                endif
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency'//prStrByVal(cummImprovement/iteration,3,2)//&
                                            '%. Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headLen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('failed')
                iterationbd%lastProcedureStatus=-1

                if (dnull) return
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Accuracy was not reached after '//prStrByVal(iteration)//' iterations'//&
                                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
                if (iterationbd%energy) then
                    write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                else
                    write (iou,'(A)') tpAdjustc('Fitness '//prStrByVal(energy(1),0,13,'exp')//' '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                endif
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(cummImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headLen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('diverged')
                iterationbd%lastProcedureStatus=-2

                if (dnull) return
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Procedure diverged and stopped after '//prStrByVal(iteration)//' iterations'//&
                                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
                if (iterationbd%energy) then
                    write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                else
                    write (iou,'(A)') tpAdjustc('Fitness '//prStrByVal(energy(1),0,13,'exp')//' '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                endif
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(cummImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headlen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('stagnated')
                iterationbd%lastProcedureStatus=-3

                if (dnull) return
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Procedure stagnated and stopped after '//prStrByVal(iteration)//' iterations'//&
                                            ' with accuracy'//prStrByVal(currentAccuracy,0,2,'exp')//'.',headLen)
                if (iterationbd%energy) then
                    write (iou,'(A)') tpAdjustc('Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                else
                    write (iou,'(A)') tpAdjustc('Fitness '//prStrByVal(energy(1),0,13,'exp')//' '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                endif
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(cummImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headlen)
                write (iou,'(A)') tpFill(headLen,'=')

            case ('interrupted')
                iterationbd%lastProcedureStatus=-4

                if (dnull) return
                fTime(1)=timeControl(fTime(2))

                write (iou,'(A)') tpFill(headLen,'=')
                write (iou,'(A)') tpAdjustc('Procedure is interrupted due to StopIteration on '//prStrByVal(iteration)//' iteration at accuracy',headLen)
                if (iterationbd%energy) then
                    write (iou,'(A)') tpAdjustc(prStrByVal(currentAccuracy,0,2,'exp')//'. '//&
                                                'Energy '//prStrByVal(energy(1),0,13,'exp')//' eV. '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                else
                    write (iou,'(A)') tpAdjustc(prStrByVal(currentAccuracy,0,2,'exp')//'. '//&
                                                'Fitness '//prStrByVal(energy(1),0,13,'exp')//' '//&
                                                'Total '//trim(adjustl(prStrByVal(fTime(1)-sTime(1),10,2)))//' seconds spent.',headLen)
                endif
                write (iou,'(A)') tpAdjustc('Avg iteration efficiency '//trim(adjustl(prStrByVal(cummImprovement/iteration,3,2)))//'%. '//&
                                            'Avg iteration time '//trim(adjustl(prStrByVal(avgTimePerIter,8,5)))//' seconds.',headlen)
                write (iou,'(A)') tpFill(headLen,'=')

        end select

        return
        end subroutine printInfo

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine iterator
