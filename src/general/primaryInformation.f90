    subroutine primaryInformation(action)

    use glob,       only: iglu,lglu,uch,date_time,mid,i4kind
    use glob,       only: getPath,convertTime,timeControl,isPosix,getPOSIXinfo
    use hdb,        only: ouWidth,ou,eu,heVersion,heDate,heCompilDate
    use hdb,        only: appPid,timeSpent,ierror,finalizeHelios,ipFailed
    use txtParser,  only: tpAdjustc,tpFill
    use argsParser, only: apGetCommandLine
    use printmod,   only: prLongText,prStrByVal

    implicit none

    character (len=*)    :: action
    type(uch)            :: commandline,apppath,login,hostname
    integer(kind=i4kind) :: memory,ppid
    logical(kind=lglu)   :: opened


    !ou=6
    select case (action)
        case ('init')
            write (ou,'(/A)') tpAdjustc('HELIOS program by Anton B. Zakharov and Vladimir V. Ivanov',ouWidth)
            write (ou,'(A/)') tpAdjustc('V.N. Karazin Kharkiv National University',ouWidth)
            call prLongText('Implemented CUE (Covalently Unbonded Ethylene molecues) approach is based on the classical'       //&
                            ' representation of pi-conjugated system as the set of signle and double bonds. The main'          //&
                            ' purpose of such approach is to decrease computational complexity of coupled cluster methods'     //&
                            ' with minor accuracy loss. For brief introduction to the features provided by the CUE basis'      //&
                            ' read "A.B. Zakharov and V.V. Ivanov, J. Struct. Chem. (Engl.Transl.) 52, 645 (2011)". For'       //&
                            ' detailed description read "Anton B. Zakharov, Vladimir V. Ivanov and Ludwik Adamowicz. '         //&
                            ' Optical Parameters of pi-Conjugated Oligomer Chains from the Semiempirical Local Coupled-Cluster'//&
                            ' Theory // Practical Aspects of Computational Chemistry IV J. Leszczynski, M. K. Shukla (Eds.).'  //&
                            ' Springer Science+Business Media, New York, 2016. Chapter 3, P. 57-102."'&
                            ,ou,'center','justified',70,ouWidth,' ')
            write (ou,*)
            call prLongText('For details about LR-cue-ccsd method read "A. B. Zakharov, V. V. Ivanov, and L. Adamowicz pi-Electron'//&
                            ' Calculations Using the Local Linear-Response Coupled-Cluster Singles and Doubles Theory // Journal '//&
                            ' of Physical Chemistry C. 2015. Vol. 119, N. 52. P. 28737-28748."'&
                            ,ou,'center','justified',70,ouWidth,' ')

            write (ou,'(/A)') tpAdjustc('Core HELIOS version 02-July-2009.',ouWidth)
            write (ou,'( A)') tpAdjustc('Current version '//heVersion//'. last modified '//heDate//', compiled '//heCompilDate//'.',ouWidth)

            commandline=apGetCommandLine()
            apppath=getPath()

            write (ou,'(/A)') tpAdjustc('Execution of HELIOS started '//date_time(),ouWidth)
            write (ou,'(/A)') tpFill(ouWidth,'/\') !' sublime enjoyes escaping quotes

            if (isPosix) then
                call getPOSIXinfo(hostname,login,memory,ppid)

                write (ou,'(/A)') 'Executed:         '//commandline%get()
                write (ou,'(A)' ) 'PID/PPID:         '//trim(prStrByVal(appPID))//'/'//trim(prStrByVal(ppid))
                write (ou,'(A)' ) 'login@host:       '//login%get()//'@'//hostname%get()
                write (ou,'(A)' ) 'Available memory: '//trim(adjustl(prStrByVal(memory,fmt='000000')))//' MB'
                write (ou,'(A/)') 'Working path:     '//apppath%get()
            else
                write (ou,'(/A)') 'Executed:     '//commandline%get()
                write (ou,'(A)' ) 'PID:          '//trim(prStrByVal(appPID))
                write (ou,'(A/)') 'Working path: '//apppath%get()
            endif

        case ('end')
            if (ipFailed) write (ou,'(/A/)') tpAdjustc('Attention! For one or more iteration procedures desired accuracy was not reached.', ouWidth)

            write (ou,'(/A)') tpFill(ouWidth,'/\') !' sublime enjoyes escaping quotes
            timeSpent(1,2)=timeControl(timeSpent(2,2))
            write (ou,'(/A)') tpAdjustc('Proper termination of HELIOS '//date_time(),ouWidth)
            write (ou,'( A)') tpAdjustc('Total time spent '//trim(convertTime(timeSpent(1,2)-timeSpent(1,1)))//&
                                        ' ('//trim(prStrByVal(timeSpent(1,2)-timeSpent(1,1),'^.00'))//' seconds )',ouWidth)
            write (ou,'( A)') tpAdjustc('Total CPU utilization '//prStrByVal(100*(timeSpent(2,2)&
                                        -timeSpent(2,1))/(timeSpent(1,2)-timeSpent(1,1)),4,2)//'%',ouWidth)
            write (ou,*)
            call finalizeHelios
            stop

        case ('error')
            ! =====> need to write pre-ambule <=====
            write (eu,'(2X,A)') 'During execution following error occured:'
            write (eu,'(2X,A)') ierror%get()

            inquire(ou, opened=opened)
            if ((ou.NE.eu).AND.(opened)) then
                write (ou,'(2X,A)') 'During execution following error occured:'
                write (ou,'(2X,A)') ierror%get()
            endif
            call finalizeHelios
            stop

    end select


    return
    end subroutine primaryInformation