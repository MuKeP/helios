	subroutine primaryInformation(action)

	use glob      , only: iglu,uch,uchGet,uchSet,date_time,mid
	use glob      , only: getPath,convertTime,timeControl
	use hdb       , only: ouWidth,ou,heVersion,heDate,heCompilDate
	use hdb       , only: appPid,timeSpent
	use txtParser , only: tpAdjustc,tpFill
	use argsParser, only: apGetCommandLine
	use printmod  , only: prLongText,prStrByVal

	implicit none

	character (len=*) :: action


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
			                ,ou,'center','justified',64,ouWidth,' ')
			write (ou,*)
			call prLongText('For details about LR-cue-ccsd method read "A. B. Zakharov, V. V. Ivanov, and L. Adamowicz pi-Electron'//&
			                ' Calculations Using the Local Linear-Response Coupled-Cluster Singles and Doubles Theory // Journal '//&
							' of Physical Chemistry C. 2015. Vol. 119, N. 52. P. 28737-28748."'&
							,ou,'center','justified',64,ouWidth,' ')

			write (ou,'(/A)') tpAdjustc('Core HELIOS version 02-July-2009.',ouWidth)
			write (ou,'( A)') tpAdjustc('Current version '//heVersion//'. last modified '//heDate//', compiled '//heCompilDate//'.',ouWidth)

			write (ou,'(/A)') tpAdjustc('Execution of HELIOS started '//date_time(),ouWidth)
			write (ou,'(/A)') tpFill(ouWidth,'/\')
			write (ou,'(/A,A,A,i<mid(appPID)>)')  'Executed:     ',uchGet(apGetCommandLine()),'with PID:',appPID
			write (ou,'(A,A/)') 'Working path: ',uchGet(getPath())

		case ('end')
			write (ou,'(/A)') tpFill(ouWidth,'/\')
			timeSpent(1,2)=timeControl(timeSpent(2,2))
			write (ou,'(/A)') tpAdjustc('Proper termination of HELIOS '//date_time(),ouWidth)
			write (ou,'( A)') tpAdjustc('Total time spent '//convertTime(timeSpent(1,2)-timeSpent(1,1)),ouWidth)
			write (ou,'( A)') tpAdjustc('Total CPU utilization '//prStrByVal(100*(timeSpent(2,2)&
			                            -timeSpent(2,1))/(timeSpent(1,2)-timeSpent(1,1)),4,2)//'%',ouWidth)
			stop

		case ('error')
			stop

	end select


	return
	end subroutine primaryInformation