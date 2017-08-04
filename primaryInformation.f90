	subroutine primaryInformation(action)

	use hdb

	implicit none

	character (len=*) :: action
	character (len=*), parameter :: months(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
	                                             'Jul','Aug','Sep','Oct','Nov','Dec'/)
	character (len=*), parameter :: days(7)   =(/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)

	integer*4                    :: values(8),p1,p2


	stop
	values=0; call date_and_time(values=values)

	!VALUE(1): The year
	!VALUE(2): The month
	!VALUE(3): The day of the month
	!VALUE(4): Time difference with UTC in minutes
	!VALUE(5): The hour of the day
	!VALUE(6): The minutes of the hour
	!VALUE(7): The seconds of the minute
	!VALUE(8): The milliseconds of the second

	select case (action)
		case ('init')
			write (ou,100) days(dayOfWeek(values(1),values(2),values(3))),&
			               values(1),months(values(2)),values(3),values(5:7)
			100 format (<ouIndent>X,'Execution of HELIOS started ',A3,1X,i4,'-',A3,'-',i2.2,2X,i2.2,':',i2.2,':',i2.2)

			write (ou,101)

			101 format (<ouIndent>X,'')


		case ('end')


		case ('error')

	end select

	!write (*,*) values

	stop




	! Date
	! Authors
	! Version/Date
	! Description

	!open (init,file=uchGet(generalbd%fname))
	!call prEchoFile(src,trgt,marker)

	return

!	contains







	end subroutine primaryInformation