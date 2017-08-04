	!fcNewID(stat)      ! get new file id
	!fcBanID(ID)        ! reserve file id
	!fcUnBanID(ID)      ! free file id
	!fcNullID(ID,[act]) ! release file id (act='delete')
	!fcCloseFiles()     ! close all files

	module fcontrol

	character (len=*), parameter :: fcVersion='2.101'
	character (len=*), parameter :: fcDate   ='2017.07.07'
	character (len=*), parameter :: fcAuthor ='Anton B. Zakharov'

	logical*1, parameter                           :: true=.true., false=.false.
	integer*4, parameter                           :: fileListLen=256

	character (len=*), parameter, dimension (-1:5) :: statname=['standard',             & !-1
                                                                'free',                 & ! 0
                                                                'reserved',             & ! 1
                                                                'reserved connected',   & ! 2
                                                                'unreserved connected', & ! 3
                                                                'banned',               & ! 4
                                                                'banned connected']       ! 5

	integer*4            :: fileList(0:fileListLen)=0
	integer*4            :: listShift
	integer*4            :: fcIO,iounit,void

	private

	public :: fcVersion,fcDate,fcAuthor

	public :: fcSetListShift,fcNewID,fcBanID,fcUnBanID,fcNullID,&
	          fcCloseFiles,fcSetIOunit,fcIO,fileList,fcFinalize

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLOBAL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		pure function lcase(str) result(outstr)
		implicit none; character (Len=*), intent(in) :: str
		character (len=len(str))      :: outstr; 
		integer*4                     :: i,ascii
		do i = 1,len(str); ascii=iachar(str(i:i))
			if ((ascii.GE.65).AND.(ascii.LE.90)) then; outstr(i:i)=char(ascii+32)
			else; outstr(i:i)=char(ascii); endif
		enddo; return
		end function lcase
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~END GLOBAL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function fcSetListShift(shift) result(ret)
		implicit none

		integer*4 :: shift

		ret=0
		if (shift.LT.0)           then; ret=-1; return; endif
		if (shift.GE.fileListLen) then; ret=-2; return; endif

		listShift=shift

		return
		end function fcSetListShift

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function fcNewID(stat) result(ret)
		implicit none

		integer*4, optional, intent(out) :: stat
		integer*4                        :: i
		logical*1                        :: isOpened


		void=checkID()
		if (present(stat)) stat=0
		do i = listShift,fileListLen
			
			if (fileList(i).EQ.0) then
				inquire(i,opened=isOpened)
				if (.NOT.isOpened) then
					fileList(i)=1; ret=i; return
				endif
			endif

		enddo

		if (present(stat)) stat=-1
		ret=-1; return
		end function fcNewID

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function fcBanID(ID,exist) result(ret)
		implicit none
		integer*4 :: ID
		logical*1, optional :: exist
		logical*1           :: dexist


		dexist=false; if (present(exist)) dexist=exist
		void=checkID()
		if (ID.LT.0)           then ; ret=-1; return; endif !<zero
		if (ID.GT.fileListLen) then ; ret=-2; return; endif !>available

		select case (fileList(ID))
			case (-1); ret=-3; return !standard
			case ( 0); ret=ID; fileList(ID)=4
			case ( 1)
				if (dexist) then
					ret=ID; fileList(ID)=4
					return 
				else
					ret=-4; return !already reserved
				endif
			case ( 2)
				if (dexist) then
					ret=ID; fileList(ID)=5
					return
				else
					ret=-5; return !already reserved
				endif
			case ( 3); ret=-6; return !already connected
			case ( 4); ret=-7; return !already banned
			case ( 5); ret=-7; return !already banned
		end select

		return
		end function fcBanID

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function fcUnBanID(ID) result(ret)
		implicit none
		integer*4 :: ID
		logical*1 :: isOpened


		void=checkID()
		if (ID.LT.0)           then ; ret=-1; return; endif !<zero
		if (ID.GT.fileListLen) then ; ret=-2; return; endif !>available

		inquire(ID,opened=isOpened)
		select case (fileList(ID))
			case (-1); continue
			case ( 0); continue
			case ( 1); continue
			case ( 2); continue
			case ( 3); continue
			case ( 4); fileList(ID)=0
			case ( 5); fileList(ID)=3
		end select

		ret=ID; return
		end function fcUnBanID

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function fcNullID(ID,act) result(ret)
		implicit none
		integer*4                   :: iost,ID
		character (len=6), optional :: act
		character (len=6)           :: dact
		logical*1                   :: isOpened


		void=checkID()
		if (ID.LE.0)           then; ret=-1; return; endif
		if (ID.EQ.5)           then; ret=-1; return; endif
		if (ID.EQ.6)           then; ret=-1; return; endif
		if (ID.GT.fileListLen) then; ret=-2; return; endif

		dact='keep'
		if (present(act)) then
			if (lcase(act).EQ.'delete') then
				dact=act
			endif
		endif

		iost=0
		select case (fileList(ID))
			case (-1); iost=-3; continue !unable to close
			case ( 0); iost=-4; continue !it is not connected
			case ( 1); iost=-5; continue !it is not connected
			case ( 2); close (ID,status=dact,iostat=iost); fileList(ID)=0
			case ( 3); iost=-6; continue !it is not reserved
			case ( 4); iost=-7; continue !it is not connected
			case ( 5); iost=-8; continue !it is banned
		end select
		inquire(ID,opened=isOpened)
		ret=iost

		return
		end function fcNullID

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function fcCloseFiles() result(ret)
		implicit none
		integer*4 :: i
		logical*1 :: isOpened


		void=checkID()

		do i = 0,fileListLen
			inquire(i,opened=isOpened); if (isOpened) void=fcNullID(i)
		enddo

		ret=0; return
		end function fcCloseFiles

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine fcSetIOunit(iunt)
		implicit none

		integer*4, intent(in) :: iunt
		integer*4             :: unt


		unt=iunt
		if (unt.EQ.5) unt=6 !stdin  (forbidden)
		if (unt.EQ.0) unt=6 !stdout (for compatibility)
		iounit=unt

		return
		end subroutine fcSetIOunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer*4 function checkID() result(ret)
		implicit none
		integer*4 :: i
		logical*1 :: isOpened


		! file status
		!	-1 - standard
		!	 0 - free
		!	 1 - reserved
		!	 2 - reserved and connected
		!	 3 - unreserved and connected
		!	 4 - banned
		!	 5 - banned and connected

		fileList(0)=-1
		fileList(5)=-1
		fileList(6)=-1

		do i = 0,fileListLen
			inquire(i,opened=isOpened)

			select case (fileList(i))
				case (-1); continue
				case ( 0); if (isOpened) fileList(i)=3
				case ( 1); if (isOpened) fileList(i)=2
				case ( 2); continue
				case ( 3); continue
				case ( 4); if (isOpened) fileList(i)=5
				case ( 5); continue
			end select
			!if (i.LT.15) write (*,*) i,fileList(i),isOpened
		enddo

		ret=0; return
		end function checkID

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine fcFinalize
		implicit none





		return
		end subroutine fcFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!
!		integer*4 function showFileIDList(fID)
!		implicit none
!		integer*4, optional    :: fID
!		integer*4, parameter   :: stdlen=20
!		integer*4              :: i
!
!		character (len=stdlen) :: acc, act, frm, pos, shr
!		character (len=512)    :: nam
!		logical*1              :: opnd
!		integer*4              :: iost
!
!
!		void=checkFileIDs()
!
!		acc   =repeat(' ',stdlen)
!		act   =repeat(' ',stdlen)
!		frm   =repeat(' ',stdlen)
!		pos   =repeat(' ',stdlen)
!		shr   =repeat(' ',stdlen)
!		nam   =repeat(' ',512)
!
!		if (present(fID)) then
!			i=fID
!
!			if (i.LT.0)           then; showFileIDList=-1; return; endif
!			if (i.GT.fileListLen) then; showFileIDList=-1; return; endif
!
!			if (fileIDs(i).EQ.0) then
!				write (iounit,101) i
!			101	format ('File ID:',1X,i<mid(i)>,1X,'is free.')
!				showFileIDList=0
!				return
!			endif
!
!			inquire(i,access=acc, action=act, opened=opnd, form=frm, &
!			&		  iostat=iost, name=nam, position=pos, share=shr)
!
!			write (iounit,100) i,trim(nam),trim(acc),trim(act),speakCorrectly(opnd),&
!			&				   trim(frm),iost,trim(pos),trim(shr),statname(fileIDs(i))
!		else
!		
!			do i = 0,maxUsed
!				if (fileIDs(i).NE.0) then
!
!					inquire(i,access=acc, action=act, opened=opnd, form=frm, &
!					&		  iostat=iost, name=nam, position=pos, share=shr)
!
!					write (iounit,100) i,trim(nam),trim(acc),trim(act),speakCorrectly(opnd),&
!					&				   trim(frm),iost,trim(pos),trim(shr),statname(fileIDs(i))
!				endif
!			enddo
!		
!		end if
!
!100		format ('File ID:       ',1X,i<mid(i)>/&
!&				'Path:          ',1X,A<len_trim(nam)>/&
!&				'Access:        ',1X,A<len_trim(acc)>/&
!&				'Action:        ',1X,A<len_trim(act)>/&
!&				'Connected:     ',1X,A<len_trim(speakCorrectly(opnd))>/&
!&				'Form:          ',1X,A<len_trim(frm)>/&
!&				'Iostat:        ',1X,i<mid(iost)>/&
!&				'Position:      ',1X,A<len_trim(pos)>/&
!&				'Share:         ',1X,A<len_trim(shr)>/&
!&				'Internal state:',1X,A/)
!
!		showFileIDList=0; return
!		end function showFileIDList
!
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module fcontrol