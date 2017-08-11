	module painter

	use glob, only: rglu,iglu,dtr,pi

	character (len=*), parameter :: tdVersion='2.100'
	character (len=*), parameter :: tdDate   ='2017.05.23'
	character (len=*), parameter :: tdAuthor ='Anton B. Zakharov'

	real(kind=rglu), allocatable :: thX(:),thY(:),thZ(:)

	contains

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine centrcord(N,X,Y,Z)

	implicit none

	integer(kind=iglu) :: N,i
	real(kind=rglu)    :: X(N),Y(N),Z(N),sumCoordinates(3)


	sumCoordinates=0
	do i = 1,N
		sumCoordinates(1)=sumCoordinates(1)+X(i)
		sumCoordinates(2)=sumCoordinates(2)+Y(i)
		sumCoordinates(3)=sumCoordinates(3)+Z(i)
	enddo; sumCoordinates=sumCoordinates/N

	do i = 1,N
		X(i)=X(i)-sumCoordinates(1)
		Y(i)=Y(i)-sumCoordinates(2)
		Z(i)=Z(i)-sumCoordinates(3)
	enddo

	return
	end subroutine centrcord

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function xyz (singleX,singleY,singleZ,Ncurr)
	implicit none

	real(kind=rglu)    :: singleX,singleY,singleZ
	integer(kind=iglu) :: Ncurr

	Ncurr=Ncurr+1
	
	if (isEveryThingFine(Ncurr).EQ.-1) then; xyz=-1; return; endif

	thX(Ncurr)=singleX; thY(Ncurr)=singleY; thZ(Ncurr)=singleZ

	xyz=Ncurr; return
	end function xyz

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function sft(refAtom,singleX,singleY,singleZ,Ncurr)
	implicit none

	real(kind=rglu)    :: singleX,singleY,singleZ
	integer(kind=iglu) :: refAtom,Ncurr

	Ncurr=Ncurr+1
	
	if (isEveryThingFine(Ncurr)  .EQ.-1) then; sft=-1; return; endif
	if (isEveryThingFine(refAtom).EQ.-1) then; sft=-1; return; endif

	thX(Ncurr)=thX(refAtom)+singleX; thY(Ncurr)=thY(refAtom)+singleY; thZ(Ncurr)=thZ(refAtom)+singleZ

	sft=Ncurr; return
	end function sft

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function plr(refAtom,alph,radi,plane,Ncurr)
	implicit none

	integer(kind=iglu) :: refAtom,Ncurr
	real(kind=rglu)    :: alph,radi,singleX,singleY,singleZ
	character          :: plane*2


	Ncurr=Ncurr+1
	
	if (isEveryThingFine(Ncurr)  .EQ.-1) then; plr=-1; return; endif
	if (isEveryThingFine(refAtom).EQ.-1) then; plr=-1; return; endif

	singleX=thX(refAtom)
	singleY=thY(refAtom)
	singleZ=thZ(refAtom)

	alph=-alph

	thX=thX-singleX; thY=thY-singleY; thZ=thZ-singleZ

	select case (plane)

		case ('xy')
		thX(Ncurr)=-radi*sin(alph*dtr); thY(Ncurr)=radi*cos(alph*dtr); thZ(Ncurr)=0

		case ('xz')
		thX(Ncurr)=-radi*sin(alph*dtr); thY(Ncurr)=0; thZ(Ncurr)=radi*cos(alph*dtr)

		case ('yz')
		thX(Ncurr)=0; thY(Ncurr)=-radi*sin(alph*dtr); thZ(Ncurr)=radi*cos(alph*dtr)

	end select

	thX=thX+singleX; thY=thY+singleY; thZ=thZ+singleZ

	plr=Ncurr; return
	end function plr

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function rem(remAtom,Ncurr)
	implicit none

	integer(kind=iglu) :: i,remAtom,Ncurr


	if (isEveryThingFine(remAtom).EQ.-1) then; rem=-1; return; endif

	do i = remAtom+1,min(size(thX),size(thY),size(thZ))
		thX(i-1)=thX(i); thY(i-1)=thY(i); thZ(i-1)=thZ(i)
	enddo
							
	Ncurr=Ncurr-1

	rem=Ncurr; return
	end function rem

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function rot (alph,plane,Ncurr)
	implicit none

	real(kind=rglu)    :: sinAlphRad,cosAlphRad,ttx,tty,alph
	character          :: plane*2
	integer(kind=iglu) :: i,Ncurr


	if (isEveryThingFine(Ncurr)   .EQ.-1) then; rot=-1; return; endif

	call centrcord (Ncurr,thX(1:Ncurr),thY(1:Ncurr),thZ(1:Ncurr))

	sinAlphRad=sin(alph*dtr)
	cosAlphRad=cos(alph*dtr)

	select case (plane)

		case ('xy')
		do i = 1,Ncurr
			ttx=thX(i)*cosAlphRad-thY(i)*sinAlphRad
			tty=thX(i)*sinAlphRad+thY(i)*cosAlphRad
			thX(i)=ttx
			thY(i)=tty
		enddo

		case ('xz')
		do i = 1,Ncurr
			ttx=thX(i)*cosAlphRad-thZ(i)*sinAlphRad
			tty=thX(i)*sinAlphRad+thZ(i)*cosAlphRad
			thX(i)=ttx
			thZ(i)=tty
		enddo

		case ('yz')
		do i = 1,Ncurr
			ttx=thY(i)*cosAlphRad-thZ(i)*sinAlphRad
			tty=thY(i)*sinAlphRad+thZ(i)*cosAlphRad
			thY(i)=ttx
			thZ(i)=tty
		enddo
		case default; rot=-1; return

	end select

	rot=Ncurr; return
	end function rot

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function rpl(refAtom,ref2Atom,radi,alphread,plane,Ncurr)
	implicit none

	integer(kind=iglu) :: pl,sig,refAtom,ref2Atom,Ncurr
	real(kind=rglu)    :: alph,alphread,radi,singleX,singleY,singleZ
	character          :: plane*2


	Ncurr=Ncurr+1

	if (isEveryThingFine(Ncurr)   .EQ.-1) then; rpl=-1; return; endif
	if (isEveryThingFine(refAtom) .EQ.-1) then; rpl=-1; return; endif
	if (isEveryThingFine(ref2Atom).EQ.-1) then; rpl=-1; return; endif


	!new   ! refAtom = aa; ref2Atom = jj
	!  \
	!    \
	!      jj -----aa
	! rotate relativly to aa around jj by angle (positive cw; negative ccw)


	singleX=thX(ref2Atom); singleY=thY(ref2Atom); singleZ=thZ(ref2Atom)
	thX=thX-singleX;         thY=thY-singleY;         thZ=thZ-singleZ

	sig=1
	select case (plane)
		case ('xy'); pl=1; if (thX(refAtom).GT.thX(ref2Atom)) sig=-1
		case ('xz'); pl=2; if (thZ(refAtom).GT.thZ(ref2Atom)) sig=-1
		case ('yz'); pl=3; if (thY(refAtom).GT.thY(ref2Atom)) sig=-1
		case default; rpl=-1; return
	end select

	alph=-alphread+dble(sig)*ang(refAtom,ref2Atom,pl)

	select case (plane)
		case ('xy'); thX(Ncurr)=-radi*sin(alph*dtr); thY(Ncurr)=radi*cos(alph*dtr); thZ(Ncurr)=0
		case ('xz'); thX(Ncurr)=-radi*sin(alph*dtr); thY(Ncurr)=0; thZ(Ncurr)=radi*cos(alph*dtr)
		case ('yz'); thX(Ncurr)=0; thY(Ncurr)=-radi*sin(alph*dtr); thZ(Ncurr)=radi*cos(alph*dtr)
		case default; rpl=-1; return
	end select

	thX=thX+singleX; thY=thY+singleY; thZ=thZ+singleZ

	rpl=Ncurr; return
	end function rpl

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function ini(polygon_n,ccbond,Ncurr)
	implicit none

	real(kind=rglu)    :: ccbond,angtorot,polygon_rad,prevX,prevY,currX,currY
	integer(kind=iglu) :: polygon_n,Ncurr,i


	thX(1:size(thX))=0
	thY(1:size(thY))=0
	thZ(1:size(thZ))=0
	Ncurr=polygon_n

	if (isEveryThingFine(Ncurr).EQ.-1) then; ini=-1; return; endif

	angtorot=dtr*real(360,rglu)/polygon_n
	polygon_rad=ccbond/sin(angtorot)*sin(dtr*(real(180,rglu)-real(360,rglu)/polygon_n)/real(2,rglu)) ! sin theorem

	thX(1)=0; thY(1)=polygon_rad; thZ(1)=0

	prevX=0
	prevY=polygon_rad

	do i = 2,polygon_n
		currX=+prevX*cos(angtorot)+prevY*sin(angtorot)
		currY=-prevX*sin(angtorot)+prevY*cos(angtorot)

		thX(i)=currX; thY(i)=currY

		prevX=currX
		prevY=currY
	enddo

	ini=Ncurr; return
	end function ini

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function isEveryThingFine(atom)
	implicit none

	integer(kind=iglu) :: atom

	if (&
	&	(atom.GT.min(size(thX),size(thY),size(thZ))).OR.&
	&	(atom.LT.1)&
	&						) then
		isEveryThingFine=-1
		return
	endif

	isEveryThingFine=0
	end function isEveryThingFine

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu) function ang(aa,jj,pl)

	! Привет, великий Рендом. Ты прибывал со мной во время написания этой хуйни.
	! п.с. сука, я не понимаю как оно работает. это сука гребаный подбор.
	! берите книгу "Как нихуя не понять".
	! ИксИгрикЙ его знает как для других плоскостей. я использую только xy.
	!
	!    zero direction: COORD(pl)=COORD(pl)+1.0d0   others COORD=COORD;
	!    i.e. delta(COORD(pl))=-1.0d0   others delta(COORD)=0.d0
	!
	! x1 y1 z1 соответствуют направлению орта, относительно которого будет
	! отсчитываться угол.
	! xy - орт j.
	! xz - орт i.
	! yz - орт k.
	! x1 y1 z1 - длины орта вдоль каждой из осей.
	!
	! Y
	! ^  zero direction
	! |  | 
	! |  |
	! | jj
	! |   \
	! |     \
	! |       aa
	!  ----------> X
	!

	implicit none

	integer(kind=iglu) :: jj,aa,pl
	real(kind=rglu)    :: x1,x2,y1,y2,z1,z2,zna,chis,cos_ang,x1s,y1s,z1s
		
	select case (pl)
		case (1); x1= 0; y1=real(-1,rglu); z1= 0 ! xy
		case (2); x1=real(-1,rglu); y1= 0; z1= 0 ! xz
		case (3); x1= 0; y1= 0; z1=real(-1,rglu) ! yz
	end select

	x2=thX(jj)-thX(aa); y2=thY(jj)-thY(aa); z2=thZ(jj)-thZ(aa)

	x1s=x1*x1; y1s=y1*y1; z1s=z1*z1

	zna  = sqrt( (x1s+y1s+z1s) * (x2*x2+y2*y2+z2*z2) )
	chis = x1*x2+y1*y2+z1*z2
	cos_ang = chis/zna
	if (abs(cos_ang) .ge. real(1,rglu))  cos_ang = sign (real(1,rglu),cos_ang)

	ang=acos ( cos_ang )/dtr

	return
	end function ang

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module painter