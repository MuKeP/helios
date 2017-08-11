	module orientation

	use glob, only: rglu,iglu,lglu,false,true,void,gluUnity
	use glob, only: definePi,pi,dtr

	character (len=*), parameter :: orVersion='1.302'
	character (len=*), parameter :: orDate   ='2015.07.04'
	character (len=*), parameter :: orAuthor ='Anton B. Zakharov'

	type cart
		real(kind=rglu) :: X,Y,Z
	end type
	type(cart)          :: prevCentre

	real(kind=rglu), allocatable :: oX(:),oY(:),oZ(:)

	integer(kind=iglu) :: N,iounit=6

	private

	public  :: orVersion,orDate,orAuthor
	public  :: shareCoords,getNewCoords,putOnPlane,putOnAxis,orSetIOunit

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine shareCoords(X,Y,Z)
		implicit none
		real(kind=rglu), intent(in) :: X(:),Y(:),Z(:)


        void=definePi()
		N=UBound(X,1); allocate ( oX(-2:N),oY(-2:N),oZ(-2:N) )

		oX=0; oX(1:N)=X(1:N)
		oY=0; oY(1:N)=Y(1:N)
		oZ=0; oZ(1:N)=Z(1:N)

		return
		end subroutine shareCoords

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine getNewCoords(X,Y,Z)
		implicit none
		real(kind=rglu), intent(out) :: X(:),Y(:),Z(:)


		X=oX(1:N); Y=oY(1:N); Z=oZ(1:N)
		deallocate (oX,oY,oZ)

		return
		end subroutine getNewCoords

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine putOnPlane(i,j,k)
		implicit none
		integer(kind=iglu), intent(in) :: i,j,k
		integer(kind=iglu)             :: acNum,kk
		real(kind=rglu)                :: a


		write (iounit,100)
		do kk = 1,N
			write (iounit,101) oX(kk),oY(kk),oZ(kk)
		enddo
		write (iounit,102)

		acNum=1
		! translate j to the O | X(j)=0; Y(j)=0; Z(j)=0
		write (iounit,150) acNum,-oX(j),-oY(j),-oZ(j)
		call recentre(j)


		acNum=acNum+1
		! put i on XZ
		oX(-2)=oX(i)               ; oY(-2)=oY(i); oZ(-2)=0 
		oX(-1)=sign(gluUnity,oX(i)); oY(-1)=0    ; oZ(-1)=0
		a=ang(-2,j,-1); call rotateMol(sign(a,oX(i)*oY(i)),'xy')
		write (iounit,151) acNum,a/dtr,'z'

		acNum=acNum+1
		! put i on Ox
		oX(-2)=oX(i)               ; oY(-2)=0; oZ(-2)=oZ(i)
		oX(-1)=sign(gluUnity,oX(i)); oY(-1)=0; oZ(-1)=0
		a=ang(-2,j,-1); call rotateMol(sign(a,-oX(i)*oZ(i)),'xz') 
		write (iounit,151) acNum,a/dtr,'y'

		acNum=acNum+1
		! put k on XY
		oX(-2)=0; oY(-2)=oY(k)               ; oZ(-2)=oZ(k)
		oX(-1)=0; oY(-1)=sign(gluUnity,oY(k)); oZ(-1)=0
		a=ang(-2,j,-1); call rotateMol(sign(a,oY(k)*oZ(k)),'yz')
		write (iounit,151) acNum,a/dtr,'x'

		if (oX(i).LT.0) then ! orient i-j to the positive direction of X.
			acNum=acNum+1
			call rotateMol(pi,'xy')
			write (iounit,151) acNum,real(180,kind=rglu),'z'
		endif

		if (oY(k).LT.0) then ! orient to set order i-j-k clockwise
			acNum=acNum+1
			call rotateMol(pi,'yz')
			write (iounit,151) acNum,real(180,kind=rglu),'x'
		endif

		acNum=acNum+1
		call centrate
		write (iounit,150) acNum,prevCentre%X,prevCentre%Y,prevCentre%Z

		write (iounit,103)
		do kk = 1,N
			write (iounit,101) oX(kk),oY(kk),oZ(kk)
		enddo
		write (iounit,104)

	100 format (//3X,'=========================ReOrientation procedure========================='//&
	&			 24X,'Coordinates before ReOrientation'/)
	101 format ( 20X,3(F11.5,2X))
	102 format (/32X,'Procedure steps:')
	103 format (/24X,'Coordinates  after ReOrientation'/)
	104 format (  3X,'========================================================================='//)
	150	format (  9X,i1,') Shift on X: ',F10.4,'; on Y: ',F10.4,'; on Z: ',F10.4)
	151 format (  9X,i1,') Rotaion on ',F7.2,' degrees around O',A1)

		return
		end subroutine putOnPlane

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!
!		integer*4 function showMOL(unt)
!		implicit none
!		integer*4 unt,i
!
!		open (unt,file='currentPOS.xyz')
!		write (unt,*) N
!		do i = 1,N
!			write (unt,100) oX(i),oY(i),oZ(i)
!		enddo
!		close (unt)
!		call system ('start Visual.exe '//'currentPOS.xyz')
!
!	100 format ( 'C',2X,3(F20.12,2X))
!
!		showMOL=0; return
!		end function showMOL
!
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine putOnAxis(i,j)
		implicit none
		integer(kind=iglu), intent(in) :: i,j
		integer(kind=iglu)             :: acNum,kk
		real(kind=rglu)                :: a


		write (iounit,100)
		do kk = 1,N
			write (iounit,101) oX(kk),oY(kk),oZ(kk)
		enddo
		write (iounit,102)

		acNum=1
		! translate centre between i and j to the O
		oX(0)=(oX(i)+oX(j))/real(2,kind=rglu)
		oY(0)=(oY(i)+oY(j))/real(2,kind=rglu)
		oZ(0)=(oZ(i)+oZ(j))/real(2,kind=rglu)

		call recentre(0)
		write (iounit,150) acNum,-oX(0),-oY(0),-oZ(0)
		oX(0)=0; oY(0)=0; oZ(0)=0

		acNum=acNum+1
		! put i on XZ
		oX(-2)=oX(i); oY(-2)=oY(i); oZ(-2)=0 
		oX(-1)=oX(i); oY(-1)=0    ; oZ(-1)=0
		a=ang(-2,0,-1); call rotateMol(sign(a,oX(i)*oY(i)),'xy')
		write (iounit,151) acNum,a/dtr,'z'

		acNum=acNum+1
		! put i on Ox
		oX(-2)=oX(i); oY(-2)=oY(i); oZ(-2)=oZ(i)
		oX(-1)=oX(i); oY(-1)=0    ; oZ(-1)=0
		a=ang(-2,0,-1); call rotateMol(sign(a,-oX(i)*oZ(i)),'xz') 
		write (iounit,151) acNum,a/dtr,'y'

		if (oX(i).LT.oX(j)) then ! orient i to the positive direction of X.
			acNum=acNum+1
			call rotateMol(pi,'xy')
			write (iounit,151) acNum,real(180,kind=rglu),'z'
		endif

		acNum=acNum+1
		call centrate
		write (iounit,150) acNum,prevCentre%X,prevCentre%Y,prevCentre%Z

		write (iounit,103)
		do kk = 1,N
			write (iounit,101) oX(kk),oY(kk),oZ(kk)
		enddo
		write (iounit,104)

	100 format (//3X,'=========================ReOrientation procedure========================='//&
	&			 24X,'Coordinates before ReOrientation'/)
	101 format ( 20X,3(F11.5,2X))
	102 format (/32X,'Procedure steps:')
	103 format (/24X,'Coordinates  after ReOrientation'/)
	104 format (  3X,'========================================================================='//)
	150	format (  9X,i1,') Shift on X: ',F10.4,'; on Y: ',F10.4,'; on Z: ',F10.4)
	151 format (  9X,i1,') Rotaion on ',F7.2,' degrees around O',A1)

		return
		end subroutine putOnAxis

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		real(kind=rglu) function ang(ii,jj,kk) result(ret)
		implicit none
		integer(kind=iglu), intent(in) :: ii,jj,kk
		real(kind=rglu)                :: x1,y1,z1,x2,y2,z2,zna,chis,cos_ang


		x1=oX(jj)-oX(ii)
		y1=oY(jj)-oY(ii)
		z1=oZ(jj)-oZ(ii)

		x2=oX(jj)-oX(kk)
		y2=oY(jj)-oY(kk)
		z2=oZ(jj)-oZ(kk)

		zna  = sqrt( (x1*x1+y1*y1+z1*z1) * (x2*x2+y2*y2+z2*z2) )
		chis = x1*x2+y1*y2+z1*z2
		if (zna.LT.real(1d-15,kind=rglu)) then
			zna=gluUnity
		endif
		cos_ang=chis/zna
		if ( abs(cos_ang).GE.gluUnity)  cos_ang=sign(gluUnity,cos_ang)
		
		ret=acos( cos_ang ); return
		end function ang

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine rotateMol(angInRad,cha)
		implicit none
		real(kind=rglu)  , intent(in) :: anginRad
		character (len=2), intent(in) :: cha
		real(kind=rglu)               :: tx,ty
		integer(kind=iglu)            :: kk


		select case (cha) ! clockwise

			case('xy')
			do kk = 1,N
				tx= oX(kk)*cos(angInRad)+oY(kk)*sin(angInRad)
				ty=-oX(kk)*sin(angInRad)+oY(kk)*cos(angInRad)
				oX(kk)=tx
				oY(kk)=ty
			enddo

			case ('xz')
			do kk = 1,N
				tx= oX(kk)*cos(angInRad)-oZ(kk)*sin(angInRad)
				ty= oX(kk)*sin(angInRad)+oZ(kk)*cos(angInRad)
				oX(kk)=tx
				oZ(kk)=ty
			enddo

			case ('yz')
			do kk = 1,N
				tx= oY(kk)*cos(angInRad)+oZ(kk)*sin(angInRad)
				ty=-oY(kk)*sin(angInRad)+oZ(kk)*cos(angInRad)
				oY(kk)=tx
				oZ(kk)=ty
			enddo

		end select

		return
		end subroutine rotateMol

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine rotateMe(anginRad,me,cha)
		implicit none
		real(kind=rglu)  , intent(in)  :: anginRad
		character (len=2), intent(in)  :: cha
		integer(kind=iglu), intent(in) :: me
		real(kind=rglu)                :: tx,ty


		select case (cha) ! clockwise

			case('xy')
			tx= oX(Me)*cos(angInRad)+oY(Me)*sin(angInRad)
			ty=-oX(Me)*sin(angInRad)+oY(Me)*cos(angInRad)
			oX(Me)=tx
			oY(Me)=ty

			case ('xz')
			tx= oX(Me)*cos(angInRad)-oZ(Me)*sin(angInRad)
			ty= oX(Me)*sin(angInRad)+oZ(Me)*cos(angInRad)
			oX(Me)=tx
			oZ(Me)=ty

			case ('yz')
			tx= oY(Me)*cos(angInRad)+oZ(Me)*sin(angInRad)
			ty=-oY(Me)*sin(angInRad)+oZ(Me)*cos(angInRad)
			oY(Me)=tx
			oZ(Me)=ty

		end select

		return
		end subroutine rotateMe

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine recentre(ii)
		implicit none
		integer(kind=iglu), intent(in) :: ii


		oX(1:N)=oX(1:N)-oX(ii)
		oY(1:N)=oY(1:N)-oY(ii)
		oZ(1:N)=oZ(1:N)-oZ(ii)

		return
		end subroutine recentre

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
		
		subroutine centrate()
		implicit none
		integer(kind=iglu) :: i


		prevCentre%X=0
		prevCentre%Y=0
		prevCentre%Z=0
		do i = 1,N
			prevCentre%X=prevCentre%X+oX(i)
			prevCentre%Y=prevCentre%Y+oY(i)
			prevCentre%Z=prevCentre%Z+oZ(i)
		enddo

		prevCentre%X=prevCentre%X/N
		prevCentre%Y=prevCentre%Y/N
		prevCentre%Z=prevCentre%Z/N

		do i = 1,N
			oX(i)=oX(i)-prevCentre%X
			oY(i)=oY(i)-prevCentre%Y
			oZ(i)=oZ(i)-prevCentre%Z
		enddo

		return
		end subroutine centrate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine orSetIOunit(iunt)
		implicit none

		integer*4, intent(in) :: iunt
		integer*4             :: unt


		unt=iunt
		if (unt.EQ.5) unt=6 !stdin  (forbidden)
		if (unt.EQ.0) unt=6 !stdout (for compatibility)
		iounit=unt

		return
		end subroutine orSetIOunit

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module orientation