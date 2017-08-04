	subroutine planarityCheck

	use hdb
	use orientation

	implicit none

	real   (kind=rglu), allocatable :: X(:),Y(:),Z(:)
	integer(kind=iglu)              :: N,k,i1,i2,i3
	real   (kind=rglu)              :: sum,A,B,C,D


	N=mol%nAtoms
	geometrybd%searchPlanar(2)=false
	geometrybd%searchLinear(2)=false

	if (geometrybd%searchPlanar(1) .OR. geometrybd%searchLinear(1)) then
		allocate (X(N),Y(N),Z(N))
		X=mol%atm(:)%coords(1); Y=mol%atm(:)%coords(2); Z=mol%atm(:)%coords(3)
	endif

	if ( geometrybd%searchLinear(1) ) then
		if (N.EQ.2) then
			i1=1; i2=2; call setLinear(i1,i2); return
		endif

		i1=rangen(1,N)
		do
			i2=rangen(1,N)
			if (i1.NE.i2) exit
		enddo

		A=(Y(i2)-Y(i1))/(X(i2)-X(i1))
		B= Y(i1)+X(i1)* (Y(i2)-Y(i1))/(X(i2)-X(i1))

		sum=0.d0
		do k = 1,N
			sum=sum+abs(Y(k)-X(k)*A-B)
		enddo

		if (sum.LT.geometrybd%searchTolerance) then
			call setLinear(i1,i2); return
		endif
	endif

	if ( geometrybd%searchPlanar(1) ) then
		if (N.EQ.3) then
			i1=1; i2=2; i3=3; call setPlanar(i1,i2,i3); return
		endif

		do
			i1=rangen(1,N)
			do
				i2=rangen(1,N)
				if (i1.NE.i2) exit
			enddo
			do
				i3=rangen(1,N)
				if ((i1.NE.i3).AND.(i2.NE.i3)) exit
			enddo

			A=(Y(i2)-Y(i1))/(X(i2)-X(i1))
			B= Y(i1)+X(i1)* (Y(i2)-Y(i1))/(X(i2)-X(i1))

			! check whether three points lie on the straight line.
			sum=abs(Y(i3)-X(i3)*A-B)

			if (sum.GT.geometrybd%searchTolerance) exit
		enddo

		!
		!A*x+B*y+C*z+D=0
		!
		!  |1  y1  z1|	  |x1  1  z1|	  |x1  y1  1|	  |x1  y1  z1|
		!A=|1  y2  z2|	B=|x2  1  z2|	C=|x2  y2  1|	D=|x2  y2  z2|
		!  |1  y3  z3|	  |x3  1  z3|	  |x3  y3  1|	  |x3  y3  z3|
		!

		A=Y(i1)*(Z(i2)-Z(i3))+Y(i2)*(Z(i3)-Z(i1))+Y(i3)*(Z(i1)-Z(i2))
		B=Z(i1)*(X(i2)-X(i3))+Z(i2)*(X(i3)-X(i1))+Z(i3)*(X(i1)-X(i2))
		C=X(i1)*(Y(i2)-Y(i3))+X(i2)*(Y(i3)-Y(i1))+X(i3)*(Y(i1)-Y(i2))

		D= X(i1)*(Y(i2)*Z(i3)-Y(i3)*Z(i2))&
		& +X(i2)*(Y(i3)*Z(i1)-Y(i1)*Z(i3))&
		& +X(i3)*(Y(i1)*Z(i2)-Y(i2)*Z(i1))

		sum=0.d0
		do k = 1,N
			sum=sum+abs(X(k)*A+Y(k)*B+Z(k)*C+D)
		enddo

		if (sum.LT.geometrybd%searchTolerance) then
			call setPlanar(i1,i2,i3); return
		endif
	endif

	return

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine setLinear(i1,i2)
		implicit none
		
		real   (kind=rglu) :: sumY,sumZ
		integer(kind=iglu) :: i1,i2,k


		sumY=0; sumZ=0
		do k = 1,N
			sumY=sumY+Y(k)**2
			sumZ=sumZ+Z(k)**2
		enddo

		if (max(sumY,sumZ).GT.geometrybd%searchTolerance) then
			call shareCoords(X,Y,Z)
			call orSetIOunit(unul)
			call putOnAxis(i1,i2)
			call getNewCoords(X,Y,Z)
		endif

		geometrybd%searchLinear(2)=true
		polarizbd%scales     =uchSet('x')
		hyperchargesbd%scales=uchSet('x')

		deallocate (X,Y,Z)

		return
		end subroutine setLinear

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine setPlanar(i1,i2,i3)
		implicit none
		real   (kind=rglu) :: sumZ
		integer(kind=iglu) :: i1,i2,i3,k


		sumZ=0
		do k = 1,N
			sumZ=sumZ+Z(k)**2
		enddo

		if (sumZ.GT.geometrybd%searchTolerance) then
			call shareCoords(X,Y,Z)
			call orSetIOunit(unul)
			call putOnPlane(i1,i2,i3)
			call getNewCoords(X,Y,Z)
		endif

		geometrybd%searchPlanar(2)=true
		polarizbd%scales     =uchSet('xy')
		hyperchargesbd%scales=uchSet('xy')

		deallocate (X,Y,Z)

		return
		end subroutine setPlanar

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end subroutine planarityCheck