	module approx

	use glob     , only: rglu,iglu,lglu,rspu,true,false
	use math     , only: tred4,gaussSLE
	use txtParser, only: tpLowerCase

	character (len=*), parameter :: axVersion='1.100'
	character (len=*), parameter :: axDate   ='2017.05.23'
	character (len=*), parameter :: axAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	! glu = real for global use  (storage and non-accuracy-demanding procedures).
	! spu = real for special use (accuracy-demanding procedures).
	!integer*4, parameter :: rglu=r8kind, rspu=r8kind
	!integer*4, parameter :: iglu=i4kind, ispu=i8kind
	!integer*4, parameter :: lglu=l1kind

	real(kind=rglu), parameter :: eValTol=real(1d-16,rglu),conditTol=real(1d+10,rglu)

	!real(kind=rglu), parameter :: gluzero=real(0,rglu)
	!real(kind=rspu), parameter :: spuzero=real(0,rspu)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	character                    :: task*3,meth*2
	logical(kind=lglu)           :: doOutput
	integer(kind=iglu)           :: iount
	integer(kind=iglu)           :: N,M

	real(kind=rglu), allocatable :: X(:),Y(:),Yth(:),kof(:), A(:,:),B(:,:),C(:,:)

	private
	public :: axVersion,axDate,axAuthor,axFinalize,approximate

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		real(kind=rglu) function approximate(gX,gY,gM,gtask,gkof,giount,gmeth) result(ret)
		implicit none

		integer(kind=iglu), intent(in)  :: gM
		real(kind=rglu)   , intent(in)  :: gX(:),gY(:)
		real(kind=rglu)   , intent(out) :: gkof(gM)

		character(len=3)  , intent(in)  :: gtask
		character(len=2)  , optional    :: gmeth
		integer(kind=iglu), optional    :: giount

		integer(kind=iglu)              :: i,j,k

		real(kind=rspu)                 :: condit,numer,sum1,sum2
		real(kind=rglu), allocatable    :: eVec(:,:),eVal(:)


		doOutput=false
		if (present(giount)) then
			iount=giount
			doOutput=true
		endif

		task=tpLowerCase(gtask)
		N=UBound(gX,1); M=gM

		meth='sd'; if (present(gmeth)) meth=tpLowerCase(gmeth)

		allocate (X(N),Y(N),Yth(N),A(M,N),B(M,M),C(M,M),kof(M))
		X=0; Y=0; Yth=0; A=0
		B=0; C=0; kof=0

		A(1,:)=1
		Y=gY; X=gX

		select case (task)
			case ('lin')
				do i = 1,N
					A(2,i)=X(i)
				enddo

			case ('lg ')
				do i = 1,N
					A(2,i)=log10(X(i))
				enddo

			case ('ln ')
				do i = 1,N
					A(2,i)=log(X(i))
				enddo

			case ('exp')
				do i = 1,N
					A(2,i)=X(i)
					Y(i)=log(Y(i))
				enddo

			case ('pol')
				do i = 1,N
					do j = 2,M
						A(j,i)=X(i)**(j-1)
					enddo
				enddo

			case ('pow')
				do i = 1,N
					A(2,i)=log(X(i))
					Y(i)=log(Y(i))
				enddo

			case default
				ret=-1; return

		end select

		do i = 1,M
		do j = 1,M
			sum1=0
			do k = 1,N
				sum1=sum1+A(i,k)*A(j,k)
			enddo
			B(i,j)=sum1
		enddo
		enddo

		allocate (eVec(M,M),eVal(M))

		select case (meth)
			case ('sd')
				call tred4 (B,eVec,eVal,M,real(1d-100,rspu),real(1d-100,rspu))
				do i = 1,M
				do j = 1,M
					sum1=0
					do k = 1,M
						if (abs(eVal(k)).LT.eValTol) cycle
						sum1=sum1+eVec(i,k)*eVec(j,k)/eVal(k)
					enddo
					C(i,j)=sum1
				enddo
				enddo
				condit=maxval(abs(eVal))/minval(abs(eVal))

			case ('jg')
				condit=gaussSLE(M,B,C)
		end select

		do i = 1,M
			sum1=0
			do j = 1,M
				sum2=0
				do k = 1,N
					sum2=sum2+Y(k)*A(j,k)
				enddo
				sum1=sum1+sum2*C(i,j)
			enddo
			kof(i)=sum1
		enddo

		select case (task)
			case ('exp')
				do i = 1,N
					Yth(i)=exp(kof(1))*exp(X(i)*kof(2))
				enddo
				do i = 1,N
					Y(i)=exp(Y(i))
				enddo

			case ('pow')
				do i = 1,N
					Yth(i)=exp(kof(1))*X(i)**kof(2)
				enddo
				do i = 1,N
					Y(i)=exp(Y(i))
				enddo

			case default
				do i = 1,N
					sum1=0
					do j = 1,M
						sum1=sum1+A(j,i)*kof(j)
					enddo
					Yth(i)=sum1
				enddo

		end select

		numer=0
		do i = 1,N
			numer=numer+(Yth(i)-Y(i))**2
		enddo
		sum1=0
		sum2=0
		do i = 1,N
			sum1=sum1+Y(i)**2
			sum2=sum2+Y(i)
		enddo
		ret=1-numer/(sum1-sum2**2/N)

		gkof=kof

		deallocate (eVec,eVal,X,Y,Yth,A,B,C,kof)

		return
		end function approximate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine axFinalize
		implicit none


		!continue

		return
		end subroutine axFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module approx
