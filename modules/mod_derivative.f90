	module derivat

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	use glob     , only: true,false,rglu,rspu,iglu,lglu,r16kind,find
	use math     , only: LagrangeDerivative,factorial,sptred4
	use txtParser, only: tpCount

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: deVersion='1.000'
	character (len=*), parameter :: deDate   ='2015.08.22'
	character (len=*), parameter :: deAuthor ='Anton B. Zakharov'

	character (len=*), parameter :: crt='xyz'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu), allocatable :: im1(:),im2(:,:),im3(:,:,:),lsmReady(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) :: iounit=6,dePnts
	real   (kind=rglu) :: deStep,lsmCondit
	logical(kind=lglu) :: lsmPrepared=false,parShared=false

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private
	public :: deVersion,deDate,deAuthor
	public :: deShareParams,deLSMShareData,deLagDeriv,deLSMDeriv,deFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function deLagDeriv(nPoints,str) result(ret)
	implicit none

	real(kind=rglu)    :: ret

	integer(kind=iglu) :: nPoints
	integer(kind=iglu) :: icart(3),lcart(3),icart2(2),lcart2(2)
	integer(kind=iglu) :: i,j,k,tp,aStart,aStop

	character (len=*)  :: str


	if (.NOT.parShared) then; write (iounit,'(4X,A)') 'call deShareParams first.'; ret=0; return; endif

	icart=0; do i = 1,3; icart(i)=tpCount(str,crt(i:i)); enddo
	if (sum(icart).GT.nPoints-1) then; ret=0; return; endif

	aStart=-(nPoints-1)/2; aStop=(nPoints-1)/2

	lcart=0; do i = 1,3; if (icart(i).NE.0) lcart(i)=1;    enddo
	tp=0   ; do i = 1,3; if (lcart(i).NE.0) tp=tp+1;       enddo

	select case (tp)
		case(1)
			do i = 1,3; if (lcart(i).NE.0) j=i; enddo

		case(2)
			do i = 1,3; if (lcart(i).EQ.0) j=i; enddo

			k=0
			do i = 1,3
				if (lcart(i).NE.0) then
					k=k+1
					lcart2(k)=lcart(i)
					icart2(k)=icart(i)
				endif
			enddo
		case(3)
			continue
	end select

	select case (tp) !
		case (1)
			select case (j)
				case (1); ret=LagrangeDerivative(nPoints,icart(1),im3(aStart:aStop,0,0),real(deStep,rspu))
				case (2); ret=LagrangeDerivative(nPoints,icart(2),im3(0,aStart:aStop,0),real(deStep,rspu))
				case (3); ret=LagrangeDerivative(nPoints,icart(3),im3(0,0,aStart:aStop),real(deStep,rspu))
			end select

		case (2)
			do i = aStart,aStop
				select case (j)
					case (1); im1(i)=LagrangeDerivative(nPoints,icart2(2),im3(0,i,aStart:aStop),real(deStep,rspu))
					case (2); im1(i)=LagrangeDerivative(nPoints,icart2(2),im3(i,0,aStart:aStop),real(deStep,rspu))
					case (3); im1(i)=LagrangeDerivative(nPoints,icart2(2),im3(i,aStart:aStop,0),real(deStep,rspu))
				end select
			enddo; ret=LagrangeDerivative(nPoints,icart2(1),im1(aStart:aStop),real(deStep,rspu))

		case (3)
			do i = aStart,aStop
				do j = aStart,aStop
					im2(i,j)=LagrangeDerivative(nPoints,icart(3),im3(i,j,aStart:aStop),real(deStep,rspu))
				enddo
			enddo
			do i = aStart,aStop
				im1(i)=LagrangeDerivative(nPoints,icart(2),im2(i,aStart:aStop),real(deStep,rspu))
			enddo; ret=LagrangeDerivative(nPoints,icart(1),im1,real(deStep,rspu))

	end select

	return
	end function deLagDeriv

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	function deLSMDeriv(str) result(ret)
	implicit none

	real   (kind=rglu) :: ret
	character  (len=*) :: str
	integer(kind=iglu) :: k,l,m,i,j,a,pshft


	!if (.NOT.lsmPrepared) then; write (iounit,'(4X,A)') 'LSM data required. call deLSMShareData.'; ret=0; return; endif
	!if (.NOT.parShared)   then; write (iounit,'(4X,A)') 'call deShareParams first.';               ret=0; return; endif

	k=tpCount(str,crt(1:1)); l=tpCount(str,crt(2:2)); m=tpCount(str,crt(3:3))
		
	pshft=(dePnts-1)/2+1
	ret=0
	do i = 1,dePnts
	do j = 1,dePnts
	do a = 1,dePnts
		ret=ret+lsmReady(k,i)*lsmReady(l,j)*lsmReady(m,a)*im3(i-pshft,j-pshft,a-pshft)
	enddo
	enddo
	enddo; ret=ret*factorial(k)*factorial(l)*factorial(m) !*transition(k+l+m)

	return
	end function deLSMDeriv

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine prepareLSM
	implicit none

	integer(kind=iglu), parameter  :: maxpower=5 ! max power of polynomial (a1+a2*x+a3*x^2...)

	integer(kind=iglu)             :: points,i,j,k
	real   (kind=rspu)             :: field,sum

	real   (kind=rspu), dimension (:,:), allocatable :: lsKof,lsWork,lsInv,lsEVec
	real   (kind=rspu), dimension   (:), allocatable :: lsRez,lsEVal


	field=deStep; points=dePnts

	allocate (lsKof(points,maxpower),lsWork(maxpower,maxpower),lsInv(maxpower,maxpower))
	allocate (lsEVec(maxpower,maxpower),lsEVal(maxpower))
	allocate (lsRez(maxpower))
	lsKof=0; lsWork=0; lsInv=0; lsEVal=0; lsEVec=0; lsRez=0

	do k = -(points-1)/2,(points-1)/2
		lsKof(k+points/2+1,2)=k*field
	enddo

	lsKof(:,1)=1
	do k = 3,maxpower
		do i = 1,points
			lsKof(i,k)=lsKof(i,2)**real(k-1,rspu) ! does matter
		enddo
	enddo

	do i = 1,maxpower
	do j = 1,maxpower
		sum=0
		do k = 1,points
			sum=sum+lsKof(k,i)*lsKof(k,j)
		enddo
		lsWork(i,j)=sum
	enddo
	enddo

	call sptred4(lsWork,lsEVec,lsEVal,maxpower,real(1.d-100,rspu),real(1.d-300,rspu))

	lsmCondit=maxval(abs(lsEVal))/minval(abs(lsEVal))
	do i = 1,maxpower
	do j = 1,maxpower
		sum=0
		do k = 1,maxpower
			sum=sum+lsEVec(i,k)*lsEVec(j,k)/lsEVal(k)
		enddo
		lsInv(i,j)=sum
	enddo
	enddo
	
	do i = 1,maxpower
		do j = 1,points
			sum=0
			do k = 1,maxpower
				sum=sum+lsInv(i,k)*lsKof(j,k)
			enddo
			lsmReady(i-1,j)=sum
		enddo
	enddo

	deallocate (lsKof,lsWork,lsInv,lsEVec,lsEVal,lsRez)
	
	return
	end subroutine prepareLSM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) function deShareParams(ppoints,pstep,valSet) result(ret)
	implicit none

	integer(kind=iglu) :: ppoints,aStart,aStop,sft
	real   (kind=rglu) :: pstep,valSet(:,:,:)


	deStep=pstep; dePnts=ppoints
	aStart=-(dePnts-1)/2; aStop=(dePnts-1)/2; sft=(UBound(valSet,1)-1)/2+1

	allocate (im1(aStart:aStop),&
	&         im2(aStart:aStop,aStart:aStop),&
	&         im3(aStart:aStop,aStart:aStop,aStart:aStop))

	im3=valSet(aStart+sft:aStop+sft,aStart+sft:aStop+sft,aStart+sft:aStop+sft)

	parShared=true

	if (rspu.EQ.r16kind) then
		allocate (lsmReady(0:4,dePnts)); call prepareLSM
	endif

	ret=0; return
	end function deShareParams

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine deFinalize
	implicit none


	if (allocated(im1))      deallocate (im1)
	if (allocated(im2))      deallocate (im2)
	if (allocated(im3))      deallocate (im3)
	if (allocated(lsmReady)) deallocate (lsmReady)

	return
	end subroutine deFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module derivat