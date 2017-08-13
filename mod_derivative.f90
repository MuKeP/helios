	module derivat

	use glob     , only: true,false,rglu,rspu,iglu,lglu,i4kind,r8kind,find
	use math     , only: LagrangeDerivative,factorial
	use txtParser, only: tpCount

	character (len=*), parameter :: deVersion='1.000'
	character (len=*), parameter :: deDate   ='2015.08.22'
	character (len=*), parameter :: deAuthor ='Anton B. Zakharov'

	character (len=*), parameter :: crt='xyz'

	private

	integer*4                         :: iounit=6
	integer(kind=i4kind)              :: dePnts
	real   (kind=r8kind)              :: deStep,lsmCondit
	real   (kind=r8kind), allocatable :: lsmReady(:,:)


	logical(kind=lglu)           :: lsmPrepared=false,parShared=false
	real(kind=rglu), allocatable :: im1(:),im2(:,:),im3(:,:,:)

	public :: deVersion,deDate,deAuthor,deShareParams,deLSMShareData,&
	&         deLagDeriv,deLSMDeriv,deFinalize

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	function deLagDeriv(nPoints,str) result(ret)
	implicit none

	real(kind=rglu)    :: ret

	integer(kind=iglu) :: nPoints
	integer(kind=iglu) :: icart(3),lcart(3),icart2(2),lcart2(2)
	integer(kind=iglu) :: sft,dPower,i,j,k,tp,aStart,aStop

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	function deLSMDeriv(str) result(ret)
	implicit none

	real(kind=rglu)    :: ret
	character (len=*)  :: str
	integer(kind=iglu) :: k,l,m,i,j,a,pshft


	if (.NOT.lsmPrepared) then; write (iounit,'(4X,A)') 'LSM data required. call deLSMShareData.'; ret=0; return; endif
	if (.NOT.parShared)   then; write (iounit,'(4X,A)') 'call deShareParams first.';               ret=0; return; endif

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer(kind=iglu) function deLSMShareData(con,ind,dat) result(ret)
	implicit none

	real(kind=r8kind)                 :: stafiel,stofiel,stefiel,a
	real(kind=r8kind), allocatable    :: fielList(:),conditional(:,:)

	integer(kind=i4kind)              :: stapnts,stopnts,stepnts,i
	integer(kind=i4kind), allocatable :: indexMatrix(:,:),pntsList(:)

	integer(kind=iglu)                :: nPoints,con,ind,dat
	integer(kind=iglu)                :: pshft,pntsCount,fielCount,curPnts,curFiel
	integer(kind=iglu)                :: lenpnts,lenfiel,tostartwith


	if (.NOT.parShared) then; write (iounit,'(4X,A)') 'call deShareParams first.'; ret=0; return; endif

	allocate (lsmReady(0:4,dePnts))
	if (eof(ind)) goto 101; read (ind,err=101) stapnts,stopnts,stepnts
	if (eof(ind)) goto 101; read (ind,err=101) stafiel,stofiel,stefiel

	lenpnts=   (stopnts-stapnts)/stepnts   +1
	lenfiel=int((stofiel-stafiel)/stefiel) +1
	
	allocate ( pntsList(lenpnts),fielList(lenfiel) )
	
	pntsCount=0; do i = stapnts,stopnts,stepnts; pntsCount=pntsCount+1; pntsList(pntsCount)=i; enddo
	fielCount=0; do a = stafiel,stofiel,stefiel; fielCount=fielCount+1; fielList(fielCount)=a; enddo
	curPnts=find(pntsList,dePnts); curFiel=find(fielList,deStep); deallocate (pntsList,fielList)
	
	if ((curPnts.LE.0).OR.(curFiel.LE.0)) goto 100

	allocate ( indexMatrix(lenpnts,lenfiel),conditional(lenpnts,lenfiel) )
		
	if (eof(ind)) goto 101; read (ind,err=101) indexMatrix; close (ind)
	if (eof(con)) goto 101; read (con,err=101) conditional; close (con)

	tostartwith=indexMatrix(curPnts,curFiel); lsmCondit=conditional(curPnts,curFiel)
	deallocate (indexMatrix,conditional)

	if (eof(dat)) goto 101; read (dat,rec=tostartwith,err=101) lsmReady; close(dat)

	lsmPrepared=true

	ret=0 ; return
100 ret=-1; return !combination not found.
101 ret=-2; return !eof during read.

	end function deLSMShareData

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

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

	ret=0; return
	end function deShareParams

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine deFinalize
	implicit none


	if (allocated(im1))      deallocate (im1)
	if (allocated(im2))      deallocate (im2)
	if (allocated(im3))      deallocate (im3)
	if (allocated(lsmReady)) deallocate (lsmReady)

	return
	end subroutine deFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module derivat