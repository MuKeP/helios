    module derivat

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,      only: true,false,rglu,rspu,iglu,lglu,r16kind,find,glControlMemory,void,i8kind
    use math,      only: LagrangeDerivative,factorial,tred4
    use txtParser, only: tpCount

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: deVersion='1.200'
    character (len=*), parameter :: deDate   ='2018.08.05'
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
    public :: deShareParams,deLSMShareData,deLagDeriv,deLSMDeriv,deFinalize,crt

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

    select case (tp)
        case (1)
            select case (j)
                case (1); ret=LagrangeDerivative(nPoints,icart(1),im3(aStart:aStop,0,0),deStep)
                case (2); ret=LagrangeDerivative(nPoints,icart(2),im3(0,aStart:aStop,0),deStep)
                case (3); ret=LagrangeDerivative(nPoints,icart(3),im3(0,0,aStart:aStop),deStep)
            end select

        case (2)
            do i = aStart,aStop
                select case (j)
                    case (1); im1(i)=LagrangeDerivative(nPoints,icart2(2),im3(0,i,aStart:aStop),deStep)
                    case (2); im1(i)=LagrangeDerivative(nPoints,icart2(2),im3(i,0,aStart:aStop),deStep)
                    case (3); im1(i)=LagrangeDerivative(nPoints,icart2(2),im3(i,aStart:aStop,0),deStep)
                end select
            enddo
            ret=LagrangeDerivative(nPoints,icart2(1),im1(aStart:aStop),deStep)

        case (3)
            do i = aStart,aStop
                do j = aStart,aStop
                    im2(i,j)=LagrangeDerivative(nPoints,icart(3),im3(i,j,aStart:aStop),deStep)
                enddo
            enddo
            do i = aStart,aStop
                im1(i)=LagrangeDerivative(nPoints,icart(2),im2(i,aStart:aStop),deStep)
            enddo
            ret=LagrangeDerivative(nPoints,icart(1),im1,deStep)

    end select

    return
    end function deLagDeriv

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    function deLSMDeriv(str) result(ret)
    implicit none

    real   (kind=rglu) :: ret
    character  (len=*) :: str
    integer(kind=iglu) :: k,l,m,i,j,a,pshft,switch


    !if (.NOT.lsmPrepared) then; write (iounit,'(4X,A)') 'LSM data required. call deLSMShareData.'; ret=0; return; endif
    !if (.NOT.parShared)   then; write (iounit,'(4X,A)') 'call deShareParams first.';               ret=0; return; endif

    k=tpCount(str,crt(1:1)); l=tpCount(str,crt(2:2)); m=tpCount(str,crt(3:3))

    pshft=(dePnts-1)/2+1

    switch=0
    if ((k.NE.0).AND.(l.EQ.0).AND.(m.EQ.0)) switch=1 ! x
    if ((k.EQ.0).AND.(l.NE.0).AND.(m.EQ.0)) switch=2 ! y
    if ((k.EQ.0).AND.(l.EQ.0).AND.(m.NE.0)) switch=3 ! z

    if ((k.NE.0).AND.(l.NE.0).AND.(m.EQ.0)) switch=4 ! xy
    if ((k.NE.0).AND.(l.EQ.0).AND.(m.NE.0)) switch=5 ! xz
    if ((k.EQ.0).AND.(l.NE.0).AND.(m.NE.0)) switch=6 ! yz
    if ((k.NE.0).AND.(l.NE.0).AND.(m.NE.0)) switch=7 ! xyz

    ret=0
    select case(switch)
        case(0)
            ret=im3(0,0,0)
        case(1)
            do i = 1,dePnts
                ret=ret+lsmReady(k,i)*im3(i-pshft,0,0)
            enddo
        case(2)
            do i = 1,dePnts
                ret=ret+lsmReady(l,i)*im3(0,i-pshft,0)
            enddo
        case(3)
            do i = 1,dePnts
                ret=ret+lsmReady(m,i)*im3(0,0,i-pshft)
            enddo
        case(4)
            do i = 1,dePnts
            do j = 1,dePnts
                ret=ret+lsmReady(k,i)*lsmReady(l,j)*im3(i-pshft,j-pshft,0)
            enddo
            enddo
        case(5)
            do i = 1,dePnts
            do a = 1,dePnts
                ret=ret+lsmReady(k,i)*lsmReady(m,a)*im3(i-pshft,0,a-pshft)
            enddo
            enddo

        case(6)
            do j = 1,dePnts
            do a = 1,dePnts
                ret=ret+lsmReady(l,j)*lsmReady(m,a)*im3(0,j-pshft,a-pshft)
            enddo
            enddo

        case(7)
            do i = 1,dePnts
            do j = 1,dePnts
            do a = 1,dePnts
                ret=ret+lsmReady(k,i)*lsmReady(l,j)*lsmReady(m,a)*im3(i-pshft,j-pshft,a-pshft)
            enddo
            enddo
            enddo

    end select

    ret=ret*factorial(k)*factorial(l)*factorial(m) !*transition(k+l+m)

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

    void=glControlMemory(int( rglu*((points+2)*maxpower+3*maxpower*maxpower) ,kind=i8kind),'tmp. LSM')
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

    call tred4(lsWork,lsEVec,lsEVal,maxpower,1.e-100_rspu,1.e-300_rspu)

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

    !call prMatrix(lsmReady,ou,'KOEF MATRIX','^.000000000',maxwidth=300)

    deallocate (lsKof,lsWork,lsInv,lsEVec,lsEVal,lsRez)
    void=glControlMemory(int( sizeof(lsKof)+sizeof(lsWork)+sizeof(lsInv)+sizeof(lsEVec)+sizeof(lsEVal)+sizeof(lsRez) ,kind=i8kind),'tmp. LSM', 'free')

    return
    end subroutine prepareLSM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function deShareParams(ppoints,pstep,valSet) result(ret)
    implicit none

    integer(kind=iglu) :: ppoints,aStart,aStop,sft
    real   (kind=rglu) :: pstep,valSet(1:,1:,1:)


    deStep=pstep; dePnts=ppoints
    aStart=-(dePnts-1)/2; aStop=(dePnts-1)/2; sft=(UBound(valSet,1)-1)/2+1

    void=glControlMemory(int( rglu*(dePnts+dePnts**2+dePnts**3) ,kind=i8kind),'tmp. LSM')
    allocate (im1(aStart:aStop),&
    &         im2(aStart:aStop,aStart:aStop),&
    &         im3(aStart:aStop,aStart:aStop,aStart:aStop))

    im3=valSet(aStart+sft:aStop+sft,aStart+sft:aStop+sft,aStart+sft:aStop+sft)

    parShared=true

    if (rspu.EQ.r16kind) then
        void=glControlMemory(int( rglu*(5*dePnts) ,kind=i8kind),'tmp. LSM')
        allocate (lsmReady(0:4,dePnts)); call prepareLSM
    endif

    ret=0; return
    end function deShareParams

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine deFinalize
    implicit none

    integer(kind=iglu) :: err


    deallocate(im1,im2,im3,lsmReady,stat=err)
    if (err.EQ.0) then
        void=glControlMemory(int( sizeof(im1)+sizeof(im2)+sizeof(im3)+sizeof(lsmReady) ,kind=i8kind),'tmp. LSM', 'free')
    endif

    return
    end subroutine deFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module derivat