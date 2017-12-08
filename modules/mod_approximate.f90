    module approx

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob     , only: rglu,iglu,lglu,true,false,void,i8kind,glControlMemory
    use math     , only: tred4,gaussSLE
    use txtParser, only: tpLowerCase

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: axVersion='1.100'
    character (len=*), parameter :: axDate   ='2017.05.23'
    character (len=*), parameter :: axAuthor ='Anton B. Zakharov'

    real  (kind=rglu), parameter :: eValTol=epsilon(eValTol)*10,conditTol=1e+10_rglu

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu), allocatable :: X(:),Y(:),Yth(:),kof(:)
    real(kind=rglu), allocatable :: A(:,:),B(:,:),C(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character                    :: task*3,meth*2
    logical(kind=lglu)           :: doOutput
    integer(kind=iglu)           :: iount
    integer(kind=iglu)           :: N,M

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: axVersion,axDate,axAuthor,axFinalize,approximate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function approximate(gX,gY,gM,gtask,gkof,giount,gmeth) result(ret)
    implicit none

    integer(kind=iglu), intent(in)  :: gM
    real(kind=rglu)   , intent(in)  :: gX(:),gY(:)
    real(kind=rglu)   , intent(out) :: gkof(gM)

    character(len=3)  , intent(in)  :: gtask
    character(len=2)  , optional    :: gmeth
    integer(kind=iglu), optional    :: giount

    integer(kind=iglu)              :: i,j,k

    real(kind=rglu)                 :: condit,numer,sum1,sum2
    real(kind=rglu), allocatable    :: eVec(:,:),eVal(:)


    doOutput=false
    if (present(giount)) then
        iount=giount
        doOutput=true
    endif

    task=tpLowerCase(gtask)
    N=UBound(gX,1); M=gM

    meth='sd'; if (present(gmeth)) meth=tpLowerCase(gmeth)

    void=glControlMemory(int( rglu*(3*N+M*N+2*M*M+M) ,kind=i8kind),'tmp. Approximate')
    allocate (X(N),Y(N),Yth(N),A(M,N),B(M,M),C(M,M),kof(M))
    X=0; Y=0; Yth=0; A=0; B=0; C=0; kof=0

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

    void=glControlMemory(int( rglu*(M*M+M) ,kind=i8kind),'tmp. Approximate')
    allocate (eVec(M,M),eVal(M))

    select case (meth)
        case ('sd')
            call tred4 (B,eVec,eVal,M,1e-100_rglu,1e-100_rglu)
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
    void=glControlMemory(int( sizeof(eVec)+sizeof(eVal)+sizeof(X)+&
                              sizeof(Y)+sizeof(Yth)+sizeof(A)+sizeof(B)+&
                              sizeof(C)+sizeof(kof) ,kind=i8kind),'tmp. Geometry displacement','free')

    return
    end function approximate

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine axFinalize
    implicit none


    !continue

    return
    end subroutine axFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module approx