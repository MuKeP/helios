    module math

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob, only: rglu,rspu,iglu,lglu,true,false,void,i8kind,glControlMemory

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: maVersion='1.100'
    character (len=*), parameter :: maDate   ='2017.12.10'
    character (len=*), parameter :: maAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    interface tred4
        module procedure gltred4, sptred4
    end interface tred4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: maVersion,maDate,maAuthor,maFinalize
    public :: tred4,gaussSLE,LagrangeDerivative,getMult,factorial,gcd,lcm,reduceFraction

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine gltred4(A,Z,D,N,AC,ATOL)
    implicit none

    integer(kind=iglu), intent(in)  :: N
    real   (kind=rglu), intent(in)  :: A(:,:)
    real   (kind=rglu), intent(in)  :: ATOL,AC

    real   (kind=rglu), intent(out) :: D(:),Z(:,:)

    integer(kind=iglu)              :: i,j,k,l,m
    real   (kind=rglu)              :: F,G,H,HH,B,P,R,C,S
    real   (kind=rglu), allocatable :: E(:)


    void=glControlMemory(int( rglu*N ,kind=i8kind),'tmp. tred4')
    allocate (E(N)); E=0
    do i=1,N
        do j=1,i
           Z(i,j)=A(i,j)
        enddo
    enddo

    if (N.GT.1) then
        do i=N,2,-1
            l=i-2
            F=Z(i,i-1)
            G=0
            if(l.NE.0) then
                do k = 1,l
                    G=G+Z(i,k)**2
                enddo
            endif
            H=G+F**2
            if (G.LE.ATOL) then
                E(i)=F
                H=0
                D(i)=H
                cycle
            endif
            l=l+1
            G=sqrt(H)
            if (F.GE.0) G=-G
            E(i)=G
            H=H-F*G
            Z(i,i-1)=F-G
            F=0
            do j=1,l
                Z(j,i)=Z(i,j)/H
                G=0
                do k=1,j
                    G=G+Z(j,k)*Z(i,k)
                enddo
                if (j+1.LE.l) then
                    do k=j+1,l
                        G=G+Z(k,j)*Z(i,k)
                    enddo
                endif
                E(j)=G/H
                F=F+G*Z(j,i)
            enddo
            HH=F/(H+H)
            do j=1,l
                F=Z(i,j)
                G=E(j)-HH*F
                E(j)=G
                do k=1,j
                    Z(j,k)=Z(j,k)-F*E(k)-G*Z(i,k)
                enddo
            enddo
            D(i)=H
        enddo
    endif

    E(1)=0
    D(1)=0
    do i=1,N
        if (D(i).NE.0) then
            do j=1,i-1
                G=0
                do k=1,i-1
                    G=G+Z(i,k)*Z(k,j)
                enddo
                do k=1,i-1
                    Z(k,j)=Z(k,j)-G*Z(k,i)
                enddo
            enddo
        endif
        D(i)=Z(i,i)
        Z(i,i)=1
        if (i-1.NE.0) then
            do j=1,i-1
                Z(i,j)=0
                Z(j,i)=0
            enddo
        endif
    enddo
    if (N.NE.1) then
        do i=2,N
            E(i-1)=E(i)
        enddo
    endif
    E(N)=0

    B=0
    F=0
    j=30*N
    do l=1,N
        H=AC*(abs(D(l))+abs(E(l)))
        if (B.LT.H) B=H
        do m=l,N
            if (abs(E(m)).LE.B) exit
        enddo
        if (m.EQ.l) then
            D(l)=D(l)+F
            cycle
        endif
1       if (j.LE.0) then
            deallocate (E)
            void=glControlMemory(int( rglu*N ,kind=i8kind),'tmp. tred4','free')
            return
        endif
        j=j-1
        G=D(l)
        H=D(l+1)-G
        if (abs(H).LT.abs(E(l))) then
            P=H*0.5_rglu/E(l)
            R=sqrt(P**2+1)
            H=P+R
            if (P.LT.0) H=P-R
            D(l)=E(l)/H
        else
            P=2._rglu*E(l)/H
            R=sqrt(P**2+1)
            D(l)=E(l)*P/(1+R)
        endif
        H=G-D(l)
        if (l+1.LE.N) then
            do i=l+1,N
                D(i)=D(i)-H
            enddo
        endif
        F=F+H
        P=D(m)
        C=1
        S=0
        do i=m-1,l,-1
            G=C*E(i)
            H=C*P
            if (abs(P).GE.abs(E(i))) then
                C=E(i)/P
                R=sqrt(C**2+1)
                E(i+1)=S*P*R
                S=C/R
                C=1/R
            else
                C=P/E(i)
                R=sqrt(C**2+1)
                E(i+1)=S*E(i)*R
                S=1/R
                C=C/R
            endif
            P=C*D(i)-S*G
            D(i+1)=H+S*(C*G+S*D(i))
            do k=1,N
                H=Z(k,i+1)
                Z(k,i+1)=S*Z(k,i)+C*H
                Z(k,i)=C*Z(k,i)-S*H
            enddo
        enddo
        E(l)=S*P
        D(l)=C*P
        if (abs(E(l)).GT.B) goto 1
        D(l)=D(l)+F
    enddo

    do i=1,N
        k=i
        P=D(i)
        if(i+1.LE.N) then
            do j=i+1,N
                if(D(j).GE.P) cycle
                k=j
                P=D(j)
            enddo
        endif
        if (k.EQ.i) cycle
        D(k)=D(i)
        D(i)=P
        do j=1,N
            P=Z(j,i)
            Z(j,i)=Z(j,k)
            Z(j,k)=P
        enddo
    enddo

    deallocate (E)
    void=glControlMemory(int( rglu*N ,kind=i8kind),'tmp. tred4','free')

    return
    end subroutine gltred4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine sptred4(A,Z,D,N,AC,ATOL)
    implicit none

    integer(kind=iglu), intent(in)  :: N
    real   (kind=rspu), intent(in)  :: A(:,:)
    real   (kind=rspu), intent(in)  :: ATOL,AC

    real   (kind=rspu), intent(out) :: D(:),Z(:,:)

    integer(kind=iglu)              :: i,j,k,l,m
    real   (kind=rspu)              :: F,G,H,HH,B,P,R,C,S

    real   (kind=rspu), allocatable :: E(:)


    void=glControlMemory(int( rspu*N ,kind=i8kind),'tmp. tred4')
    allocate (E(N)); E=0
    do i=1,N
        do j=1,i
           Z(i,j)=A(i,j)
        enddo
    enddo

    if (N.GT.1) then
        do i=N,2,-1
            l=i-2
            F=Z(i,i-1)
            G=0
            if(l.NE.0) then
                do k = 1,l
                    G=G+Z(i,k)**2
                enddo
            endif
            H=G+F**2
            if (G.LE.ATOL) then
                E(i)=F
                H=0
                D(i)=H
                cycle
            endif
            l=l+1
            G=sqrt(H)
            if (F.GE.0) G=-G
            E(i)=G
            H=H-F*G
            Z(i,i-1)=F-G
            F=0
            do j=1,l
                Z(j,i)=Z(i,j)/H
                G=0
                do k=1,j
                    G=G+Z(j,k)*Z(i,k)
                enddo
                if (j+1.LE.l) then
                    do k=j+1,l
                        G=G+Z(k,j)*Z(i,k)
                    enddo
                endif
                E(j)=G/H
                F=F+G*Z(j,i)
            enddo
            HH=F/(H+H)
            do j=1,l
                F=Z(i,j)
                G=E(j)-HH*F
                E(j)=G
                do k=1,j
                    Z(j,k)=Z(j,k)-F*E(k)-G*Z(i,k)
                enddo
            enddo
            D(i)=H
        enddo
    endif

    E(1)=0
    D(1)=0
    do i=1,N
        if (D(i).NE.0) then
            do j=1,i-1
                G=0
                do k=1,i-1
                    G=G+Z(i,k)*Z(k,j)
                enddo
                do k=1,i-1
                    Z(k,j)=Z(k,j)-G*Z(k,i)
                enddo
            enddo
        endif
        D(i)=Z(i,i)
        Z(i,i)=1
        if (i-1.NE.0) then
            do j=1,i-1
                Z(i,j)=0
                Z(j,i)=0
            enddo
        endif
    enddo
    if (N.NE.1) then
        do i=2,N
            E(i-1)=E(i)
        enddo
    endif
    E(N)=0

    B=0
    F=0
    j=30*N
    do l=1,N
        H=AC*(abs(D(l))+abs(E(l)))
        if (B.LT.H) B=H
        do m=l,N
            if (abs(E(m)).LE.B) exit
        enddo
        if (m.EQ.l) then
            D(l)=D(l)+F
            cycle
        endif
1       if (j.LE.0) then
            deallocate (E)
            void=glControlMemory(int( rspu*N ,kind=i8kind),'tmp. tred4','free')
            return
        endif
        j=j-1
        G=D(l)
        H=D(l+1)-G
        if (abs(H).LT.abs(E(l))) then
            P=H*0.5_rspu/E(l)
            R=sqrt(P**2+1)
            H=P+R
            if (P.LT.0) H=P-R
            D(l)=E(l)/H
        else
            P=2._rspu*E(l)/H
            R=sqrt(P**2+1)
            D(l)=E(l)*P/(1+R)
        endif
        H=G-D(l)
        if (l+1.LE.N) then
            do i=l+1,N
                D(i)=D(i)-H
            enddo
        endif
        F=F+H
        P=D(m)
        C=1
        S=0
        do i=m-1,l,-1
            G=C*E(i)
            H=C*P
            if (abs(P).GE.abs(E(i))) then
                C=E(i)/P
                R=sqrt(C**2+1)
                E(i+1)=S*P*R
                S=C/R
                C=1/R
            else
                C=P/E(i)
                R=sqrt(C**2+1)
                E(i+1)=S*E(i)*R
                S=1/R
                C=C/R
            endif
            P=C*D(i)-S*G
            D(i+1)=H+S*(C*G+S*D(i))
            do k=1,N
                H=Z(k,i+1)
                Z(k,i+1)=S*Z(k,i)+C*H
                Z(k,i)=C*Z(k,i)-S*H
            enddo
        enddo
        E(l)=S*P
        D(l)=C*P
        if (abs(E(l)).GT.B) goto 1
        D(l)=D(l)+F
    enddo

    do i=1,N
        k=i
        P=D(i)
        if(i+1.LE.N) then
            do j=i+1,N
                if(D(j).GE.P) cycle
                k=j
                P=D(j)
            enddo
        endif
        if (k.EQ.i) cycle
        D(k)=D(i)
        D(i)=P
        do j=1,N
            P=Z(j,i)
            Z(j,i)=Z(j,k)
            Z(j,k)=P
        enddo
    enddo

    deallocate (E)
    void=glControlMemory(int( rspu*N ,kind=i8kind),'tmp. tred4','free')

    return
    end subroutine sptred4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function gaussSLE(N,matrix,inverted) result(ret)

    implicit none

    integer(kind=iglu), intent(in) :: N
    integer(kind=iglu)             :: k,l

    real(kind=rglu), intent(in)    :: matrix(N,N)
    real(kind=rglu), intent(out)   :: inverted(N,N)

    real(kind=rglu)                :: tol
    real(kind=rglu), allocatable   :: dmatrix(:,:),dinverted(:,:)


    tol=epsilon(tol)*10; ret=1
    void=glControlMemory(int( 2*rglu*N*N ,kind=i8kind),'tmp. gaussSLE')
    allocate (dmatrix(N,N),dinverted(N,N))
    dmatrix=matrix

    dinverted=0
    do k = 1,N
        dinverted(k,k)=1
    enddo

    do k = 1,N
        if (abs(dmatrix(k,k)).LT.tol) then
            do l = 1,N
                if (dmatrix(l,k).GT.tol) exit
            enddo

            dinverted(k,:)=dinverted(k,:)-dinverted(l,:)
            dmatrix(k,:)  =dmatrix(k,:)  -dmatrix(l,:)
        endif
    enddo

    do k = 1,N

        ret=ret*dmatrix(k,k)
        dinverted(k,:)=dinverted(k,:)/dmatrix(k,k)
        dmatrix(k,:)  =dmatrix(k,:)  /dmatrix(k,k)

        do l = k+1,N
            if (abs(dmatrix(l,k)).LT.tol) cycle
            dinverted(l,:)=dinverted(l,:)-dmatrix(l,k)*dinverted(k,:)
            dmatrix(l,:)  =dmatrix(l,:)  -dmatrix(l,k)*dmatrix(k,:)
        enddo

    enddo

    do k = N,1,-1
        do l = k-1,1,-1
            if (abs(dmatrix(l,k)).LT.tol) cycle
            dinverted(l,:)=dinverted(l,:)-dmatrix(l,k)*dinverted(k,:)
            dmatrix(l,:)  =dmatrix(l,:)  -dmatrix(l,k)*dmatrix(k,:)
        enddo
    enddo

    inverted=dinverted

    deallocate (dmatrix,dinverted)
    void=glControlMemory(int( 2*rglu*N*N ,kind=i8kind),'tmp. gaussSLE','free')

    return
    end function gaussSLE

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rspu) function LagrangeDerivative(nump,power,enGrid,ff) result(ret)

    implicit none

    integer(kind=iglu), intent(in)  :: nump,power
    real(kind=rglu)   , intent(in)  :: ff,enGrid( -(nump-1)/2 : (nump-1)/2 )

    integer(kind=iglu)              :: shif,maxpower,i,j,k
    integer(kind=iglu), allocatable :: lagWork(:)

    real(kind=rspu)                 :: sum,rez
    real(kind=rspu)   , allocatable :: lagZnam(:),lagChis(:,:),lagDerivat(:),lagFr(:,:)


    maxpower=nump-1; shif=(nump-1)/2
    allocate ( lagZnam(-shif:shif),lagChis(-shif:shif,maxpower),lagWork(maxpower),&
    &           lagDerivat(-shif:shif),lagFr(-shif:shif,2) )

    do i = -shif,shif
        lagZnam(i)=1
        do j = -shif,shif
            if (i.EQ.j) cycle
            lagZnam(i)=lagZnam(i)*real(i-j,rspu)
        enddo
    enddo

    do i = -shif,shif
        k=0
        do j = -shif,shif
            if (i.EQ.j) cycle
            k=k+1; lagWork(k)=j
        enddo
        do j = 2,maxpower
            call getMult(lagWork,maxpower,maxpower-j,rez)
            lagChis(i,j-1)=rez
        enddo
    enddo

    lagChis(:,maxpower)=1

    do j = -shif,shif
        lagDerivat(j)=(-1.)**power*lagChis(j,power)*factorial(power) !change
    enddo

    do j = -shif,shif
        lagFr(j,1)=lagDerivat(j)
        lagFr(j,2)=lagZnam(j)
    enddo

    do j = -shif,shif
        lagFr(j,1)=lagFr(j,1)/lagFr(j,2)
    enddo

    sum=0
    do i = -shif,shif
        sum=sum+lagFr(i,1)*enGrid(i)
    enddo

    deAllocate ( lagZnam,lagChis,lagWork,lagDerivat,lagFr )

    ret=sum/(ff**power)

    return
    end function LagrangeDerivative

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine getMult(arr,i,incl,rez)
    implicit none
    integer(kind=iglu), intent(in)  :: arr(:),incl
    integer(kind=iglu)              :: i,k
    real   (kind=rspu)              :: bb
    real   (kind=rspu), intent(out) :: rez


    if (incl.EQ.0) then
        rez=real(sum(arr(:i)),rspu)
    else
        bb=0
        do k = i,1,-1
            call getMult(arr,k-1,incl-1,rez)
            bb=bb+real(arr(k),rspu)*rez
        enddo
        rez=bb
    endif

    return
    end subroutine getMult

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure real(kind=rglu) function factorial(x) result(ret)
    implicit none
    integer(kind=iglu), intent(in) :: x
    integer(kind=iglu)             :: i


    ret=1
    do i = 1,x
        ret=ret*real(i,rglu)
    enddo

    return
    end function factorial

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure real(kind=rspu) function gcd(i,j) result(ret)
    implicit none
    real(kind=rspu), intent(in) :: i,j
    real(kind=rspu)             :: a,b


    a=i; b=j
    a=abs(a); b=abs(b)
    do
        if (floor(a*b).EQ.0) exit

        if (a.GE.b) then
            a=a-floor(a/b)*b
        else
            b=b-floor(b/a)*a
        endif
    enddo
    ret=a+b

    return
    end function gcd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    pure real(kind=rspu) function lcm(i,j) result(ret)
    implicit none
    real(kind=rspu), intent(in) :: i,j
    real(kind=rspu)             :: k,l,a,b


    a=i; b=j; k=a; l=b
    k=abs(k); l=abs(l)
    do
        if (floor(k*l).EQ.0) exit
        if (k.GE.l) then
            k=k-floor(k/l)*l
        else
            l=l-floor(l/k)*k
        endif
    enddo
    ret=floor(a*b)/floor(k+l)

    return
    end function lcm

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine reduceFraction(i,j)
    implicit none
    real(kind=rspu) :: devisor,i,j


    devisor=gcd(i,j); i=i/devisor; j=j/devisor

    return
    end subroutine reduceFraction

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine maFinalize
    implicit none


    !continue

    return
    end subroutine maFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module math