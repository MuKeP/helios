    module mbpt

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob     , only: assignment(=)
    use glob     , only: uch,iglu,rglu,lglu,true,false,void,i8kind,glControlMemory
    use printmod , only: prMatrix,prEigenProblem
    use scf      , only: setSCFParameters,initSCF,iterationSCF,getSCFResult
    use scf      , only: energySCF,finalizeSCF,printSCFSolution
    use hdb      , only: mol,systembd,scfbd,ou,ouWidth

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: ptVersion='1.000'
    character (len=*), parameter :: ptDate   ='2017.08.11'
    character (len=*), parameter :: ptAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: density(:,:),V(:,:),Vs(:,:),F(:,:),R(:,:,:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type(uch)          :: umethod
    integer(kind=iglu) :: N,Nel,No,Nocc,Nth
    real   (kind=rglu) :: accuracy(5),refeEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    public :: prVersion,ptDate,ptAuthor
    public :: setMBPTParameters,initMBPT,energyMBPT,finalizeMBPT

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setMBPTParameters(method)
    implicit none

    character (len=*), intent(in) :: method


    select case (method)
        case ('mp2','mp3')
            umethod=method

        case default
            stop 'MBPT: Unknown method'

    end select

    N=mol%nAtoms; Nel=mol%nEls; Nth=systembd%nNodes; No=2*N; Nocc=Nel/2
    call controlMemoryMBPT('general','allocate')
    call setSCFParameters

    return
    end subroutine setMBPTParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine initMBPT
    implicit none

    integer(kind=iglu) :: i,j,a,b


    call initSCF
    call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true)
    call getSCFResult(vectors=V)

    call printSCFSolution
    call prepareDensity(V)
    call prepareFock

    do i = 1,N
    do j = 1,N
        Vs(i,2*j-1)=V(i,j)
        Vs(i,2*j  )=V(i,j)
    enddo
    enddo

    !$omp parallel default(shared) private(i,j,a,b)
    !$omp do
    do i = 1,No
    do j = 1,No
    do a = 1,No
    do b = 1,No
        R(i,j,a,b)=spin_hf_int(i,j,a,b)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel

    return
    end subroutine initMBPT

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine energyMBPT(energy)

    implicit none

    real   (kind=rglu), intent(out) :: energy(5)
    real   (kind=rglu)              :: Ax,sum0,sum1,sum2,sum3,sum4
    integer(kind=iglu)              :: i,j,a,b,k,c,l,d


    energy=0
    select case (umethod%get())
        case ('mp2')
            sum1=0
            !$omp parallel default(shared) private(i,j,a,b) reduction(+:sum1)
            !$omp do
            do i = 1,Nel-1
            do j = i+1,Nel
                do a = Nel+1,No-1
                do b = a+1,No
                    sum1=sum1+R(i,a,j,b)**2/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
                enddo
                enddo
            enddo
            enddo
            !$omp end parallel
            energy(1)=sum1+refeEnergy
            energy(2)=refeEnergy
            energy(3)=sum1

        case ('mp3')
            sum0=0
            !$omp parallel default(shared) private(i,j,a,b) reduction(+:sum0)
            !$omp do
            do i = 1,Nel-1
            do j = i+1,Nel
                do a = Nel+1,No-1
                do b = a+1,No
                    sum0=sum0+R(i,a,j,b)**2/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
                enddo
                enddo
            enddo
            enddo
            !$omp end parallel

            sum4=0
            do i = 1,Nel
            do j = 1,Nel
            do a = Nel+1,No
            do b = Nel+1,No

                Ax=R(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))

                sum1=0
                !$omp parallel default(shared) private(c,d) reduction(+:sum1)
                !$omp do
                do c = Nel+1,No
                do d = Nel+1,No
                    sum1=sum1+R(a,c,b,d)*R(c,i,d,j)/(F(i,i)+F(j,j)-F(c,c)-F(d,d))
                enddo
                enddo
                !$omp end parallel

                sum2=0
                !$omp parallel default(shared) private(k,l) reduction(+:sum2)
                !$omp do
                do k = 1,Nel
                do l = 1,Nel
                    sum2=sum2+R(i,k,j,l)*R(a,k,b,l)/(F(k,k)+F(l,l)-F(a,a)-F(b,b))
                enddo
                enddo
                !$omp end parallel

                sum3=0
                !$omp parallel default(shared) private(k,c) reduction(+:sum3)
                !$omp do
                do k = 1,Nel
                do c = Nel+1,No
                    sum3=sum3+R(i,a,c,k)*R(j,b,k,c)/(F(j,j)+F(k,k)-F(b,b)-F(c,c))
                enddo
                enddo
                !$omp end parallel

                sum4=sum4+Ax*( (sum1+sum2)/8+sum3 )
            enddo
            enddo
            enddo
            enddo
            energy(1)=sum0+sum4+refeEnergy
            energy(2)=refeEnergy
            energy(3)=sum0+sum4
            energy(4)=sum0

    end select


    return
    end subroutine energyMBPT

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeMBPT

    implicit none


    call controlMemoryMBPT('general','deallocate')
    call finalizeSCF

    return
    end subroutine finalizeMBPT

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareFock

    implicit none

    real   (kind=rglu), allocatable :: X(:,:)
    real   (kind=rglu)              :: sum
    integer(kind=iglu)              :: mu,nu,i,j


    allocate (X(N,N)); X=0

    !$omp parallel default(shared) private(mu,nu,sum)
    !$omp do
    do mu = 1,N
        sum=0
        do nu = 1,N
            sum=sum+mol%G(mu,nu)*density(nu,nu)
        enddo
        do nu = 1,N
            X(mu,nu)=mol%core(mu,nu)-mol%G(mu,nu)*density(mu,nu)
        enddo
        X(mu,mu)=X(mu,mu)+2*sum
    enddo
    !$omp end parallel

    call referenceEnergy(X)

    F=0
    !$omp parallel default(shared) private(i,j,mu,nu,sum)
    !$omp do
    do i = 1,N
        do j = i,N
            sum=0
            do mu = 1,N
            do nu = 1,N
                sum=sum+V(mu,i)*V(nu,j)*X(mu,nu)
            enddo
            enddo

            F(2*i-1,2*j-1)=sum; F(2*i,2*j)=sum
            F(2*j-1,2*i-1)=sum; F(2*j,2*i)=sum
        enddo
    enddo
    !$omp end parallel

    deallocate (X)

    return
    end subroutine prepareFock

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine referenceEnergy(X)
    implicit none

    real   (kind=rglu) :: Eel,Enuc,X(:,:)
    integer(kind=iglu) :: mu,nu


    Eel=0
    !$omp parallel default(shared) private(mu,nu) reduction(+:Eel)
    !$omp do
    do mu = 1,N
        do nu = 1,N
            Eel=Eel+density(mu,nu)*( mol%core(mu,nu)+X(mu,nu) )
        enddo
    enddo
    !$omp end parallel

    Enuc=0
    !$omp parallel default(shared) private(mu,nu) reduction(+:Enuc)
    !$omp do
    do mu = 1,N
        do nu = 1,N
            if (mu.EQ.nu) cycle
            Enuc=Enuc+mol%G(mu,nu)/2
         enddo
    enddo
    !$omp end parallel
    refeEnergy=Eel+Enuc

    return
    end subroutine referenceEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareDensity(V)

    implicit none

    real(kind=rglu)    :: V(N,N)
    real(kind=rglu)    :: sum
    integer(kind=iglu) :: i,mu,nu

    !$omp parallel default(shared) private(mu,nu,i,sum)
    !$omp do
    do mu = 1,N
        do nu = 1,N
            sum=0
            do i = 1,Nocc
                sum=sum+V(mu,i)*V(nu,i)
            enddo
            density(mu,nu)=sum
        enddo
    enddo
    !$omp end parallel

    end subroutine prepareDensity

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function spin_hf_int(a,b,c,d) result(ret) ! [ab||cd]=[ab|cd]-[ad|bc]
    implicit none
    integer(kind=iglu), intent(in) :: a,b,c,d
    real   (kind=rglu)             :: sum1,sum2
    integer(kind=iglu)             :: mu,nu


    sum1=0
    if ((.NOT.btest(a+b,0)).AND.(.NOT.btest(c+d,0))) then ! [ab|cd]
        do mu = 1,N
        do nu = 1,N
            sum1=sum1+Vs(mu,a)*Vs(mu,b)*Vs(nu,c)*Vs(nu,d)*mol%G(mu,nu)
        enddo
        enddo
    endif

    sum2=0
    if ((.NOT.btest(a+d,0)).AND.(.NOT.btest(b+c,0))) then ! [ad|bc]
        do mu = 1,N
        do nu = 1,N
            sum2=sum2+Vs(mu,a)*Vs(mu,d)*Vs(nu,b)*Vs(nu,c)*mol%G(mu,nu)
        enddo
        enddo
    endif

    ret=sum1-sum2; return
    end function spin_hf_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine controlMemoryMBPT(section,action)
    implicit none

    character  (len=*) :: section,action
    integer(kind=iglu) :: err


    select case (section)
        case ('general')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(N*No+N*N+No*No*No*No+No*No) ,kind=i8kind),'MBPT module')
                    allocate (Vs(N,No),V(N,N),density(N,N),R(No,No,No,No),F(No,No))
                    Vs=0; V=0; density=0; R=0; F=0
                case ('deallocate')
                    deallocate (Vs,V,density,R,F, stat=err)
                    void=glControlMemory(int( rglu*(N*No+N*N+No*No*No*No+No*No) ,kind=i8kind),'MBPT module','free')
            end select
    end select

    return
    end subroutine controlMemoryMBPT

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module mbpt