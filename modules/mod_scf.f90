    module scf

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob      , only: iglu,rglu,void,true,false,uch,i8kind,glControlMemory
    use hdb       , only: mol,scfbd,ou,ouWidth
    use txtParser , only: tpFill
    use math      , only: tred4
    use printmod  , only: prEigenProblem

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: hfVersion='1.000'
    character (len=*), parameter :: hfDate   ='2017.08.11'
    character (len=*), parameter :: hfAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: F(:,:),V(:,:),E(:),D(:,:)
    real   (kind=rglu), allocatable :: Fmin(:,:),Comut(:,:),Refl(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) :: N,Nocc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: setSCFParameters,initSCF,iterationSCF,energySCF,getSCFResult,&
              finalizeSCF,printSCFSolution

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setSCFParameters
    implicit none


    N=mol%nAtoms; Nocc=mol%nEls/2

    call controlMemorySCF('general','allocate')

    return
    end subroutine setSCFParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine initSCF(vecs)
    implicit none

    real(kind=rglu), optional :: vecs(:,:)


    if (present(vecs)) then
        call guessSCF(vecs)
    else
        call guessSCF
    endif

    call prepareFockian

    return
    end subroutine initSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine iterationSCF(iteration,epsilon,accuracy)
    implicit none

    integer(kind=iglu), intent(in)  :: iteration
    real   (kind=rglu), intent(in)  :: epsilon
    real   (kind=rglu), intent(out) :: accuracy(5)

    real   (kind=rglu)              :: eps
    integer(kind=iglu)              :: mu,nu


    call prepareFmin; eps=maxval(abs(Fmin))

    !$omp parallel default(shared) private(mu,nu)
    !$omp do
    do mu = 1,N
        do nu = 1,N
            D(mu,nu)=D(mu,nu)-scfbd%iterStep*Fmin(mu,nu)
        enddo
    enddo
    !$omp end parallel

    call prepareFockian
    call tred4(F,V,E,N,1.e-100_rglu,1.e-300_rglu)
    call prepareDensity(V)
    call prepareFockian

    accuracy=-1; accuracy(1)=eps

    return
    end subroutine iterationSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine energySCF(energy)
    implicit none

    real   (kind=rglu) :: energy(5)
    real   (kind=rglu) :: Eel,Enuc
    integer(kind=iglu) :: mu,nu


    Eel=0
    !$omp parallel default(shared) private(mu,nu) reduction(+:Eel)
    !$omp do
    do mu = 1,N
        do nu = 1,N
            Eel=Eel+D(mu,nu)*(mol%core(mu,nu)+F(mu,nu))
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

    energy=0
    energy(1)=Eel+Enuc
    energy(2)=Eel
    energy(3)=Enuc

    return
    end subroutine energySCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getSCFResult(fockian,vectors,energies)
    implicit none

    real   (kind=rglu), optional :: fockian(:,:),vectors(:,:),energies(:)
    integer(kind=iglu)           :: i,j


    if (present(fockian)) then
        if (UBound(fockian,1).NE.N) stop 'Incorrect size of Fockian matrix.'

        do i = 1,N
            do j = 1,N
                fockian(i,j)=F(i,j)
            enddo
        enddo
    endif

    if (present(vectors)) then
        if (UBound(vectors,1).NE.N) stop 'Incorrect size of Fockian eigenvectors.'

        do i = 1,N
            do j = 1,N
                vectors(i,j)=V(i,j)
            enddo
        enddo
    endif

    if (present(energies)) then
        if (UBound(energies,1).NE.N) stop 'Incorrect size of Fockian eigenvalues.'

        do i = 1,N
            energies(i)=E(i)
        enddo
    endif

    return
    end subroutine getSCFResult

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine printSCFSolution
    implicit none


    call prEigenProblem(V,E,ou,'SCF procedure solution','^.0000',maxwidth=ouWidth)

    return
    end subroutine printSCFSolution

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeSCF
    implicit none


    call controlMemorySCF('general','deallocate')

    return
    end subroutine finalizeSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine guessSCF(vecs)
    implicit none

    character (len=10)           :: guess
    integer(kind=iglu)           :: mu
    real   (kind=rglu), optional :: vecs(N,N)


    guess=tpFill(guess); guess=scfbd%guess%get()
    if (present(vecs)) then
        guess=tpFill(guess); guess='defined'
    endif

    select case (guess)
        case ('huckel')
            call tred4(mol%huckelCore,V,E,N,1.e-100_rglu,1.e-300_rglu)
            call prepareDensity(V)

        case ('unitmatrix')
            D=0
            do mu = 1,N
                D(mu,mu)=1
            enddo

        case ('defined')
            call prepareDensity(vecs)

    end select

    return
    end subroutine guessSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareDensity(vec)
    implicit none

    real(kind=rglu)    :: vec(N,N)
    real(kind=rglu)    :: sum
    integer(kind=iglu) :: i,mu,nu


    !$omp parallel default(shared) private(mu,nu,sum,i)
    !$omp do
    do mu = 1,N
        do nu = 1,N
            sum=0
            do i = 1,Nocc
                sum=sum+vec(mu,i)*vec(nu,i)
            enddo
            D(mu,nu)=sum
        enddo
    enddo
    !$omp end parallel

    return
    end subroutine prepareDensity

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareFockian
    implicit none

    real(kind=rglu)    :: sum
    integer(kind=iglu) :: mu,nu,i


    !$omp parallel default(shared) private(mu,nu,sum,i)
    !$omp do
    do mu = 1,N
        sum=0
        do i = 1,N
            sum=sum+mol%G(mu,i)*D(i,i)
        enddo

        do nu = 1,N
            F(mu,nu)=mol%core(mu,nu)-mol%G(mu,nu)*D(mu,nu)
        enddo

        F(mu,mu)=F(mu,mu)+2*sum
    enddo
    !$omp end parallel

    return
    end subroutine prepareFockian

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareFmin
    implicit none

    real   (kind=rglu) :: sum1,sum2
    integer(kind=iglu) :: i,j,mu


    !$omp parallel default(shared) private(i,j)
    !$omp do
    do i = 1,N
        do j = 1,N
            Refl(i,j)=2*D(i,j)
        enddo
        Refl(i,i)=Refl(i,i)-1
    enddo
    !$omp end parallel

    !$omp parallel default(shared) private(i,j,sum1,sum2)
    !$omp do
    do i = 1,N
        do j = 1,N
            sum1=0; sum2=0
            do mu = 1,N
                sum1=sum1+F(i,mu)*Refl(mu,j)
                sum2=sum2+F(mu,j)*Refl(i,mu)
            enddo
            Comut(i,j)=sum1-sum2
        enddo
    enddo
    !$omp end parallel

    !$omp parallel default(shared) private(i,j,sum1)
    !$omp do
    do i = 1,N
        do j = 1,N
            sum1=0
            do mu = 1,N
                sum1=sum1+Comut(i,mu)*Refl(mu,j)
            enddo
            Fmin(i,j)=sum1
        enddo
    enddo
    !$omp end parallel

    return
    end subroutine prepareFmin

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine controlMemorySCF(section,action)
    implicit none

    character (len=*)  :: section,action
    integer(kind=iglu) :: err


    select case (section)
        case ('general')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int(rglu*(6*N*N+N),kind=i8kind),'SCF module')
                    allocate (F(N,N),V(N,N),E(N),D(N,N))
                    allocate (Fmin(N,N),Comut(N,N),Refl(N,N))
                    F=0; V=0; E=0; D=0; Fmin=0; Comut=0; Refl=0

                case ('deallocate')
                    deallocate (F,V,E,D,Fmin,Comut,Refl, stat=err)
                    void=glControlMemory(int(rglu*(6*N*N+N),kind=i8kind),'SCF module','free')

            end select
    end select

    return
    end subroutine controlMemorySCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module scf