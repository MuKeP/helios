    module localizer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,     only: iglu,rglu,lglu,gluCompare,void,true,false,nullsub
    use hdb,      only: ou,ouwidth,localbd,mol,resetIterationState,setIterationHeader
    use printmod, only: prEigenProblem

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: lcVersion='1.000'
    character (len=*), parameter :: lcDate   ='2019.07.09'
    character (len=*), parameter :: lcAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu), allocatable :: lvecs(:,:),dmo(:),mocentroids(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) :: N,Nocc
    logical(kind=lglu) :: converged,block_converged(2)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: localize,transformationMatrix

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine localize(initvecs,localvecs,silent,fitness,centroids)
    implicit none

    real(kind=rglu),    intent(in)            :: initvecs(:,:)
    logical(kind=lglu), intent(in)            :: silent
    real(kind=rglu),    intent(out)           :: localvecs(:,:)
    real(kind=rglu),    intent(out), optional :: fitness(:),centroids(:,:)

    real(kind=rglu) :: energy(5)


    N=mol%nAtoms; Nocc=mol%nEls/2
    allocate(lvecs(N,N),dmo(N),mocentroids(N,3))

    lvecs=initvecs

    if (.NOT.silent) then
        if (localbd%method%get().EQ.'pipek-mezey') then
            call setIterationHeader(' Localization procedure: Pipek-Mezey ')
        else
            ! not implemented
        endif
    endif

    block_converged=false

    ! if verbose
    call energyLC(energy)

    if (.NOT.silent) then
        call prEigenProblem(lvecs,dmo,ou,'Localization: Init vectors (with localization numbers)','^.0000',ouWidth)
    endif

    call resetIterationState(false)
    call iterator(iterationLC,energyLC,localbd%maxiters,localbd%accuracy,true,nullsub,false,converged)

    call checkOrthonormality

    ! if verbose
    if (.NOT.silent) then
        call prEigenProblem(lvecs,dmo,ou,'Localization: Resulting vectors (with localization numbers)','^.0000',ouWidth)
    endif

    localvecs=lvecs

    if (present(fitness)) then
        fitness=dmo
    endif

    if (present(centroids)) then
        call computeCentroids
        centroids=mocentroids
    endif

    deallocate(lvecs,dmo,mocentroids)
    if (.NOT.silent) then
        call setIterationHeader
    endif

    return
    end subroutine localize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine transformationMatrix(v1, v2, transformation)
    implicit none

    real(kind=rglu), intent(in)  :: v1(:,:), v2(:,:)
    real(kind=rglu), intent(out) :: transformation(:,:)

    integer(kind=iglu)           :: N,i,j,mu


    N=UBound(v1,1)
    do i = 1,N
        do j = 1,N
            transformation(i,j)=0
            do mu = 1,N
                transformation(i,j)=transformation(i,j)+v1(mu,i)*v2(mu,j)
            enddo
        enddo
    enddo

    return
    end subroutine transformationMatrix

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine iterationLC(iteration,epsilon,accuracy)
    implicit none

    integer(kind=iglu), intent(in)  :: iteration
    real   (kind=rglu), intent(in)  :: epsilon
    real   (kind=rglu), intent(out) :: accuracy(5)

    real   (kind=rglu)              :: accur(2)


    accuracy=0
    if (localbd%method%get().EQ.'pipek-mezey') then
        call pipek_mezey(epsilon,accur); accuracy(1)=accur(1); accuracy(2)=accur(2)
    else
        ! not implemented
    endif

    return
    end subroutine iterationLC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine energyLC(energy)
    implicit none

    real   (kind=rglu) :: energy(5),asum
    integer(kind=iglu) :: i,mu


    do i = 1,N
        asum=0
        !$omp parallel default(shared) private(mu) reduction(+:asum)
        !$omp do
        do mu = 1,N
            asum=asum+lvecs(mu,i)**4
        enddo
        !$omp end parallel
        dmO(i)=1/asum
    enddo

    asum=0
    !$omp parallel default(shared) private(i) reduction(+:asum)
    !$omp do
    do i = 1,N
        asum=asum+1/dmO(i)
    enddo
    !$omp end parallel

    energy=1/(asum/N)

    return
    end subroutine energyLC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine computeCentroids
    implicit none

    integer(kind=iglu) :: i,k,mu
    real   (kind=rglu) :: asum


    do i = 1,N
        do k = 1,3
            asum=0
            !$omp parallel default(shared) private(mu) reduction(+:asum)
            !$omp do
            do mu = 1,N
                asum=asum+mol%atm(mu)%coords(k)*(lvecs(mu,i)**2)
            enddo
            !$omp end parallel
            mocentroids(i,k)=asum
        enddo
    enddo

    ! if verbose
    write(ou,'(4X,A/)') 'Localized MOs centroids (X, Y, Z):'

    write(ou,'(2X,A)') 'occupied'
    do i = 1,Nocc
        write(ou, '(4X,i3,3(2X,F11.8))') i, (mocentroids(i,k), k=1,3)
    enddo
    write(ou,'(/4X,A)') 'vacant'
    do i = Nocc+1,N
        write(ou, '(4X,i3,3(2X,F11.8))') i, (mocentroids(i,k), k=1,3)
    enddo
    write(ou,*)

    return
    end subroutine computeCentroids

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine checkOrthonormality
    implicit none

    integer(kind=iglu) :: i,j,mu
    real   (kind=rglu) :: norm,scalar
    logical(kind=lglu) :: notnormal,notorthogonal


    notnormal=false; notorthogonal=false
    do i = 1,N
        norm=0
        do mu = 1,N
            norm=norm+lvecs(mu,i)**2
        enddo
        if (norm-1.GT.10*localbd%accuracy) notnormal=true
        do j = 1,N
            if (i.EQ.j) cycle
            scalar=0
            do mu = 1,N
                scalar=scalar+lvecs(mu,i)*lvecs(mu,j)
            enddo
            if (scalar.GT.10*localbd%accuracy) notorthogonal=true
        enddo
    enddo

    if (notnormal.OR.notorthogonal) then
        write(ou,*)
        if (notnormal) then
            write(ou,'(4X,A)') 'Warning! One or more orbitals are not normalized.'
        endif
        if (notorthogonal) then
            write(ou,'(4X,A)') 'Warning! One or more pairs of orbitals are not orthgonal.'
        endif
        write(ou,*)
    endif

    return
    end subroutine checkOrthonormality

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine pipek_mezey(epsilon,accuracy)
    implicit none

    real   (kind=rglu), intent(in)  :: epsilon
    real   (kind=rglu), intent(out) :: accuracy(2)
    integer(kind=iglu)              :: k,i,j,mu,blocks(2,2)
    real   (kind=rglu)              :: teta,Pij,Pii,Pjj,tcos,tsin,v1,v2,amax,val
    real   (kind=rglu)              :: accur,norm,asum(2)


    ! occupied
    blocks(1,:)=[1,Nocc]
    ! vacant
    blocks(2,:)=[Nocc+1,N]

    norm=(Nocc**2-Nocc)/2
    accuracy=0

    do k = 1,2
        accur=0
        if (.NOT.block_converged(k)) then
            do i = blocks(k,1),blocks(k,2)-1
                do j = i+1,blocks(k,2)

                    asum=0
                    !$omp parallel default(shared) private(mu,Pij,Pii,Pjj) reduction(+:asum)
                    !$omp do
                    do mu = 1,N
                        Pij=lvecs(mu,i)*lvecs(mu,j)
                        Pii=lvecs(mu,i)*lvecs(mu,i)
                        Pjj=lvecs(mu,j)*lvecs(mu,j)

                        asum(1)=asum(1)+Pij**2-((Pii-Pjj)**2)/4
                        asum(2)=asum(2)+Pij*(Pii-Pjj)
                    enddo
                    !$omp end parallel

                    val=-asum(1)/sqrt(asum(1)**2+asum(2)**2)
                    teta=sign(acos(val)/4,asum(2))
                    tcos=cos(teta); tsin=sin(teta)

                    !$omp parallel default(shared) private(mu,v1,v2)
                    !$omp do
                    do mu = 1,N
                        v1= lvecs(mu,i)*tcos+lvecs(mu,j)*tsin
                        v2=-lvecs(mu,i)*tsin+lvecs(mu,j)*tcos
                        lvecs(mu,i)=v1
                        lvecs(mu,j)=v2
                    enddo
                    !$omp end parallel

                    accur=accur+teta**2
                enddo
            enddo
            accuracy(k)=sqrt(accur/norm)
            if (accuracy(k).LT.epsilon) block_converged(k)=true
        endif
    enddo

    return
    end subroutine pipek_mezey

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module localizer