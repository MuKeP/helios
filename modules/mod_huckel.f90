    module huckel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob    , only: iglu,rglu,lglu
    use hdb     , only: mol,ou,ouWidth
    use math    , only: tred4
    use printmod, only: prEigenProblem

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: V(:,:),E(:),D(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu)              :: N,Nocc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: getHuckelResult

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getHuckelResult(energy,vectors,energies,density,output)
    implicit none

    real   (kind=rglu), optional, intent(out) :: energy(:)
    real   (kind=rglu), optional, intent(out) :: vectors(:,:),density(:,:),energies(:)
    logical(kind=lglu), optional, intent(in)  :: output

    integer(kind=iglu) :: i,mu,nu
    real   (kind=rglu) :: sum,Enuc,Eel


    N=mol%nAtoms; Nocc=mol%nEls/2
    allocate(V(N,N),E(N),D(N,N)); V=0; E=0; D=0
    call tred4(mol%huckelCore,V,E,N,1.e-100_rglu,1.e-300_rglu)

    if (present(output)) then
        if (output) call prEigenProblem(V,E,ou,'Huckel solution','^.0000',maxwidth=ouWidth)
    endif

    Eel=0
    !$omp parallel default(shared) private(mu) reduction(+:Eel)
    !$omp do
    do mu = 1,Nocc
        Eel=Eel+2*E(mu)
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

    !$omp parallel default(shared) private(mu,nu,i,sum)
    !$omp do
    do mu = 1,N
        do nu = 1,N
            sum=0
            do i = 1,Nocc
                sum=sum+V(mu,i)*V(nu,i)
            enddo
            D(mu,nu)=sum
        enddo
    enddo
    !$omp end parallel

    if (present(vectors))  vectors=V
    if (present(energies)) energies=E
    if (present(density))  density=D

    if (present(energy)) then
        energy=0; energy(1)=Eel
    endif

    deallocate (V,E,D)

    return
    end subroutine getHuckelResult

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module huckel