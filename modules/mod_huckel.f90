    module huckel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob    , only: iglu,rglu,lglu,void,i8kind,glControlMemory
    use hdb     , only: mol,ou,ouWidth
    use math    , only: tred4
    use printmod, only: prEigenProblem

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: huVersion='1.100'
    character (len=*), parameter :: huDate   ='2019.05.24'
    character (len=*), parameter :: huAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: V(:,:),E(:),D(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu)              :: N,Nocc
    real(kind=rglu)                 :: Eel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: setHuckelParameters,energyHuckel,getHuckRDMElement,getHuckelResult,finalizeHuck
    public :: getHuckelCoulsonPolarizability

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setHuckelParameters
    implicit none


    N=mol%nAtoms; Nocc=mol%nEls/2
    call controlMemoryHuck('general','allocate')

    return
    end subroutine setHuckelParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine energyHuckel(energy)
    implicit none

    real   (kind=rglu) :: energy(5)
    real   (kind=rglu) :: sum
    integer(kind=iglu) :: mu,nu,i


    call tred4(mol%huckelCore,V,E,N,1.e-100_rglu,1.e-300_rglu)

    ! write (*,*) E

    Eel=0
    !$omp parallel default(shared) private(mu) reduction(+:Eel)
    !$omp do
    do mu = 1,Nocc
        Eel=Eel+2*E(mu)
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

    energy=0; energy(1)=Eel

    return
    end subroutine energyHuckel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function getHuckRDMElement(i,j) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: i,j


    ret=2*D(i,j); return
    end function getHuckRDMElement

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function getHuckelCoulsonPolarizability(mu,nu) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: mu,nu
    integer(kind=iglu)             :: i,j
    real(kind=rglu)                :: Ax


    ret=0
    do i = 1,Nocc
        Ax=V(mu,i)*V(nu,i)
        do j = Nocc+1,N
            ret=ret+Ax*V(mu,j)*V(nu,j)/(E(i)-E(j))
        enddo
    enddo
    ret=4*ret

    return
    end function getHuckelCoulsonPolarizability

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getHuckelResult(energy,vectors,energies,density,output)
    implicit none

    real   (kind=rglu), optional, intent(out) :: energy(:)
    real   (kind=rglu), optional, intent(out) :: vectors(:,:),density(:,:),energies(:)
    logical(kind=lglu), optional, intent(in)  :: output


    if (present(output)) then
        if (output) call prEigenProblem(V,E,ou,'Huckel solution','^.0000',maxwidth=ouWidth)
    endif

    if (present(vectors))  vectors=V
    if (present(energies)) energies=E
    if (present(density))  density=2*D

    if (present(energy)) then
        energy=0; energy(1)=Eel
    endif

    return
    end subroutine getHuckelResult

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeHuck
    implicit none


    call controlMemoryHuck('general','deallocate')

    return
    end subroutine finalizeHuck

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine controlMemoryHuck(section,action)
    implicit none

    character (len=*)  :: section,action
    integer(kind=iglu) :: err


    select case (section)
        case ('general')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(N*N+N+N*N) ,kind=i8kind),'Huckel')
                    allocate(V(N,N),E(N),D(N,N)); V=0; E=0; D=0

                    ! void=glControlMemory(int(rglu*(6*N*N+N),kind=i8kind),'SCF module')
                    ! allocate (F(N,N),V(N,N),E(N),D(N,N))
                    ! allocate (Fmin(N,N),Comut(N,N),Refl(N,N))
                    ! F=0; V=0; E=0; D=0; Fmin=0; Comut=0; Refl=0

                case ('deallocate')
                    deallocate (V,E,D)
                    void=glControlMemory(int( sizeof(V)+sizeof(E)+sizeof(D) ,kind=i8kind),'Huckel','free')

                    ! deallocate (F,V,E,D,Fmin,Comut,Refl, stat=err)
                    ! void=glControlMemory(int(rglu*(6*N*N+N),kind=i8kind),'SCF module','free')

            end select
    end select

    return
    end subroutine controlMemoryHuck

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine huckDensity2

    use math,     only: LagrangeDerivative
    use printmod, only: prMatrix
    use hdb,      only: perturbate,fieldbd,densitybd

    implicit none

    integer(kind=iglu) :: k,N,fatom,satom,field,sta,sto,derivPoints
    real(kind=rglu)    :: derivStep,zpEnergy,energy(5)
    real(kind=rglu), allocatable :: D(:,:),derivEnergies(:),E(:)


    call setHuckelParameters

    N=mol%nAtoms
    derivPoints=densitybd%nPoints
    derivStep=densitybd%derivStep; sta=-derivPoints/2; sto= derivPoints/2


    allocate (D(N,N),E(N),derivEnergies(sta:sto)); D=0; derivEnergies=0

    call setCore(1,1,0)
    call prMatrix(mol%connect,6,'Connectivity','^.00000',maxwidth=78)
    call prMatrix(mol%huckelCore,6,'Huckel core','^.00000',maxwidth=78)
    call energyHuckel(energy)
    zpEnergy = energy(1)

    do fatom = 1,N
    do satom = fatom,N
        derivEnergies(0)=zpEnergy

        do field = sta,sto
            if (field.EQ.0) cycle

            call setCore(fatom,satom,field)
            call energyHuckel(energy); derivEnergies(field)=energy(1)
        enddo

        write (*,'(1X,i2,1X,i2,<derivPoints>(1X,F15.10))') fatom, satom, derivEnergies

        call getDeriv(fatom,satom)
    enddo
    enddo

    call prMatrix(D,6,'Huckel numerical','^.00000',maxwidth=78)

    call setCore(1,1,0)

    call prMatrix(mol%connect,6,'Connectivity__','^.00000',maxwidth=78)
    call prMatrix(mol%huckelCore,6,'Huckel core__','^.00000',maxwidth=78)

    call energyHuckel(energy)

    call getHuckelResult(energies=E)
    write (*,*) E

    do fatom = 1,N
    do satom = 1,N
        D(fatom,satom)=getHuckRDMElement(fatom,satom)
    enddo
    enddo

    call prMatrix(D,6,'Huckel analytical','^.00000',maxwidth=78)

    call finalizeHuck

    return

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(ip1,ip2,fshift)
        implicit none
        integer(kind=iglu), intent(in) :: ip1,ip2,fshift


        ! set core perturbation
        mol%perturbation=0
        call applyField

        if (ip1.EQ.ip2) then
            mol%perturbation(ip1,ip2)=mol%perturbation(ip1,ip2)+fShift*derivStep
        else
            mol%perturbation(ip1,ip2)=fShift*derivStep
            mol%perturbation(ip2,ip1)=fShift*derivStep
        endif
        call perturbate

        return
        end subroutine setCore

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine applyField
        implicit none

        integer(kind=iglu) :: k,mu


        do k = 1,3
            if (abs(fieldbd%strength(k)).GT.1D-10) then
                do mu = 1,mol%nAtoms
                    mol%perturbation(mu,mu)=fieldbd%strength(k)*mol%atm(mu)%coords(k)
                enddo
            endif
        enddo

        return
        end subroutine applyField

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine getDeriv(ip1,ip2)
        implicit none
        integer(kind=iglu), intent(in) :: ip1,ip2


        ! get derivatives
        D(ip1,ip2)=LagrangeDerivative(derivPoints,1,derivEnergies,derivStep)
        if (ip1.EQ.ip2) then
            D(ip2,ip1)=D(ip1,ip2)
        else
            D(ip1,ip2)=D(ip1,ip2)/2
            D(ip2,ip1)=D(ip1,ip2)
        endif

        return
        end subroutine getDeriv

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine huckDensity2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine huckDensity

    use math,     only: LagrangeDerivative
    use printmod, only: prMatrix
    use hdb,      only: perturbate,fieldbd,densitybd

    implicit none

    integer(kind=iglu) :: k,N,fatom,satom,field,sta,sto,derivPoints
    real(kind=rglu)    :: derivStep,zpEnergy,energy(5)
    real(kind=rglu), allocatable :: D(:,:),derivEnergies(:),E(:),H(:,:),V(:,:)


    N=mol%nAtoms
    derivPoints=densitybd%nPoints
    derivStep=densitybd%derivStep; sta=-derivPoints/2; sto= derivPoints/2


    allocate (D(N,N),E(N),derivEnergies(sta:sto)); D=0; derivEnergies=0
    allocate (H(N,N),V(N,N)); H=0; V=0

    call setCore(1,1,0)

    zpEnergy = huckEnergy()
    write (*,*) E

    do fatom = 1,N
    do satom = fatom,N
        derivEnergies(0)=zpEnergy

        do field = sta,sto
            if (field.EQ.0) cycle

            call setCore(fatom,satom,field)
            derivEnergies(field)=huckEnergy()
        enddo

        write (*,'(1X,i2,1X,i2,<derivPoints>(1X,F15.10))') fatom, satom, derivEnergies

        call getDeriv(fatom,satom)
    enddo
    enddo

    call prMatrix(D,6,'Huckel numerical','^.00000',maxwidth=78)

    return

    call setCore(1,1,0)

    call energyHuckel(energy)
    do fatom = 1,N
    do satom = 1,N
        D(fatom,satom)=getHuckRDMElement(fatom,satom)
    enddo
    enddo

    call prMatrix(D,6,'Huckel analytical','^.00000',maxwidth=78)

    call finalizeHuck

    return

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(ip1,ip2,fshift)
        implicit none
        integer(kind=iglu), intent(in) :: ip1,ip2,fshift


        H=mol%connect
        call applyField

        if (ip1.EQ.ip2) then
            H(ip1,ip2)=H(ip1,ip2)+fShift*derivStep
        else
            H(ip1,ip2)=H(ip1,ip2)+fShift*derivStep
            H(ip2,ip1)=H(ip2,ip1)+fShift*derivStep
        endif

        return
        end subroutine setCore

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine applyField
        implicit none

        integer(kind=iglu) :: k,mu


        do k = 1,3
            do mu = 1,N
                H(mu,mu)=H(mu,mu)+fieldbd%strength(k)*mol%atm(mu)%coords(k)
            enddo
        enddo

        return
        end subroutine applyField

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine getDeriv(ip1,ip2)
        implicit none
        integer(kind=iglu), intent(in) :: ip1,ip2


        ! get derivatives
        D(ip1,ip2)=LagrangeDerivative(derivPoints,1,derivEnergies,derivStep)
        if (ip1.EQ.ip2) then
            D(ip2,ip1)=D(ip1,ip2)
        else
            D(ip1,ip2)=D(ip1,ip2)/2
            D(ip2,ip1)=D(ip1,ip2)
        endif

        return
        end subroutine getDeriv

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        real(kind=rglu) function huckEnergy() result(ret)
        implicit none

        integer(kind=iglu) :: mu


        call tred4(H,V,E,N,1.e-100_rglu,1.e-300_rglu)

        ret=0
        do mu = 1,N/2
            ret=ret+2*E(mu)
        enddo

        return
        end function huckEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine huckDensity

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module huckel