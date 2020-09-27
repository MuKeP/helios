    module scf

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,      only: iglu,rglu,lglu,void,true,false,uch,i8kind,glControlMemory
    use hdb,       only: mol,scfbd,ou,ouWidth,iterationbd,ipFailed,dipoleToDeby,throughbd
    use hdb,       only: singleSession,localbd
    use localizer, only: localize
    use txtParser, only: tpFill,tpAdjustc
    use math,      only: tred4
    use printmod,  only: prEigenProblem,prMatrix,prStrByVal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: hfVersion='1.100'
    character (len=*), parameter :: hfDate   ='2017.12.10'
    character (len=*), parameter :: hfAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: F(:,:),V(:,:),E(:),D(:,:)
    real   (kind=rglu), allocatable :: Fmin(:,:),Comut(:,:),Refl(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) :: N,Nocc
    logical(kind=lglu) :: hasGuess,saveKeepStatus,localized
    real   (kind=rglu) :: iterstep,variation,Eel,Enuc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: setSCFParameters,initSCF,iterationSCF,energySCF,getSCFResult,&
              finalizeSCF,printSCFSolution,getSCFRDMElement,callbackSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setSCFParameters
    implicit none


    N=mol%nAtoms; Nocc=mol%nEls/2

    call controlMemorySCF('general','allocate')
    hasGuess=false

    variation=0
    scfbd%exceeded=false

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

    hasGuess=true
    localized=false

    if ((iterationbd%lastProcedureStatus.GE.0).OR.(scfbd%exceeded)) then
        variation=0
    endif

    call prepareFockian

    return
    end subroutine initSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!    subroutine showGuess
!    implicit none

!    real(kind=rglu), allocatable :: V(:,:),E(:)


!    allocate (V(N,N),E(N)); V=0; E=0

!    call tred4(mol%huckelCore,V,E,N,1.e-100_rglu,1.e-300_rglu)
!    call prepareDensity(V)

!    call prMatrix(V,100,'MO LKAO','^.0000',maxwidth=200)
!    call prMatrix(D,100,'Density','^.0000',maxwidth=200)

!    deallocate(V,E)

!    return
!    end subroutine showGuess

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
            D(mu,nu)=D(mu,nu)-(scfbd%iterStep+variation)*Fmin(mu,nu)
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

    subroutine callbackSCF
    implicit none


    variation=variation+scfbd%iterStepChange
    if (scfbd%achieveSolution) then
        if (variation.LT.scfbd%iterStepVariation) then
            scfbd%exceeded=false
            write(ou,'(A,F5.3,A)') '**** Restarting SCF procedure. Current step: ',&
                                   scfbd%iterStep+variation,' ****'
            saveKeepStatus=scfbd%keep
            scfbd%keep=false
            call initSCF
            scfbd%keep=saveKeepStatus
            return
        endif
    endif

    iterationbd%doRestart=false
    scfbd%exceeded=true

    return
    end subroutine callbackSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine energySCF(energy)
    implicit none

    real   (kind=rglu) :: energy(5)
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

    subroutine getSCFResult(fockian,vectors,energies,density)
    implicit none

    real   (kind=rglu), optional :: fockian(:,:),vectors(:,:),energies(:),density(:,:)
    integer(kind=iglu)           :: i,j


    if (localbd%enabled) then
        if (.NOT.localized) then
            call localize(V,V,true)
            localized=true
        endif
    endif

    if (present(fockian)) then
        if (UBound(fockian,1).NE.N) then
            stop 'Internal error (scf:: getSCFResult): Incorrect size of Fockian matrix.'
        endif

        !$omp parallel default(shared) private(i,j)
        !$omp do
        do i = 1,N
            do j = 1,N
                fockian(i,j)=F(i,j)
            enddo
        enddo
        !$omp end parallel
    endif

    if (present(vectors)) then
        if (UBound(vectors,1).NE.N) then
            stop 'Internal error (scf:: getSCFResult): Incorrect size of Fockian eigenvectors.'
        endif

        !$omp parallel default(shared) private(i,j)
        !$omp do
        do i = 1,N
            do j = 1,N
                vectors(i,j)=V(i,j)
            enddo
        enddo
        !$omp end parallel
    endif

    if (present(energies)) then
        if (UBound(energies,1).NE.N) then
            stop 'Internal error (scf:: getSCFResult): Incorrect size of Fockian eigenvalues.'
        endif

        !$omp parallel default(shared) private(i)
        !$omp do
        do i = 1,N
            energies(i)=E(i)
        enddo
        !$omp end parallel
    endif

    if (present(density)) then
        if (UBound(density,1).NE.N) then
            stop 'Internal error (scf:: getSCFResult): Incorrect size of Fockian eigenvalues.'
        endif

        !$omp parallel default(shared) private(i,j)
        !$omp do
        do i = 1,N
            do j = 1,N
                density(i,j)=2*D(i,j)
            enddo
        enddo
        !$omp end parallel
    endif

    return
    end subroutine getSCFResult

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function getSCFRDMElement(i,j) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: i,j


    ret=2*D(i,j); return
    end function getSCFRDMElement

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine printSCFSolution
    implicit none

    integer(kind=iglu) :: c,paccur
    real   (kind=rglu) :: dx,dy,dz,dmod


    paccur=5

    write (ou,'(/A)') tpAdjustc(' RHF results ',ouWidth,'=')
    call prEigenProblem(V,E,ou,'RHF solution','^.0000',maxwidth=ouWidth)
    write (ou,99) Eel,Enuc,Eel+Enuc,E(Nocc+1)-E(Nocc)
    call prMatrix(2*D,ou,'RHF RDM1', '^.0000',maxwidth=ouWidth)

    ! write(*,*) throughbd%enabled(1), throughbd%property%get()
    if (throughbd%enabled(1) .AND. (throughbd%property%get().EQ.'gap')) then
         call singleSession('  '//prStrByVal(E(Nocc+1)-E(Nocc),'__.0000')//'  ')
    endif

    write (ou,100) 'Atom',tpAdjustc('Density',paccur+3),tpAdjustc('Charge',paccur+4)
    dx=0; dy=0; dz=0
    do c = 1,N
        write (ou,101) c,2*D(c,c),1-2*D(c,c)
        dx=dx+mol%atm(c)%coords(1)*(1-2*D(c,c))
        dy=dy+mol%atm(c)%coords(2)*(1-2*D(c,c))
        dz=dz+mol%atm(c)%coords(3)*(1-2*D(c,c))
    enddo
    dmod=sqrt(dx**2+dy**2+dz**2)
    write (ou,102) 'X',dx*dipoleToDeby,'Y',dy*dipoleToDeby,&
                   'Z',dz*dipoleToDeby,'|D|',dmod*dipoleToDeby

    write (ou,'(/A)') tpAdjustc(' End of RHF results ',ouWidth,'=')

 99 format(4X,'Electron energy:',1X,F16.8,1X,'eV'/&
           4X,'Nuclear energy: ',1X,F16.8,1X,'eV'/&
           4X,'Total energy:   ',1X,F16.8,1X,'eV'/&
           4X,'HOMO-LUMO gap:  ',5X,F8.4,5X,'eV')
100 format(/2X,A,2X,A,3X,A)
101 format(2X,i4,2X,F<3+paccur>.<paccur>,2X,F<4+paccur>.<paccur>)
102 format(<17+2*paccur>('_')//3(<10+paccur>X,A,' =',F<4+paccur>.<paccur>,1X,'D'/),&
                                  <8+paccur>X,A,' =',F<4+paccur>.<paccur>,1X,'D'/)

    return
    end subroutine printSCFSolution

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeSCF
    implicit none


    call controlMemorySCF('general','deallocate')
    hasGuess=false

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

    if ((.NOT.scfbd%keep).OR.(.NOT.hasGuess)) then
        select case (guess)
            case ('huckel')
                call tred4(mol%huckelCore,V,E,N,1.e-100_rglu,1.e-300_rglu)
                call prepareDensity(V)

            case ('unitmatrix')
                D=0
                !$omp parallel default(shared) private(mu)
                !$omp do
                do mu = 1,N
                    D(mu,mu)=0.5_rglu
                enddo
                !$omp end parallel

            case ('defined')
                call prepareDensity(vecs)

        end select
    endif

    return
    end subroutine guessSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareDensity(vec)
    implicit none

    real(kind=rglu)    :: vec(N,N)
    real(kind=rglu)    :: sum
    integer(kind=iglu) :: i,mu,nu


    D=0
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


    F=0
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


    Refl=0
    !$omp parallel default(shared) private(i,j)
    !$omp do
    do i = 1,N
        do j = 1,N
            Refl(i,j)=2*D(i,j)
        enddo
        Refl(i,i)=Refl(i,i)-1
    enddo
    !$omp end parallel

    Comut=0
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

    Fmin=0
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