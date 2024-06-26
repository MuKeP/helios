    module lrccsdModule

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,      only: assignment(=)
    use glob,      only: uch,rglu,iglu,lglu,true,false,void,i8kind,glControlMemory,nullsub
    use glob,      only: ivarVector,rvarVector
    use hdb,       only: mol,statesbd,lrbd,lrst,ccbd
    use hdb,       only: scfbd,ou,ouWidth,cueConstant1
    use scf,       only: setSCFParameters,initSCF,iterationSCF,getSCFResult,callbackSCF
    use scf,       only: energySCF,finalizeSCF,printSCFSolution
    use txtParser, only: tpFill,operator(.in.)
    use printmod,  only: prEigenProblem,prMatrix
    use math,      only: tred4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: lrVersion='3.100'       !5
    character (len=*), parameter :: lrDate   ='2017.12.10'  !10
    character (len=*), parameter :: lrAuthor ='Vladimir V. Ivanov' !18

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! return accuracy holder
    real(kind=rglu) :: accuracy(5)

    ! amplitudes
    real(kind=rglu), allocatable    :: r1(:,:),r2(:,:,:,:)

    ! deltas
    real(kind=rglu), allocatable    :: d1(:,:),d2(:,:,:,:)

    ! CC solution and integrals
    real(kind=rglu), allocatable    :: t1(:,:),t2(:,:,:,:),F(:,:),R(:,:,:,:)

    ! arrays for integrals calculation
    real   (kind=rglu), allocatable :: G(:,:),hVs(:,:)
    integer   (kind=1), allocatable :: V(:,:)
    integer(kind=iglu), allocatable :: cueIndex(:,:)

    ! orbital accordance for cue
    integer(kind=iglu), allocatable :: iapairs(:)

    ! cue orbitals for integrals

    ! intermediates
    real(kind=rglu), dimension (:,:,:,:), allocatable :: ai1,ai2,ai3,ai4,ai5,ai6,ai7
    real(kind=rglu), dimension     (:,:), allocatable :: Fab,Fij,Fia,at1,at2,at3,at4,at5,at6

    ! diis storage
    real   (kind=rglu), allocatable :: st1(:,:),st2(:,:)
    real   (kind=rglu), allocatable :: sd1(:,:),sd2(:,:)

    ! diis arrays
    real   (kind=rglu), allocatable :: diisVectors(:,:),diisValues(:)
    real   (kind=rglu), allocatable :: diisMatrix(:,:),diisCoefficients(:)

    ! guess preparation
    integer(kind=iglu), allocatable :: excSet(:,:),confs(:)
    real   (kind=rglu), allocatable :: hfV(:,:),cueV(:,:),hfE(:)
    real   (kind=rglu), allocatable :: HH(:,:),Vectors(:,:),Values(:),coefs(:)
    real   (kind=rglu), allocatable :: basisTransformation(:,:)
    type  (ivarVector), allocatable :: guessConfigurations(:)
    type  (rvarVector), allocatable :: guessCoefficients(:)

    ! storage arrays
    real(kind=rglu),    allocatable :: lrHoldStateVectorR2(:,:,:,:,:)
    real(kind=rglu),    allocatable :: lrHoldStateVectorR1(:,:,:)
    real(kind=rglu),    allocatable :: lrHoldStateEnergy(:)
    real(kind=rglu),    allocatable :: lrHoldStateProperties(:,:)
    real(kind=rglu),    allocatable :: lrOrthogonality(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type(uch)          :: umethod
    integer(kind=iglu) :: N,Nocc,Nel,No,Nd,Ne,currentState
    real   (kind=rglu) :: omega,currentEnergy,AEL,RSCon
    logical(kind=lglu) :: converged,guessReady

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    public :: lrVersion,lrDate,lrAuthor
    public :: setLRParameters,initLR,energyLR,finalizeLR,analizewfLR

    ! to harvest the results
    public :: lrHoldStateEnergy,lrHoldStateProperties

    ! access for projection subroutines
    public :: Nel,No
    public :: r1,r2,d1,d2,t1,t2,R,F,iapairs
    public :: Fab,Fij,Fia
    public :: ai1,ai2,ai3,ai4,ai5,ai6,ai7
    public :: at1,at2,at3,at4,at5,at6

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setLRParameters(method)
    implicit none

    character (len=*)  :: method


    if (method .in. ['r-ccsd','u-ccsd','cue-ccsd','spin-cue-ccsd','spin-u-ccsd','spin-r-ccsd'] ) then
        umethod='lr-'//method
    else
        stop 'Internal error (lrccsd::setLRParameters): Unknown method.'
    endif

    N=mol%nAtoms; Nel=mol%nEls; No=2*N; Nocc=Nel/2; Ne=Nocc*(N-Nocc); Nd=lrbd%diisSteps
    guessReady=false

    call controlMemoryLR('general','allocate')

    return
    end subroutine setLRParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine initLR
    implicit none

    integer(kind=iglu) :: k
    logical(kind=lglu) :: convergence


    if (.NOT.guessReady) call prepareGuess

    call getCCSolution

    if (lrbd%diisEnabled) call controlMemoryLR('diis','allocate')
    call controlMemoryLR('intermediates','allocate')

    call projection_lrccsd_spin_hf_intermediates_1

    !write (*,*) 'DONE'
    !stop

    convergence=true

    lrHoldStateVectorR1=0; lrHoldStateVectorR2=0

    lrHoldStateEnergy=0
    do k = 1,statesbd%nStates
        currentState=k

        call guessState(k)

        call iterator(iterationLR,energylr,lrbd%maxiters,lrbd%accuracy*real(2**(k-1),rglu),false,nullsub,false,converged)

        ! convergence=convergence.AND.converged

        lrHoldStateVectorR1  (    :,:,k)=r1
        lrHoldStateVectorR2  (:,:,:,:,k)=r2
        lrHoldStateEnergy    (        k)=omega
        lrHoldStateProperties(      :,k)=[AEL,RSCon]
    enddo

    if (convergence) then
        call dumpSolution
    endif
    call controlMemoryLR('intermediates','deallocate')
    call controlMemoryLR('diis','deallocate')
    call resortStates

    return
    end subroutine initLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine iterationLR(iteration,epsilon,saccuracy)
    implicit none

    integer(kind=iglu) :: iteration
    real   (kind=rglu) :: epsilon,saccuracy(5)


    call projection_lrccsd_spin_hf_intermediates_2
    select case (umethod%get())
        case ('lr-spin-cue-ccsd','lr-cue-ccsd')
            call projection_lrccsd_singles_spin_hf
            call projection_lrccsd_doubles_spin_hf
            !call projection_lrccsd_singles_spin_cue
            !call projection_lrccsd_doubles_spin_cue

        case ('lr-spin-u-ccsd','lr-spin-r-ccsd','lr-u-ccsd','lr-r-ccsd')
            call projection_lrccsd_singles_spin_hf
            call projection_lrccsd_doubles_spin_hf

    end select

    call stepLR; saccuracy=accuracy

    if (maxval(accuracy).LT.epsilon) return

    if (lrbd%diisEnabled) then
        call pushLRVectors

        if (iteration.GT.Nd) then
            call newLRVectors(Nd)
        else
            call newLRVectors(iteration)
        endif
    endif

    return
    end subroutine iterationLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine stepLR
    implicit none

    integer(kind=iglu) :: i,a,j,b
    real(kind=rglu)    :: Ax,rez,s1,s2


    accuracy=-1

    s1=0; s2=0
    !$omp parallel default(shared) private (i,a) reduction(+:s1,s2)
    !$omp do
    do i=1,Nel
    do a=Nel+1,No
        s1=s1+d1(i,a)*r1(i,a)
        s2=s2+r1(i,a)*r1(i,a)
    enddo
    enddo
    !$omp end parallel

    !$omp parallel default(shared) private (i,j,a,b) reduction(+:s1,s2)
    !$omp do
    do i=1,Nel-1
    do j=i+1,Nel
    do a=Nel+1,No-1
    do b=a+1,No
        s1=s1+d2(i,j,a,b)*r2(i,j,a,b)
        s2=s2+r2(i,j,a,b)*r2(i,j,a,b)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel

    omega=s1/s2; currentEnergy=omega

    !$omp parallel default(shared) private (i,a)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        d1(i,a)=d1(i,a)-omega*r1(i,a)
    enddo
    enddo
    !$omp end parallel

    !$omp parallel default(shared) private (i,a)
    !$omp do
    do i=1,Nel
    do a=Nel+1,No
        r1(i,a)=r1(i,a)-lrbd%iterStep(1)*d1(i,a)
    enddo
    enddo
    !$omp end parallel
    accuracy(1)=maxval(abs(d1))

    !$omp parallel default(shared) private (i,j,a,b)
    !$omp do
    do i=1,Nel
    do a=Nel+1,No
    do j=1,Nel
    do b=Nel+1,No
        d2(i,j,a,b)=d2(i,j,a,b)-omega*r2(i,j,a,b)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel

    !$omp parallel default(shared) private (i,j,a,b)
    !$omp do
    do i=1,Nel
    do a=Nel+1,No
    do j=1,Nel
    do b=Nel+1,No
        r2(i,j,a,b)=r2(i,j,a,b)-lrbd%iterStep(2)*d2(i,j,a,b)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel
    accuracy(2)=maxval(abs(d2))

    if(statesbd%spin.EQ.0) then
        !$omp parallel default(shared) private(i,a,rez)
        !$omp do
        do i=1,Nel,2
        do a=Nel+1,No,2
            rez=( r1(i,a)+r1(i+1,a+1) )/2
            r1(i,a)    = rez
            r1(i+1,a+1)= rez
        enddo
        enddo
        !$omp end parallel
    else
        !$omp parallel default(shared) private(i,a,rez)
        !$omp do
        do i=1,Nel,2
        do a=Nel+1,No,2
            rez=( r1(i,a)-r1(i+1,a+1) )/2
            r1(i,a)    = rez
            r1(i+1,a+1)=-rez
        enddo
        enddo
        !$omp end parallel
    endif

    call normalize
    if (lrbd%orthogonalize) call orthogonalization(currentState)

    return
    end subroutine stepLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine energyLR(energy)
    implicit none

    real   (kind=rglu), intent(out) :: energy(5)
    real   (kind=rglu)              :: r0,avgExcitationLvl,sum1,sum2
    integer(kind=iglu)              :: i,j,a,b


    energy=0

    sum1=0
    !$omp parallel default(shared) private(i,a) reduction(+:sum1)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        sum1=sum1+r1(i,a)**2
    enddo
    enddo
    !$omp end parallel

    sum2=0
    !$omp parallel default(shared) private(i,j,a,b) reduction(+:sum2)
    !$omp do
    do i = 1,Nel-1
    do j = i+1,Nel
    do a = Nel+1,No-1
    do b = a+1,No
        sum2=sum2+r2(i,j,a,b)**2
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel
    avgExcitationLvl=(sum1+2*sum2)/(sum1+sum2)

    sum1=0; sum2=0
    !$omp parallel default(shared) private(i,j,a,b) reduction(+:sum1,sum2)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        sum1=sum1+Fia(i,a)*r1(i,a)
        do j = 1,Nel
        do b = Nel+1,No
            sum2=sum2+R(j,b,i,a)*r2(i,j,a,b)
        enddo
        enddo
    enddo
    enddo
    !$omp end parallel
    r0=(sum1+sum2/4)/omega

    energy(1)=currentEnergy
    energy(2)=r0
    energy(3)=avgExcitationLvl

    AEL=avgExcitationLvl
    RSCon=r0

    return
    end subroutine energyLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine normalize
    implicit none

    integer(kind=iglu) :: a,b,i,j
    real   (kind=rglu) :: sum1,sum2


    sum1=0
    !$omp parallel default(shared) private (i,a) reduction(+:sum1)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        sum1=sum1+r1(i,a)**2
    enddo
    enddo
    !$omp end parallel

    sum2=0
    !$omp parallel default(shared) private (i,j,a,b) reduction(+:sum2)
    !$omp do
    do i = 1,Nel-1
    do j = i+1,Nel
    do a = Nel+1,No-1
    do b = a+1,No
        sum2=sum2+r2(i,j,a,b)**2
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel

    sum1=1/sqrt(sum1+sum2)

    !$omp parallel default(shared) private (i,a)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        r1(i,a)=r1(i,a)*sum1
    enddo
    enddo
    !$omp end parallel

    !$omp parallel default(shared) private (i,j,a,b,sum2)
    !$omp do
    do i = 1,Nel-1
    do j = i+1,Nel
    do a = Nel+1,No-1
    do b = a+1,No
        sum2=r2(i,j,a,b)*sum1
        r2(i,j,a,b)= sum2
        r2(j,i,a,b)=-sum2
        r2(i,j,b,a)=-sum2
        r2(j,i,b,a)= sum2
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel

    return
    end subroutine normalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine orthogonalization(cstate)
    implicit none

    integer(kind=iglu) :: i,j,a,b,k,cstate,nstates
    real(kind=rglu)    :: sum1,sum2,normm


    nstates=cstate-1; if (nstates.LE.0) return

    lrOrthogonality=0
    do k = 1,nstates

        sum1=0
        !$omp parallel default(shared) private (i,a) reduction(+:sum1)
        !$omp do
        do i = 1,Nel
        do a = Nel+1,No
            sum1=sum1+r1(i,a)*lrHoldStateVectorR1(i,a,k)
        enddo
        enddo
        !$omp end parallel

        sum2=0
        !$omp parallel default(shared) private (i,j,a,b) reduction(+:sum2)
        !$omp do
        do i = 1,Nel
        do j = 1,Nel
        do a = Nel+1,No
        do b = Nel+1,No
            sum2=sum2+r2(i,j,a,b)*lrHoldStateVectorR2(i,j,a,b,k)
        enddo
        enddo
        enddo
        enddo
        !$omp end parallel

        lrOrthogonality(k)=sum1+sum2
    enddo

    !$omp parallel default(shared) private (i,a,k,sum1)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        sum1=0
        do k = 1,nstates
            sum1=sum1+lrHoldStateVectorR1(i,a,k)*lrOrthogonality(k)
        enddo
        r1(i,a)=r1(i,a)-sum1
    enddo
    enddo
    !$omp end parallel

    !$omp parallel default(shared) private (i,a,j,b,k,sum1)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        do j = 1,Nel
        do b = Nel+1,No
            sum1=0
            do k = 1,nstates
                sum1=sum1+lrHoldStateVectorR2(i,j,a,b,k)*lrOrthogonality(k)
            enddo
            r2(i,j,a,b)=r2(i,j,a,b)-sum1
        enddo
        enddo
    enddo
    enddo
    !$omp end parallel

    return
    end subroutine orthogonalization

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine resortStates
    implicit none

    integer(kind=iglu) :: glcount,swcount,i,j,a,b,k,l
    real(kind=rglu)    :: swEnergy,swVal,swRez(8)


    do
        swcount=0
        do k = 1,statesbd%nStates-1
            l=k+1

            if (lrHoldStateEnergy(l).LT.lrHoldStateEnergy(k)) then
                swcount=swcount+1

                swEnergy=lrHoldStateEnergy(l)
                lrHoldStateEnergy(l)=lrHoldStateEnergy(k)
                lrHoldStateEnergy(k)=swEnergy

                !$omp parallel default(shared) private (i,a,j,b,swVal)
                !$omp do
                do i = 1,Nel
                do a = Nel+1,No
                    swVal=lrHoldStateVectorR1(i,a,l)
                    lrHoldStateVectorR1(i,a,l)=lrHoldStateVectorR1(i,a,k)
                    lrHoldStateVectorR1(i,a,k)=swVal
                    do j = 1,Nel
                    do b = Nel+1,No
                        swVal=lrHoldStateVectorR2(i,j,a,b,l)
                        lrHoldStateVectorR2(i,j,a,b,l)=lrHoldStateVectorR2(i,j,a,b,k)
                        lrHoldStateVectorR2(i,j,a,b,k)=swVal
                    enddo
                    enddo
                enddo
                enddo
                !$omp end parallel

            endif
        enddo
        if (swcount.EQ.0) exit
    enddo

    return
    end subroutine resortStates

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeLR
    implicit none


    call controlMemoryLR('general','deallocate')
    call controlMemoryLR('intermediates','deallocate')
    call controlMemoryLR('diis','deallocate')

    if (umethod%get() .in. ['lr-spin-cue-ccsd','lr-cue-ccsd']) call finalizeSCF

    return
    end subroutine finalizeLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine dumpSolution
    implicit none

    integer(kind=iglu) :: i,j,a,b,k,nonzero(2)


    if (.NOT.lrbd%storeSolution(1)) return

    rewind(lrst)
    do k = 1,statesbd%nStates
        nonzero(1)=0
        do i = 1,Nel
        do a = Nel+1,No
            if (abs(lrHoldStateVectorR1(i,a,k)).GT.lrbd%storeSolutionThreshold) then
                nonzero(1)=nonzero(1)+1
            endif
        enddo
        enddo

        nonzero(2)=0
        do i = 1,Nel-1
        do j = i+1,Nel
        do a = Nel+1,No-1
        do b = a+1,No
            if (abs(lrHoldStateVectorR2(i,j,a,b,k)).GT.lrbd%storeSolutionThreshold) then
                nonzero(2)=nonzero(2)+1
            endif
        enddo
        enddo
        enddo
        enddo

        ! write(ou,'(A,1X,i1,1X,A,i5,1X,i5)') ' >>> state',k,'elements to be dumped',nonzero

        write(lrst) k
        write(lrst) nonzero

        select case (lrbd%storeSolutionMode%get())

            case('r1')
                do i = 1,Nel
                do a = Nel+1,No
                    if (abs(lrHoldStateVectorR1(i,a,k)).GT.lrbd%storeSolutionThreshold) then
                        write(lrst) i,a,lrHoldStateVectorR1(i,a,k)
                    endif
                enddo
                enddo

            case('r2')
                do i = 1,Nel-1
                do j = i+1,Nel
                do a = Nel+1,No-1
                do b = a+1,No
                    if (abs(lrHoldStateVectorR2(i,j,a,b,k)).GT.lrbd%storeSolutionThreshold) then
                        write(lrst) i,j,a,b,lrHoldStateVectorR2(i,j,a,b,k)
                    endif
                enddo
                enddo
                enddo
                enddo

            case('r1r2')
                do i = 1,Nel
                do a = Nel+1,No
                    if (abs(lrHoldStateVectorR1(i,a,k)).GT.lrbd%storeSolutionThreshold) then
                        write(lrst) i,a,lrHoldStateVectorR1(i,a,k)
                    endif
                enddo
                enddo

                do i = 1,Nel-1
                do j = i+1,Nel
                do a = Nel+1,No-1
                do b = a+1,No
                    if (abs(lrHoldStateVectorR2(i,j,a,b,k)).GT.lrbd%storeSolutionThreshold) then
                        write(lrst) i,j,a,b,lrHoldStateVectorR2(i,j,a,b,k)
                    endif
                enddo
                enddo
                enddo
                enddo

        end select
    enddo

    lrbd%storeSolution(2)=true

    return
    end subroutine dumpSolution

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine loadSolution(istate)
    implicit none

    integer(kind=iglu), intent(in) :: istate
    integer(kind=iglu)             :: i,j,a,b,k,c,dmp,nonzero(2)
    real(kind=rglu)                :: value


    if (istate.EQ.1) then
        rewind(lrst)
    endif

    read(lrst) dmp
    read(lrst) nonzero

    select case (lrbd%storeSolutionMode%get())

        case('r1')
            do c = 1,nonzero(1)
                read(lrst) i,a,value
                r1(i,a)=value
            enddo

        case('r2')
            do c = 1,nonzero(2)
                read(lrst) i,j,a,b,value
                r2(i,j,a,b)= value
                r2(j,i,a,b)=-value
                r2(i,j,b,a)=-value
                r2(j,i,b,a)= value
            enddo

        case('r1r2')
            do c = 1,nonzero(1)
                read(lrst) i,a,value
                r1(i,a)=value
            enddo

            do c = 1,nonzero(2)
                read(lrst) i,j,a,b,value
                r2(i,j,a,b)= value
                r2(j,i,a,b)=-value
                r2(i,j,b,a)=-value
                r2(j,i,b,a)= value
            enddo

    end select

    return
    end subroutine loadSolution

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getCCSolution

    use coupledCluster, only: ccmethod => umethod
    use coupledCluster, only: ccNe=> Ne, ccNocc=> Nocc
    use coupledCluster, only: cct1=> t1, cct2=> t2
    use coupledCluster, only: ccF=> F, ccR=> R
    use coupledCluster, only: cciapairs=> iapairs, ccexcSet=> excSet
    use coupledCluster, only: cchV=> hV
    use coupledCluster, only: convertSpatialToSpinCC

    implicit none

    integer(kind=iglu) :: mm,nn,i,j,a,b,k,l
    real   (kind=rglu) :: Ax


    call controlMemoryLR('spin_transformation','allocate')

    if (convertSpatialToSpinCC(ccNocc, Nel, cct1, cct2, t1, t2, 'sd').EQ.-1) then
        stop 'Internal error (lrccsd::getCCSolution): Error while converting spatial to spin orbitals'
    endif

    select case ( ccmethod%get() ) !fockian
        case ('cue-ccsd','u-ccsd','r-ccsd')
            do i = 1,N
            do j = 1,N
                F(2*i-1,2*j-1)=ccF(i,j)
                F(2*i-1,2*j  )=0
                F(2*i  ,2*j-1)=0
                F(2*i  ,2*j  )=ccF(i,j)
            enddo
            enddo

        case ('spin-cue-ccsd','spin-u-ccsd','spin-r-ccsd')
            F=ccF

    end select

    select case ( ccmethod%get() ) !for integrals
        case ('cue-ccsd','spin-cue-ccsd')
            G=mol%G
            ! correction to the calculation of two-electron integrals
            ! in case of the atoms with unshared electron pairs.
            do k = 1,N
            do l = k,N
                if ((mol%atm(k)%nels.EQ.2).AND.(mol%atm(l)%nels.EQ.2)) then
                    G(k,l)=4*G(k,l)
                    G(l,k)=G(k,l)
                    cycle
                endif

                if ((mol%atm(k)%nels.EQ.2).OR. (mol%atm(l)%nels.EQ.2)) then
                    G(k,l)=2*G(k,l)
                    G(l,k)=G(k,l)
                endif
            enddo
            enddo
        case ('r-ccsd','u-ccsd','spin-r-ccsd','spin-u-ccsd')
            G=mol%G

    end select

    select case ( ccmethod%get() ) !integrals
        case ('cue-ccsd')
            do k = 1,ccNocc
                i=mol%orb(k)%atoms(1)
                j=mol%orb(k)%atoms(2)
                l=mol%orb(k)%ova

                if (mol%orb(k)%nels.EQ.2) then
                    V(i,2*k-1)=1; V(i,2*k)=1

                    cueIndex(1,2*k-1)=i; cueIndex(1,2*k)=i
                    cueIndex(2,2*k-1)=j; cueIndex(2,2*k)=j

                    iapairs(2*k-1)=2*l-1; iapairs(2*k)=2*l
                    cycle
                endif

                V(i,2*k-1)=1; V(i,2*l-1)=-1
                V(i,2*k  )=1; V(i,2*l  )=-1
                V(j,2*k-1)=1; V(j,2*l-1)= 1
                V(j,2*k  )=1; V(j,2*l  )= 1

                cueIndex(1,2*k-1)=i; cueIndex(1,2*l-1)=i
                cueIndex(1,2*k  )=i; cueIndex(1,2*l  )=i
                cueIndex(2,2*k-1)=j; cueIndex(2,2*l-1)=j
                cueIndex(2,2*k  )=j; cueIndex(2,2*l  )=j

                iapairs(2*k-1)=2*l-1; iapairs(2*k)=2*l
                iapairs(2*l-1)=2*k-1; iapairs(2*l)=2*k
            enddo

            !$omp parallel default(shared) private(i,j,a,b)
            !$omp do
            do i = 1,No
            do j = 1,No
            do a = 1,No
            do b = 1,No
                R(i,j,a,b)=spin_cue_int(i,j,a,b)/4
            enddo
            enddo
            enddo
            enddo
            !$omp end parallel

        case ('u-ccsd','r-ccsd')
            do i = 1,N
            do j = 1,N
                hVs(i,2*j-1)=cchV(i,j)
                hVs(i,2*j  )=cchV(i,j)
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

        case ('spin-cue-ccsd')
            R=ccR; iapairs=cciapairs

        case ('spin-u-ccsd','spin-r-ccsd')
            R=ccR

    end select

    call controlMemoryLR('spin_transformation','deallocate')

    return
    end subroutine getCCSolution

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine analizewfLR(state)
    implicit none

    integer(kind=iglu) :: state


    r1=lrHoldStateVectorR1(:,:,state)
    r2=lrHoldStateVectorR2(:,:,:,:,state)

    call wfAnalize(umethod%get(), ccbd%wfSwitches, true)

    return
    end subroutine analizewfLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine controlMemoryLR(section,action)
    implicit none

    character  (len=*) :: section,action
    integer(kind=iglu) :: i,j,a,b,vv,err

    !write (*,*) '====> Control memory: '//section//' '//action

    select case (section)

        case ('general')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(3*Nel*(No-Nel)+3*Nel*Nel*(No-Nel)*(No-Nel)+No*No+No*No*No*No+2*N*N+N)+&
                                              iglu*(2*Ne+No)+&
                                              1*(0) ,kind=i8kind),'Linear response module')
                    allocate (r1(Nel,Nel+1:No),r2(Nel,Nel,Nel+1:No,Nel+1:No)); r1=0; r2=0
                    allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No)); t1=0; t2=0
                    allocate (d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No)); d1=0; d2=0
                    allocate (F(No,No),R(No,No,No,No)); F=0; R=0

                    allocate (hfV(N,N),cueV(N,N),hfE(N),excSet(2,Ne),iapairs(No))
                    hfV=0; cueV=0; hfE=0; excSet=0; iapairs=0

                    void=glControlMemory(int( rglu*(N*N+statesbd%nStates*(Nel*(No-Nel)+Nel*Nel*(No-Nel)*(No-Nel))+2*statesbd%nStates)+&
                                              iglu*(0)+&
                                              1*(0) ,kind=i8kind),'Linear response module')

                    allocate (basisTransformation(N,N)); basisTransformation=0

                    allocate (lrHoldStateVectorR1(Nel,Nel+1:No,statesbd%nStates),&
                              lrHoldStateVectorR2(Nel,Nel,Nel+1:No,Nel+1:No,statesbd%nStates),&
                              lrHoldStateEnergy(statesbd%nStates),&
                              lrOrthogonality(statesbd%nStates),&
                              lrHoldStateProperties(2,statesbd%nStates))

                    lrHoldStateVectorR1=0; lrHoldStateVectorR2=0; lrHoldStateEnergy=0; lrOrthogonality=0

                    ! hard to estimate
                    allocate (guessConfigurations(statesbd%nStates),guessCoefficients(statesbd%nStates))

                case ('deallocate')

                    void=glControlMemory(int( sizeof(lrHoldStateEnergy)+sizeof(lrHoldStateVectorR1)+&
                                              sizeof(lrHoldStateVectorR2)+sizeof(lrOrthogonality)+&
                                              sizeof(hfV)+sizeof(cueV)+sizeof(hfE)+sizeof(excSet)+&
                                              sizeof(iapairs)+sizeof(lrHoldStateProperties)&
                                               ,kind=i8kind),'Linear response module', 'free')

                    void=glControlMemory(int( sizeof(r1)+sizeof(r2)+sizeof(t1)+sizeof(t2)+sizeof(d1)+&
                                              sizeof(d2)+sizeof(F)+sizeof(R)+sizeof(basisTransformation)&
                                               ,kind=i8kind),'Linear response module', 'free')

                    deallocate (lrHoldStateEnergy,lrHoldStateVectorR1,lrHoldStateVectorR2,&
                                lrOrthogonality,lrHoldStateProperties,hfV,cueV,hfE,excSet,iapairs,stat=err)

                    deallocate (r1,r2,t1,t2,d1,d2,F,R, stat=err)
                    deallocate (basisTransformation,stat=err)

                    call guessConfigurations%del()
                    call guessCoefficients%del()

                    deallocate (guessConfigurations,guessCoefficients,stat=err)

            end select

        case ('intermediates')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*( (No-Nel)**2+Nel**2+Nel*(No-Nel)+&
                                                     (No-Nel)**2*(Nel**2)+(No-Nel)**1*(Nel**3)+&
                                                     (No-Nel)**3*(Nel**1)+(No-Nel)**0*(Nel**4)+&
                                                     (No-Nel)**4*(Nel**0)+(No-Nel)**1*(Nel**3)+&
                                                     (No-Nel)**3*(Nel**1)+&
                                                     3*Nel**2+3*(No-Nel)**2)&
                                              ,kind=i8kind),'Linear response module. intermediates')

                    allocate (Fab(Nel+1:No,Nel+1:No),Fij(Nel,Nel),Fia(Nel,Nel+1:No))
                    Fab=0; Fij=0; Fia=0

                    allocate (ai1(Nel,Nel+1:No,Nel,Nel+1:No),ai2(Nel,Nel,Nel,Nel+1:No))
                    allocate (ai3(Nel,Nel+1:No,Nel+1:No,Nel+1:No),ai4(Nel,Nel,Nel,Nel))
                    allocate (ai5(Nel+1:No,Nel+1:No,Nel+1:No,Nel+1:No),ai6(Nel,Nel,Nel,Nel+1:No))
                    allocate (ai7(Nel,Nel+1:No,Nel+1:No,Nel+1:No))
                    ai1=0; ai2=0; ai3=0; ai4=0; ai5=0; ai6=0; ai7=0

                    allocate (at1(Nel,Nel),at2(Nel+1:No,Nel+1:No),at3(Nel,Nel))
                    allocate (at4(Nel+1:No,Nel+1:No),at5(Nel,Nel),at6(Nel+1:No,Nel+1:No))
                    at1=0; at2=0; at3=0; at4=0; at5=0; at6=0

                case ('deallocate')
                    deallocate (Fab,Fij,Fia,ai1,ai2,ai3,ai4,ai5,ai6,ai7,at1,at2,at3,at4,at5,at6, stat=err)
                    void=glControlMemory(int( sizeof(Fab)+sizeof(Fij)+sizeof(ai1)+sizeof(ai2)+sizeof(ai3)+&
                                              sizeof(ai4)+sizeof(ai5)+sizeof(ai6)+sizeof(ai7)+sizeof(at1)+&
                                              sizeof(at2)+sizeof(at3)+sizeof(at4)+sizeof(at5)+sizeof(at6)&
                                              ,kind=i8kind),'Linear response module. intermediates', 'free')

            end select

        case ('diis')
            select case (action)
                case ('allocate')
                    vv=0
                    !$omp parallel default(shared) private (i,a) reduction(+:vv)
                    !$omp do
                    do i = 1,Nel
                    do a = Nel+1,No
                        if (btest(a,0).NE.btest(i,0)) cycle
                        vv=vv+1
                    enddo
                    enddo
                    !$omp end parallel

                    void=glControlMemory(int( rglu*(2*vv*Nd) ,kind=i8kind),'Linear response module. diis')
                    allocate (st1(vv,Nd),sd1(vv,Nd))

                    vv=0
                    !$omp parallel default(shared) private (i,j,a,b) reduction(+:vv)
                    !$omp do
                    do i = 1,Nel-1
                    do a = Nel+1,No-1
                    do j = i+1,Nel
                    do b = a+1,No
                        if (btest(a+b,0).NE.btest(i+j,0)) cycle
                        vv=vv+1
                    enddo
                    enddo
                    enddo
                    enddo
                    !$omp end parallel

                    void=glControlMemory(int( rglu*(2*vv*Nd) ,kind=i8kind),'Linear response module. diis')
                    allocate (st2(vv,Nd),sd2(vv,Nd))

                    void=glControlMemory(int( rglu*( 2*(Nd+1)**2 + Nd+1 + Nd ) ,kind=i8kind),'Linear response module. diis')
                    allocate (diisVectors(Nd+1,Nd+1),diisValues(Nd+1))
                    allocate (diisMatrix (Nd+1,Nd+1),diisCoefficients(Nd))

                case ('deallocate')
                    deallocate (diisVectors,diisValues,diisMatrix,diisCoefficients, stat=err)
                    deallocate (st1,sd1,st2,sd2, stat=err)
                    void=glControlMemory(int( sizeof(diisVectors)+sizeof(diisValues)+&
                                              sizeof(diisMatrix)+sizeof(diisCoefficients)+&
                                              sizeof(st1)+sizeof(sd1)+sizeof(st2)+sizeof(sd2)&
                                              ,kind=i8kind),'Linear response module. diis', 'free')

            end select

        case('guess')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(2*(2*Ne)**2+3*(2*Ne)) ,kind=i8kind),'Linear response module. guess')
                    allocate (HH(2*Ne,2*Ne),Vectors(2*Ne,2*Ne),Values(2*Ne),coefs(2*Ne),confs(2*Ne))
                    HH=0; Vectors=0; Values=0; coefs=0; confs=0

                case ('deallocate')
                    deallocate (HH,Vectors,Values,coefs,confs, stat=err)
                    void=glControlMemory(int( sizeof(HH)+sizeof(Vectors)+sizeof(Values)+&
                                              sizeof(coefs)+sizeof(confs)&
                                              ,kind=i8kind),'Linear response module. guess', 'free')

            end select

        case ('spin_transformation')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(N*N+N*No)+iglu*(2*No)+N*No ,kind=i8kind),'Linear response module')
                    allocate (G(N,N),hVs(N,No),cueIndex(2,No),V(N,No))
                    G=0; hVs=0; cueIndex=0; V=0
                case ('deallocate')
                    deallocate (G,hVs,cueIndex,V)
                    void=glControlMemory(int( sizeof(G)+sizeof(hVs)+sizeof(cueIndex)+sizeof(V)&
                                              ,kind=i8kind),'Linear response module', 'free')

            end select

    end select

    return
    end subroutine controlMemoryLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine pushLRVectors
    implicit none

    integer(kind=iglu) :: i,j,k,a,b,vv


    do k = Nd,2,-1
        do vv = 1,UBound(sd1,1)
            sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
        enddo

        do vv = 1,UBound(sd2,1)
            sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
        enddo
    enddo

    vv=0
    do i = 1,Nel
    do a = Nel+1,No
        if (btest(i,0).NE.btest(a,0)) cycle
        vv=vv+1; sd1(vv,1)=d1(i,a); st1(vv,1)=r1(i,a)
    enddo
    enddo

    vv=0
    do i = 1,Nel-1
    do j = i+1,Nel
    do a = Nel+1,No-1
    do b = a+1,No
        if (btest(i+j,0).NE.btest(a+b,0)) cycle
        vv=vv+1; sd2(vv,1)=d2(i,j,a,b); st2(vv,1)=r2(i,j,a,b)
    enddo
    enddo
    enddo
    enddo

    return
    end subroutine pushLRVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine newLRVectors(iteration)
    implicit none

    integer(kind=iglu) :: iteration,k,l,vv,pp,i,j,a,b,c,mm,nn,ub
    real   (kind=rglu) :: sum


    do k = iteration,2,-1
    do l = iteration,2,-1
        diisMatrix(k,l)=diisMatrix(k-1,l-1)
    enddo
    enddo

    do k = 1,iteration
        call diisDotProduct(k,1)
    enddo

    if (iteration.EQ.1) return

    do k = 1,iteration+1
        diisMatrix(k,iteration+1)=1
        diisMatrix(iteration+1,k)=1
    enddo
    diisMatrix(iteration+1,iteration+1)=0

    ub=iteration+1
    call tred4(diisMatrix(1:ub,1:ub),diisVectors(1:ub,1:ub),diisValues(1:ub),ub,1.e-100_rglu,1.e-300_rglu)

    do k = 1,iteration
        diisCoefficients(k)=0
        do l = 1,iteration+1
            diisCoefficients(k)=diisCoefficients(k)+diisVectors(k,l)*diisVectors(iteration+1,l)/diisValues(l)
        enddo
    enddo

!    sum=0
!    do k = 1,iteration
!        sum=sum+diisCoefficients(k)
!    enddo
!    write (*,*) 'sum of diis coefficients',sum

    vv=0
    do i = 1,Nel
    do a = Nel+1,No
        if (btest(i,0).NE.btest(a,0)) cycle
        vv=vv+1; sum=0
        do pp = 1,iteration
            sum=sum+diisCoefficients(pp)*st1(vv,pp)
        enddo
        r1(i,a)=sum
    enddo
    enddo

    vv=0
    do i = 1,Nel-1
    do j = i+1,Nel
    do a = Nel+1,No-1
    do b = a+1,No
        if (btest(i+j,0).NE.btest(a+b,0)) cycle
        vv=vv+1; sum=0
        do pp = 1,iteration
            sum=sum+diisCoefficients(pp)*st2(vv,pp)
        enddo
        r2(i,j,a,b)= sum
        r2(j,i,a,b)=-sum
        r2(i,j,b,a)=-sum
        r2(j,i,b,a)= sum
    enddo
    enddo
    enddo
    enddo

    return
    end subroutine newLRVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine diisDotProduct(k1,k2)
    implicit none

    integer(kind=iglu) :: k1,k2,vv
    real   (kind=rglu) :: sum1,sum2


    sum1=0 ! singles
    !$omp parallel default(shared) private(vv) reduction(+:sum1)
    !$omp do
    do vv = 1,UBound(sd1,1)
        sum1=sum1+sd1(vv,k1)*sd1(vv,k2)
    enddo
    !$omp end parallel

    sum2=0 ! doubles
    !$omp parallel default(shared) private(vv) reduction(+:sum2)
    !$omp do
    do vv = 1,UBound(sd2,1)
        sum2=sum2+sd2(vv,k1)*sd2(vv,k2)
    enddo
    !$omp end parallel

    diisMatrix(k1,k2)=sum1+sum2; diisMatrix(k2,k1)=diisMatrix(k1,k2)

    return
    end subroutine diisDotProduct

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareGuess
    implicit none

    character (len=3)  :: guess
    integer(kind=iglu) :: i,j,a,b,k,l,state,vv
    real   (kind=rglu) :: Ax,Bx,sum


    guessReady=true

    call controlMemoryLR('guess','allocate')

    guess=lrbd%guess%get()
    select case (umethod%get())
        case ('lr-spin-cue-ccsd','lr-cue-ccsd')
            call setSCFParameters
            call initSCF
            call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true,callbackSCF,true,converged)
            call getSCFResult(vectors=hfV,energies=hfE)

        case ('lr-spin-u-ccsd','lr-spin-r-ccsd','lr-u-ccsd','lr-r-ccsd')
            call getSCFResult(vectors=hfV,energies=hfE)
    end select

    k=0
    do i = 1,Nocc
    do j = Nocc+1,N
        k=k+1; excSet(1,k)=i; excSet(2,k)=j
    enddo
    enddo

    !$omp parallel default(shared) private(k,i,a,l,j,b,Ax,Bx)
    !$omp do
    do k = 1,Ne
        i=excSet(1,k); a=excSet(2,k)
        do l = 1,Ne
            j=excSet(1,l); b=excSet(2,l)

            Ax=spat_hf_int(a,i,j,b)

            if ((i.EQ.j).AND.(a.EQ.b)) Ax=Ax+hfE(a)-hfE(i)

            select case (guess)
                case ('cis'); Bx=0
                case ('rpa'); Bx=spat_hf_int(a,i,b,j)
            end select

            HH(k   ,l)=Ax; HH(k+Ne,l+Ne)=Ax
            HH(k+Ne,l)=Bx; HH(k   ,l+Ne)=Bx
        enddo
    enddo
    !$omp end parallel

    call tred4(HH,Vectors,Values,2*Ne,1.e-100_rglu,1.e-300_rglu)

    ! call prMatrix(hfV,6,'HF orbitals','^.00000',maxwidth=79)
    ! call prEigenProblem(Vectors(1:Ne,1:Ne),Values(1:Ne),6,guess//' solution','^.00000',maxwidth=79)

    select case (guess)
        case ('cis')
            do state = 1,statesbd%nStates
                vv=0
                do k = 1,2*Ne
                    Ax=Vectors(k,2*state-1) !first state (by spin)

                    if (abs(Ax).GE.lrbd%guessThreshold) then
                        vv=vv+1
                        coefs(vv)=Ax
                        confs(vv)=k
                    endif
                enddo
                guessConfigurations(state)=confs(1:vv)
                guessCoefficients  (state)=coefs(1:vv)
            enddo

        case ('rpa')
            do state = 1,statesbd%nStates
                vv=0
                do k = 1,Ne
                    Ax=Vectors(k,state)

                    if (abs(Ax).GE.lrbd%guessThreshold) then
                        vv=vv+1
                        coefs(vv)=Ax
                        confs(vv)=k
                    endif
                enddo
                guessConfigurations(state)=confs(1:vv)
                guessCoefficients  (state)=coefs(1:vv)
            enddo

    end select

    do k = 1,Nocc
        l=mol%orb(k)%ova

        i=mol%orb(k)%atoms(1)
        j=mol%orb(k)%atoms(2)

        cueV(i,k)= cueConstant1
        cueV(j,k)= cueConstant1

        i=mol%orb(l)%atoms(1)
        j=mol%orb(l)%atoms(2)

        cueV(i,l)=-cueConstant1
        cueV(j,l)= cueConstant1
    enddo

    basisTransformation=0
    !$omp parallel default(shared) private(i,j,k,sum)
    !$omp do
    do i = 1,N
        do j = 1,N
            sum=0
            do k = 1,N
                sum=sum+cueV(k,i)*hfV(k,j)
            enddo
            basisTransformation(i,j)=sum
        enddo
    enddo
    !$omp end parallel

    call controlMemoryLR('guess','deallocate')

    return
    end subroutine prepareGuess

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine guessState(istate)
    implicit none

    integer(kind=iglu) :: istate,k,l,i,j,a,b
    real(kind=rglu)    :: Ax,Bx


    r1=0; r2=0

!    do i = 1,Nocc
!        r1(2*i-1,2*i-1+Nel)=0.5_rglu
!        if (statesbd%spin.EQ.0) then
!            r1(2*i,2*i+Nel)=0.5_rglu
!        else
!            r1(2*i,2*i+Nel)=-0.5_rglu
!        endif
!    enddo
!
!    return

    if (lrbd%storeSolution(2)) then
        call loadSolution(istate)
        return
    endif

    do k = 1,guessConfigurations(istate)%ln
        l=guessConfigurations(istate)%v(k)
        if (l.GT.Ne) l=l-Ne
        i =excSet(1,l); a =excSet(2,l)
        Ax=guessCoefficients(istate)%v(k)

        !write (70,*) i,a,Ax

        select case (umethod%get())
            case ('lr-spin-cue-ccsd','lr-cue-ccsd')
                do j = 1,Nocc
                do b = Nocc+1,N
                    Bx=r1(2*j-1,2*b-1)+Ax*basisTransformation(j,i)*basisTransformation(b,a)

                    r1(2*j-1,2*b-1)=Bx
                    if (statesbd%spin.EQ.0) then
                        r1(2*j,2*b)= Bx
                    else
                        r1(2*j,2*b)=-Bx
                    endif
                enddo
                enddo

            case ('lr-spin-r-ccsd','lr-spin-u-ccsd','lr-r-ccsd','lr-u-ccsd')
                r1(2*i-1,2*a-1)= Ax
                if (statesbd%spin.EQ.0) then
                    r1(2*i,2*a)= Ax
                else
                    r1(2*i,2*a)=-Ax
                endif

        end select
    enddo

    return
    end subroutine guessState

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function spat_hf_int(a,b,c,d) result(ret) ! 2*[ab|cd]-[ad|cb]
    implicit none

    integer(kind=iglu), intent(in) :: a,b,c,d
    real   (kind=rglu)             :: sum,irt,ir1,ir2
    integer(kind=iglu)             :: mu,nu


    sum=0
    !$omp parallel default(shared) private(mu,ir1,ir2,nu,irt) reduction(+:sum)
    !$omp do
    do mu = 1,N
        ir1=2*hfV(mu,a)*hfV(mu,b)
        ir2=  hfV(mu,a)*hfV(mu,d)
        do nu = 1,N
            irt=( ir1*hfV(nu,d)-ir2*hfV(nu,b) )*hfV(nu,c)
            sum=sum+irt*mol%G(mu,nu)
        enddo
    enddo
    !$omp end parallel

    ret=sum; return
    end function spat_hf_int

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
            sum1=sum1+hVs(mu,a)*hVs(mu,b)*hVs(nu,c)*hVs(nu,d)*G(mu,nu)
        enddo
        enddo
    endif

    sum2=0
    if ((.NOT.btest(a+d,0)).AND.(.NOT.btest(b+c,0))) then ! [ad|bc]
        do mu = 1,N
        do nu = 1,N
            sum2=sum2+hVs(mu,a)*hVs(mu,d)*hVs(nu,b)*hVs(nu,c)*G(mu,nu)
        enddo
        enddo
    endif

    ret=sum1-sum2; return
    end function spin_hf_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function spin_cue_int(a,b,c,d) result(ret) ! [ab|cd]-[ad|cb]
    implicit none

    integer(kind=iglu), intent(in) :: a,b,c,d
    real   (kind=rglu)             :: sum1,sum2
    integer(kind=1)                :: irt,ir1
    integer(kind=iglu)             :: mu,nu

    sum1=0
    if ((.NOT.btest(a+b,0)).AND.(.NOT.btest(c+d,0))) then ! [ab|cd]
        mu=cueIndex(1,a); ir1=    V(mu,a)*V(mu,b)
        nu=cueIndex(1,c); irt=ir1*V(nu,d)*V(nu,c); sum1=     irt*G(mu,nu)
        nu=cueIndex(2,c); irt=ir1*V(nu,d)*V(nu,c); sum1=sum1+irt*G(mu,nu)

        mu=cueIndex(2,a); ir1=    V(mu,a)*V(mu,b)
        nu=cueIndex(1,c); irt=ir1*V(nu,d)*V(nu,c); sum1=sum1+irt*G(mu,nu)
        nu=cueIndex(2,c); irt=ir1*V(nu,d)*V(nu,c); sum1=sum1+irt*G(mu,nu)
    endif

    sum2=0
    if ((.NOT.btest(a+d,0)).AND.(.NOT.btest(c+b,0))) then ! [ad|cb]
        mu=cueIndex(1,a); ir1=    V(mu,a)*V(mu,d)
        nu=cueIndex(1,c); irt=ir1*V(nu,b)*V(nu,c); sum2=     irt*G(mu,nu)
        nu=cueIndex(2,c); irt=ir1*V(nu,b)*V(nu,c); sum2=sum2+irt*G(mu,nu)

        mu=cueIndex(2,a); ir1=    V(mu,a)*V(mu,d)
        nu=cueIndex(1,c); irt=ir1*V(nu,b)*V(nu,c); sum2=sum2+irt*G(mu,nu)
        nu=cueIndex(2,c); irt=ir1*V(nu,b)*V(nu,c); sum2=sum2+irt*G(mu,nu)
    endif

    ret=sum1-sum2; return
    end function spin_cue_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module lrccsdModule
