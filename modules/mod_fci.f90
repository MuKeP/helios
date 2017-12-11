    module fci

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob    , only: true,false,rglu,iglu,lglu,void,glControlMemory,i8kind
    use math    , only: tred4
    use fcontrol, only: fcNewID,fcNullID
    use hdb     , only: mol,fcibd,statesbd,densitybd,ou

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: ciVersion='4.500'       !5
    character (len=*), parameter :: ciDate   ='2017.12.10'  !10
    character (len=*), parameter :: ciAuthor ='Vladimir V. Ivanov'

    integer(kind=iglu), parameter :: maxAtoms=16,maxSteps=8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: fciHoldStateEnergy(:),fciHoldStateOrthog(:)
    real   (kind=rglu), allocatable :: fciHoldStateVector(:,:),fciHoldStateRDM(:,:,:)

    real   (kind=rglu), allocatable :: H(:,:),G(:,:)    ! Hamiltonian, Gamma
    real   (kind=rglu), allocatable :: X(:),Z(:),GAM(:) ! D,GRD,integrals

    real   (kind=rglu), allocatable :: TAU(:,:)

    ! arrays for Knuth sparse algorythm
    real   (kind=rglu), allocatable :: AN(:)
    integer(kind=iglu), allocatable :: JPA(:),NEXTR(:),JR(:),IA(:,:),IND(:)

    ! arrays for iteration procedure
    real   (kind=rglu), allocatable :: eiMas(:,:),eiVec(:,:),eiVal(:)
    real   (kind=rglu), allocatable :: OV(:),suppMass(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) :: fciNSteps,fciNStates,iount

    real   (kind=rglu) :: currentEnergy,coreEnergy
    integer(kind=iglu) :: N,Nocca,Noccb,Nel,KF
    integer(kind=iglu) :: NN,NN1,currentState,GIJ

    logical(kind=lglu) :: doOffDiagonal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    public :: ciVersion,ciDate,ciAuthor
    public :: setFCIParameters,initFCI,energyFCI,finalizeFCI
    public :: getFCIRDM,getFCIRDMElement
    public :: fciHoldStateEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine setFCIParameters
    implicit none

    integer (kind=iglu) :: i,j


    iount=ou
    N=mol%nAtoms; Nocca=mol%naEls; Noccb=mol%nbEls
    fciNSteps=fcibd%nSteps
    fciNStates=statesbd%nStates

    KF=(NOCCA+NOCCB)/2

    if (N.GT.maxAtoms) then; write (iount,'(A)') ' Too big molecule.'; stop; endif

    if (N.EQ.maxAtoms) then
        if (fciNSteps.GT.maxSteps) then
            write (iount,100) fciNSteps,maxSteps
            fciNSteps=maxSteps
        endif
    endif
100 format ('Too much memory is necessary for ',i2,' steps. It will be set ',i2,'.')

    call controlMemoryFCI('general','allocate')

    do i = 1,N
    do j = 1,N
        G(i,j)=mol%G(i,j)
    enddo
    enddo

    return
    end subroutine setFCIParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine initFCI
    implicit none

    integer(kind=iglu) :: I,J,IJ,II
    integer(kind=iglu) :: istate
    real   (kind=rglu) :: AN0


    TAU=0; eiMas=0; eiVec=0; eiVal=0
    fciHoldStateEnergy=0
    fciHoldStateOrthog=0

    H=0
    do i = 1,N
    do j = 1,N
        H(i,j)=mol%core(i,j)
    enddo
    enddo

    call MULTI

    IJ=0
    do I=1,NN
    do J=1,I
        IJ=IJ+1
    enddo
    enddo
    GIJ=IJ

    call controlMemoryFCI('init','allocate')

    coreEnergy=0
    do i= 1,N-1
    do j= i+1,N
        coreEnergy=coreEnergy+G(I,J)
    enddo
    enddo

    do istate = 0,fciNStates
        currentState=istate
        call SUPER
        call FSPARS

        IJ=0
        do I=1,NN
        do J=1,I
           IJ=IJ+1
           X(IJ)=0
           if(I.EQ.J) X(IJ)=1
        enddo
        enddo

        if (istate.GE.2) then
            do I=1,NN
                II=index(i,i)
                X(II)=(-1)**i
            enddo
        endif

        if (istate.EQ.0) then
            if (fciNStates.NE.0) write (iount,100)
        else
            write (iount,101) istate
        endif
        AN0=ANOR(X,NN,NN1)

        call iterator(iteractionFCI,energyFCI,fcibd%maxiters,fcibd%accuracy*real(2**istate,rglu),false)

        !write (*,*) 'fci state',istate,currentEnergy

        fciHoldStateEnergy(  istate)=currentEnergy
        fciHoldStateVector(:,istate)=X

        call controlMemoryFCI('iter','deallocate')
    enddo

    if (fciNStates.GT.2) call resortStates

    call controlMemoryFCI('init','deallocate')

    call prepareRDM

100 format (/18X,'   ==========Ground state==========    '/)
101 format (/18X,'==========Excited state # ',i1,'=========='/)

    return
    end subroutine initFCI

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine iteractionFCI(iteration,epsilon,accuracy)
    implicit none

    integer(kind=iglu), intent(in)  :: iteration
    real   (kind=rglu), intent(in)  :: epsilon
    real   (kind=rglu), intent(out) :: accuracy(5)
    real   (kind=rglu)              :: accur,energy


    accuracy=-1
    call operg(X,Z,energy,accur); accuracy(1)=accur; currentEnergy=energy+coreEnergy
    if(accur.LT.epsilon*real(2**currentState,rglu)) return
    call stepms(X,Z,energy); currentEnergy=energy+coreEnergy
    call orthogonalization(currentState)

    return
    end subroutine iteractionFCI

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine energyFCI(energy)
    implicit none

    real   (kind=rglu) :: energy(5)
    integer(kind=iglu) :: k


    energy=0
    !do k = 0,min(UBound(fciHoldStateEnergy,1),UBound(energy,1))
    !    energy(k+1)=fciHoldStateEnergy(k)
    !enddo
    energy(1)=currentEnergy

    !do k = 1,UBound(fciHoldStateEnergy,1)+1
    !    write (*,*) 'fci energy ready',k,energy(k)
    !enddo

    return
    end subroutine energyFCI

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine prepareRDM

    implicit none

    integer(kind=iglu) :: k


    select case ( densitybd%dtype%get() )
        case ('none')
            fciHoldStateRDM=0

        case ('all','orders')
            doOffDiagonal=true
            do k = 0,fciNStates
                call density(k)
            enddo

        case ('charges')
            doOffDiagonal=false
            do k = 0,fciNStates
                call density(k)
            enddo

    end select

    return
    end subroutine prepareRDM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getFCIRDM(RDM,state)

    implicit none

    integer(kind=iglu), intent(in)  :: state
    real   (kind=rglu), intent(out) :: RDM(:,:)


    RDM=fciHoldStateRDM(:,:,state)

    return
    end subroutine getFCIRDM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function getFCIRDMElement(i,j,state) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: i,j,state


    ret=fciHoldStateRDM(i,j,state)
    return
    end function getFCIRDMElement

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine finalizeFCI
    implicit none


    call controlMemoryFCI('general','deallocate')

    return
    end subroutine finalizeFCI

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine resortStates

    implicit none

    integer(kind=iglu) :: K
    integer(kind=iglu) :: i,j,l,d,swcount,glcount
    real   (kind=rglu) :: swEnergy,swRDM


    glcount=0
    do
        swcount=0
        do i = 1,fciNStates-1
            j=i+1
            if (fciHoldStateEnergy(j).LT.fciHoldStateEnergy(i)) then
                swcount=swcount+1

                swEnergy=fciHoldStateEnergy(j)
                fciHoldStateEnergy(j)=fciHoldStateEnergy(i)
                fciHoldStateEnergy(i)=swEnergy

                do d = 1,N
                do l = 1,N
                    swRDM=fciHoldStateRDM(d,l,j)
                    fciHoldStateRDM(d,l,j)=fciHoldStateRDM(d,l,i)
                    fciHoldStateRDM(d,l,i)=swRDM
                enddo
                enddo

                do K = 1,UBound(X,1)
                    X(K)=fciHoldStateVector(K,j)
                    fciHoldStateVector(K,j)=fciHoldStateVector(K,i)
                    fciHoldStateVector(K,i)=X(K)
                enddo
            endif
        enddo
        glcount=glcount+swcount
        if (swcount.EQ.0) exit
    enddo

    if (glcount.NE.0) write (iount,100)

100 format (4X,'Incorrect order of excited states was fixed.')

    return
    end subroutine resortStates

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine orthogonalization(state)
    implicit none

    integer(kind=iglu) :: I,J,IJ
    integer(kind=iglu) :: state,k
    real   (kind=rglu) :: spr,sen


    if (state.EQ.0) return

    do k = 0,state
        IJ=0; spr=0
        do I=1,NN
        do J=1,I
            IJ=IJ+1
            sen=X(IJ)*fciHoldStateVector(IJ,k)
            if (I.NE.J) sen=sen+sen
            spr=spr+sen
        enddo
        enddo
        fciHoldStateOrthog(k)=spr
    enddo

!    spr=maxval(abs(fciHoldStateOrthog(1:state))); if (spr.LT.fciTol) return
!    write (iount,"(5X,'orthogonality=',1X,ES12.4)") spr

    do k = 0,state
        do I=1,NN1
           X(I)=X(I)-fciHoldStateOrthog(k)*fciHoldStateVector(I,k)
        enddo
    enddo

    return
    end subroutine orthogonalization

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine density(state)
    implicit none

    integer(kind=iglu) :: I,J
    integer(kind=iglu) :: state,kk,ll
    real   (kind=rglu) :: T


    fciHoldStateRDM(:,:,state)=0
    write (iount,'(/4X,A,L/)') 'RDM preparation. Do off-diagonal elemnts:',doOffDiagonal
    do I=1,NN
        if (mod(I,200).EQ.0) write(iount,'(1X,F5.1\)') 100.*I/NN
        do J=I,NN
            T=VSPS(I,J)
            call TAY(I,J)
            if (doOffDiagonal) then
                do kk=1,N
                do ll = kk,N
                    fciHoldStateRDM(kk,ll,state)=fciHoldStateRDM(kk,ll,state)+2*T*TAU(kk,ll)
                enddo
                enddo
                do kk = 1,N
                do ll = kk,N
                    fciHoldStateRDM(ll,kk,state)=fciHoldStateRDM(kk,ll,state)
                enddo
                enddo
            else
                do kk=1,N
                    fciHoldStateRDM(kk,kk,state)=fciHoldStateRDM(kk,kk,state)+2*T*TAU(kk,kk)
                enddo
            endif
        enddo
    enddo
    write (iount,*)

    return
    end subroutine density

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine stepms(D,GRD,ALAM)

    implicit none

    integer(kind=iglu) :: I,K
    integer(kind=iglu) :: ii,inom
    real   (kind=rglu) :: D(:), GRD(:),ALAM, TT,AMA,AN0


    suppMass=0; suppMass(:,1)=D

    eiMas=0
    eiMas(1,1)=ALAM
    eiMas(2,1)=ANOR( GRD, NN, NN1)
    eiMas(1,2)=eiMas(2,1)

    do K = 1,UBound(GRD,1)
        suppMass(K,2)=GRD(K)
    enddo
    do ii=2, fciNSteps
        inom=ii-1
        do K = 1,UBound(GRD,1)
            GRD(K)=suppMass(K,ii)
        enddo
        call OPERG(GRD,OV,TT,AMA)
        eiMas(ii,ii)=TT

        do K = 1,UBound(GRD,1)
            GRD(K)=suppMass(K,inom)
        enddo

        do I=1,NN1
            OV(I)=OV(I)-eiMas(inom,ii)*GRD(I)
        enddo

        if(ii.EQ.fciNSteps) exit

        eiMas(ii,ii+1)=ANOR(OV,NN,NN1)
        eiMas(ii+1,ii)=eiMas(ii,ii+1)
        do K = 1,UBound(OV,1)
            suppMass(K,ii+1)=OV(K)
        enddo
    enddo

    call tred4(eiMas,eiVec,eiVal,fciNSteps,1.e-100_rglu,1.e-300_rglu)

    D=0
    do ii = 1,fciNSteps
        OV=suppMass(:,ii)
        do I=1,NN1
            D(I)=D(I)+eiVec(ii,1)*OV(I)
        enddo
    enddo
    AN0=ANOR(D,NN,NN1)

    return
    end subroutine stepms

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine operg(D,GRD,EEE,AMA)
    implicit none

    integer(kind=iglu) :: I,J,IJ,IE,NX,L,LJ,IL
    real   (kind=rglu) :: D(:),GRD(:)
    real   (kind=rglu) :: EEE,AMA,SUM1,SUM2,SEN


    IJ=0; EEE=0
    do I=1,NN
        !$omp parallel default(shared) private(J,IJ,IE,SUM1,L,LJ,NX,SUM2,IL,SEN) reduction(+:EEE)
        !$omp do
        do J=1,I
            IJ=INDEX(I,J)

            IE=JR(I)
            SUM1=0
            do
                L=JPA(IE)
                LJ=INDEX(L,J)
                SUM1=SUM1+AN(IE)*D(LJ)
                NX=NEXTR(IE)
                if(NX.EQ.0) exit
                IE=NX
            enddo

            IE=JR(J)
            SUM2=0
            do
                L=JPA(IE)
                IL=INDEX(I,L)
                SUM2=SUM2+AN(IE)*D(IL)
                NX=NEXTR(IE)
                if(NX.EQ.0) exit
                IE=NX
            enddo

            GRD(IJ)=SUM1+SUM2+D(IJ)*GAM(IJ)

            SEN=GRD(IJ)*D(IJ)
            if(I.NE.J) SEN=SEN+SEN
            EEE=EEE+SEN
        enddo
        !$omp end parallel
    enddo

    AMA=0
    do i=1,NN1
        GRD(i)=GRD(i)-EEE*D(i)
        if (abs(GRD(i)).GT.AMA) AMA=abs(GRD(i))
    enddo

    return
    end subroutine operg

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine TAY(I,J)
    implicit none

    integer(kind=iglu) :: K,L,I,J,ISI,ISJ,IP,JP,IAA,IBB,IIA,IIB,KOL,M


    TAU=0; KOL=0
    do K=1,KF
    do L=1,KF
        if(IA(I,K).EQ.IA(J,L)) KOL=KOL+1
    enddo
    enddo

    if((KF-KOL).EQ.1) then
        ISI=0; ISJ=0
        IP=1;  JP=1
        do K=1,KF
            ISI=ISI+IA(I,K)
            ISJ=ISJ+IA(J,K)
            IP=IP*IA(I,K)
            JP=JP*IA(J,K)
        enddo
        IAA=(IP*(ISI-ISJ))/(IP-JP)
        IBB=(JP*(ISI-ISJ))/(IP-JP)
        do K=1,KF
            if(IAA.EQ.IA(I,K)) IIA=K
            if(IBB.EQ.IA(J,K)) IIB=K
        enddo
        TAU(IAA,IBB)=(-1)**(IIA+IIB)
    endif
    if(KOL.EQ.KF) then
         do M=1,KF
            K=IA(I,M)
            TAU(K,K)=1
        enddo
    endif

    return
    end subroutine TAY

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function ANOR(D,NN,NN1)
    implicit none

    integer(kind=iglu) :: IJ,I,J,NN,NN1
    real   (kind=rglu) :: D(:)
    real   (kind=rglu) :: SS,SEN


    SS=0; IJ=0
    do I=1,NN
        do J=1,I
            IJ=IJ+1
            SEN=D(IJ)*D(IJ)
            if(I.NE.J) SEN=SEN+SEN
            SS=SS+SEN
        enddo
    enddo

    ANOR=sqrt(SS)
    SS=1/ANOR

    do I=1,NN1
        D(I)=SS*D(I)
    enddo

    return
    end function ANOR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine fspars
    implicit none

    integer(kind=iglu) :: I,J,IJ,NJ


    do I=1,NN
        NJ=0
        do J=1,NN
            IJ=INDEX(I,J)
            if( abs(Z(IJ)).LT.fcibd%zeroThreshold) cycle
            NEL=NEL+1
        enddo
    enddo

    call controlMemoryFCI('iter','allocate')

    AN=0; JPA=0; NEXTR=0; JR=0

    NEL=0
    do I=1,NN
        NJ=0
        do J=1,NN
            IJ=INDEX(I,J)
            if( abs(Z(IJ)).LT.fcibd%zeroThreshold) cycle
            NEL=NEL+1
            NJ=NJ+1
            AN(NEL)=Z(IJ)

            JPA(NEL)=J
            if(NJ.NE.1) then
                NEXTR(NEL-1)=NEL
                cycle
            endif
            JR(I)=NEL
        enddo
    enddo

    return
    end subroutine fspars

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine super
    implicit none

    integer(kind=iglu) :: ILA,ILB,IP,JP,ISS,JSS,LAA,LBB,LL
    integer(kind=iglu) :: I,J,K,L,M,IJ,K1,L1,KK,M1,KOL
    real   (kind=rglu) :: GAM1,GAM2


    GAM=0; X=0; Z=0; K1=KF-1; IJ=0
    do I=1,NN
    do J=1,I
        IJ=IJ+1; Z(IJ)=0; KOL=0
        do M=1,KF
        do L=1,KF
            if (IA(I,M).EQ.IA(J,L)) KOL=KOL+1
        enddo
        enddo

        if (KOL+1.EQ.KF) then
            ISS=0; JSS=0; IP=1; JP=1

            do L=1,KF
                ISS=ISS+IA(I,L)
                JSS=JSS+IA(J,L)
                IP=IP*IA(I,L)
                JP=JP*IA(J,L)
            enddo
            LAA=int(IP,kind=8)*(ISS-JSS)/(IP-JP)
            LBB=int(JP,kind=8)*(ISS-JSS)/(IP-JP)

            if (dabs(H(LAA,LBB)).GE.fcibd%zeroThreshold) then
                do L=1,KF
                    if(LAA.eq.IA(I,L)) ILA=L
                    if(LBB.eq.IA(J,L)) ILB=L
                enddo
                LL=ILA+ILB
                Z(IJ)=((-1)**LL)*H(LAA,LBB)
            endif
        endif

        if(KOL.EQ.KF) then
            do K=1,KF
                KK=IA(I,K)
                Z(IJ)=Z(IJ)+H(KK,KK)
            enddo
        endif

        GAM1=0
        if(K1.NE.0) then
            do M=1,K1
            do L=M+1,KF
                GAM1=GAM1+G(IA(I,M),IA(I,L))+G(IA(J,M),IA(J,L))
            enddo
            enddo
        endif

        GAM2=0
        do M=1,KF
        do L=1,KF
           M1=IA(I,M)
           L1=IA(J,L)
           GAM2=GAM2+G(M1,L1)
        enddo
        enddo

        GAM(IJ)= GAM1+GAM2
    enddo
    enddo

    return
    end subroutine super

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine multi
    implicit none

    integer(kind=iglu)              :: I,L,P
    integer(kind=iglu)              :: fileid


    fileid=fcNewID(); open (fileid,status='scratch',form='binary')

    NN=1
    do I=1,KF
        IND(I)=I
        write (fileid) NN,I,I
    enddo

    P=KF
    do
        if(P.LT.1) exit
        if(IND(KF).EQ.N) then
            P=P-1
        else
            P=KF
        endif

        if(P.LT.1) exit

        I=KF+1
        do
            I=I-1
            if(I.LT.P) exit
            IND(I)=IND(P)+I-P+1
        enddo
        NN=NN+1
        do L=1,KF
            write (fileid) NN,L,IND(L) !IA(NN,L)=IND(L)
        enddo
    enddo

    if (allocated(IA)) then
        deallocate (IA)
        void=glControlMemory(int( sizeof(IA) ,kind=i8kind),'tmp. FCI module', 'free')
    endif
    void=glControlMemory(int( iglu*NN*KF ,kind=i8kind),'tmp. FCI module')
    allocate (IA(NN,KF))
    rewind(fileid)

    do
        if (eof(fileid)) exit
        read (fileid) NN,L,P; IA(NN,L)=P
    enddo

    close (fileid); void=fcNullID(fileid)

    NN1=NN*(NN+1)/2

    return
    end subroutine multi

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function index(I,J)
    implicit none

    integer(kind=iglu) :: I,J


    if(I.GT.J) then
        index=I*(I-1)/2+J
    else
        index=J*(J-1)/2+I
    endif

    return
    end function index

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function VSPS(I,J)
    implicit none

    integer(kind=iglu) :: I,J,K,IK,JK
    real   (kind=rglu) :: T


    T=0
    do K=1,NN
        IK=INDEX(I,K)
        JK=INDEX(J,K)
        T=T+X(IK)*X(JK)
    enddo

    VSPS=T; return
    end function VSPS

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine controlMemoryFCI(section,action)
    implicit none

    character (len=*)   :: section,action
    integer (kind=iglu) :: err


    select case (section)
        case ('general')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(N*N*(3+fciNStates+1)+fciNsteps*(fciNSteps+1))+iglu*(KF)+2*(fciNStates+1) ,kind=i8kind),'FCI module')
                    allocate(G(N,N),H(N,N),TAU(N,N),IND(KF),fciHoldStateRDM(N,N,0:fciNStates))
                    allocate(eiMas(fciNSteps,fciNSteps),eiVec(fciNSteps,fciNSteps),eiVal(fciNSteps))
                    allocate(fciHoldStateEnergy(0:fciNStates),fciHoldStateOrthog(0:fciNStates))

                    G=0; H=0; TAU=0; IND=0
                    eiMas=0; eiVec=0; eiVal=0
                    fciHoldStateEnergy=0; fciHoldStateOrthog=0; fciHoldStateRDM=0

                case ('deallocate')
                    deallocate(G,H,TAU,eiMas,eiVec,eiVal,IND,stat=err)
                    deallocate(fciHoldStateEnergy,fciHoldStateOrthog,fciHoldStateRDM,stat=err)
                    void=glControlMemory(int( sizeof(G)+sizeof(H)+sizeof(TAU)+sizeof(eiMas)+&
                                              sizeof(eiVec)+sizeof(eiVal)+sizeof(IND)+&
                                              sizeof(fciHoldStateEnergy)+sizeof(fciHoldStateOrthog)+&
                                              sizeof(fciHoldStateRDM) ,kind=i8kind),'FCI module', 'free')
            end select

        case ('init')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(3*GIJ+NN1*fciNSteps+NN1+GIJ*(fciNStates+1)) ,kind=i8kind),'FCI module')
                    allocate (GAM(GIJ),X(GIJ),Z(GIJ)); GAM=0; X=0; Z=0
                    allocate (suppMass(NN1,fciNSteps),OV(NN1)); suppMass=0; OV=0
                    allocate (fciHoldStateVector(GIJ,0:fciNStates))

                    fciHoldStateVector=0

                case ('deallocate')
                    deallocate (GAM,X,Z,suppMass,OV,fciHoldStateVector)
                    void=glControlMemory(int( sizeof(GAM)+sizeof(X)+sizeof(Z)+sizeof(OV)+&
                                              sizeof(suppMass)+sizeof(fciHoldStateVector) ,kind=i8kind),'FCI module', 'free')
            end select

        case ('iter')
            select case (action)
                case ('allocate')
                    void=glControlMemory(int( rglu*(NEL)+iglu*(NEL+NEL+NN) ,kind=i8kind),'FCI module')
                    allocate(AN(NEL),JPA(NEL),NEXTR(NEL),JR(NN))

                case ('deallocate')
                    deallocate (AN,JPA,NEXTR,JR)
                    void=glControlMemory(int( sizeof(AN)+sizeof(JPA)+sizeof(NEXTR)+sizeof(JR) ,kind=i8kind),'FCI module', 'free')
            end select
    end select

    return
    end subroutine controlMemoryFCI

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module fci