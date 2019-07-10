    module property

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,           only: assignment(=),glGetIOunit,glSetIOunit,purifyValues,glControlMemory
    use glob,           only: rglu,iglu,lglu,true,false,uch,void,gluCompare,timeControlCheckpoint
    use glob,           only: i8kind,i4kind,r8kind,nullsub
    use txtParser,      only: tpSplit,tpSplitLen,tpSplitHold,tpFill,tpAdjustc,operator(.in.)
    use txtParser,      only: tpIntByStr,tpIntArrayByStr,tpCount
    use printmod,       only: prStrByVal,prMatrix,prEigenProblem,prLongText,prJoin
    use math,           only: LagrangeDerivative,tred4
    use hdb,            only: mol,scfbd,ccbd,statesbd,generalbd,polarizbd,ou,ouWidth,gEnergyHolder
    use hdb,            only: pointToCalc,pointSet,GlEt,perturbate,perturbationID,getMethodNumber
    use hdb,            only: noPerturbation,densitybd,hyperchargesbd,atomEqu,bondEqu
    use hdb,            only: HartreeEnergy,BohrRadius,dipoleToDeby,coulsonbd,systembd,fieldbd
    use hdb,            only: singleSession,setIterationHeader
    use huckel,         only: getHuckelResult,getHuckRDMElement,setHuckelParameters,energyHuckel
    use huckel,         only: finalizeHuck
    use scf,            only: setSCFParameters,initSCF,iterationSCF,energySCF,getSCFRDMElement
    use scf,            only: getSCFResult,finalizeSCF,printSCFSolution,callbackSCF
    use mbpt,           only: setMBPTParameters,initMBPT,energyMBPT,finalizeMBPT
    use coupledCluster, only: setCCParameters,initCC,iterationCC,energyCC,finalizeCC,analizewfCC
    use coupledCluster, only: getCCResults
    use fci,            only: setFCIParameters,initFCI,energyFCI,finalizeFCI,fciHoldStateEnergy
    use fci,            only: getFCIRDM,getFCIRDMElement
    use excitedStates,  only: getExcitationEnergy
    use lrccsdModule,   only: setLRParameters,initLR,energyLR,finalizeLR
    use lrccsdModule,   only: lrHoldStateEnergy,analizewfLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type(uch), allocatable  :: methodSet(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical(kind=lglu)              :: firstRun=true,excitedReady,converged
    type(uch)                       :: lastMethod,iterationHeader
    integer(kind=iglu)              :: lastPerturbationID=-1
    real   (kind=rglu), allocatable :: lastEnergyHolder(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: getEnergy,getPolarizability,getRDM,getWaveFunctionAnalize,getHypercharges,getCoulson
    public :: getResults,getOnlyEnergies

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function getEnergy(method,state) result(ret)
    implicit none

    character  (len=*), intent(in)            :: method
    integer(kind=iglu), intent(in), optional  :: state
    integer(kind=iglu)                        :: dstate

    integer(kind=iglu)                        :: N,i

    real   (kind=rglu)                        :: excEnergy
    real   (kind=rglu)                        :: holdEnergy(5)
    real   (kind=rglu), allocatable           :: V(:,:),E(:)


    ! write (*,*) '*** Method before '//lastMethod%get()
    ! write (*,*) '*** New method '//method
    ! write (*,*) '*** First run? ',firstRun
    ! write (*,*) '*** PerturbationID',perturbationID

    if (.NOT.allocated(lastEnergyHolder)) then
        allocate (lastEnergyHolder(0:statesbd%nStates))
        lastEnergyHolder=0
    endif

    dstate=0; if (present(state)) dstate=state

    if (dstate.GT.statesbd%nStates) then
        ret=0
        return
    endif

    if ((lastPerturbationID.EQ.perturbationID).AND.(method.EQ.lastMethod%get())) then
        if (dstate.EQ.0) then
            ret=lastEnergyHolder(0)
            return
        elseif ((dstate.GT.0).AND.(excitedReady)) then
            ret=lastEnergyHolder(dstate)
            return
        else
            continue
        endif
    endif

    N=mol%nAtoms; allocate (V(N,N),E(N)); V=0; E=0
    do
        if (firstRun) then
            lastMethod=method
            lastEnergyHolder=0

            select case (method)

                case ('huckel') !done
                    call setHuckelParameters
                    call energyHuckel(holdEnergy)

                    lastEnergyHolder(0)=holdEnergy(1)

                    if (dstate.GT.0) then
                        ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
                        excitedReady=true
                    else
                        ret=holdEnergy(1)
                        excitedReady=false
                    endif
                    lastPerturbationID=perturbationID

                case ('rhf') !done
                    call setSCFParameters
                    call initSCF
                    call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,false,callbackSCF,true,converged)
                    call energySCF(holdEnergy)
                    call getSCFResult(vectors=V,energies=E)

                    !call printSCFSolution

                    lastEnergyHolder(0)=holdEnergy(1)

                    if (dstate.GT.0) then
                        ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
                        excitedReady=true
                    else
                        ret=holdEnergy(1)
                        excitedReady=false
                    endif
                    lastPerturbationID=perturbationID

                case ('mp2','mp3') !done
                    call setMBPTParameters(method)
                    call initMBPT
                    call energyMBPT(holdEnergy)

                    lastEnergyHolder(0)=holdEnergy(1)

                    if (dstate.GT.0) then
                        ret=0; lastEnergyHolder(1:statesbd%nStates)=0
                        excitedReady=true
                    else
                        ret=holdEnergy(1)
                        excitedReady=false
                    endif
                    lastPerturbationID=perturbationID

                case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd',&
                      'r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt') !done
                    call setCCParameters(method)
                    call initCC
                    call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false,nullsub,false,converged)
                    call energyCC(holdEnergy)

                    lastEnergyHolder(0)=holdEnergy(1)

                    if ((method .in. ['cue-ccsd','r-ccsd','u-ccsd']).AND.(dstate.GT.0)) then
                        call setLRParameters(method)
                        call initLR

                        ret=lrHoldStateEnergy(dstate)
                        lastEnergyHolder(1:statesbd%nStates)=lrHoldStateEnergy(1:statesbd%nStates)
                        excitedReady=true
                    else
                        ret=holdEnergy(1)
                        excitedReady=false
                    endif
                    lastPerturbationID=perturbationID

                case ('fci') !done
                    call setFCIParameters
                    call initFCI
                    call energyFCI(holdEnergy)

                    !lastEnergyHolder(0)=holdEnergy(1)
                    lastEnergyHolder(0)=fciHoldStateEnergy(0)

                    excitedReady=true
                    lastEnergyHolder(1:statesbd%nStates)=fciHoldStateEnergy(1:statesbd%nStates)-&
                                                         fciHoldStateEnergy(0)
                    ret=lastEnergyHolder(dstate)
                    lastPerturbationID=perturbationID

                case default
                    ret=0
                    lastEnergyHolder=0
                    excitedReady=false
                    lastPerturbationID=perturbationID

            end select
            !void=putMethodEnergy(ret, method, dstate)
            firstRun=false
            exit
        else
            if (method.EQ.lastMethod%get()) then
                select case (method)

                    case ('huckel') !done
                        call energyHuckel(holdEnergy)

                        lastEnergyHolder(0)=holdEnergy(1)

                        if (dstate.GT.0) then
                            ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
                            excitedReady=true
                        else
                            ret=holdEnergy(1)
                            excitedReady=false
                        endif
                        lastPerturbationID=perturbationID

                    case ('rhf') !done
                        call initSCF
                        call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,false,callbackSCF,true,converged)
                        call energySCF(holdEnergy)
                        call getSCFResult(vectors=V,energies=E)

                        lastEnergyHolder(0)=holdEnergy(1)

                        if (dstate.GT.0) then
                            ret=getExcitationEnergy(V,E,dstate,lastEnergyHolder(1:statesbd%nStates))
                            excitedReady=true
                        else
                            ret=holdEnergy(1)
                            excitedReady=false
                        endif
                        lastPerturbationID=perturbationID

                    case ('mp2','mp3') !done
                        call initMBPT
                        call energyMBPT(holdEnergy)

                        lastEnergyHolder(0)=holdEnergy(1)

                        if (dstate.GT.0) then
                            ret=0; lastEnergyHolder(1:statesbd%nStates)=0
                            excitedReady=true
                        else
                            ret=holdEnergy(1)
                            excitedReady=false
                        endif
                        lastPerturbationID=perturbationID

                    case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd',&
                          'r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt') !done
                        call initCC
                        call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false,nullsub,false,converged)
                        call energyCC(holdEnergy)

                        lastEnergyHolder(0)=holdEnergy(1)

                        if ((method .in. ['cue-ccsd','r-ccsd','u-ccsd']).AND.(dstate.GT.0)) then
                            call initLR

                            ret=lrHoldStateEnergy(dstate)
                            lastEnergyHolder(1:statesbd%nStates)=lrHoldStateEnergy(1:statesbd%nStates)
                            excitedReady=true
                        else
                            ret=holdEnergy(1)
                            excitedReady=false
                        endif
                        lastPerturbationID=perturbationID

                    case ('fci') !done
                        call initFCI
                        call energyFCI(holdEnergy)

                        lastEnergyHolder(0)=holdEnergy(1)

                        excitedReady=true
                        if (dstate.GT.0) then
                            lastEnergyHolder(1:statesbd%nStates)=fciHoldStateEnergy(1:statesbd%nStates)-fciHoldStateEnergy(0)
                            ret=lastEnergyHolder(dstate)
                        else
                            ret=lastEnergyHolder(0)
                        endif
                        lastPerturbationID=perturbationID

                    case default
                        ret=0
                        lastEnergyHolder=0
                        excitedReady=false
                        lastPerturbationID=perturbationID

                end select
                !void=putMethodEnergy(ret, method, dstate)
                exit
            else
                select case (lastMethod%get())

                    case ('huckel')
                        call finalizeHuck

                    case ('rhf')
                        call finalizeSCF

                    case ('mp2','mp3')
                        call finalizeMBPT

                    case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd',&
                          'r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt')
                        call finalizeCC
                        if (excitedReady) call finalizeLR

                    case ('fci')
                        call finalizeFCI

                    case default
                        continue

                end select
                firstRun=true
                cycle
            endif
        endif
    enddo
    deallocate (V,E)

    return
    end function getEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getOnlyEnergies
    implicit none


    ! TODO

    return
    end subroutine getOnlyEnergies

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getResults
    implicit none


    ! TODO

    select case (lastMethod%get())

        case ('huckel')

        case ('rhf')
            call printSCFSolution

        case ('mp2','mp3')
            call printSCFSolution

        case ('cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd','u-ccsd',&
              'r-ccsd(t)','cue-ccsdt','u-ccsdt','r-ccsdt')
            call getCCResults

        case ('fci')

        case default

    end select

    return
    end subroutine getResults

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getPolarizability
    implicit none

    integer(kind=iglu)             :: pp,pX,pY,pZ,i,k,state,nmethods,meth,pcurr,svio
    real   (kind=rglu)             :: rez
    character (len=:), allocatable :: cmethod,sendString


    void=tpSplit(generalbd%methods%get(),'+'); nmethods=tpSplitLen

    allocate (methodSet(nmethods))

    ! setting output stream for timecontrol
    svio=glGetIOunit(); call glSetIOunit(ou)

    do i = 1,nmethods
        methodSet(i)=tpSplitHold(i)
    enddo

    do meth = 1,nmethods
        cmethod=methodSet(meth)%get()

        pcurr=0; k=getMethodNumber(cmethod)
        do pp = 1,pointToCalc
            pX=pointSet(1,pp); pY=pointSet(2,pp); pZ=pointSet(3,pp)

            if (cmethod .in. ['mp2','mp3']) then
                void=timeControlCheckpoint('Polarizability: '//cmethod//&
                                           ' '//prStrByVal(pp)//'/'//&
                                           prStrByVal(pointToCalc),raw=true)
            endif

            call setCore(pX,pY,pZ); call updateHeader(pX,pY,pZ,statesbd%nStates)

            rez=getEnergy(cmethod,statesbd%nStates)
            do state = 0,statesbd%nStates
                gEnergyHolder(pX,pY,pZ,state,k)=getEnergy(cmethod,state)
            enddo
        enddo

        do state = 0,statesbd%nStates
            GlEt(:,:,:)=gEnergyHolder(:,:,:,state,k); call putPoint

            call showPolarizability(cmethod,state,ou,ouWidth)
        enddo
    enddo

    deallocate(methodSet)

    ! recover global output stream
    call glSetIOunit(svio)

    return

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(pX,pY,pZ)
        implicit none

        integer(kind=iglu), intent(in) :: pX,pY,pZ
        integer(kind=iglu)             :: mu


        mol%perturbation=0
        void=applyField()

        do mu = 1,mol%nAtoms
            mol%perturbation(mu,mu)=mol%perturbation(mu,mu)&
                                    +pX*polarizbd%derivStep*mol%atm(mu)%coords(1)&
                                    +pY*polarizbd%derivStep*mol%atm(mu)%coords(2)&
                                    +pZ*polarizbd%derivStep*mol%atm(mu)%coords(3)
        enddo
        call perturbate

        return
        end subroutine setCore

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine updateHeader(pX,pY,pZ,states)
        implicit none

        integer(kind=iglu), intent(in) :: pX,pY,pZ,states


        pcurr=pcurr+1
        call setIterationHeader(' Polarizability: '//trim(cmethod)//', States '//prStrByVal(states+1)//&
                                ', Field '//prStrByVal(pX)//', '//prStrByVal(pY)//', '//prStrByVal(pZ)//&
                                ' ['//prStrByVal(pcurr)//'/'//prStrByVal(pointToCalc)//'] ')

        return
        end subroutine updateHeader

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine getPolarizability

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getRDM

    implicit none

    integer(kind=iglu)              :: N,k,l,a,b,field,state,fatom,satom,derivPoints,sta,sto,pair(2)
    integer(kind=iglu)              :: pcount,pcurr,meth,nmethods,ncharges,norders,ln,swap,i,j
    real   (kind=rglu)              :: rez,derivStep,zpEnergy(0:statesbd%nStates)

    character (len=:), allocatable  :: cmethod

    real   (kind=rglu), allocatable :: D(:,:,:),derivEnergies(:,:)
    real   (kind=rglu), allocatable :: dmEVec(:,:),dmEVal(:)
    integer(kind=iglu), allocatable :: compGrid(:,:)
    integer(kind=iglu), allocatable :: orderSet(:,:),prCharSet(:),prOrderSet(:,:)


    ! parse method string
    void=tpSplit(generalbd%methods%get(),'+'); nmethods=tpSplitLen

    allocate (methodSet(nmethods))
    do meth = 1,nmethods
        methodSet(meth)=tpSplitHold(meth)
    enddo

    ! extracting parameters
    N=mol%nAtoms
    derivPoints=densitybd%nPoints
    derivStep=densitybd%derivStep; sta=-derivPoints/2; sto= derivPoints/2

    ! memory allocation
    void=glControlMemory(int( rglu*( N*N*(statesbd%nStates+1) + derivPoints*(statesbd%nStates+1) +&
                                     N*N+N )+iglu*N*N ,kind=i8kind),&
                         'RDM calculation')

    allocate (D(N,N,0:statesbd%nStates),compGrid(N,N))
    allocate (derivEnergies(sta:sto,0:statesbd%nStates))
    allocate (dmEVec(N,N),dmEVal(N))

    ! select RDM elements to be computed
    compGrid=0
    if (densitybd%dtype%get().EQ.'scharges') then

        ! parse selected charges
        void=tpSplit(densitybd%scharges%get(),';'); ncharges=tpSplitLen
        allocate (prCharSet(N)); prCharSet=0
        do k = 1,ncharges
            fatom=tpIntByStr(tpSplitHold(k)%get())
            if (atomEqu(fatom,fatom).LT.0) then
                satom=abs(atomEqu(fatom,fatom))
                compGrid(satom,satom)=1
            else
                compGrid(fatom,fatom)=1
            endif
        enddo

        ! prepare set of output charges
        do k = 1,N
            if (compGrid(k,k).EQ.1) then
                prCharSet(k)=1
                i=atomEqu(k,k)
                do j = 1,N
                    if (abs(atomEqu(j,j)).EQ.i) then
                        prCharSet(j)=1
                    endif
                enddo
            endif
        enddo

    elseif (densitybd%dtype%get().EQ.'sorders') then

        ! parese selected orders
        void=tpSplit(densitybd%sorders%get(),';'); norders=tpSplitLen
        allocate (orderSet(norders,2))
        allocate (prOrderSet(N,N)); prOrderSet=0
        do k = 1,norders
            pair=tpIntArrayByStr(tpSplitHold(k)%get(),2)
            i=pair(1); j=pair(2)
            if (i.GT.j) then
                swap=i; i=j; j=swap
            endif
            orderSet(k,1)=i
            orderSet(k,2)=j

            if (atomEqu(i,j).LT.0) then
                fatom=max(int(abs(atomEqu(i,j))/N),1)
                satom=mod(abs(atomEqu(i,j)),N)
                if (satom.EQ.0) then
                    satom=N
                    fatom=fatom-1
                endif
            else
                fatom=i; satom=j
            endif
            compGrid(fatom,satom)=1
        enddo

        ! prepare set of output orders
        do k = 1,N
        do l = k,N
            if (compGrid(k,l).EQ.1) then
                prOrderSet(k,l)=1
                if (k.EQ.l) then
                    a=k
                else
                    a=k*N+l
                endif
                do i = 1,N
                do j = i,N
                    if (abs(atomEqu(i,j)).EQ.a) then
                        prOrderSet(i,j)=1
                    endif
                enddo
                enddo
            endif
        enddo
        enddo
    else

        ! normal mode
        do fatom = 1,N
            if (densitybd%dtype%get().EQ.'charges') then
                if (atomEqu(fatom,fatom).GT.0) compGrid(fatom,fatom)=1
            elseif (densitybd%dtype%get().EQ.'all') then
                do satom = fatom,N
                    if (atomEqu(fatom,satom).GT.0) compGrid(fatom,satom)=1
                enddo
            else
                do satom = fatom,N
                    if (satom.EQ.fatom) cycle
                    if (atomEqu(fatom,satom).GT.0) compGrid(fatom,satom)=1
                enddo
            endif
        enddo
    endif

    ! count number of points to be computed
    pcount=0
    do fatom = 1,N
    do satom = 1,N
        if (compGrid(fatom,satom).EQ.1) then
            pcount=pcount+1
        endif
    enddo
    enddo
    pcount=pcount*(derivPoints-1)+1

    ! computaion of RDM1 elements
    do meth = 1,nmethods
        cmethod=methodSet(meth)%get()
        pcurr=0; state=0
        D=0; derivEnergies=0

        ! methods that require numerical RDM elements computation
        if (densitybd%forceNumerical.OR.(.NOT.(cmethod .in. ['huckel','rhf','fci']))) then
            call setCore(1,1,0); call updateHeader('zero point')
            rez=getEnergy(cmethod,statesbd%nStates)
            call getResults
            do state = 0,statesbd%nStates
                zpEnergy(state)=getEnergy(cmethod,state)
            enddo

            do fatom = 1,N
            do satom = 1,N
                if (compGrid(fatom,satom).EQ.0) cycle

                derivEnergies(0,:)=zpEnergy

                do field = sta,sto
                    if (field.EQ.0) cycle

                    call setCore(fatom,satom,field)
                    call updateHeader('Atoms: '//prStrByVal(fatom)//', '//prStrByVal(satom)//&
                                      ', Field:'//prStrByVal(field))

                    if (cmethod .in. ['cue-ccsd','u-ccsd','r-ccsd']) then
                        rez=getEnergy(cmethod,statesbd%nStates)
                        do state = 0, statesbd%nStates
                            if (state.NE.0) then
                                derivEnergies(field,state)=getEnergy(cmethod,0)+&
                                                           getEnergy(cmethod,state)
                            else
                                derivEnergies(field,state)=getEnergy(cmethod,0)
                            endif
                        enddo
                    else
                        derivEnergies(field,0)=getEnergy(cmethod)
                    endif
                enddo
                if (cmethod .in. ['cue-ccsd','u-ccsd','r-ccsd']) then
                    do state = 0,statesbd%nStates
                        call getDeriv(fatom,satom,state)
                    enddo
                else
                    call getDeriv(fatom,satom,0)
                endif
            enddo
            enddo
        else ! methods with analytical RDM computation
            call setCore(1,1,0)
            call updateHeader(cmethod,1,1)
            select case(cmethod)
                case('huckel')
                    rez=getEnergy(cmethod,state)
                    do fatom = 1,N
                        do satom = 1,N
                            D(fatom,satom,0)=getHuckRDMElement(fatom,satom)
                        enddo
                    enddo
                    call getResults
                case('rhf')
                    do state = 0,statesbd%nStates
                        rez=getEnergy(cmethod,state)
                        call getResults
                        call getSCFResult(density=D(:,:,state))
                    enddo
                case('fci')
                    rez=getEnergy(cmethod,statesbd%nStates)
                    call getResults
                    do state = 0,statesbd%nStates
                        call getFCIRDM(D(:,:,state),state)
                    enddo
            end select
        endif
        call printRDM
    enddo

    ! free allocated memory
    void=glControlMemory(int( sizeof(D)+sizeof(compGrid)+sizeof(derivEnergies)+sizeof(dmEVec)+&
                              sizeof(dmEVal) ,kind=i8kind),&
                        'RDM1 calculation', 'free')

    deallocate (D,compGrid,derivEnergies,dmEVec,dmEVal)
    if (densitybd%dtype%get().EQ.'scharges') then
        continue
    elseif (densitybd%dtype%get().EQ.'sorders') then
        deallocate(orderSet,prOrderSet)
    endif
    deallocate (methodSet)

    call setIterationHeader

    return

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine getDeriv(ip1,ip2,state)
        implicit none
        integer(kind=iglu), intent(in) :: ip1,ip2,state


        ! get derivatives
        D(ip1,ip2,state)=LagrangeDerivative(derivPoints,1,derivEnergies(:,state),derivStep)
        if (ip1.EQ.ip2) then
            D(ip2,ip1,state)=D(ip1,ip2,state)
        else
            D(ip1,ip2,state)=D(ip1,ip2,state)/2
            D(ip2,ip1,state)=D(ip1,ip2,state)
        endif

        return
        end subroutine getDeriv

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(ip1,ip2,fshift)
        implicit none
        integer(kind=iglu), intent(in) :: ip1,ip2,fshift


        ! set core perturbation
        mol%perturbation=0
        void=applyField()

        if (ip1.EQ.ip2) then
            mol%perturbation(ip1,ip2)=mol%perturbation(ip1,ip2)+fShift*derivStep
        else
            mol%perturbation(ip1,ip2)=mol%perturbation(ip1,ip2)+fShift*derivStep
            mol%perturbation(ip2,ip1)=mol%perturbation(ip2,ip1)+fShift*derivStep
        endif
        call perturbate

        return
        end subroutine setCore

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine printRDM
        implicit none
        integer(kind=iglu) :: c,b,k,l,paccur,ngroups,count
        real   (kind=rglu) :: swap,dx,dy,dz,dmod,dsum,csum,inds(3),shell_check
        type(uch)          :: chState

        integer(kind=iglu), allocatable :: array(:)


        paccur=densitybd%prntAccuracy

        write (ou,'(/A)') tpAdjustc(' RDM1 results for '//cmethod//' ', ouWidth, '#')

        ! put points by symmetry
        do state = 0,statesbd%nStates
            call putAtom(D(:,:,state))
        enddo
        do state = 0,statesbd%nStates
            if (state.EQ.0) then
                chState=', ground state'
            else
                chState=', excited state #'//prStrByVal(state)
            endif

            select case (densitybd%dtype%get())
                case ('charges') ! normal mode
                    write (ou,102) 'Atom',tpAdjustc('Density',paccur+3),tpAdjustc('Charge',paccur+4)

                    dx=0; dy=0; dz=0
                    do c = 1,N
                        write (ou,100) c,D(c,c,state),1-D(c,c,state)
                        dx=dx+mol%atm(c)%coords(1)*(1-D(c,c,state))
                        dy=dy+mol%atm(c)%coords(2)*(1-D(c,c,state))
                        dz=dz+mol%atm(c)%coords(3)*(1-D(c,c,state))
                    enddo
                    write (ou,103) 'X',dx*dipoleToDeby,'Y',dy*dipoleToDeby,'Z',dz*dipoleToDeby

                    ! if group set is defined
                    if (densitybd%gcharges%get().NE.'none') then
                        void=tpSplit(densitybd%gcharges%get(),';'); ngroups=tpSplitLen

                        do k = 1,ngroups

                            count=tpCount(tpSplitHold(k)%get(),',')+1
                            allocate (array(count)); array=tpIntArrayByStr(tpSplitHold(k)%get(),count)

                            dx=0; dy=0; dz=0; dsum=0; csum=0
                            write (ou,105) 'Sum for group '//tpSplitHold(k)%get()//':'
                            write (ou,102) 'Atom',tpAdjustc('Density',paccur+3),&
                                                  tpAdjustc('Charge',paccur+4)
                            do l = 1,count
                                c=array(l)

                                write (ou,100) c,D(c,c,state),1-D(c,c,state)

                                dsum=dsum+D(c,c,state)
                                csum=csum+(1-D(c,c,state))
                                dx=dx+mol%atm(c)%coords(1)*(1-D(c,c,state))
                                dy=dy+mol%atm(c)%coords(2)*(1-D(c,c,state))
                                dz=dz+mol%atm(c)%coords(3)*(1-D(c,c,state))

                            enddo
                            dmod=sqrt(dx**2+dy**2+dz**2)
                            write (ou,106) dsum,csum,'X',dx*dipoleToDeby,&
                                                     'Y',dy*dipoleToDeby,&
                                                     'Z',dz*dipoleToDeby,&
                                                   '|D|',dmod*dipoleToDeby
                            deallocate (array)
                        enddo
                    endif

                ! selected charges
                case ('scharges')
                    write (ou,102) 'Atom',tpAdjustc('Density',paccur+3),tpAdjustc('Charge',paccur+4)
                    do c = 1,N
                        if (prCharSet(c).EQ.1) then
                            write (ou,100) c,D(c,c,state),1-D(c,c,state)
                        endif
                    enddo

                ! selected orders
                case ('sorders')
                    write (ou,'(/2X,A,2X,A)') tpAdjustc('Atoms',10),&
                                              tpAdjustc('Bond order',max(3+paccur,10))
                    do c = 1,N-1
                        do b = c,N
                            if (prOrderSet(c,b).EQ.1) then
                                write (ou,101) c,b,D(c,b,state)
                            endif
                        enddo
                    enddo
                    write (ou,*)

                case default

                    ! if whole RDM1 is computed
                    if (densitybd%dtype%get().EQ.'all') then
                        write (ou,102) 'Atom',tpAdjustc('Density',paccur+3),&
                                              tpAdjustc('Charge' ,paccur+4)

                        dx=0; dy=0; dz=0
                        do c = 1,N
                            write (ou,100) c,D(c,c,state),1-D(c,c,state)
                            dx=dx+mol%atm(c)%coords(1)*(1-D(c,c,state))
                            dy=dy+mol%atm(c)%coords(2)*(1-D(c,c,state))
                            dz=dz+mol%atm(c)%coords(3)*(1-D(c,c,state))
                        enddo
                        write (ou,103) 'X',dx*dipoleToDeby,'Y',dy*dipoleToDeby,'Z',dz*dipoleToDeby
                    endif

                    ! if group set is defined
                    if (densitybd%gcharges%get().NE.'none') then
                        void=tpSplit(densitybd%gcharges%get(),';'); ngroups=tpSplitLen

                        do k = 1,ngroups

                            count=tpCount(tpSplitHold(k)%get(),',')+1
                            allocate (array(count)); array=tpIntArrayByStr(tpSplitHold(k)%get(),count)

                            dx=0; dy=0; dz=0; dsum=0; csum=0
                            write (ou,105) 'Sum for group '//tpSplitHold(k)%get()//':'
                            write (ou,102) 'Atom',tpAdjustc('Density',paccur+3),&
                                                  tpAdjustc('Charge' ,paccur+4)
                            do l = 1,count
                                c=array(l)

                                write (ou,100) c,D(c,c,state),1-D(c,c,state)

                                dsum=dsum+D(c,c,state)
                                csum=csum+(1-D(c,c,state))
                                dx=dx+mol%atm(c)%coords(1)*(1-D(c,c,state))
                                dy=dy+mol%atm(c)%coords(2)*(1-D(c,c,state))
                                dz=dz+mol%atm(c)%coords(3)*(1-D(c,c,state))

                            enddo
                            dmod=sqrt(dx**2+dy**2+dz**2)
                            write (ou,106) dsum,csum,'X',dx*dipoleToDeby,&
                                                     'Y',dy*dipoleToDeby,&
                                                     'Z',dz*dipoleToDeby,&
                                                   '|D|',dmod*dipoleToDeby
                            deallocate (array)
                        enddo
                    endif

                    ! output whole RDM1
                    call prMatrix(D(:,:,state),ou,&
                                  'RDM1('//cmethod//')'//chState%get(),&
                                  '^.'//tpFill(paccur,'0'),maxwidth=ouWidth)

            end select

            ! if natural orbitals required
            if (densitybd%NOAnalize) then
                call tred4(D(:,:,state),dmEVec,dmEVal,N,1.e-100_rglu,1.e-300_rglu)

                do c = 1,N/2
                    do k = 1,N
                        swap=dmEVec(k,c); dmEVec(k,c)=dmEVec(k,N-c+1); dmEVec(k,N-c+1)=swap
                    enddo
                    swap=dmEVal(c); dmEVal(c)=dmEVal(N-c+1); dmEVal(N-c+1)=swap
                enddo

                ! compute NO-occupancy correlation indices
                inds=0
                do c = 1,N/2
                    inds(1)=inds(1)+(2-dmEVal(c))*dmEVal(c)
                    inds(2)=inds(2)+min(2-dmEVal(c),dmEVal(c))
                enddo

                do c = N/2+1,N
                    inds(3)=inds(3)+dmEVal(c)
                enddo
                inds(3)=2*inds(3)

                shell_check=0
                do c = 1,N
                    shell_check=shell_check+dmEVal(c)
                enddo

                ! print eigenvalues and eigenvectors of RDM1
                call prEigenProblem(dmEVec,dmEVal,ou,cmethod//'. Natural orbitals'//&
                                    chState%get(),'^.'//tpFill(paccur,'0'),maxwidth=ouWidth)
                write (ou,104) inds
                write (ou,107) shell_check
            endif
        enddo

        write (ou,'(A/)') tpFill(ouWidth, '#')

100     format(2X,i4,2X,F<3+paccur>.<paccur>,2X,F<4+paccur>.<paccur>)
101     format(2X,i4,2X,i4,3X,F<3+paccur>.<paccur>)
102     format(/2X,A,2X,A,2X,A)
103     format(<17+2*paccur>('_')//3(<10+paccur>X,A,' =',F<4+paccur>.<paccur>/))
104     format(2X,'NO-occupancy-based indices:',3(1X,F<4+paccur>.<paccur>)/)
105     format(/2X,A)
106     format(<17+2*paccur>('_')/&
               8X,F<3+paccur>.<paccur>,2X,F<4+paccur>.<paccur>//&
               3(<10+paccur>X,A,' =',F<4+paccur>.<paccur>,1X,'D'/)&
                  <8+paccur>X,A,' =',F<4+paccur>.<paccur>,1X,'D')
107     format(2X,'Check the sum of occupancies:',1X,F10.5,1X,'electrons.'/)

        return
        end subroutine printRDM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine updateHeader(message,p1,p2)
        implicit none

        character (len=*)  :: message
        integer(kind=iglu), optional :: p1,p2
        integer(kind=iglu) :: pr1,pr2


        pcurr=pcurr+1
        pr1=pcurr; pr2=pcount

        if (present(p1)) pr1=p1
        if (present(p1)) pr2=p2

        ! header for iteration procedure
        call setIterationHeader(' RDM: '//trim(cmethod)//', '//message//&
                                ' ['//prStrByVal(pr1)//'/'//prStrByVal(pr2)//'] ')

        return
        end subroutine updateHeader

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine getRDM

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getHypercharges
    implicit none

    integer(kind=iglu)              :: N,nScales,naPoints,nfPoints,asta,asto,fsta,fsto,nmethods
    integer(kind=iglu)              :: i,j,k,l,c,d,sc,pcount,pcurr,aCount,svio,meth,isc
    real   (kind=rglu)              :: fStep,aStep,coeff,zpEnergy,rez,hc(8),hm(8)
    character(len=1)                :: scl,transform
    character (len=:), allocatable  :: cmethod

    real   (kind=rglu), allocatable :: fgrid(:,:,:),forZeroAlpha(:,:),hypercharges(:,:),alphas(:),&
                                       forZeroField(:,:)
    integer(kind=iglu), allocatable :: compGrid(:),prCharSet(:)

    real   (kind=rglu)              :: transition(0:8)
    real   (kind=rglu), parameter   :: derivThreshold(0:4)=[1E-5_rglu,1E-4_rglu,&
                                                            1E-3_rglu,1.E-2_rglu,1E-1_rglu]


    ! prepare transition constants to atom units
    transition(0)=1; transition(1)=-dipoleToDeby
    do i = 2,UBound(transition,1)
        transition(i)=-HartreeEnergy**(i-1)/BohrRadius**(i)
    enddo

    ! define bounds for derivative arrays
    N=mol%nAtoms
    naPoints=hyperchargesbd%naPoints; asta=-(naPoints-1)/2; asto=+(naPoints-1)/2; aStep=hyperchargesbd%derivaStep
    nfPoints=hyperchargesbd%nfPoints; fsta=-(nfPoints-1)/2; fsto=+(nfPoints-1)/2; fStep=hyperchargesbd%derivfStep

    ! define dimensionality
    nScales = hyperchargesbd%scales%ln

    ! set rdm calculation to charges only mode
    densitybd%dtype='charges'

    ! memory allocation
    void=glControlMemory(int( iglu*N ,kind=i8kind),'Hypercharges calculation')
    allocate (compGrid(N)); compGrid=0

    void=glControlMemory(int( rglu*(nfPoints*(3*N+3)+10*N+naPoints*(N+1) ) ,kind=i8kind),&
                         'Hypercharges calculation')
    allocate (forZeroAlpha(fsta:fsto,3),fgrid(N,fsta:fsto,3),hypercharges(N,0:9),alphas(asta:asto),&
              forZeroField(asta:asto,N))

    forZeroField=0; forZeroAlpha=0; hypercharges=0; fgrid=0

    !################ TODO: accout hypercharges (as for RDM1 elements)
    ! collect all unique atoms
    aCount=0
    do i = 1,N
        if (atomEqu(i,i).GT.0) then
            compGrid(i)=1
            aCount=aCount+1
        endif
    enddo

    ! setting output stream for timecontrol
    svio=glGetIOunit(); call glSetIOunit(ou)

    ! prepare methods
    void=tpSplit(generalbd%methods%get(),'+'); nmethods=tpSplitLen
    allocate (methodSet(nmethods))
    do i = 1,nmethods
        methodSet(i)=tpSplitHold(i)
    enddo

    do meth = 1,nmethods
        cmethod=methodSet(meth)%get()

        ! define number of points to be calculated
        if (cmethod .in. ['huckel', 'rhf', 'fci']) then
            pcount=1 + nScales*(nfPoints-1)
        else
            !pcount=nScales*aCount*nfPoints*naPoints
            pcount=1+nScales*(nfPoints-1)+aCount*(naPoints-1)+nScales*aCount*(nfPoints-1)*(naPoints-1)
        endif

        ! zero point calculation
        pcurr=0; call updateHeader('zero point energy')
        call setCore(0,0,0,0,0); zpEnergy=getEnergy(cmethod); call getResults

        ! methods with "simple" RDM calculation
        if (cmethod .in. ['huckel', 'rhf', 'fci']) then

            ! calculation of zero field charges
            do i = 1,N
                select case (cmethod)
                    case('huckel')
                        fgrid(i,0,1)=getHuckRDMElement(i,i)
                        fgrid(i,0,2)=getHuckRDMElement(i,i)
                        fgrid(i,0,3)=getHuckRDMElement(i,i)
                    case('rhf')
                        fgrid(i,0,1)=getSCFRDMElement(i,i)
                        fgrid(i,0,2)=getSCFRDMElement(i,i)
                        fgrid(i,0,3)=getSCFRDMElement(i,i)
                    case('fci')
                        fgrid(i,0,1)=getFCIRDMElement(i,i,0)
                        fgrid(i,0,2)=getFCIRDMElement(i,i,0)
                        fgrid(i,0,3)=getFCIRDMElement(i,i,0)
                end select
            enddo

            ! calculation of in-field charges
            do sc = 1,nScales
                scl=hyperchargesbd%scales%get(sc,sc)

                do k = fsta,fsto
                    if (k.EQ.0) cycle
                    if (cmethod.EQ.'rhf') then
                        void=timeControlCheckpoint('Hypercharges: scale '//scl//', Field '//prStrByVal(k),raw=true)
                    endif
                    call updateHeader('Scale '//scl//', Field '//prStrByVal(k))

                    ! setting core and calculation of in-field RDM
                    select case (scl)
                        case ('x'); call setCore(0,0,k,0,0); rez=getEnergy(cmethod); isc=1
                        case ('y'); call setCore(0,0,0,k,0); rez=getEnergy(cmethod); isc=2
                        case ('z'); call setCore(0,0,0,0,k); rez=getEnergy(cmethod); isc=3
                    end select

                    ! collecting in-filed charges
                    do i = 1,N
                        select case (cmethod)
                            case('huckel'); fgrid(i,k,isc)=getHuckRDMElement(i,i)
                            case('rhf')   ; fgrid(i,k,isc)=getSCFRDMElement(i,i)
                            case('fci')   ; fgrid(i,k,isc)=getFCIRDMElement(i,i,0)
                        end select
                    enddo
                enddo
            enddo

        ! methods that require differentiation for every RDM element
        else
            forZeroAlpha(0,:)=zpEnergy

            ! calculation of in-field energies without core perturbation
            do k = fsta,fsto
                if (k.EQ.0) cycle
                do sc = 1,nScales
                    scl=hyperchargesbd%scales%get(sc,sc)
                    if (cmethod .in. ['mp2', 'mp3']) then
                        void=timeControlCheckpoint('Hypercharges: zero alpha, scale '//scl//', field '//prStrByVal(k),raw=true)
                    endif

                    call updateHeader('zero alpha preparation')
                    select case (scl)
                        case('x'); call setCore(0,0,k,0,0); forZeroAlpha(k,1)=getEnergy(cmethod)
                        case('y'); call setCore(0,0,0,k,0); forZeroAlpha(k,2)=getEnergy(cmethod)
                        case('z'); call setCore(0,0,0,0,k); forZeroAlpha(k,3)=getEnergy(cmethod)
                    end select
                enddo
            enddo

            ! calculation of zero-field charges
            forZeroField(0,:)=zpEnergy
            do i = 1,N
                if (compGrid(i).EQ.0) cycle

                do l = asta,asto
                    if (l.EQ.0) cycle
                    call updateHeader('zero field preparation')
                    call setCore(i,l,0,0,0); forZeroField(l,i)=getEnergy(cmethod)
                enddo
            enddo

            ! put non-unique atoms
            do l = asta,asto
                call putAtomDiagonal(forZeroField(l,:),'n')
            enddo

            ! calculation of charges in zero-field
            do i = 1,N
                rez=LagrangeDerivative(naPoints,1,forZeroField(:,i),aStep)
                fgrid(i,0,1)=rez
                fgrid(i,0,2)=rez
                fgrid(i,0,3)=rez
            enddo

            ! calculation of in-field energies
            do sc = 1,nScales
                do i = 1,N
                    scl=hyperchargesbd%scales%get(sc,sc)
                    if (cmethod .in. ['mp2', 'mp3']) then
                        void=timeControlCheckpoint('Hypercharges: scale '//scl//' Atom '//prStrByVal(i),raw=true)
                    endif
                    if (compGrid(i).EQ.0) cycle

                    do k = fsta,fsto
                        if (k.EQ.0) cycle
                        alphas=0
                        do l = asta,asto
                            if (l.EQ.0) then
                                alphas(l)=forZeroAlpha(k,1)
                            else
                                call updateHeader('Scale '//scl//', Atom '//prStrByVal(i))
                                select case(scl)
                                    case('x'); call setCore(i,l,k,0,0); alphas(l)=getEnergy(cmethod); isc=1
                                    case('y'); call setCore(i,l,0,k,0); alphas(l)=getEnergy(cmethod); isc=2
                                    case('z'); call setCore(i,l,0,0,k); alphas(l)=getEnergy(cmethod); isc=3
                                end select
                            endif
                        enddo
                        fgrid(i,k,isc)=LagrangeDerivative(naPoints,1,alphas,aStep)
                    enddo
                enddo
            enddo
        endif

        hypercharges=0
        do i = 1,N
            if (.NOT.(cmethod .in. ['huckel', 'rhf', 'fci'])) then
                if (compGrid(i).EQ.0) cycle
            endif

            ! calculation of hypercharges
            do sc = 1,nScales
                select case (hyperchargesbd%scales%get(sc,sc))
                    case('x')
                        hypercharges(i,0)=fgrid(i,0,1)
                        hypercharges(i,1)=LagrangeDerivative(nfPoints,1,fgrid(i,:,1),fStep)
                        hypercharges(i,4)=LagrangeDerivative(nfPoints,2,fgrid(i,:,1),fStep)
                        hypercharges(i,7)=LagrangeDerivative(nfPoints,3,fgrid(i,:,1),fStep)
                    case('y')
                        hypercharges(i,0)=fgrid(i,0,2)
                        hypercharges(i,2)=LagrangeDerivative(nfPoints,1,fgrid(i,:,2),fStep)
                        hypercharges(i,5)=LagrangeDerivative(nfPoints,2,fgrid(i,:,2),fStep)
                        hypercharges(i,8)=LagrangeDerivative(nfPoints,3,fgrid(i,:,2),fStep)
                    case('z')
                        hypercharges(i,0)=fgrid(i,0,3)
                        hypercharges(i,3)=LagrangeDerivative(nfPoints,1,fgrid(i,:,3),fStep)
                        hypercharges(i,6)=LagrangeDerivative(nfPoints,2,fgrid(i,:,3),fStep)
                        hypercharges(i,9)=LagrangeDerivative(nfPoints,3,fgrid(i,:,3),fStep)
                end select
            enddo
        enddo

        if (.NOT.(cmethod .in. ['huckel', 'rhf', 'fci'])) then
            ! put hypercharges with account of symmetry
            do j = 0,9
                transform='n'
                select case (j)
                    case (1,7); transform='x'
                    case (2,8); transform='y'
                    case (3,9); transform='z'
                end select
                call putAtomDiagonal(hypercharges(:,j),transform)
            enddo
        endif

        write (ou,'(/A)') tpAdjustc(' Hypercharges values for '//cmethod//' ', 85, '*')

        ! output hypercharges
        do sc = 1,nScales
            scl=hyperchargesbd%scales%get(sc,sc)
            write (ou,100) '',           tpFill(1,scl),&
                           tpFill(1,scl),tpFill(2,scl),&
                           tpFill(2,scl),tpFill(3,scl),&
                           tpFill(3,scl),tpFill(4,scl)

            hc(5:8)=0; hm(5:8)=0
            do i = 1,N
                select case (scl)
                    case('x'); l=1
                    case('y'); l=2
                    case('z'); l=3
                end select
                hc(1)=hypercharges(i,0  ); hm(1)=getAtomMoment(i,0  ,scl,transition(1))
                hc(2)=hypercharges(i,l  ); hm(2)=getAtomMoment(i,l  ,scl,transition(2))
                hc(3)=hypercharges(i,l+3); hm(3)=getAtomMoment(i,l+3,scl,transition(3))
                hc(4)=hypercharges(i,l+6); hm(4)=getAtomMoment(i,l+6,scl,transition(4))

                void=purifyValues(derivThreshold(0),hc(1),hc(2),hc(3),hc(4))
                do k = 1,4
                    void=purifyValues(derivThreshold(k),hm(k))
                enddo

                write (ou,101) i, (hc(k), hm(k), k=1,4)

                hc(5)=hc(5)+hc(1); hm(5)=hm(5)+hm(1)
                hc(6)=hc(6)+hc(2); hm(6)=hm(6)+hm(2)
                hc(7)=hc(7)+hc(3); hm(7)=hm(7)+hm(3)
                hc(8)=hc(8)+hc(4); hm(8)=hm(8)+hm(4)
            enddo
            void=purifyValues(derivThreshold(0),hc(5),hc(6),hc(7),hc(8))
            do k = 5,8
                void=purifyValues(derivThreshold(k-4),hm(k))
            enddo
            write (ou,102) repeat('_',80),(hc(k), hm(k), k=5,8)
        enddo

        write (ou,'(/A/)') tpFill(85, '*')
    enddo

    call prLongText('Given values are calculated using Real-Space Finite-Field'//&
                    ' approach proposed by Bredas et. al.'&
                    ,ou,'center','justified',65,80,' ')

    call prLongText('** Geskin V.M., Bredas J.L. Evolution of the third-order molecular'      //&
                    ' polarizability in polyenes: A local view from atomic charge derivatives'//&
                    ' / J. Chem. Phys. - 1998. - Vol. 109, No. 14. - P. 6163.'&
                    ,ou,'center','justified',65,80,' ')

    call prLongText('** Geskin V.M., Lambert C., Bredas J.L. Origin of high second- and'       //&
                    ' third-order nonlinear optical response in ammonio/borato diphenylpolyene'//&
                    ' zwitterions: the remarkable role of polarized aromatic'                 //&
                    ' groups / J. Am. Chem. Soc. - 2003. - Vol. 125, No. 50. - P. 15651.'&
                    ,ou,'center','justified',65,80,' ')

    void=glControlMemory(int( sizeof(compGrid)+sizeof(forZeroAlpha)+sizeof(alphas)+&
                              sizeof(fgrid)+sizeof(hypercharges)+sizeof(forZeroField) ,kind=i8kind),&
                        'Hypercharges calculation', 'free')

    deallocate(compGrid,forZeroAlpha,alphas,fgrid,hypercharges,forZeroField)

    deallocate (methodSet)

    ! recover global output stream
    call glSetIOunit(svio)
    call setIterationHeader

100 format (/1X,'Atom',4X,'Q',A,9X,'M',A,8X,'Q',A,7X,'M',A,8X,'Q',A,6X,'M',A,7X,'Q',A,5X,'M',A)
101 format ( 1X,i3    ,1X,4(1X,F7.4,1X,ES11.4))
102 format ( 5X,A/5X,F8.4,1X,ES11.4,3(1X,F7.4,1X,ES11.4))

    return
    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(anumber,ashift,xshift,yshift,zshift)
        implicit none

        integer(kind=iglu), intent(in) :: anumber,ashift,xshift,yshift,zshift
        integer(kind=iglu)             :: atom


        mol%perturbation=0
        void=applyField()

        do atom = 1,N
            mol%perturbation(atom,atom)=mol%perturbation(atom,atom)&
                                        +xshift*fStep*mol%atm(atom)%coords(1)&
                                        +yshift*fStep*mol%atm(atom)%coords(2)&
                                        +zshift*fStep*mol%atm(atom)%coords(3)
        enddo
        if (anumber.NE.0) then
            mol%perturbation(anumber,anumber)=mol%perturbation(anumber,anumber)+ashift*aStep
        endif
        call perturbate

        return
        end subroutine setCore

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        real(kind=rglu) function getAtomMoment(atom,charge,coord,multiplier) result(ret)
        implicit none

        integer(kind=iglu) :: charge,atom
        character (len=*)  :: coord
        real   (kind=rglu) :: multiplier,coordinate


        select case (coord)
            case ('x'); ret=hypercharges(atom,charge)*mol%atm(atom)%coords(1)*multiplier
            case ('y'); ret=hypercharges(atom,charge)*mol%atm(atom)%coords(2)*multiplier
            case ('z'); ret=hypercharges(atom,charge)*mol%atm(atom)%coords(3)*multiplier
        end select
        return
        end function getAtomMoment

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine updateHeader(message)
        implicit none

        character (len=*)  :: message


        pcurr=pcurr+1
        call setIterationHeader(' hypercharges: '//trim(cmethod)//', '//message//&
                                ' ['//prStrByVal(pcurr)//'/'//prStrByVal(pcount)//'] ')

        return
        end subroutine updateHeader

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine getHypercharges

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getCoulson

    implicit none

    integer(kind=iglu) :: N,M,UBnd(2),felem,selem,pcount,pcurr,dpoints,dpower,paccur
    integer(kind=iglu) :: ii,jj,sta,sto,svio,meth,nmethods,ndiagonals,noffdiagonals
    integer(kind=iglu) :: i,j,k,l,a,b,c,pair(2),swap

    real(kind=rglu)    :: zpEnergy,statist(2)

    character (len=:) , allocatable :: cmethod

    integer(kind=iglu), allocatable :: compGrid(:,:),prDiagSet(:),prOffDiagSet(:,:),through(:,:)
    real(kind=rglu)   , allocatable :: d1Energies(:,:),d2Energies(:),coPolariz(:,:)


    N=mol%nAtoms; M=mol%nBonds
    select case (coulsonbd%ctype%get())
        case ('atom-atom'); UBnd(1)=N; UBnd(2)=N
        case ('atom-bond'); UBnd(1)=N; UBnd(2)=M
        case ('bond-bond'); UBnd(1)=M; UBnd(2)=M
    end select

    dpoints=coulsonbd%nPoints; sta=-dpoints/2; sto=dpoints/2
    dpower =coulsonbd%derivPower
    paccur =coulsonbd%prntAccuracy

    void=glControlMemory(int( rglu*(dpower**2+dpower) ,kind=i8kind),'Coulson calculation')
    allocate ( d1Energies(sta:sto,sta:sto),d2Energies(sta:sto) )

    void=glControlMemory(int( UBnd(1)*UBnd(2)*(rglu+iglu) ,kind=i8kind),'Coulson calculation')
    allocate (compGrid(UBnd(1),UBnd(2)),coPolariz(UBnd(1),UBnd(2)))
    allocate (through(UBnd(1),UBnd(2))); through=0

    ! select Coulson Polarizability matrix elements to be computed
    compGrid=0
    if (coulsonbd%selected%get().EQ.'diagonals') then

        ! parse selected diagonals
        void=tpSplit(coulsonbd%sdiagonals%get(),';'); ndiagonals=tpSplitLen
        allocate(prDiagSet(UBnd(1))); prDiagSet=0
        do k = 1,ndiagonals
            felem=tpIntByStr(tpSplitHold(k)%get())
            through(felem,felem)=1

            select case (coulsonbd%ctype%get())
                case('atom-atom')
                    if (atomEqu(felem,felem).LT.0) then
                        selem=abs(atomEqu(felem,felem))
                        compGrid(selem,selem)=1
                    else
                        compGrid(felem,felem)=1
                    endif

                case('bond-bond')
                    if (bondEqu(felem,felem).LT.0) then
                        selem=abs(bondEqu(felem,felem))
                        compGrid(selem,selem)=1
                    else
                        compGrid(felem,felem)=1
                    endif

            end select
        enddo

        ! prepare set of output diagonals
        do k = 1,UBnd(1)
            if (compGrid(k,k).EQ.1) then
                prDiagSet(k)=1

                select case (coulsonbd%ctype%get())
                    case('atom-atom')
                        i=atomEqu(k,k)
                        do j = 1,UBnd(1)
                            if (abs(atomEqu(j,j)).EQ.i) then
                                prDiagSet(j)=1
                            endif
                        enddo

                    case('bond-bond')
                        i=bondEqu(k,k)
                        do j = 1,UBnd(1)
                            if (abs(bondEqu(j,j)).EQ.i) then
                                prDiagSet(j)=1
                            endif
                        enddo

                end select
            endif
        enddo

    elseif (coulsonbd%selected%get().EQ.'offdiagonals') then

        ! parse selected offdiagonals
        void=tpSplit(coulsonbd%soffdiagonals%get(),';'); noffdiagonals=tpSplitLen
        allocate (prOffDiagSet(UBnd(1),UBnd(1))); prOffDiagSet=0
        do k = 1,noffdiagonals
            pair=tpIntArrayByStr(tpSplitHold(k)%get(),2)
            i=pair(1); j=pair(2)
            if (i.GT.j) then
                swap=i; i=j; j=swap
            endif
            through(i,j)=1

            select case (coulsonbd%ctype%get())
                case('atom-atom')
                    if (atomEqu(i,j).LT.0) then
                        felem=max(int(abs(atomEqu(i,j))/UBnd(1)), 1)
                        selem=mod(abs(atomEqu(i,j)),UBnd(1))
                        if (selem.EQ.0) then
                            selem=UBnd(1)
                            felem=felem-1
                        endif
                    else
                        felem=i; selem=j
                    endif

                case('bond-bond')
                    if (bondEqu(i,j).LT.0) then
                        felem=max(int(abs(bondEqu(i,j))/UBnd(1)),1)
                        selem=mod(abs(bondEqu(i,j)),UBnd(1))
                        if (selem.EQ.0) then
                            selem=UBnd(1)
                            felem=felem-1
                        endif
                    else
                        felem=i; selem=j
                    endif

                case ('atom-bond')
                    felem=i; selem=j
            end select

            compGrid(felem,selem)=1
        enddo

        ! ==========> TEST ON ATOM-BOND TYPE <==========
        ! prepare set of output offdiagonals
        do k = 1,UBnd(1)
        do l = k,UBnd(1)
            if (compGrid(k,l).EQ.1) then
                prOffDiagSet(k,l)=1
                if (k.EQ.l) then
                    a=k
                else
                    a=k*UBnd(1)+l
                endif
                do i = 1,UBnd(1)
                do j = i,UBnd(1)

                    select case (coulsonbd%ctype%get())
                        case('atom-atom')
                            if (abs(atomEqu(i,j)).EQ.a) prOffDiagSet(i,j)=1

                        case('bond-bond')
                            if (abs(bondEqu(i,j)).EQ.a) prOffDiagSet(i,j)=1
                    end select

                enddo
                enddo
            endif
        enddo
        enddo
    else

        ! normal mode
        compGrid=0
        do felem = 1,UBnd(1)
            do selem = felem,UBnd(2)
                select case (coulsonbd%ctype%get())
                    case ('atom-atom')
                        if (atomEqu(felem,selem).GT.0) compGrid(felem,selem)=1

                    case ('atom-bond')
                        compGrid(felem,selem)=1
                        compGrid(selem,felem)=1

                    case ('bond-bond')
                        if (bondEqu(felem,selem).GT.0) compGrid(felem,selem)=1
                end select
            enddo
        enddo
    endif

    ! count number of points to be computed
    pcount=0
    do felem = 1,UBnd(1)
    do selem = 1,UBnd(2)
        if (compGrid(felem,selem).EQ.1) then
            pcount=pcount+1
        endif
    enddo
    enddo
    pcount=(dpoints**2-1)*pcount+1

    ! setting output stream for timecontrol
    svio=glGetIOunit(); call glSetIOunit(ou)

    ! prepare methods
    void=tpSplit(generalbd%methods%get(),'+'); nmethods=tpSplitLen
    allocate (methodSet(nmethods))
    do ii = 1,nmethods
        methodSet(ii)=tpSplitHold(ii)
    enddo

    ! call singleSession('  '//prJoin([(prStrByVal(fieldbd%strength(k), '____.000'), k=1,3)], '  '))
    call singleSession(systembd%throughPrefix%get())

    do meth = 1,nmethods
        cmethod=methodSet(meth)%get()

        pcurr=0; call updateHeader('zero point')
        call setCore(0,0,0,0); zpEnergy=getEnergy(cmethod)

        coPolariz=0
        do felem = 1,UBnd(1)
        do selem = 1,UBnd(2)
            if (compGrid(felem,selem).EQ.0) cycle

            do ii = sta,sto
            do jj = sta,sto
                if ((ii.EQ.0).AND.(jj.EQ.0)) then
                    d1Energies(ii,jj)=zpEnergy
                    cycle
                endif

                if (cmethod .in. ['mp2', 'mp3']) then
                    void=timeControlCheckpoint('Coulson: Elements '//prStrByVal(felem)//' '//prStrByVal(selem)//&
                                               ', Grid '//prStrByVal(ii)//' '//prStrByVal(jj),raw=true)
                endif

                call updateHeader('Elements: '//prStrByVal(felem)//' '//prStrByVal(selem)//&
                                  ', Grid:'//prStrByVal(ii)//' '//prStrByVal(jj))

                call setCore(felem,selem,ii,jj); d1Energies(ii,jj)=getEnergy(cmethod)
            enddo
            enddo
            call getDeriv(felem,selem)
        enddo
        enddo

        select case (coulsonbd%ctype%get())
            case ('atom-atom'); call putAtom(coPolariz)
            case ('bond-bond'); call putBond(coPolariz); coPolariz=coPolariz/2
        end select

        if (coulsonbd%selected%get().NE.'none') then

            select case(coulsonbd%selected%get())
                ! selected diagonals
                case ('diagonals')
                    select case (coulsonbd%ctype%get())
                        case('atom-atom'); write (ou,102) 'Atom',tpAdjustc('Polariz',paccur+3)
                        case('bond-bond'); write (ou,102) 'Bond',tpAdjustc('Polariz',paccur+3)
                    end select

                    do c = 1,UBnd(1)
                        if (prDiagSet(c).EQ.1) then
                            write (ou,100) c,coPolariz(c,c)
                        endif
                    enddo

                ! selected offdiagonals
                case ('offdiagonals')
                    select case (coulsonbd%ctype%get())
                        case('atom-atom'); write (ou,102) tpAdjustc('Atoms',10),' '//tpAdjustc('Polariz',paccur+3)
                        case('bond-bond'); write (ou,102) tpAdjustc('Bonds',10),' '//tpAdjustc('Polariz',paccur+3)
                    end select

                    do c = 1,UBnd(1)
                        do b = c,UBnd(1)
                            if (prOffDiagSet(c,b).EQ.1) then
                                write (ou,101) c,b,coPolariz(c,b)
                            endif
                        enddo
                    enddo
                    write (ou,*)
            end select
            do c = 1,UBnd(1)
                do b = 1,UBnd(2)
                    if (through(c,b).EQ.1) then
                        call singleSession('  '//prStrByVal(coPolariz(c,b), '____.00000'))
                    endif
                enddo
            enddo
        else
            call prMatrix(coPolariz,ou,'Coulson '//coulsonbd%ctype%get()//' polarizabilities ('//cmethod//')',&
                          '^.'//tpFill(paccur, '0'),maxwidth=ouWidth)
            if (coulsonbd%ctype%get().EQ.'atom-atom') then
                statist=Cstatistics(coPolariz)
                write(ou,'(3X,2(1X,A,2X,F<3+paccur>.<paccur>))') 'Rindex',statist(1),'Dispersion',statist(2)
            endif
        endif
    enddo

    void=glControlMemory(int( sizeof(d1Energies)+sizeof(d2Energies)+&
                              sizeof(compGrid)  +sizeof(coPolariz)   ,kind=i8kind),&
                        'Coulson calculation', 'free')

    deallocate ( d1Energies,d2Energies,compGrid,coPolariz,through )

    ! recover global output stream
    call glSetIOunit(svio)
    call setIterationHeader

100 format(2X,i4,2X,F<3+paccur>.<paccur>,2X,F<4+paccur>.<paccur>)
101 format(2X,i4,2X,i4,3X,F<3+paccur>.<paccur>)
102 format(/2X,A,2X,A,2X,A)

    return

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(fenumber,senumber,fShift,sShift)
        implicit none

        integer(kind=iglu) :: fenumber,senumber,fShift,sShift
        integer(kind=iglu) :: mu,nu,rh,si


        mol%perturbation=0

        void=applyField()

        if ((fenumber.EQ.0).AND.(senumber.EQ.0)) then
            call perturbate
            return
        endif

        select case(coulsonbd%ctype%get())
            case('atom-atom')
                mu=fenumber; nu=senumber
                mol%perturbation(mu,mu)=mol%perturbation(mu,mu)+fShift*coulsonbd%derivStep
                mol%perturbation(nu,nu)=mol%perturbation(nu,nu)+sShift*coulsonbd%derivStep

            case('atom-bond')
                mu=fenumber
                rh=mol%bnd(senumber)%atoms(1)
                si=mol%bnd(senumber)%atoms(2)

                mol%perturbation(mu,mu)=mol%perturbation(mu,mu)+fShift*coulsonbd%derivStep

                mol%perturbation(rh,si)=mol%perturbation(rh,si)+sShift*coulsonbd%derivStep
                mol%perturbation(si,rh)=mol%perturbation(si,rh)+sShift*coulsonbd%derivStep

            case('bond-bond')
                mu=mol%bnd(fenumber)%atoms(1)
                nu=mol%bnd(fenumber)%atoms(2)
                rh=mol%bnd(senumber)%atoms(1)
                si=mol%bnd(senumber)%atoms(2)

                mol%perturbation(mu,nu)=mol%perturbation(mu,nu)+fShift*coulsonbd%derivStep
                mol%perturbation(nu,mu)=mol%perturbation(nu,mu)+fShift*coulsonbd%derivStep

                mol%perturbation(rh,si)=mol%perturbation(rh,si)+sShift*coulsonbd%derivStep
                mol%perturbation(si,rh)=mol%perturbation(si,rh)+sShift*coulsonbd%derivStep

        end select
        call perturbate

        return
        end subroutine setCore

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine getDeriv(felem,selem)
        implicit none

        integer(kind=iglu), intent(in) :: felem,selem
        integer(kind=iglu)             :: ii


        do ii = sta,sto
            d2Energies(ii)=LagrangeDerivative(dpoints,dpower/2,d1Energies(ii,:),coulsonbd%derivStep)
        enddo
        coPolariz(felem,selem)=LagrangeDerivative(dpoints,dpower/2,d2Energies,coulsonbd%derivStep)

        if (coulsonbd%ctype%get() .in. ['atom-atom', 'bond-bond']) then
            coPolariz(selem,felem)=coPolariz(felem,selem)
        endif

        return
        end subroutine getDeriv

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        function Cstatistics(array) result(ret) ! atom-atom only
        implicit none

        real(kind=rglu)    :: ret(2)
        real(kind=rglu)    :: array(N,N)
        real(kind=rglu)    :: Ax,Bx,Cx
        integer(kind=iglu) :: i,j,a


        Ax=0; Bx=0
        do i = 1,N-1
        do j = i+1,N
            Ax=Ax+mol%atmdist(i,j)*array(i,j)**2
            Bx=Bx+array(i,j)**2
        enddo
        enddo
        ret(1)=Ax/Bx

        a=0; Cx=0
        do i = 1,N-1
        do j = i+1,N
            a=a+1
            Cx=Cx+(mol%atmdist(i,j)*array(i,j)**2/Bx-ret(1))**2
        enddo
        enddo
        ret(2)=Cx/a

        return
        end function Cstatistics

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine updateHeader(message)
        implicit none

        character (len=*)  :: message


        pcurr=pcurr+1
        call setIterationHeader(' coulson: '//trim(cmethod)//', '//message//&
                                ' ['//prStrByVal(pcurr)//'/'//prStrByVal(pcount)//'] ')

        return
        end subroutine updateHeader

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine getCoulson

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine getWaveFunctionAnalize

    implicit none

    integer(kind=iglu)      :: i,k,state,methods,nmethods
    real   (kind=rglu)      :: rez


    void=tpSplit(generalbd%methods%get(),'+'); nmethods=tpSplitLen

    allocate (methodSet(nmethods))

    do i = 1,nmethods
        methodSet(i)=tpSplitHold(i)
    enddo

    do i = 1,UBound(methodSet,1)
        if (.NOT.(methodSet(i)%get() .in. ['cue-ccs','r-ccd','u-ccd','cue-ccsd','r-ccsd',&
                                           'u-ccsd','cue-ccsdt','u-ccsdt','r-ccsdt'])) cycle

        call setIterationHeader(' Wave-function analysis: '//methodSet(i)%get()//' ')
        rez=getEnergy(methodSet(i)%get(),statesbd%nStates)
        call analizewfCC

        if (statesbd%nStates.GT.0) then
            do state = 1,statesbd%nStates
                call analizewfLR(state)
            enddo
        endif
    enddo

    deallocate(methodSet)

    return
    end subroutine getWaveFunctionAnalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function putMethodEnergy(energy, method, state, px, py, pz) result(ret)
    implicit none

    real(kind=rglu),  intent(in) :: energy
    character(len=*), intent(in) :: method
    integer(kind=iglu), optional :: state,px,py,pz
    integer(kind=iglu)           :: istate,ipx,ipy,ipz,imethod


    ret=0

    istate=0; if (present(state)) istate=state
    ipx=0; if (present(px)) ipx=px
    ipy=0; if (present(py)) ipy=py
    ipz=0; if (present(pz)) ipz=pz


    imethod=getMethodNumber(method)

    if (imethod.LE.0) then
        ret=-1
        return
    endif

    if (noPerturbation) then
        gEnergyHolder(ipx,ipy,ipz,istate,imethod)=energy
    endif

    return
    end function putMethodEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function getMethodEnergy(method, state, px, py, pz, calc) result(ret)
    implicit none

    character(len=*), intent(in) :: method
    integer(kind=iglu), optional :: state,px,py,pz
    logical(kind=lglu), optional :: calc
    integer(kind=iglu)           :: istate,ipx,ipy,ipz,imethod,mu
    logical(kind=lglu)           :: icalc
    real(kind=rglu)              :: energy,field


    ret=0

    istate=0; if (present(state)) istate=state
    ipx=0; if (present(px)) ipx=px
    ipy=0; if (present(py)) ipy=py
    ipz=0; if (present(pz)) ipz=pz
    icalc=false; if (present(calc)) icalc=calc


    imethod=getMethodNumber(method)

    if (imethod.LE.0) return

    energy=gEnergyHolder(ipx,ipy,ipz,istate,imethod)

    if (abs(energy).GT.gluCompare) then
        ret=energy
        return
    endif

    if (icalc) then
        mol%perturbation=0
        field=polarizbd%derivStep
        do mu = 1,mol%nAtoms
            mol%perturbation(mu,mu)=+ipX*field*mol%atm(mu)%coords(1)&
                                    +ipY*field*mol%atm(mu)%coords(2)&
                                    +ipZ*field*mol%atm(mu)%coords(3)
        enddo
        call perturbate
        gEnergyHolder(ipx,ipy,ipz,istate,imethod)=getEnergy(method,istate)
        energy=gEnergyHolder(ipx,ipy,ipz,istate,imethod)
    endif

    return
    end function getMethodEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function applyField() result(ret)
    implicit none

    integer(kind=iglu) :: k,mu


    ret=0

    do k = 1,3
        if (abs(fieldbd%strength(k)).GT.1D-10) then
            do mu = 1,mol%nAtoms
                mol%perturbation(mu,mu)=fieldbd%strength(k)*mol%atm(mu)%coords(k)
            enddo
        endif
    enddo

    return
    end function applyField

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module property