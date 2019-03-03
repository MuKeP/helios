    subroutine hyperdensity

    use glob,           only: assignment(=),glGetIOunit,glSetIOunit,purifyValues,glControlMemory
    use glob,           only: rglu,iglu,lglu,true,false,uch,void,gluCompare,timeControlCheckpoint
    use glob,           only: i8kind
    use txtParser,      only: tpSplit,tpSplitLen,tpSplitHold,tpFill,tpAdjustc,operator(.in.)
    use txtParser,      only: tpIntByStr,tpIntArrayByStr,tpCount
    use printmod,       only: prStrByVal,prMatrix,prEigenProblem,prLongText,prJoin
    use math,           only: LagrangeDerivative,tred4
    use hdb,            only: mol,scfbd,ccbd,statesbd,generalbd,polarizbd,ou,ouWidth,gEnergyHolder
    use hdb,            only: pointToCalc,pointSet,GlEt,perturbate,perturbationID,getMethodNumber
    use hdb,            only: noPerturbation,iheader,densitybd,hyperchargesbd,atomEqu,bondEqu
    use hdb,            only: HartreeEnergy,BohrRadius,dipoleToDeby,coulsonbd
    use huckel,         only: getHuckelResult,getHuckRDMElement,setHuckelParameters,energyHuckel
    use huckel,         only: finalizeHuck
    use scf,            only: setSCFParameters,initSCF,iterationSCF,energySCF,getSCFRDMElement
    use scf,            only: getSCFResult,finalizeSCF,printSCFSolution
    use mbpt,           only: setMBPTParameters,initMBPT,energyMBPT,finalizeMBPT
    use coupledCluster, only: setCCParameters,initCC,iterationCC,energyCC,finalizeCC,analizewfCC
    use coupledCluster, only: getCCResults
    use fci,            only: setFCIParameters,initFCI,energyFCI,finalizeFCI,fciHoldStateEnergy
    use fci,            only: getFCIRDM,getFCIRDMElement
    use excitedStates,  only: getExcitationEnergy
    use lrccsdModule,   only: setLRParameters,initLR,energyLR,finalizeLR
    use lrccsdModule,   only: lrHoldStateEnergy,analizewfLR
    use property,       only: getEnergy

    implicit none

    integer(kind=iglu)              :: N,nScales,naPoints,nfPoints,asta,asto,fsta,fsto,nmethods
    integer(kind=iglu)              :: i,j,k,l,a,b,c,sc,pcount,pcurr,aCount,svio,meth,isc
    real   (kind=rglu)              :: fStep,aStep,coeff,zpEnergy,rez,swap
    character(len=1)                :: scl,transform
    character (len=:), allocatable  :: cmethod,cstr

    real   (kind=rglu), allocatable :: fgrid(:,:,:,:,:),D(:,:,:),alphas(:),interm(:)
    real   (kind=rglu), allocatable :: dmEVal(:),dmEVec(:,:)

    real   (kind=rglu)              :: transition(0:8)
    real   (kind=rglu), parameter   :: derivThreshold(0:4)=[1E-5_rglu,1E-4_rglu,&
                                                            1E-3_rglu,1.E-2_rglu,1E-1_rglu]


    svio=glGetIOunit(); call glSetIOunit(ou)

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

    allocate (fgrid(N,N,fsta:fsto,fsta:fsto,fsta:fsto),D(N,N,0:25),alphas(asta:asto),interm(fsta:fsto))

    allocate (dmEVal(N),dmEVec(N,N))

    ! prepare methods
    cmethod='cue-ccsd'

    fgrid=0
    pcurr=0
    do i = 1,N
    do j = i,N
        do a = fsta,fsto
        do b = fsta,fsto
        do c = fsta,fsto
            if ((a.NE.0).AND.(b.NE.0).AND.(c.NE.0)) cycle
            alphas=0
            do l = asta,asto
                pcurr=pcurr+1
                cstr=prJoin(['Hypercharges',&
                             prStrByVal(i),&
                             prStrByVal(j),&
                             prStrByVal(a),&
                             prStrByVal(b),&
                             prStrByVal(c),&
                             prStrByVal(l),&
                             prStrByVal(pcurr)])

                void=timeControlCheckpoint(cstr,raw=true)

                call setCore(i,j,l,a,b,c); alphas(l)=getEnergy(cmethod)
            enddo
            fgrid(i,j,a,b,c)=LagrangeDerivative(naPoints,1,alphas,aStep)
        enddo
        enddo
        enddo
    enddo
    enddo

    D=0
    D(:,:,0)=fgrid(:,:,0,0,0)
    do i = 1,N
    do j = 1,N
        D(i,j, 1)=LagrangeDerivative(nfPoints,1,fgrid(i,j,:,0,0),fStep) ! x
        D(i,j, 2)=LagrangeDerivative(nfPoints,1,fgrid(i,j,0,:,0),fStep) ! y
        D(i,j, 3)=LagrangeDerivative(nfPoints,1,fgrid(i,j,0,0,:),fStep) ! z

        D(i,j, 4)=LagrangeDerivative(nfPoints,2,fgrid(i,j,:,0,0),fStep) ! xx
        D(i,j, 5)=LagrangeDerivative(nfPoints,2,fgrid(i,j,0,:,0),fStep) ! yy
        D(i,j, 6)=LagrangeDerivative(nfPoints,2,fgrid(i,j,0,0,:),fStep) ! zz

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,1,fgrid(i,j,:,l,0),fStep)
        enddo
        D(i,j, 7)=LagrangeDerivative(nfPoints,1,interm,fStep) ! xy

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,1,fgrid(i,j,l,:,0),fStep) ! useless
        enddo
        D(i,j, 8)=LagrangeDerivative(nfPoints,1,interm,fStep) ! yx

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,1,fgrid(i,j,:,0,l),fStep)
        enddo
        D(i,j, 9)=LagrangeDerivative(nfPoints,1,interm,fStep) ! xz

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,1,fgrid(i,j,l,0,:),fStep) ! useless
        enddo
        D(i,j,10)=LagrangeDerivative(nfPoints,1,interm,fStep) ! zx

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,1,fgrid(i,j,0,:,l),fStep)
        enddo
        D(i,j,11)=LagrangeDerivative(nfPoints,1,interm,fStep) ! yz

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,1,fgrid(i,j,0,l,:),fStep) ! useless
        enddo
        D(i,j,12)=LagrangeDerivative(nfPoints,1,interm,fStep) ! zy

        D(i,j,13)=LagrangeDerivative(nfPoints,3,fgrid(i,j,:,0,0),fStep) ! xxx
        D(i,j,14)=LagrangeDerivative(nfPoints,3,fgrid(i,j,0,:,0),fStep) ! yyy
        D(i,j,15)=LagrangeDerivative(nfPoints,3,fgrid(i,j,0,0,:),fStep) ! zzz

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,2,fgrid(i,j,:,l,0),fStep)
        enddo
        D(i,j,16)=LagrangeDerivative(nfPoints,1,interm,fStep) ! xxy

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,2,fgrid(i,j,l,:,0),fStep)
        enddo
        D(i,j,17)=LagrangeDerivative(nfPoints,1,interm,fStep) ! yyx

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,2,fgrid(i,j,:,0,l),fStep)
        enddo
        D(i,j,18)=LagrangeDerivative(nfPoints,1,interm,fStep) ! xxz

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,2,fgrid(i,j,l,0,:),fStep)
        enddo
        D(i,j,19)=LagrangeDerivative(nfPoints,1,interm,fStep) ! zzx

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,2,fgrid(i,j,0,:,l),fStep)
        enddo
        D(i,j,20)=LagrangeDerivative(nfPoints,1,interm,fStep) ! yyz

        do l = fsta,fsto
            interm(l)=LagrangeDerivative(nfPoints,2,fgrid(i,j,0,l,:),fStep)
        enddo
        D(i,j,21)=LagrangeDerivative(nfPoints,1,interm,fStep) ! zzy
    enddo
    enddo

    do k = 0,21
        do i = 1,N-1
        do j = i+1,N
            D(i,j,k)=D(i,j,k)/2
            D(j,i,k)=D(i,j,k)
        enddo
        enddo
    enddo

    write (ou,'(/A)') tpAdjustc(' Hypercharges values for '//cmethod//' ', 85, '*')

    do k = 0,21
        call prMatrix(D(:,:,k),ou,'HDensity '//prStrByVal(k),'^.00000',maxwidth=ouWidth)

        call tred4(D(:,:,k),dmEVec,dmEVal,N,1.e-100_rglu,1.e-300_rglu)

        do c = 1,N/2
            do b = 1,N
                swap=dmEVec(b,c); dmEVec(b,c)=dmEVec(b,N-c+1); dmEVec(b,N-c+1)=swap
            enddo
            swap=dmEVal(c); dmEVal(c)=dmEVal(N-c+1); dmEVal(N-c+1)=swap
        enddo

        ! print eigenvalues and eigenvectors of RDM1
        call prEigenProblem(dmEVec,dmEVal,ou,cmethod//'. Natural orbitals HDensity '//prStrByVal(k),'^.00000',maxwidth=ouWidth)

    enddo

!     ! output hypercharges
!     do sc = 1,nScales
!         scl=hyperchargesbd%scales%get(sc,sc)
!         write (ou,100) '',           tpFill(1,scl),&
!                        tpFill(1,scl),tpFill(2,scl),&
!                        tpFill(2,scl),tpFill(3,scl),&
!                        tpFill(3,scl),tpFill(4,scl)

!         hc(5:8)=0; hm(5:8)=0
!         do i = 1,N
!             select case (scl)
!                 case('x'); l=1
!                 case('y'); l=2
!                 case('z'); l=3
!             end select
!             hc(1)=hypercharges(i,0  ); hm(1)=getAtomMoment(i,0  ,scl,transition(1))
!             hc(2)=hypercharges(i,l  ); hm(2)=getAtomMoment(i,l  ,scl,transition(2))
!             hc(3)=hypercharges(i,l+3); hm(3)=getAtomMoment(i,l+3,scl,transition(3))
!             hc(4)=hypercharges(i,l+6); hm(4)=getAtomMoment(i,l+6,scl,transition(4))

!             void=purifyValues(derivThreshold(0),hc(1),hc(2),hc(3),hc(4))
!             do k = 1,4
!                 void=purifyValues(derivThreshold(k),hm(k))
!             enddo

!             write (ou,101) i, (hc(k), hm(k), k=1,4)

!             hc(5)=hc(5)+hc(1); hm(5)=hm(5)+hm(1)
!             hc(6)=hc(6)+hc(2); hm(6)=hm(6)+hm(2)
!             hc(7)=hc(7)+hc(3); hm(7)=hm(7)+hm(3)
!             hc(8)=hc(8)+hc(4); hm(8)=hm(8)+hm(4)
!         enddo
!         void=purifyValues(derivThreshold(0),hc(5),hc(6),hc(7),hc(8))
!         do k = 5,8
!             void=purifyValues(derivThreshold(k-4),hm(k))
!         enddo
!         write (ou,102) repeat('_',80),(hc(k), hm(k), k=5,8)
!     enddo

!     write (ou,'(/A/)') tpFill(85, '*')

!     ! recover global output stream
!     call glSetIOunit(svio)

!     iheader=''

! 100 format (/1X,'Atom',4X,'Q',A,9X,'M',A,8X,'Q',A,7X,'M',A,8X,'Q',A,6X,'M',A,7X,'Q',A,5X,'M',A)
! 101 format ( 1X,i3    ,1X,4(1X,F7.4,1X,ES11.4))
! 102 format ( 5X,A/5X,F8.4,1X,ES11.4,3(1X,F7.4,1X,ES11.4))

    stop
    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(fnumber,snumber,ashift,xshift,yshift,zshift)
        implicit none

        integer(kind=iglu), intent(in) :: fnumber,snumber,ashift,xshift,yshift,zshift
        integer(kind=iglu)             :: atom


        mol%perturbation=0
        do atom = 1,N
            mol%perturbation(atom,atom)=+xshift*fStep*mol%atm(atom)%coords(1)&
                                        +yshift*fStep*mol%atm(atom)%coords(2)&
                                        +zshift*fStep*mol%atm(atom)%coords(3)
        enddo
        if (fnumber.NE.0) then
            if (fnumber.EQ.snumber) then
                mol%perturbation(fnumber,snumber)=mol%perturbation(fnumber,snumber)+ashift*aStep
            else
                mol%perturbation(fnumber,snumber)=mol%perturbation(fnumber,snumber)+ashift*aStep
                mol%perturbation(snumber,fnumber)=mol%perturbation(snumber,fnumber)+ashift*aStep
            endif
        endif
        call perturbate

        return
        end subroutine setCore

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        ! real(kind=rglu) function getAtomMoment(atom,charge,coord,multiplier) result(ret)
        ! implicit none

        ! integer(kind=iglu) :: charge,atom
        ! character (len=*)  :: coord
        ! real   (kind=rglu) :: multiplier,coordinate


        ! select case (coord)
        !     case ('x'); ret=hypercharges(atom,charge)*mol%atm(atom)%coords(1)*multiplier
        !     case ('y'); ret=hypercharges(atom,charge)*mol%atm(atom)%coords(2)*multiplier
        !     case ('z'); ret=hypercharges(atom,charge)*mol%atm(atom)%coords(3)*multiplier
        ! end select
        ! return
        ! end function getAtomMoment

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        ! subroutine updateHeader(message)
        ! implicit none

        ! character (len=*)  :: message


        ! pcurr=pcurr+1
        ! iheader=' hypercharges: '//cmethod//', '//message//' ['//prStrByVal(pcurr)//'/'//prStrByVal(pcount)//'] '

        ! return
        ! end subroutine updateHeader

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
    end subroutine hyperdensity