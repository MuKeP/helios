    subroutine ivv_dens_hyper

    use glob,           only: iglu,rglu,system,purifyValues,NaN
    use glob,           only: lglu,true,false,uch,void,gluCompare,timeControlCheckpoint
    use printmod,       only: prStrByVal,prMatrix,prEigenProblem,prLongText
    use hdb,            only: mol,hyperchargesbd,HartreeEnergy,BohrRadius,dipoleToDeby
    use hdb,            only: perturbate,ou,ccbd,setIterationHeader
    use coupledCluster, only: R,F,Nel,No,t1,t2,hV,putCUEMOs
    use property,       only: getEnergy
    use math,           only: LagrangeDerivative
    use txtParser,      only: tpSplit,tpSplitLen,tpSplitHold,tpFill,tpAdjustc,operator(.in.)

    implicit none

    real(kind=rglu)  , allocatable  :: allDens(:,:)
    integer(kind=iglu)              :: N,nScales,naPoints,nfPoints,asta,asto,fsta,fsto
    integer(kind=iglu)              :: i,j,k,l,c,d,sc,pcount,pcurr,aCount,svio,meth,isc
    real   (kind=rglu)              :: fStep,aStep,zpEnergy,rez,hc(8),hm(8)
    character(len=1)                :: scl,transform
    character (len=:), allocatable  :: cmethod

    real   (kind=rglu), allocatable :: fgrid(:,:,:,:),hypercharges(:,:,:)

    real   (kind=rglu)              :: transition(0:8)
    real   (kind=rglu), parameter   :: derivThreshold(0:4)=[1E-5_rglu,1E-4_rglu,1E-3_rglu,0.01_rglu,0.1_rglu]


    !open(50,file='alldens')
    ! prepare transition constants to atom units
    transition=NaN
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

    write (100,*) N,nScales,hyperchargesbd%scales%get()
    write (100,*) naPoints,nfPoints
    write (100,*) aStep,fStep

    allocate (fgrid(N,fsta:fsto,3,4),hypercharges(N,0:9,4))
    hypercharges=NaN; fgrid=NaN

    allocate(allDens(N,4)); allDens=NaN

    cmethod='r-ccsd'

    pcount=1 + nScales*(nfPoints-1)
    pcurr=0; void=updateHeader('zero point energy')
    call setCore(0,0,0,0,0); zpEnergy=getEnergy(cmethod)

    call communicate

    ! calculation of zero field charges
    do k = 1,4
        fgrid(:,0,1,k)=allDens(:,k)
        fgrid(:,0,2,k)=allDens(:,k)
        fgrid(:,0,3,k)=allDens(:,k)
    enddo

    do i = 1,N
        write(ou,'(2X,i2,4(2X,F7.4))') i,allDens(i,:)
    enddo
    stop
    ! calculation of in-field charges
    do sc = 1,nScales
        scl=hyperchargesbd%scales%get(sc,sc)

        do k = fsta,fsto
            if (k.EQ.0) cycle
            void=updateHeader('Scale '//scl//', Field '//prStrByVal(k))

            ! setting core and calculation of in-field RDM
            select case (scl)
                case ('x'); call setCore(0,0,k,0,0); rez=getEnergy(cmethod); isc=1
                case ('y'); call setCore(0,0,0,k,0); rez=getEnergy(cmethod); isc=2
                case ('z'); call setCore(0,0,0,0,k); rez=getEnergy(cmethod); isc=3
            end select

            call communicate

            ! collecting in-filed charges
            do l = 1,4
                fgrid(:,k,isc,l)=allDens(:,l)
            enddo

        enddo
    enddo

    hypercharges=NaN
    do k = 1,4
        do i = 1,N

            ! calculation of hypercharges
            do sc = 1,nScales
                scl=hyperchargesbd%scales%get(sc,sc)
                select case (scl)
                    case ('x'); isc=1
                    case ('y'); isc=2
                    case ('z'); isc=3
                end select

                hypercharges(i,0,k)=fgrid(i,0,isc,k)
                hypercharges(i,isc+0,k)=LagrangeDerivative(nfPoints,1,fgrid(i,:,isc,k),fStep)
                hypercharges(i,isc+3,k)=LagrangeDerivative(nfPoints,2,fgrid(i,:,isc,k),fStep)
                hypercharges(i,isc+6,k)=LagrangeDerivative(nfPoints,3,fgrid(i,:,isc,k),fStep)
            enddo

        enddo
    enddo


    ! output hypercharges
    do d = 1,4
        do sc = 1,nScales
            hc=0; hm=0
            scl=hyperchargesbd%scales%get(sc,sc)
            write (ou,100) '',           tpFill(1,scl),&
                           tpFill(1,scl),tpFill(2,scl),&
                           tpFill(2,scl),tpFill(3,scl),&
                           tpFill(3,scl),tpFill(4,scl)

            !hc=NaN; hm=NaN
            do i = 1,N
                select case (scl)
                    case('x'); l=1
                    case('y'); l=2
                    case('z'); l=3
                end select
                hc(1)=hypercharges(i,0  ,d); hm(1)=getAtomMoment(i,0  ,scl,transition(1),d)
                hc(2)=hypercharges(i,l  ,d); hm(2)=getAtomMoment(i,l  ,scl,transition(2),d)
                hc(3)=hypercharges(i,l+3,d); hm(3)=getAtomMoment(i,l+3,scl,transition(3),d)
                hc(4)=hypercharges(i,l+6,d); hm(4)=getAtomMoment(i,l+6,scl,transition(4),d)

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

100 format (/1X,'Atom',4X,'Q',A,9X,'M',A,8X,'Q',A,7X,'M',A,8X,'Q',A,6X,'M',A,7X,'Q',A,5X,'M',A)
101 format ( 1X,i3    ,1X,4(1X,F7.4,1X,ES11.4))
102 format ( 5X,A/5X,F8.4,1X,ES11.4,3(1X,F7.4,1X,ES11.4))

    stop

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine setCore(anumber,ashift,xshift,yshift,zshift)
        implicit none

        integer(kind=iglu), intent(in) :: anumber,ashift,xshift,yshift,zshift
        integer(kind=iglu)             :: atom


        mol%perturbation=0
        do atom = 1,N
            mol%perturbation(atom,atom)=+xshift*fStep*mol%atm(atom)%coords(1)&
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

        real(kind=rglu) function getAtomMoment(atom,charge,coord,multiplier,type) result(ret)
        implicit none

        integer(kind=iglu) :: charge,atom,type
        character (len=*)  :: coord
        real   (kind=rglu) :: multiplier,coordinate


        select case (coord)
            case ('x'); ret=hypercharges(atom,charge,type)*mol%atm(atom)%coords(1)*multiplier
            case ('y'); ret=hypercharges(atom,charge,type)*mol%atm(atom)%coords(2)*multiplier
            case ('z'); ret=hypercharges(atom,charge,type)*mol%atm(atom)%coords(3)*multiplier
        end select
        return
        end function getAtomMoment

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        integer(kind=iglu) function updateHeader(message) result(ret)
        implicit none

        character (len=*)  :: message


        pcurr=pcurr+1
        call setIterationHeader(' hypercharges: '//cmethod//', '//message//' ['//prStrByVal(pcurr)//'/'//prStrByVal(pcount)//'] ')

        ret=0; return
        end function updateHeader

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine communicate
        implicit none

        integer(kind=iglu) :: i,j,a,b,k,l,N,err
        real(kind=rglu)    :: threshold !,val


        threshold=1D-12
        N=mol%nAtoms

        ! prepare input
        open(80, file='input.inp')
        write (80,'(A)') mol%name%get()
        write (80,*) mol%nAtoms
        write (80,*) mol%nBonds
        write (80,*) mol%naEls
        write (80,*) mol%nbEls
        close (80)

        ! prepare fort.88
        open(80, file='fort.88')
        do i = 1,No
        do j = 1,No
        do k = 1,No
        do l = 1,No
            if (abs(R(i,k,j,l)).GT.threshold) then
                write (80,*) i,j,k,l,R(i,k,j,l)
            endif
        enddo
        enddo
        enddo
        enddo
        close (80)

        ! prepare fort.89
        open(80, file='fort.89')
        do i = 1,No
        do j = 1,No
            if (abs(F(i,j)).GT.threshold) then
                write (80,*) i,j,F(i,j)
            endif
        enddo
        enddo
        close(80)

        ! prepare t1.bin
        open(80, file='t1.bin', form='unformatted')
        write (80) Nel,Nel
        write (80) ((t1(i,a), i=1,Nel), a=Nel+1,No)
        close(80)

        ! prepare t2.bin
        open(80, file='t2.bin', form='unformatted')
        write (80) Nel,Nel
        write (80) ((((t2(i,j,a,b), i=1,Nel), j=1,Nel), a=Nel+1,No), b=Nel+1,No)
        close(80)

        ! prepare MO.dat
        call putCUEMOs
        open(80, file='MO.dat')
        write (80,*) N
        do i = 1,N
            write (80,'(<N>F20.15)') hV(i,:)
        enddo
        close(80)
        call putCUEMOs(dealloc=true)

        ! calling ivv density
        err=system('./lrdens.v1.x')

        ! getting densities
        open(80, file='density.dat')
        do k = 1,4
            do i = 1,N
                read(80,*) allDens(i,k)
            enddo
        enddo
        close(80)

        err=system('rm fort.88 fort.89 t1.bin t2.bin MO.dat out_LR-CCSD-dens l01.bin l02.bin density.dat input.inp')
        return
        end subroutine communicate


    end subroutine ivv_dens_hyper