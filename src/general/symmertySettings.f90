    subroutine symmetrySettings

    use glob,      only: assignment (=)
    use glob,      only: rglu,iglu,lglu,true,false,mid,void,i8kind,glControlMemory
    use glob,      only: uch,mid
    use hdb,       only: mol,geometrybd,polarizbd,generalbd,fieldbd,ou,ouWidth
    use hdb,       only: gEnergyHolder,GlEt,Et,MEt,MMEt,MethodListLen,cuebd,ierror
    use hdb,       only: pointAccordance,pointSet,pointToPut,pointToCalc,atomEqu,bondEqu
    use txtParser, only: tpAdjustc,tpGetSplit,operator(.in.)
    use printmod,  only: prMatrix
    use math,      only: tred4

    implicit none

    real   (kind=rglu) :: symmetryTolerance,gridRange,sum,ul,dl,diff
    integer(kind=iglu) :: N,M,Np,Ncue,dSize,i,j,k,l,a,b,c,sta,sto,cElements,ubndDistr,&
                          ipos,jpos,zpos,uniqueAtoms,writeMypnts,scaleIndices(3)

    real   (kind=rglu), allocatable :: coords(:,:),distance(:,:),eigenVectors(:,:),eigenValues(:),&
                                       d1gridEquivalence(:),&
                                       d2gridEquivalence(:,:),evecDifference(:,:),&
                                       d3gridEquivalence(:,:,:)

    integer(kind=iglu), allocatable :: d1grid(:),d2grid(:,:),d3grid(:,:,:),pntDistribution(:)
    logical(kind=lglu)              :: dcue


    symmetryTolerance=geometrybd%symmetryTolerance; N=mol%nAtoms; M=mol%nBonds; Ncue=mol%nEls/2

    ! Symmetry analysis is switched off
    if (.NOT.geometrybd%symmetryAccount) then
        select case( generalbd%task%get() )

            case ('polarizability')
                Np=polarizbd%nPoints; sta=-(Np-1)/2; sto=-sta
                void=glControlMemory(int( rglu*(Np+Np**2+Np**3*(2+6*MethodListLen)) ,kind=i8kind),'Symmetry settings')
                allocate (GlEt(sta:sto,sta:sto,sta:sto))
                allocate (gEnergyHolder(sta:sto,sta:sto,sta:sto,0:5,MethodListLen))
                allocate (Et(sta:sto),MEt(sta:sto,sta:sto),MMEt(sta:sto,sta:sto,sta:sto))

                allocate (pointSet       (3,Np**3  ))
                allocate (pointAccordance(Np**3,2,3))

                pointAccordance=0; pointSet=0; pointToCalc=1
                select case (polarizbd%scales%get())

                    ! single scale mode
                    case ('x','y','z')

                        scaleIndices=getScaleIndices()
                        ipos=scaleIndices(1)

                        pointToPut=0
                        do i = sto,sta,-1
                            ! define accordance of grid points (child == parent)
                            pointToPut=pointToPut+1
                            pointAccordance(pointToPut,1,ipos)=i
                            pointAccordance(pointToPut,2,ipos)=i
                        enddo

                        ! select all points
                        do i = sta,sto
                            ! zero point is already included
                            if (i.EQ.0) cycle
                            pointToCalc=pointToCalc+1
                            pointSet(ipos,pointToCalc)=i
                        enddo

                    ! double scale mode
                    case ('xy','yx','xz','zx','yz','zy')

                        scaleIndices=getScaleIndices()
                        ipos=scaleIndices(1)
                        jpos=scaleIndices(2)
                        zpos=scaleIndices(3)

                        pointToPut=0
                        do i = sto,sta,-1
                        do j = sto,sta,-1
                            ! define accordance of grid points (child == parent)
                            pointToPut=pointToPut+1
                            pointAccordance(pointToPut,1,ipos)=i
                            pointAccordance(pointToPut,1,jpos)=j

                            pointAccordance(pointToPut,2,ipos)=i
                            pointAccordance(pointToPut,2,jpos)=j
                        enddo
                        enddo

                        ! select all points
                        do i = sta,sto
                        do j = sta,sto
                            ! zero point is already included
                            if ((i.EQ.0).AND.(j.EQ.0)) cycle
                            pointToCalc=pointToCalc+1
                            pointSet(ipos,pointToCalc)=i
                            pointSet(jpos,pointToCalc)=j
                        enddo
                        enddo

                    ! full scale mode
                    case ('xyz')

                        scaleIndices=getScaleIndices()

                        pointToPut=0
                        do i = sto,sta,-1
                        do j = sto,sta,-1
                        do k = sto,sta,-1
                            ! define accordance of grid points (child == parent)
                            pointToPut=pointToPut+1
                            pointAccordance(pointToPut,1,1)=i
                            pointAccordance(pointToPut,1,2)=j
                            pointAccordance(pointToPut,1,3)=k

                            pointAccordance(pointToPut,2,1)=i
                            pointAccordance(pointToPut,2,2)=j
                            pointAccordance(pointToPut,2,3)=k
                        enddo
                        enddo
                        enddo

                        ! select all points
                        do i = sta,sto
                        do j = sta,sto
                        do k = sta,sto
                            ! such combinations do not exist
                            if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle

                            ! zero point is already included
                            if ((i.EQ.0).AND.(j.EQ.0).AND.(k.EQ.0)) cycle
                            pointToCalc=pointToCalc+1
                            pointSet(1,pointToCalc)=i
                            pointSet(2,pointToCalc)=j
                            pointSet(3,pointToCalc)=k
                        enddo
                        enddo
                        enddo

                end select

            case ('density','coulson','hypercharges')
                allocate (atomEqu(N,N),bondEqu(M,M)); atomEqu=1; bondEqu=1

                do i = 1,N
                    atomEqu(i,i)=i

                    if (i.EQ.N) exit
                    do j = i+1,N
                        atomEqu(i,j)=i*N+j
                    enddo
                enddo

                do i = 1,M
                    bondEqu(i,i)=i

                    if (i.EQ.M) exit
                    do j = i+1,M
                        bondEqu(i,j)=i*N+j
                    enddo
                enddo

        end select
        return
    endif

    dcue = ['cue-ccs', 'cue-ccsd', 'cue-ccsdt'] .in. tpGetSplit(generalbd%methods%get(), '+')
    if (dcue .AND. cuebd%nonset) then
        ierror='If some cue calculation is required, all bonds kind must be "single" or "double", not "nonset".'
        call primaryInformation('error')
    endif

    ! Symmetry non-demanding cases
    if (generalbd%task%get() .in. ['energy','wf-analysis']) then
        void=glControlMemory(int( MethodListLen*2*2*2*6 ,kind=i8kind),'Symmetry settings')
        allocate(gEnergyHolder(0:1,0:1,0:1,0:5,MethodListLen))
        return
    endif

    ! Start of symmetry analysis
    write (ou,'(/A/)') tpAdjustc('Symmetry analysis',ouWidth,'=')

    ! N atoms, 3 for field points
    dSize=N+3

    write (ou,120) symmetryTolerance
    select case( generalbd%task%get() )
        case ('polarizability')
            Np=polarizbd%nPoints; sta=-(Np-1)/2; sto=-sta
            ! Np*Np*Np: for every of three grid points.
            dSize=dSize+Np**3
            if (dcue) dSize=dSize+Ncue

        case ('density','coulson','hypercharges')
            ! N*(N-1)/2: centroids for all atom pairs.
            dSize=dSize+N*(N-1)/2
            if (dcue) dSize=dSize+Ncue

        case default
            stop 'Internal error (general::symmetrySettings): Unexpected general:jobtype value.'

    end select

    ! temporary array for analysis
    void=glControlMemory(int( rglu*(4*dSize+2*dSize*dSize) ,kind=i8kind),'tmp. Symmetry settings')
    allocate (coords(3,dSize),distance(dSize,dSize),eigenVectors(dSize,dSize),eigenValues(dSize))

    coords=0

    ! put atoms
    cElements=0
    do i = 1,N
        cElements=cElements+1
        do k = 1,3
            coords(k,cElements)=mol%atm(i)%coords(k)
        enddo
    enddo

    ! put cue orbitals
    if (dcue) then
        do i = 1,Ncue
            cElements=cElements+1
            do k = 1,3
                coords(k,cElements)=mol%orb(i)%coords(k)
            enddo
        enddo
    endif

    ! shift for field
    do k = 1,3
        cElements=cElements+1
        if (abs(fieldbd%strength(k)).GT.1D-10) then
            coords(k,cElements)=sign(maxval(abs(mol%atm(:)%coords(k))),fieldbd%strength(k))
        endif
    enddo

    ! set distance matrix structure
    select case( generalbd%task%get() )

        ! add grid to the set of distance matrix elements
        case ('polarizability')

            void=glControlMemory(int( 9*iglu*Np**3 ,kind=i8kind),'Symmetry settings')
            allocate (pointSet       (3,Np**3  ))
            allocate (pointAccordance(Np**3,2,3))

            gridRange=max(maxval(abs(mol%atm(:)%coords(1))),&
                          maxval(abs(mol%atm(:)%coords(2))),&
                          maxval(abs(mol%atm(:)%coords(3))) )/(Np-1)

            write (ou,121) gridRange

            do i = sta,sto
            do j = sta,sto
            do k = sta,sto
                cElements=cElements+1
                coords(1,cElements)=gridRange*i
                coords(2,cElements)=gridRange*j
                coords(3,cElements)=gridRange*k
            enddo
            enddo
            enddo

        ! add atom pairs centroids to the set of distance matrix elements
        case ('density','coulson','hypercharges')
            do i = 1,N-1
                do j = i+1,N
                    cElements=cElements+1
                    do k = 1,3
                        coords(k,cElements)=( mol%atm(i)%coords(k)+mol%atm(j)%coords(k) )/2
                    enddo
                enddo
            enddo

    end select

    ! proceed distance matrix
    !$omp parallel default(shared) private(i,j,k,sum)
    !$omp do
    do i = 1,dSize
    do j = i,dSize
        sum=0
        do k = 1,3
            sum=sum+(coords(k,i)-coords(k,j))**2
        enddo
        sum=sqrt(sum); distance(i,j)=sum; distance(j,i)=sum
    enddo
    enddo
    !$omp end parallel
    call tred4(distance,eigenVectors,eigenValues,dSize,1e-100_rglu,1e-100_rglu)

    select case( generalbd%task%get() )

        case ('polarizability')
            void=glControlMemory(int( rglu*(Np+Np**2+Np**3) ,kind=i8kind),'Symmetry settings')
            allocate (d1gridEquivalence(sta:sto),&
            &         d2gridEquivalence(sta:sto,sta:sto),&
            &         d3gridEquivalence(sta:sto,sta:sto,sta:sto))

            void=glControlMemory(int( rglu*(Np+Np**2+Np**3*(2+6*MethodListLen)) ,kind=i8kind),'Symmetry settings')
            allocate (GlEt(sta:sto,sta:sto,sta:sto))
            allocate (gEnergyHolder(sta:sto,sta:sto,sta:sto,0:5,MethodListLen))
            allocate (Et(sta:sto),MEt(sta:sto,sta:sto),MMEt(sta:sto,sta:sto,sta:sto))

            ! collect grid equivalence
            d3gridEquivalence=0
            do i = dSize-Np**3+1,dSize
                d3gridEquivalence( int(coords(1,i)/gridRange),int(coords(2,i)/gridRange),int(coords(3,i)/gridRange) )=eigenVectors(i,dSize)
            enddo

            ! proceed all possible points and get distribution
            ubndDistr=16; allocate (pntDistribution(0:ubndDistr))
            do l = 0,ubndDistr
                pntDistribution(l)=0; ul=10._rglu**(-l); dl=10._rglu**(-l-1)

                if (l.EQ.ubndDistr) dl=-1

                do i = sta,sto
                do j = sta,sto
                do k = sta,sto

                    ! such combinations do not exist
                    if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle

                    do a = sta,sto
                    do b = sta,sto
                    do c = sta,sto

                        ! such combinations do not exist
                        if ((a.NE.0).AND.(b.NE.0).AND.(c.NE.0)) cycle

                        ! full match
                        if ((i.EQ.a).AND.(j.EQ.b).AND.(k.EQ.c)) cycle

                        diff=abs(d3gridEquivalence(i,j,k)-d3gridEquivalence(a,b,c))

                        if ( (diff.GT.dl) .AND. (diff.LE.ul) ) then
                            pntDistribution(l)=pntDistribution(l)+1
                        endif
                    enddo
                    enddo
                    enddo
                enddo
                enddo
                enddo
            enddo

            call prMatrix(d3gridEquivalence(:,:,0),ou,'Matrix of xy grid equivalence','0.00000E00')
            call prMatrix(d3gridEquivalence(:,0,:),ou,'Matrix of xz grid equivalence','0.00000E00')
            call prMatrix(d3gridEquivalence(0,:,:),ou,'Matrix of yz grid equivalence','0.00000E00')

            void=glControlMemory(int( rglu*(Np+Np**2+Np**3) ,kind=i8kind),'tmp. Symmetry settings')
            allocate (d1grid(sta:sto),d2grid(sta:sto,sta:sto),d3grid(sta:sto,sta:sto,sta:sto))

            pointAccordance=0; pointSet=0; pointToCalc=1
            select case (polarizbd%scales%get())

                ! single scale mode
                case ('x','y','z')

                    scaleIndices=getScaleIndices()
                    ipos=scaleIndices(1)

                    select case (polarizbd%scales%get())
                        case ('x'); d1gridEquivalence=d3gridEquivalence(:,0,0)
                        case ('y'); d1gridEquivalence=d3gridEquivalence(0,:,0)
                        case ('z'); d1gridEquivalence=d3gridEquivalence(0,0,:)
                    end select

                    d1grid=0; pointToPut=0
                    do i = sto,sta,-1

                        ! already proceeded
                        if (d1grid(i).NE.0) cycle

                        do a = sto,sta,-1

                            ! already proceeded
                            if (d1grid(a).NE.0) cycle

                            diff=abs(d1gridEquivalence(i)-d1gridEquivalence(a))

                            ! define accordance of grid points (child == parent)
                            if (diff.LT.symmetryTolerance) then
                                pointToPut=pointToPut+1
                                pointAccordance(pointToPut,1,ipos)=i
                                pointAccordance(pointToPut,2,ipos)=a

                                ! label as child
                                d1grid(a)=-1
                            endif
                        enddo

                        ! label as parent
                        d1grid(i)=1
                    enddo

                    ! select points to be computed
                    do i = sta,sto
                        if (d1grid(i).EQ.1) then
                            if (i.EQ.0) cycle
                            pointToCalc=pointToCalc+1
                            pointSet(ipos,pointToCalc)=i
                        endif
                    enddo

                ! double scale mode
                case ('xy','yx','xz','zx','yz','zy')

                    scaleIndices=getScaleIndices()
                    ipos=scaleIndices(1)
                    jpos=scaleIndices(2)
                    zpos=scaleIndices(3)

                    select case (zpos)
                        case (1); d2gridEquivalence=d3gridEquivalence(0,:,:)
                        case (2); d2gridEquivalence=d3gridEquivalence(:,0,:)
                        case (3); d2gridEquivalence=d3gridEquivalence(:,:,0)
                    end select

                    d2grid=0; pointToPut=0
                    do i = sto,sta,-1
                    do j = sto,sta,-1

                        ! already proceeded
                        if (d2grid(i,j).NE.0) cycle

                        do a = sto,sta,-1
                        do b = sto,sta,-1

                            ! already proceeded
                            if (d2grid(a,b).NE.0) cycle

                            diff=abs(d2gridEquivalence(i,j)-d2gridEquivalence(a,b))

                            ! define accordance of grid points (child == parent)
                            if (diff.LT.symmetryTolerance) then
                                pointToPut=pointToPut+1
                                pointAccordance(pointToPut,1,ipos)=i
                                pointAccordance(pointToPut,1,jpos)=j

                                pointAccordance(pointToPut,2,ipos)=a
                                pointAccordance(pointToPut,2,jpos)=b

                                ! label as child
                                d2grid(a,b)=-1
                            endif
                        enddo
                        enddo

                        ! label as parent
                        d2grid(i,j)=1
                    enddo
                    enddo

                    ! select points to be computed
                    do i = sta,sto
                    do j = sta,sto
                        if (d2grid(i,j).EQ.1) then

                            ! zero point is already included
                            if ((i.EQ.0).AND.(j.EQ.0)) cycle
                            pointToCalc=pointToCalc+1
                            pointSet(ipos,pointToCalc)=i
                            pointSet(jpos,pointToCalc)=j
                        endif
                    enddo
                    enddo

                ! full scale mode
                case ('xyz')

                    scaleIndices=getScaleIndices()

                    d3grid=0; pointToPut=0
                    do i = sto,sta,-1
                    do j = sto,sta,-1
                    do k = sto,sta,-1

                        ! such combinations do not exist
                        if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle

                        ! already proceeded
                        if (d3grid(i,j,k).NE.0) cycle

                        do a = sto,sta,-1
                        do b = sto,sta,-1
                        do c = sto,sta,-1

                            ! already proceeded
                            if (d3grid(a,b,c).NE.0) cycle

                            diff=abs(d3gridEquivalence(i,j,k)-d3gridEquivalence(a,b,c))

                            ! define accordance of grid points (child == parent)
                            if (diff.LT.symmetryTolerance) then
                                pointToPut=pointToPut+1
                                pointAccordance(pointToPut,1,1)=i
                                pointAccordance(pointToPut,1,2)=j
                                pointAccordance(pointToPut,1,3)=k

                                pointAccordance(pointToPut,2,1)=a
                                pointAccordance(pointToPut,2,2)=b
                                pointAccordance(pointToPut,2,3)=c

                                ! label as child
                                d3grid(a,b,c)=-1
                            endif
                        enddo
                        enddo
                        enddo

                        ! label as parent
                        d3grid(i,j,k)=1
                    enddo
                    enddo
                    enddo

                    ! select points to be computed
                    do i = sta,sto
                    do j = sta,sto
                    do k = sta,sto

                        ! such combinations do not exist
                        if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle

                        if (d3grid(i,j,k).EQ.1) then

                            ! zero point is already included
                            if ((i.EQ.0).AND.(j.EQ.0).AND.(k.EQ.0)) cycle
                            pointToCalc=pointToCalc+1
                            pointSet(1,pointToCalc)=i
                            pointSet(2,pointToCalc)=j
                            pointSet(3,pointToCalc)=k
                        endif
                    enddo
                    enddo
                    enddo

            end select

            ! free temporary arrays
            deallocate ( d1grid,d2grid,d3grid )
            deallocate ( d1gridEquivalence,d2gridEquivalence,d3gridEquivalence)
            void=glControlMemory(int( sizeof(d1grid)+&
                                      sizeof(d2grid)+&
                                      sizeof(d3grid) ,kind=i8kind),&
                                      'tmp. Symmetry settings','free')

            void=glControlMemory(int( sizeof(d1gridEquivalence)+&
                                      sizeof(d2gridEquivalence)+&
                                      sizeof(d3gridEquivalence) ,kind=i8kind),&
                                      'tmp. Symmetry settings','free')

            ! full set of points
            select case (polarizbd%scales%ln)
                case (1); writeMypnts=Np
                case (2); writeMypnts=Np**2
                case (3); writeMypnts=3*(Np**2-2*Np+1)+3*(Np-1)+1 !unique+cross{x00,0y0,00z}+zero{000}
            end select

            ! show grid points accordance
            write (ou,100)
            do i = 1,pointToPut
                write (ou,101) (pointAccordance(i,1,k), k=1,3),(pointAccordance(i,2,k), k=1,3)
            enddo
            write (ou,110) pointToCalc,writeMypnts
            write (ou,*)

            ! show distribution
            write (ou,130) (3*(Np**2-2*Np+1)+3*(Np-1)+1) * (3*(Np**2-2*Np+1)+3*(Np-1))
            do l = 0,ubndDistr
                ul=10._rglu**(-l); dl=10._rglu**(-l-1)
                if (l.EQ.ubndDistr) dl=0
                write (ou,131) ul,dl,pntDistribution(l)
            enddo
            write (ou,*)

            deallocate (pntDistribution)

        case ('density','coulson','hypercharges')
            void=glControlMemory(int( iglu*(N*N+M*M)+rglu*dSize*dSize ,kind=i8kind),'Symmetry settings')
            allocate (atomEqu(N,N),bondEqu(M,M)); atomEqu=0; bondEqu=0
            allocate (evecDifference(dSize,dSize))

            ! prepare storage
            do i = 1,dSize
            do j = i,dSize
                evecDifference(i,j)=abs(eigenVectors(i,dSize)-eigenVectors(j,dSize))
                evecDifference(j,i)=evecDifference(i,j)
            enddo
            enddo

            ! get symmetry information
            call defAtomsEquivalence
            call defBondsEquivalence

            ubndDistr=16
            allocate (pntDistribution(0:ubndDistr))

            uniqueAtoms=0
            do i = 1,N
                do j = i,N
                    if (atomEqu(i,j).GT.0) uniqueAtoms=uniqueAtoms+1
                enddo
            enddo
            write (ou,133) uniqueAtoms

            ! compute distribution
            do l = 0,ubndDistr
                pntDistribution(l)=0; ul=10._rglu**(-l); dl=10._rglu**(-l-1)

                if (l.EQ.ubndDistr) dl=-1

                do i = 1,N
                do j = 1,N

                    if ( (evecDifference(i,j).GT.dl) .AND. (evecDifference(i,j).LE.ul) ) then
                        pntDistribution(l)=pntDistribution(l)+1
                    endif
                enddo
                enddo
            enddo

            ! show distribution
            write (ou,132) N*N
            do l = 0,ubndDistr
                ul=10._rglu**(-l); dl=10._rglu**(-l-1)
                if (l.EQ.ubndDistr) dl=0
                write (ou,131) ul,dl,pntDistribution(l)
            enddo
            write (ou,*)

            ! call prMatrix(evecDifference(1:N,1:N),ou,'Atom equivalence','0.00E00',maxwidth=ouWidth)

            ! free temporary arrays
            deallocate (pntDistribution,evecDifference)
            void=glControlMemory(int( sizeof(evecDifference) ,kind=i8kind),'tmp. Symmetry settings','free')

        case ('wf-analysis')
            continue

    end select

    write (ou,'(/A/)') tpAdjustc('The start of computation',ouWidth,'=')

100 format ( 4X,'Point accordance as a consequence of symmetry transformation.')
110 format ( 4X,'Points to be calculated: ',i<mid(pointToCalc)>,' from ',i<mid(writeMypnts)>)
101 format ( 4X,'(',i2,',',i2,',',i2,')--->(',i2,',',i2,',',i2,')')
120 format ( 4X,'Tollerance in grid symmetry definition:',1X,ES9.3)
121 format ( 4X,'Range for grid points:',1X,ES9.3)
130 format ( 4X,'Combinations in 3d grid',1X,i6,'. Distribution:')
131 format ( 4X,ES9.3,1X,'<-->',1X,ES9.3,1X,i5,1X,'combinations')
132 format ( 4X,'Atom-atom combinations ',1X,i6,'. Distribution:')
133 format (/4X,'Unique atom-atom pairs:',1X,i5)

    return

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine defAtomsEquivalence
        implicit none

        logical(kind=lglu) :: once
        integer(kind=iglu) :: i,j,k,l,a,b,shft


        shft=dSize-N*(N-1)/2
        do i = 1,N
            if (atomEqu(i,i).NE.0) cycle
            atomEqu(i,i)=i

            if (i+1.GT.N) exit

            do j = i+1,N
                if (evecDifference(i,j).LE.symmetryTolerance) then
                    atomEqu(j,j)=-i
                endif
            enddo
        enddo

        k=0
        do i = 1,N-1
        do j = i+1,N
            k=k+1; if (atomEqu(i,j).NE.0) cycle
            atomEqu(i,j)=i*N+j

            l=0
            do a = 1,N-1
            do b = a+1,N
                l=l+1; if (atomEqu(a,b).NE.0) cycle

                if (evecDifference(shft+k,shft+l).LE.symmetryTolerance) then
                    if ((abs(distance(i,j)-distance(a,b)).LE.symmetryTolerance)&
                    &  .AND.( abs(atomEqu(i,i)).EQ.abs(atomEqu(a,a)) )&
                    &  .AND.( abs(atomEqu(j,j)).EQ.abs(atomEqu(b,b)) )) atomEqu(a,b)=-atomEqu(i,j)

                    if ((abs(distance(i,j)-distance(a,b)).LE.symmetryTolerance)&
                    &  .AND.( abs(atomEqu(i,i)).EQ.abs(atomEqu(b,b)) )&
                    &  .AND.( abs(atomEqu(j,j)).EQ.abs(atomEqu(a,a)) )) atomEqu(a,b)=-atomEqu(i,j)
                endif
            enddo
            enddo
        enddo
        enddo
        call prMatrix(atomEqu,ou,'Atom equivalence','^')

        k=0
        write (ou,'(4X,A/)') 'Unique atoms:'
        do i = 1,N
            if (atomEqu(i,i).LE.0) cycle
            k=k+1
            once=false
            write (ou,'(4X,i<mid(N)>,")",1X,i<mid(N)>\)') k,i
            do j = 1,N
                if ((atomEqu(j,j).LT.0).AND.(abs(atomEqu(j,j)).EQ.atomEqu(i,i))) then
                    if (once) then
                        write (ou,'(",",1X,i<mid(N)>\)') j
                    else
                        write (ou,'(1X,"- equal:",1X,i<mid(N)>\)') j
                        once=true
                    endif
                endif
            enddo
            write (ou,*)
        enddo
        write (ou,*)

        call defInverse

        return
        end subroutine defAtomsEquivalence

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine defBondsEquivalence
        implicit none

        logical(kind=lglu), allocatable :: d(:,:)
        integer(kind=iglu), allocatable :: bond(:,:)
        real   (kind=rglu), allocatable :: X(:),Y(:),Z(:)
        integer(kind=iglu)              :: i,j,a,b,k(4),l(4),aa,bb,cc,dd
        real   (kind=rglu)              :: centr(2)
        logical(kind=lglu)              :: once


        allocate (bond(2,M))
        do b = 1,M
            bond(1,b)=mol%bnd(b)%atoms(1)
            bond(2,b)=mol%bnd(b)%atoms(2)
        enddo

        allocate (X(N),Y(N),Z(N))
        do i = 1,N
            X(i)=mol%atm(i)%coords(1)
            Y(i)=mol%atm(i)%coords(2)
            Z(i)=mol%atm(i)%coords(3)
        enddo

        allocate (d(N,N))

        !$omp parallel default(shared) private(j)
        !$omp do
        do i = 1,N
            do j = 1,N
                d(i,j)=abs(atomEqu(i,i)).EQ.abs(atomEqu(j,j))
            enddo
        enddo
        !$omp end parallel

        write (ou,'(4X,A/)') 'Bonds:'
        do i = 1,M
            write (ou,'(4X,i<mid(M)>,")",2X,i<mid(N)>," <--> ",i<mid(N)>)') i,bond(1,i),bond(2,i)
        enddo
        write (ou,*)

        do i = 1,M
            if (bondEqu(i,i).NE.0) cycle
            bondEqu(i,i)=i

            if (i+1.GT.M) exit

            aa=bond(1,i); bb=bond(2,i)

            do j = i+1,M
                cc=bond(1,j); dd=bond(2,j)

                if ((abs(distance(aa,bb)-distance(cc,dd)).LE.symmetryTolerance)&
                &  .AND.( d(aa,cc) )&
                &  .AND.( d(bb,dd) )) bondEqu(j,j)=-i

                if ((abs(distance(aa,bb)-distance(cc,dd)).LE.symmetryTolerance)&
                &  .AND.( d(aa,dd) )&
                &  .AND.( d(bb,cc) )) bondEqu(j,j)=-i

            enddo
        enddo

        do i = 1,M-1
            k(1)=bond(1,i); k(2)=bond(2,i)
            do j = i+1,M
                if (bondEqu(i,j).NE.0) cycle
                bondEqu(i,j)=i*M+j

                k(3)=bond(1,j); k(4)=bond(2,j)

                centr(1)=0
                centr(1)=centr(1)+( X(k(1))+X(k(2))-X(k(3))-X(k(4)) )**2
                centr(1)=centr(1)+( Y(k(1))+Y(k(2))-Y(k(3))-Y(k(4)) )**2
                centr(1)=centr(1)+( Z(k(1))+Z(k(2))-Z(k(3))-Z(k(4)) )**2
                centr(1)=sqrt(centr(1))

                do a = 1,M-1
                    l(1)=bond(1,a); l(2)=bond(2,a)
                    do b = a+1,M
                        if (bondEqu(a,b).NE.0) cycle


                        l(3)=bond(1,b); l(4)=bond(2,b)

                        centr(2)=0
                        centr(2)=centr(2)+( X(l(1))+X(l(2))-X(l(3))-X(l(4)) )**2
                        centr(2)=centr(2)+( Y(l(1))+Y(l(2))-Y(l(3))-Y(l(4)) )**2
                        centr(2)=centr(2)+( Z(l(1))+Z(l(2))-Z(l(3))-Z(l(4)) )**2
                        centr(2)=sqrt(centr(2))

                        if ( abs(centr(1)-centr(2)).LE.symmetryTolerance) then
                            if (bEqu(k,l,d)) bondEqu(a,b)=-bondEqu(i,j)
                        endif
                    enddo
                enddo
            enddo
        enddo
        call prMatrix(bondEqu,ou,'Bond equivalence','^')

        a=0
        write (ou,'(4X,A/)') 'Unique bonds:'
        do i = 1,M
            if (bondEqu(i,i).LE.0) cycle
            a=a+1
            once=false
            write (ou,'(4X,i<mid(M)>,")",1X,i<mid(M)>\)') a,i
            do j = 1,M
                if ((bondEqu(j,j).LT.0).AND.(abs(bondEqu(j,j)).EQ.bondEqu(i,i))) then
                    if (once) then
                        write (ou,'(",",1X,i<mid(M)>\)') j
                    else
                        write (ou,'(1X,"- equal:",1X,i<mid(M)>\)') j
                        once=true
                    endif
                endif
            enddo
            write (ou,*)
        enddo
        write (ou,*)

        deallocate (d,bond,X,Y,Z)

        return
        end subroutine defBondsEquivalence

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        logical*1 function bEqu(k1,k2,d) result(ret)
        implicit none
        integer(kind=iglu), intent(in) :: k1(4),k2(4)
        logical(kind=lglu), intent(in) :: d(:,:)
        integer(kind=iglu)             :: i,j,a,b


        ret=false
        do i = 1,4
            if (.NOT.d(k1(1),k2(i))) cycle
            do j = 1,4
                if (i.EQ.j) cycle; if (.NOT.d(k1(2),k2(j))) cycle
                do a = 1,4
                    if ((a.EQ.i).OR.(a.EQ.j)) cycle; if (.NOT.d(k1(3),k2(a))) cycle
                    do b = 1,4
                        ret=(b.NE.i).AND.(b.NE.j).AND.(b.NE.a).AND.d(k1(4),k2(b))
                        if (ret) return
                    enddo
                enddo
            enddo
        enddo

        return
        end function bEqu

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        subroutine defInverse
        implicit none

        integer(kind=iglu) :: i,j
        real   (kind=rglu) :: dx,dy,dz


        ! to be finished
        mol%inverse=0
        do i = 1,N-1
            do j = i+1,N
                dx=abs(mol%atm(i)%coords(1)+mol%atm(j)%coords(1))
                dy=abs(mol%atm(i)%coords(2)+mol%atm(j)%coords(2))
                dz=abs(mol%atm(i)%coords(3)+mol%atm(j)%coords(3))

                if ((max(dx,dy,dz).LT.symmetryTolerance).AND.(abs(atomEqu(j,j)).EQ.abs(atomEqu(i,i)))) then
                    mol%inverse(i)=j
                    mol%inverse(j)=i
                endif
            enddo
        enddo

        return
        end subroutine defInverse

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

        function getScaleIndices() result(ret)
        implicit none

        integer(kind=iglu) :: ipos,jpos,kpos,ret(3)


        select case (polarizbd%scales%ln)

            case(1)
                select case (polarizbd%scales%get())
                    case ('x'); ipos=1
                    case ('y'); ipos=2
                    case ('z'); ipos=3
                end select

                ret=[ipos,0,0]

            case(2)
                ! normalize scale label
                if (polarizbd%scales%get().EQ.'yx') then
                    polarizbd%scales='xy'
                elseif(polarizbd%scales%get().EQ.'zx') then
                    polarizbd%scales='xz'
                elseif(polarizbd%scales%get().EQ.'zy') then
                    polarizbd%scales='yz'
                endif

                ! determine appropriate case
                select case (polarizbd%scales%get(1,1))
                    case ('x'); ipos=1
                    case ('y'); ipos=2
                    case ('z'); ipos=3
                end select

                select case (polarizbd%scales%get(2,2))
                    case ('x'); jpos=1
                    case ('y'); jpos=2
                    case ('z'); jpos=3
                end select

                if (ipos.GT.jpos) then
                    i=jpos; jpos=ipos; ipos=i
                endif

                if (ipos.EQ.1) then
                    if (jpos.EQ.2) zpos=3
                    if (jpos.EQ.3) zpos=2
                else;              zpos=1
                endif

                ret=[ipos,jpos,zpos]

            case(3)
                ret=[1,2,3]

        end select

        return
        end function getScaleIndices

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine symmetrySettings
