    subroutine readMoleculeInformation

    use glob,      only: assignment(=)
    use glob,      only: iglu,rglu,lglu,true,false,void,i8kind,glControlMemory
    use glob,      only: uch
    use glob,      only: mid,rangen
    use hdb,       only: in,ou,ouWidth,ierror
    use hdb,       only: BohrRadius,HartreeEnergy
    use hdb,       only: mol,generalbd,geometrybd,parametrizationbd,cuebd
    use txtParser, only: operator(.in.)
    use printmod,  only: prMatrix

    implicit none

    integer(kind=iglu)  :: N,M,Nel,No,i,j,k,l
    real(kind=rglu)     :: coords(3),val
    character (len=256) :: rname


    read (in,*,err=666) rname; mol%name=trim(rname)
    read (in,*,err=666) mol%nAtoms
    read (in,*,err=666) mol%nBonds
    read (in,*,err=666) mol%naEls
    read (in,*,err=666) mol%nbEls

    call checkInput

    N=mol%nAtoms; M=mol%nBonds

    mol%nEls=mol%naEls+mol%nbEls
    mol%nOrbs=mol%nAtoms
    mol%nEthylenes=mol%nEls/2

    Nel=mol%nEls; No=mol%nOrbs

    void=glControlMemory(int( 4*rglu*N*N + N*(3*iglu+6*rglu+2) + M*(3*iglu+5*rglu+6) + N*(7*iglu+3*rglu) ,kind=i8kind),'Molecule information')
    allocate (mol%atm(N),mol%bnd(M),mol%orb(N))
    allocate (mol%g(N,N),mol%connect(N,N),mol%atmDist(N,N),mol%cueDist(N,N))
    allocate (mol%inverse(N))

    mol%connect=0; mol%atmDist=0; mol%cueDist=0; mol%inverse=0

    do i = 1,N
        read (in,*,err=667) mol%atm(i)%symbol,&
                            mol%atm(i)%coords(1),mol%atm(i)%coords(2),mol%atm(i)%coords(3),&
                            mol%atm(i)%alpha    ,mol%atm(i)%gamma    ,mol%atm(i)%nels

        mol%atm(i)%number=i

        if (mol%atm(i)%symbol.EQ.'C ') then
            mol%atm(i)%atype=1
        elseif ( (mol%atm(i)%symbol.EQ.'N ').AND.(mol%atm(i)%nels.EQ.1) ) then
            mol%atm(i)%atype=2
        elseif ( (mol%atm(i)%symbol.EQ.'O ').AND.(mol%atm(i)%nels.EQ.2) ) then
            mol%atm(i)%atype=3
        elseif ( (mol%atm(i)%symbol.EQ.'N ').AND.(mol%atm(i)%nels.EQ.2) ) then
            mol%atm(i)%atype=4
        else
            mol%atm(i)%atype=0
        endif
    enddo

    call geometryDisplacement
    call centrateCoordinates
    call planarityCheck

    !$omp parallel default(shared) private(i,j,k,val)
    !$omp do
    do i = 1,N-1
        do j = i+1,N
            val=0
            do k = 1,3
                val=val+( mol%atm(i)%coords(k)-mol%atm(j)%coords(k) )**2
            enddo
            mol%atmDist(i,j)=sqrt(val)
            mol%atmDist(j,i)=sqrt(val)
        enddo
    enddo
    !$omp end parallel

    call gamma

    do l = 1,3
        mol%clearance(l)=maxval(mol%atm(:)%coords(l))-minval(mol%atm(:)%coords(l))
    enddo

    do k = 1,M
        read (in,*,err=668) mol%bnd(k)%atoms(1),mol%bnd(k)%atoms(2),&
                            mol%bnd(k)%beta    ,mol%bnd(k)%kind

        if (mol%bnd(k)%kind.EQ.'nonset') then
            cuebd%nonset=true
        endif

        i=mol%bnd(k)%atoms(1); j=mol%bnd(k)%atoms(2)

        do l = 1,3
            coords(l)=(mol%atm(i)%coords(l)+mol%atm(j)%coords(l))/2
        enddo
        mol%bnd(k)%coords=coords

        do l = 1,3
            coords(l)=(mol%atm(i)%coords(l)-mol%atm(j)%coords(l))**2
        enddo
        mol%bnd(k)%length=sqrt(sum(coords))
    enddo

    do i = 1,N
        mol%connect(i,i)=mol%atm(i)%alpha
    enddo

    do k = 1,M
        i=mol%bnd(k)%atoms(1); j=mol%bnd(k)%atoms(2)

        mol%connect(i,j)=mol%bnd(k)%beta
        mol%connect(j,i)=mol%bnd(k)%beta
    enddo

    call checkUniqueness

    if (parametrizationbd%bondsAlternated) then
        do k = 1,M
            i=mol%bnd(k)%atoms(1); j=mol%bnd(k)%atoms(2)

            select case (mol%bnd(k)%kind)
                case ('single')
                    mol%connect(i,j)=mol%bnd(k)%beta*(1-parametrizationbd%alternation)
                    mol%connect(j,i)=mol%bnd(k)%beta*(1-parametrizationbd%alternation)
                case ('double')
                    mol%connect(i,j)=mol%bnd(k)%beta*(1+parametrizationbd%alternation)
                    mol%connect(j,i)=mol%bnd(k)%beta*(1+parametrizationbd%alternation)
                case ('nonset')
                    exit
                case default
                    goto 668
            endselect
        enddo
    endif

    void=glControlMemory(int( 4*rglu*N*N ,kind=i8kind),'Molecule information')
    allocate (mol%core(N,N),mol%coreImage(N,N),mol%huckelCore(N,N),mol%perturbation(N,N))
    mol%perturbation=0; mol%huckelCore=0
    do i = 1,N
    do j = 1,N
        mol%core(i,j)=mol%connect(i,j)
        mol%huckelCore(i,j)=mol%connect(i,j)
        mol%coreImage(i,j)=mol%connect(i,j)
    enddo
    enddo

    do k = 1,N
        val=0
        do l = 1,N
            if (k.EQ.l) cycle
            val=val+mol%atm(l)%nels*mol%G(k,l)
        enddo
        mol%core(k,k)=mol%core(k,k)-val
        mol%coreImage(k,k)=mol%core(k,k)
    enddo

    !ou=6

    write (ou, 99) mol%name%get(),mol%nAtoms,mol%nBonds,mol%nEls
    write (ou,100)
    do i = 1,N
        write (ou,101) mol%atm(i)%symbol,mol%atm(i)%coords(1),mol%atm(i)%coords(2),&
                       mol%atm(i)%coords(3),mol%atm(i)%alpha,int(mol%atm(i)%nels)
    enddo
    call prMatrix(mol%atmdist,ou,'Distance matrix','^.0000',maxwidth=ouWidth)

    write (ou,102) mol%clearance(:),mol%clearance(:)/mol%nEls

    write (ou,103)
    do k = 1,M
        i=mol%bnd(k)%atoms(1)
        j=mol%bnd(k)%atoms(2)

        write (ou,104) k,i,j,mol%atmdist(i,j),mol%bnd(k)%beta,mol%bnd(k)%kind
    enddo

    call prMatrix(mol%G,ou,'Gammas: '//parametrizationbd%coulombType%get(),'^.0000',maxwidth=ouWidth)

    if (.NOT.cuebd%nonset) then
        call cueOrbitals
    endif
    call symmetrySettings

    return

 99 format (/2X,'Molecule:',1X,A/&
             2X,'Atoms: ',i<mid(mol%nAtoms)>,' Bonds: ',i<mid(mol%nBonds)>,' Electrons: ',i<mid(mol%nEls)>//)
100 format ( 2X,'Type',23X,'Cartesian Coordinates',22X,'IP',3X,'e'/)
101 format ( 4X,A,1X,3(1X,F19.12),1X,F8.3,1X,i1)

102 format ( 2X,'Molecule overall dimensions:  X:',F11.5,3X,'Y:',F11.5,3X,'Z:',F11.5/&
             2X,'    specific (per electron):  X:',F11.5,3X,'Y:',F11.5,3X,'Z:',F11.5)

103 format (/2X,'bond',8X,'mu',10X,'nu',4X,'distance',8X,'beta',8X,'kind'/)
104 format ( 2X,i4,6X,i4,1X,'----->',1X,i4,2X,F10.5,2X,F10.5,1X,'eV',4X,A)

666 ierror='Molecule block. FileFormatError. Incorrect header.'
    call primaryInformation('error')

667 ierror='Molecule block. FileFormatError. Incorrect atoms information.'
    call primaryInformation('error')

668 ierror='Molecule block. FileFormatError. Incorrect bonds information.'
    call primaryInformation('error')

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

        subroutine gamma
        implicit none

        integer(kind=iglu) :: N,mu,nu
        real   (kind=rglu) :: med,trConst


        trConst=1/(BohrRadius*HartreeEnergy) !real coefficient
        !write(*,*) trConst
        !trConst=0.0694589_rglu               !old value for compatibility
        !write(*,*) trConst
        !stop

        N=mol%nAtoms; mol%G=0
        select case (parametrizationbd%coulombType%get())

            case ('ohno-klopman')
                !$omp parallel default(shared) private(mu,nu,med)
                !$omp do
                do mu = 1,N
                    do nu = mu,N
                        med=( mol%atm(mu)%gamma+mol%atm(nu)%gamma )/real(2,rglu)
                        mol%G(mu,nu)=med/sqrt(1+( trConst*med*mol%atmDist(mu,nu) )**2)
                        mol%G(nu,mu)=mol%G(mu,nu)
                    enddo
                enddo
                !$omp end parallel

            case ('mataga-nishimoto')
                !$omp parallel default(shared) private(mu,nu,med)
                !$omp do
                do mu = 1,N
                    do nu = mu,N
                        med=( mol%atm(mu)%gamma+mol%atm(nu)%gamma )/real(2,rglu)
                        mol%G(mu,nu)=med/(1+trConst*med*mol%atmDist(mu,nu))
                        mol%G(nu,mu)=mol%G(mu,nu)
                    enddo
                enddo
                !$omp end parallel

            case ('hubbard')
                do mu = 1,N
                    mol%G(mu,mu)=mol%atm(mu)%gamma
                enddo

        end select

        return
        end subroutine gamma

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

        subroutine geometryDisplacement
        implicit none

        real   (kind=rglu), allocatable :: X(:),Y(:),Z(:)
        real   (kind=rglu)              :: dRange,shift,sum
        integer(kind=iglu)              :: N,k


        N=mol%nAtoms
        void=glControlMemory(int( 3*rglu*N ,kind=i8kind),'tmp. Geometry displacement')
        allocate (X(N),Y(N),Z(N))

        X=mol%atm(:)%coords(1)
        Y=mol%atm(:)%coords(2)
        Z=mol%atm(:)%coords(3)

        dRange=geometrybd%randomDisplacement(1)
        !$omp parallel default(shared) private(k,shift)
        !$omp do
        do k = 1,N
            shift=rangen(-dRange,dRange); X(k)=X(k)+shift
            shift=rangen(-dRange,dRange); Y(k)=Y(k)+shift
            shift=rangen(-dRange,dRange); Z(k)=Z(k)+shift
        enddo
        !$omp end parallel

        sum=0
        !$omp parallel default(shared) private(k) reduction(+:sum)
        !$omp do
        do k = 1,N
            sum=sum+(X(k)-mol%atm(k)%coords(1))**2
            sum=sum+(Y(k)-mol%atm(k)%coords(2))**2
            sum=sum+(Z(k)-mol%atm(k)%coords(3))**2
        enddo
        !$omp end parallel
        geometrybd%randomDisplacement(2)=sqrt(sum/3*N)

        deallocate (X,Y,Z)

        void=glControlMemory(int( sizeof(X)+sizeof(Y)+sizeof(Z) ,kind=i8kind),'tmp. Geometry displacement','free')

        return
        end subroutine geometryDisplacement

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

        subroutine checkUniqueness
        implicit none

        integer(kind=iglu), allocatable :: connect(:)
        logical(kind=lglu), allocatable :: uniqueAtom(:)
        integer(kind=iglu)              :: N,M,Nbcc,Nbch,Na,i,j,k,c


        N=mol%nAtoms; M=mol%nBonds

        void=glControlMemory(int( iglu*N+lglu*N ,kind=i8kind),'tmp. Check uniqueness')
        allocate (connect(N),uniqueAtom(N)); uniqueAtom=true

        Na=N
        do i = 1,N
            connect(i)=mol%atm(i)%nels
        enddo

        do i = 1,N-1
            do j = i+1,N

                c=0
                do k = 1,3
                    if ( abs( mol%atm(i)%coords(k)-mol%atm(j)%coords(k) ) .LT.geometrybd%symmetryTolerance) c=c+1
                enddo

                ! find and count unique atoms
                if (c.EQ.3) then
                    Na=Na-1; uniqueAtom(j)=false
                endif
            enddo
        enddo
        mol%uniqueAtoms=Na

        Nbcc=M
        do k = 1,M
            i=mol%bnd(k)%atoms(1)
            j=mol%bnd(k)%atoms(2)

            ! find and count unique bonds
            if ((.NOT.uniqueAtom(i)).OR.(.NOT.uniqueAtom(j))) Nbcc=Nbcc-1

            ! find and count unique connectivity for sigma-core
            if (uniqueAtom(i).AND.uniqueAtom(j)) then
                connect(i)=connect(i)+1; connect(j)=connect(j)+1
            endif
        enddo
        mol%uniqueBonds=Nbcc

        connect=4-connect

        Nbch=0
        do i = 1,N
            if (uniqueAtom(i)) Nbch=Nbch+connect(i)
        enddo
        Nbch=Nbch-(N-Na)

        mol%chBonds=Nbch

        deAllocate (connect,uniqueAtom)
        void=glControlMemory(int( sizeof(connect)+sizeof(uniqueAtom) ,kind=i8kind),'tmp. Check uniqueness','free')

        return
        end subroutine checkUniqueness

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

        subroutine centrateCoordinates
        implicit none

        real(kind=rglu) :: sumCoords(3)


        sumCoords(1)=sum(mol%atm(:)%coords(1))
        sumCoords(2)=sum(mol%atm(:)%coords(2))
        sumCoords(3)=sum(mol%atm(:)%coords(3))
        sumCoords=sumCoords/mol%nAtoms

        mol%atm(:)%coords(1)=mol%atm(:)%coords(1)-sumCoords(1)
        mol%atm(:)%coords(2)=mol%atm(:)%coords(2)-sumCoords(2)
        mol%atm(:)%coords(3)=mol%atm(:)%coords(3)-sumCoords(3)

        return
        end subroutine centrateCoordinates

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    end subroutine readMoleculeInformation