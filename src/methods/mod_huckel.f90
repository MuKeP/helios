    module huckel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,      only: iglu,rglu,lglu,void,i8kind,glControlMemory
    use hdb,       only: mol,ou,ouWidth,dipoleToDeby,throughbd,singleSession
    use math,      only: tred4
    use printmod,  only: prEigenProblem,prMatrix,prStrByVal
    use txtparser, only: tpAdjustc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: huVersion='1.200'
    character (len=*), parameter :: huDate   ='2019.07.27'
    character (len=*), parameter :: huAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real   (kind=rglu), allocatable :: V(:,:),E(:),D(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu)              :: N,Nocc
    real(kind=rglu)                 :: Eel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private
    public :: setHuckelParameters,energyHuckel,getHuckRDMElement,getHuckelResult,finalizeHuck
    public :: getHuckelCoulsonPolarizability,printHuckelSolution

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
    ! real(kind=rglu)                :: Bx,Cx,Dx,Ex,Fx
    ! integer(kind=iglu)             :: k


    ret=0
    do i = 1,Nocc
        Ax=V(mu,i)*V(nu,i)
        do j = Nocc+1,N
            ret=ret+Ax*V(mu,j)*V(nu,j)/(E(i)-E(j))
        enddo
    enddo
    ret=4*ret

    ! write(*,'(1X,"pi(",i2.2,",",i2.2,")=",F8.5)') mu,nu,ret
    ! write(*,*)
    !
    ! k=0
    ! Fx = 0
    ! Ex = 5
    ! do i = 1,Nocc
    !     Ax=V(mu,i)*V(nu,i)
    !     do j = Nocc+1,N
    !         Cx = 4*Ax*V(mu,j)*V(nu,j)
    !         Dx = (E(i)-E(j))
    !         Bx = (100*Cx/Dx)/ret
    !         if (abs(Bx).GT.Ex) then
    !             Fx = Fx + Bx
    !             k=k+1
    !             write(*,'(1X,A,1X,i2,1X,i2,F12.6," (",F8.2,"%)",1X,E14.6,1X,F12.6)') 'Contribution', i, j, Cx/Dx, Bx, Cx, Dx
    !         endif
    !     enddo
    ! enddo
    ! write(*,'(1X,"Total contributions over ",F4.1,"%:",1X,i2)') Ex,k
    ! write(*,'(1X,A,1X,F7.2)') 'Contribution from all greater than threshold:',Fx

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

    subroutine printHuckelSolution
    implicit none

    integer(kind=iglu) :: c,paccur
    real   (kind=rglu) :: dx,dy,dz,dmod


    paccur=5

    write (ou,'(/A)') tpAdjustc(' Huckel results ',ouWidth,'=')
    call prEigenProblem(V,E,ou,'Huckel solution','^.0000',maxwidth=ouWidth)
    write (ou,99) Eel,E(Nocc+1)-E(Nocc)
    call prMatrix(2*D,ou,'Huckel RDM1', '^.0000',maxwidth=ouWidth)

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

    write (ou,'(/A)') tpAdjustc(' End of Huckel results ',ouWidth,'=')

 99 format(4X,'Electron energy:',1X,F16.8,1X,'eV'/&
           4X,'HOMO-LUMO gap:  ',5X,F8.4,5X,'eV')
100 format(/2X,A,2X,A,3X,A)
101 format(2X,i4,2X,F<3+paccur>.<paccur>,2X,F<4+paccur>.<paccur>)
102 format(<17+2*paccur>('_')//3(<10+paccur>X,A,' =',F<4+paccur>.<paccur>,1X,'D'/),&
                                  <8+paccur>X,A,' =',F<4+paccur>.<paccur>,1X,'D'/)

    return
    end subroutine printHuckelSolution

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

                case ('deallocate')
                    deallocate (V,E,D)
                    void=glControlMemory(int( sizeof(V)+sizeof(E)+sizeof(D) ,kind=i8kind),'Huckel','free')

            end select
    end select

    return
    end subroutine controlMemoryHuck

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!     subroutine huckDensity2

!     use math,     only: LagrangeDerivative
!     use printmod, only: prMatrix
!     use hdb,      only: perturbate,fieldbd,densitybd

!     implicit none

!     integer(kind=iglu) :: k,N,fatom,satom,field,sta,sto,derivPoints
!     real(kind=rglu)    :: derivStep,zpEnergy,energy(5)
!     real(kind=rglu), allocatable :: D(:,:),derivEnergies(:),E(:)


!     call setHuckelParameters

!     N=mol%nAtoms
!     derivPoints=densitybd%nPoints
!     derivStep=densitybd%derivStep; sta=-derivPoints/2; sto= derivPoints/2


!     allocate (D(N,N),E(N),derivEnergies(sta:sto)); D=0; derivEnergies=0

!     call setCore(1,1,0)
!     call prMatrix(mol%connect,6,'Connectivity','^.00000',maxwidth=78)
!     call prMatrix(mol%huckelCore,6,'Huckel core','^.00000',maxwidth=78)
!     call energyHuckel(energy)
!     zpEnergy = energy(1)

!     do fatom = 1,N
!     do satom = fatom,N
!         derivEnergies(0)=zpEnergy

!         do field = sta,sto
!             if (field.EQ.0) cycle

!             call setCore(fatom,satom,field)
!             call energyHuckel(energy); derivEnergies(field)=energy(1)
!         enddo

!         write (*,'(1X,i2,1X,i2,<derivPoints>(1X,F15.10))') fatom, satom, derivEnergies

!         call getDeriv(fatom,satom)
!     enddo
!     enddo

!     call prMatrix(D,6,'Huckel numerical','^.00000',maxwidth=78)

!     call setCore(1,1,0)

!     call prMatrix(mol%connect,6,'Connectivity__','^.00000',maxwidth=78)
!     call prMatrix(mol%huckelCore,6,'Huckel core__','^.00000',maxwidth=78)

!     call energyHuckel(energy)

!     call getHuckelResult(energies=E)
!     write (*,*) E

!     do fatom = 1,N
!     do satom = 1,N
!         D(fatom,satom)=getHuckRDMElement(fatom,satom)
!     enddo
!     enddo

!     call prMatrix(D,6,'Huckel analytical','^.00000',maxwidth=78)

!     call finalizeHuck

!     return

!     contains

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!         subroutine setCore(ip1,ip2,fshift)
!         implicit none
!         integer(kind=iglu), intent(in) :: ip1,ip2,fshift


!         ! set core perturbation
!         mol%perturbation=0
!         call applyField

!         if (ip1.EQ.ip2) then
!             mol%perturbation(ip1,ip2)=mol%perturbation(ip1,ip2)+fShift*derivStep
!         else
!             mol%perturbation(ip1,ip2)=fShift*derivStep
!             mol%perturbation(ip2,ip1)=fShift*derivStep
!         endif
!         call perturbate

!         return
!         end subroutine setCore

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!         subroutine applyField
!         implicit none

!         integer(kind=iglu) :: k,mu


!         do k = 1,3
!             if (abs(fieldbd%strength(k)).GT.1D-10) then
!                 do mu = 1,mol%nAtoms
!                     mol%perturbation(mu,mu)=fieldbd%strength(k)*mol%atm(mu)%coords(k)
!                 enddo
!             endif
!         enddo

!         return
!         end subroutine applyField

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!         subroutine getDeriv(ip1,ip2)
!         implicit none
!         integer(kind=iglu), intent(in) :: ip1,ip2


!         ! get derivatives
!         D(ip1,ip2)=LagrangeDerivative(derivPoints,1,derivEnergies,derivStep)
!         if (ip1.EQ.ip2) then
!             D(ip2,ip1)=D(ip1,ip2)
!         else
!             D(ip1,ip2)=D(ip1,ip2)/2
!             D(ip2,ip1)=D(ip1,ip2)
!         endif

!         return
!         end subroutine getDeriv

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!     end subroutine huckDensity2

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!     subroutine huckDensity

!     use math,     only: LagrangeDerivative
!     use printmod, only: prMatrix
!     use hdb,      only: perturbate,fieldbd,densitybd

!     implicit none

!     integer(kind=iglu) :: k,N,fatom,satom,field,sta,sto,derivPoints
!     real(kind=rglu)    :: derivStep,zpEnergy,energy(5)
!     real(kind=rglu), allocatable :: D(:,:),derivEnergies(:),E(:),H(:,:),V(:,:)


!     N=mol%nAtoms
!     derivPoints=densitybd%nPoints
!     derivStep=densitybd%derivStep; sta=-derivPoints/2; sto= derivPoints/2


!     allocate (D(N,N),E(N),derivEnergies(sta:sto)); D=0; derivEnergies=0
!     allocate (H(N,N),V(N,N)); H=0; V=0

!     call setCore(1,1,0)

!     zpEnergy = huckEnergy()
!     write (*,*) E

!     do fatom = 1,N
!     do satom = fatom,N
!         derivEnergies(0)=zpEnergy

!         do field = sta,sto
!             if (field.EQ.0) cycle

!             call setCore(fatom,satom,field)
!             derivEnergies(field)=huckEnergy()
!         enddo

!         write (*,'(1X,i2,1X,i2,<derivPoints>(1X,F15.10))') fatom, satom, derivEnergies

!         call getDeriv(fatom,satom)
!     enddo
!     enddo

!     call prMatrix(D,6,'Huckel numerical','^.00000',maxwidth=78)

!     return

!     call setCore(1,1,0)

!     call energyHuckel(energy)
!     do fatom = 1,N
!     do satom = 1,N
!         D(fatom,satom)=getHuckRDMElement(fatom,satom)
!     enddo
!     enddo

!     call prMatrix(D,6,'Huckel analytical','^.00000',maxwidth=78)

!     call finalizeHuck

!     return

!     contains

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!         subroutine setCore(ip1,ip2,fshift)
!         implicit none
!         integer(kind=iglu), intent(in) :: ip1,ip2,fshift


!         H=mol%connect
!         call applyField

!         if (ip1.EQ.ip2) then
!             H(ip1,ip2)=H(ip1,ip2)+fShift*derivStep
!         else
!             H(ip1,ip2)=H(ip1,ip2)+fShift*derivStep
!             H(ip2,ip1)=H(ip2,ip1)+fShift*derivStep
!         endif

!         return
!         end subroutine setCore

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!         subroutine applyField
!         implicit none

!         integer(kind=iglu) :: k,mu


!         do k = 1,3
!             do mu = 1,N
!                 H(mu,mu)=H(mu,mu)+fieldbd%strength(k)*mol%atm(mu)%coords(k)
!             enddo
!         enddo

!         return
!         end subroutine applyField

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!         subroutine getDeriv(ip1,ip2)
!         implicit none
!         integer(kind=iglu), intent(in) :: ip1,ip2


!         ! get derivatives
!         D(ip1,ip2)=LagrangeDerivative(derivPoints,1,derivEnergies,derivStep)
!         if (ip1.EQ.ip2) then
!             D(ip2,ip1)=D(ip1,ip2)
!         else
!             D(ip1,ip2)=D(ip1,ip2)/2
!             D(ip2,ip1)=D(ip1,ip2)
!         endif

!         return
!         end subroutine getDeriv

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!         real(kind=rglu) function huckEnergy() result(ret)
!         implicit none

!         integer(kind=iglu) :: mu


!         call tred4(H,V,E,N,1.e-100_rglu,1.e-300_rglu)

!         ret=0
!         do mu = 1,N/2
!             ret=ret+2*E(mu)
!         enddo

!         return
!         end function huckEnergy

! !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!     end subroutine huckDensity

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module huckel
