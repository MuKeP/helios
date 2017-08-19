	module huckel

	use glob    , only: iglu,rglu,rspu
	use math    , only: gltred4
	use hdb     , only: mol,ou,ouWidth
	use printmod, only: prEigenProblem

	real   (kind=rglu), allocatable :: V(:,:),E(:),D(:,:)
	integer(kind=iglu)              :: N,Nocc


	private

	public :: energyHuckel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine energyHuckel(energy)
	implicit none

	integer(kind=iglu) :: mu,nu
	real   (kind=rglu) :: energy(5),sum,Enuc


	N=mol%nAtoms; Nocc=mol%nEls/2
	allocate(V(N,N),E(N)); V=0; E=0
	call gltred4(mol%core,V,E,N,real(1.d-100,rglu),real(1.d-300,rglu))

	call prEigenProblem(V,E,ou,'Huckel solution','^.0000',maxwidth=ouWidth)

	energy=0; sum=0
	do mu = 1,Nocc
		sum=sum+2*E(mu)
	enddo

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

	energy(1)=sum+Enuc

	deallocate (V,E)

	return
	end subroutine energyHuckel

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module huckel