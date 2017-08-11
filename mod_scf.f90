	module scf

	use hdb , only: iglu,rglu,rspu,void,true,false,uch,uchGet
	use hdb , only: mol,scfbd
	use hdb , only: tpFill
	use math, only: tred4

	integer(kind=iglu) :: N,Nocc
	real   (kind=rglu), allocatable :: F(:,:),V(:,:),E(:),D(:,:)
	real   (kind=rglu), allocatable :: Fmin(:,:),Comut(:,:),Refl(:,:)

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine setSCFParameters

		implicit none


		N=mol%nAtoms; Nocc=mol%nEls/2

		allocate (F(N,N),V(N,N),E(N),D(N,N)); F=0; V=0; E=0; D=0
		allocate (Fmin(N,N),Comut(N,N),Refl(N,N)); Fmin=0; Comut=0; Refl=0

		return
		end subroutine setSCFParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine initSCF(vecs)

		implicit none

		real(kind=rglu), optional :: vecs(:,:)


		if (present(vecs)) then
			call scfGuess(vecs)
		else
			call scfGuess
		endif

		call prepareFockian

		return
		end subroutine initSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine iterationSCF(iteration,epsilon,accuracy)

		implicit none

		integer(kind=iglu), intent(in)  :: iteration
		real   (kind=rglu), intent(in)  :: epsilon
		real   (kind=rglu), intent(out) :: accuracy(5)

		real   (kind=rglu)              :: eps
		integer(kind=iglu)              :: mu,nu



		call prepareFmin; eps=maxval(abs(Fmin))
		do mu = 1,N
			do nu = 1,N
				D(mu,nu)=D(mu,nu)-scfbd%iterStep*Fmin(mu,nu)
			enddo
		enddo
		call prepareFockian
		call TRED4 (F,V,E,N,real(1.d-100,rspu),real(1.d-300,rspu))
		call prepareDensity(V)
		call prepareFockian

		accuracy=-1; accuracy(1)=eps

		return
		end subroutine iterationSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine energySCF(energy)

		implicit none

		real   (kind=rglu) :: energy(5)
		real   (kind=rglu) :: Eel,Enuc
		integer(kind=iglu) :: mu,nu


		Eel=0
		do mu = 1,N
			do nu = 1,N
				Eel=Eel+D(mu,nu)*(mol%core(mu,nu)+F(mu,nu))
			enddo
		enddo

		Enuc=0
		do mu = 1,N
			do nu = 1,N
				if (mu.EQ.nu) cycle
				Enuc=Enuc+mol%G(mu,nu)/2
	 		enddo
		enddo

		energy=0; energy(1)=Eel+Enuc

		return
		end subroutine energySCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine getSCFResult(fockian,vectors,energies)

		implicit none

		real   (kind=rglu), optional :: fockian(:,:),vectors(:,:),energies(:)
		integer(kind=iglu) :: i,j


		if (present(fockian)) then
			if (UBound(fockian,1).NE.N) stop 'Incorrect size of Fockian matrix.'

			do i = 1,N
				do j = 1,N
					fockian(i,j)=F(i,j)
				enddo
			enddo
		endif

		if (present(vectors)) then
			if (UBound(vectors,1).NE.N) stop 'Incorrect size of Fockian eigenvectors.'

			do i = 1,N
				do j = 1,N
					vectors(i,j)=V(i,j)
				enddo
			enddo
		endif

		if (present(energies)) then
			if (UBound(energies,1).NE.N) stop 'Incorrect size of Fockian eigenvalues.'

			do i = 1,N
				energies(i)=E(i)
			enddo
		endif

		return
		end subroutine getSCFResult

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine finalizeSCF

		implicit none


		deallocate (F,V,E,D,Fmin,Comut,Refl)

		return
		end subroutine finalizeSCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine scfGuess(vecs)

		implicit none

		character (len=10) :: guess
		integer(kind=iglu) :: mu
		real   (kind=rglu), optional :: vecs(N,N)


		guess=tpFill(guess); guess=uchGet(scfbd%guess)
		if (present(vecs)) then
			guess=tpFill(guess); guess='defined'
		endif

		select case (guess)
			case ('huckel')
				call tred4(mol%connect,V,E,N,real(1.d-100,rspu),real(1.d-300,rspu))
				call prepareDensity(V)

			case ('unitmatrix')
				D=0
				do mu = 1,N
					D(mu,mu)=1
				enddo

			case ('defined')
				call prepareDensity(vecs)					

		end select

		return
		end subroutine scfGuess

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prepareDensity(vec)

		implicit none

		real(kind=rglu)    :: vec(N,N)
		real(kind=rglu)    :: sum
		integer(kind=iglu) :: i,mu,nu


		do mu = 1,N
			do nu = 1,N

				sum=0
				do i = 1,Nocc
					sum=sum+vec(mu,i)*vec(nu,i)
				enddo
				D(mu,nu)=sum
			enddo
		enddo

		return
		end subroutine prepareDensity

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prepareFockian

		implicit none

		real(kind=rglu)    :: sum
		integer(kind=iglu) :: mu,nu,i


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

		return
		end subroutine prepareFockian

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prepareFmin

		implicit none

		real   (kind=rglu) :: sum1,sum2
		integer(kind=iglu) :: i,j,mu


		do i = 1,N
			do j = 1,N
				Refl(i,j)=2*D(i,j)
			enddo
			Refl(i,i)=Refl(i,i)-1
		enddo 

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

		do i = 1,N
			do j = 1,N
				sum1=0
				do mu = 1,N
					sum1=sum1+Comut(i,mu)*Refl(mu,j)
				enddo 
				Fmin(i,j)=sum1
			enddo 
		enddo 

		return
		end subroutine prepareFmin

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module scf