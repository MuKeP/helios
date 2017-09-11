	module excitedStates

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	use glob, only: iglu,rglu,uchGet,true
	use math, only: gltred4
	use hdb , only: mol,statesbd

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: esVersion='1.000'
	character (len=*), parameter :: esDate   ='2017.08.28'
	character (len=*), parameter :: esAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu), allocatable :: excSet(:,:),confs(:)
	real   (kind=rglu), allocatable :: hfV(:,:),hfE(:),HH(:,:),Vectors(:,:),Values(:),coefs(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu) :: N,Nocc,Ne

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private
	public :: getExcitationEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu) function getExcitationEnergy(V,E,state,allStates) result(ret)
	implicit none

	real   (kind=rglu), intent(in)            :: V(:,:),E(:)
	integer(kind=iglu), intent(in) , optional :: state
	real   (kind=rglu), intent(out), optional :: allStates(:)
	real   (kind=rglu)                        :: energy(1:statesbd%nStates)

	character (len=3)  :: method
	integer(kind=iglu) :: i,j,a,b,k,l
	real   (kind=rglu) :: Ax,Bx


	N=mol%nAtoms; Nocc=mol%nEls/2; Ne=Nocc*(N-Nocc)
	call controlMemory('general','allocate')

	method=uchGet(statesbd%hftype)

	hfV=V; hfE=E

	k=0
	do i = 1,Nocc; do j = Nocc+1,N
		k=k+1; excSet(1,k)=i; excSet(2,k)=j
	enddo; enddo

	!$omp parallel default(shared) private(k,i,a,l,j,b,Ax,Bx)
	!$omp do
	do k = 1,Ne
		i=excSet(1,k); a=excSet(2,k)
		do l = 1,Ne
			j=excSet(1,l); b=excSet(2,l)

			Ax=spat_hf_int(a,i,j,b)

			if ((i.EQ.j).AND.(a.EQ.b)) Ax=Ax+hfE(a)-hfE(i)

			select case (method)
				case ('cis'); Bx=0
				case ('rpa'); Bx=spat_hf_int(a,i,b,j)
			end select

			HH(k   ,l)=Ax; HH(k+Ne,l+Ne)=Ax
			HH(k+Ne,l)=Bx; HH(k   ,l+Ne)=Bx
		enddo
	enddo
	!$omp end parallel

	call gltred4(HH,Vectors,Values,2*Ne,real(1.d-100,rglu),real(1.d-300,rglu))

	select case (method)
		case ('cis')
			do k = 1,statesbd%nStates
				energy(k)=Values(2*k-1)
			enddo

		case ('rpa')
			do k = 1,statesbd%nStates
				energy(k)=Values(k)
			enddo

	end select

	call controlMemory('general','deallocate')

	if (present(allStates)) then
		allStates=energy(1:UBound(allStates,1))
	endif

	ret=energy(1)
	if (present(state)) then
		ret=energy(state)
	endif

	return
	end function getExcitationEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine controlMemory(section,action)
	implicit none

	character  (len=*) :: section,action
	integer(kind=iglu) :: err


	select case (section)
		case ('general')
			select case (action)
				case ('allocate')
					allocate (HH(2*Ne,2*Ne),Vectors(2*Ne,2*Ne),Values(2*Ne),coefs(2*Ne),confs(2*Ne))
					HH=0; Vectors=0; Values=0; coefs=0; confs=0

					allocate (hfV(N,N),hfE(N),excSet(2,Ne))
					hfV=0; hfE=0; excSet=0

				case ('deallocate')
					deallocate (HH,Vectors,Values,coefs,confs,hfV,hfE,excSet, stat=err)

			end select

		!case ()

	end select

	return
	end subroutine controlMemory

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

	end module excitedStates