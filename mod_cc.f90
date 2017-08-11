	module coupledCluster

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	use glob     , only: uch,uchGet,uchSet,iglu,rglu,lglu,rspu,true,false
	use txtParser, only: operator(.in.)
	use printmod , only: prMatrix
	use scf      , only: setSCFParameters,initSCF,iterationSCF,getSCFResult,energySCF
	use hdb      , only: mol,ccbd,diisbd,cuebd,systembd,scfbd,ou,gluCompare

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: ccVersion='1.000'
	character (len=*), parameter :: ccDate   ='2017.08.11'
	character (len=*), parameter :: ccAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=iglu), allocatable :: excSet(:,:)

	integer(kind=1)   , allocatable :: V(:,:)
	real   (kind=rglu), allocatable :: hV(:,:),hVs(:,:)
	integer(kind=iglu), allocatable :: cueIndex(:,:),iapairs(:)
	real   (kind=rglu), allocatable :: density(:,:),cueDistance(:,:),G(:,:)

	! General CC arrays
	real   (kind=rglu), allocatable :: F(:,:),R(:,:,:,:)
	real   (kind=rglu), allocatable :: t1(:,:),t2(:,:,:,:),t3(:,:,:,:,:,:)
	real   (kind=rglu), allocatable :: d1(:,:),d2(:,:,:,:),d3(:,:,:,:,:,:)
	integer(kind=iglu), allocatable :: Fnz(:,:)

	! DIIS storage
	real   (kind=rglu), allocatable :: st1(:,:),st2(:,:),st3(:,:)
	real   (kind=rglu), allocatable :: sd1(:,:),sd2(:,:),sd3(:,:)

	! DIIS arrays
	real   (kind=rglu), allocatable :: diisVectors(:,:),diisValues(:)
	real   (kind=rglu), allocatable :: diisMatrix(:,:),diisCoefficients(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type(uch)          :: umethod
	integer(kind=iglu) :: N,Nel,No,Nocc,Ne,NFnz,Nth,Nd
	real   (kind=rglu) :: accuracy(5),refeEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	public

	private :: umethod,accuracy,Nd
	private :: diisVectors,diisValues,diisMatrix,diisCoefficients
	private :: st1,st2,st3,sd1,sd2,sd3
	private :: V,hV,hVs,cueIndex,density,G
	private :: guessCC,referenceEnergy,pushCCVectors,newCCVectors,prepareDIIS 
	private :: diisDotProduct,prepareDensity,prepareFock,prepareFockCUE,stepCC

	private :: uch,uchGet,uchSet,iglu,rglu,lglu,rspu,true,false,ou,gluCompare
	private :: operator(.in.),prMatrix,mol,ccbd,diisbd,cuebd,systembd,scfbd
	private :: setSCFParameters,initSCF,iterationSCF,getSCFResult,energySCF

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine setCCParameters(method)
	implicit none

	character (len=*), intent(in) :: method
	integer(kind=iglu)            :: i,j,k,l


	select case (method)
		case ('cue-ccs','cue-ccsd','cue-ccsdt')
			do
				if (cuebd%sparse .AND. (method.EQ.'cue-ccsd')) then
					umethod=uchSet('spare-'//method)
					exit
				endif

				if (ccbd%forceSpin) then
					umethod=uchSet('spin-'//method)
					exit
				endif

				if (method.EQ.'cue-ccsdt') then
					umethod=uchSet('spin-'//method)
					exit
				endif

				umethod=uchSet(method)
				exit
			enddo

		case ('u-ccd','u-ccsd','u-ccsdt','r-ccd','r-ccsd','r-ccsdt','r-ccsd(t)')
			do
				if (ccbd%forceSpin) then
					umethod=uchSet('spin-'//method)
					exit
				endif

				if (method .in. ['u-ccsdt','r-ccsdt','r-ccsd(t)'] ) then
					umethod=uchSet('spin-'//method)
					exit
				endif

				umethod=uchSet(method)
				exit
			enddo

		case default
			stop 'Unknown method'

	end select

	N=mol%nAtoms; Nel=mol%nEls; Nth=systembd%nNodes; Nd=diisbd%steps

	do i = 1,N
	do j = 1,N
		mol%core(i,j)=mol%holdCore(i,j)
	enddo
	enddo

	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			No=2*N; Nocc=Nel/2
			allocate (cueIndex(2,No),V(0:N,No),iapairs(No))
			cueIndex=0; V=0; iapairs=0
			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)
				l=mol%orb(k)%ova

				if (mol%orb(k)%nels.EQ.2) then
					V(i,2*k-1)=int(1,kind=1); V(i,2*k)=int(1,kind=1)

					cueIndex(1,2*k-1)=i; cueIndex(1,2*k)=i
					cueIndex(2,2*k-1)=j; cueIndex(2,2*k)=j

					iapairs(2*k-1)=2*l-1; iapairs(2*k)=2*l
					cycle
				endif

				V(i,2*k-1)=1; V(i,2*l-1)=-1
				V(i,2*k  )=1; V(i,2*l  )=-1
				V(j,2*k-1)=1; V(j,2*l-1)= 1
				V(j,2*k  )=1; V(j,2*l  )= 1

				cueIndex(1,2*k-1)=i; cueIndex(1,2*l-1)=i
				cueIndex(1,2*k  )=i; cueIndex(1,2*l  )=i
				cueIndex(2,2*k-1)=j; cueIndex(2,2*l-1)=j
				cueIndex(2,2*k  )=j; cueIndex(2,2*l  )=j

				iapairs(2*k-1)=2*l-1; iapairs(2*k)=2*l
				iapairs(2*l-1)=2*k-1; iapairs(2*l)=2*k
			enddo

			allocate (G(N,N),density(N,N)); G=mol%G; density=0

			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)

				if (mol%orb(k)%nels.EQ.2) then
					density(i,i)=1
					cycle
				endif

				density(i,j)=real(0.5,rglu); density(j,i)=real(0.5,rglu)
				density(i,i)=real(0.5,rglu); density(j,j)=real(0.5,rglu)
			enddo
			
			allocate ( F(No,No) ); F=0

			allocate (cueDistance(No,No)); cueDistance=0

			do i = 1,N
				do j = i,N
					cueDistance(2*i-1,2*j-1)=mol%cuedist(i,j)
					cueDistance(2*i  ,2*j-1)=mol%cuedist(i,j)
					cueDistance(2*i-1,2*j  )=mol%cuedist(i,j)
					cueDistance(2*i  ,2*j  )=mol%cuedist(i,j)

					cueDistance(2*j-1,2*i-1)=mol%cuedist(i,j)
					cueDistance(2*j  ,2*i-1)=mol%cuedist(i,j)
					cueDistance(2*j-1,2*i  )=mol%cuedist(i,j)
					cueDistance(2*j  ,2*i  )=mol%cuedist(i,j)
				enddo
			enddo

			call prepareSparseIndexInformation
			deallocate (cueDistance)

		case ('cue-ccs','cue-ccsd')
			No=N; Nocc=Nel/2; Ne=(N-Nocc)*Nocc

			allocate (cueIndex(2,N),V(N,N),iapairs(N))
			cueIndex=0; V=0; iapairs=0

			! index information
			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)
				l=mol%orb(k)%ova

				if (mol%orb(k)%nels.EQ.2) then
					V(i,k)=int(1,kind=1)
					cueIndex(1,k)=i; cueIndex(2,k)=j
					iapairs(k)=l
					cycle
				endif

				V(i,k)=1; V(i,l)=-1
				V(j,k)=1; V(j,l)= 1

				cueIndex(1,k)=i; cueIndex(1,l)=i
				cueIndex(2,k)=j; cueIndex(2,l)=j

				iapairs(k)=l   ; iapairs(l)=k
			enddo

			allocate ( excSet(Ne,2) )
			k=0
			do i = 1,Nocc
			do j = Nocc+1,N
				k=k+1; excSet(k,1)=i; excSet(k,2)=j
			enddo
			enddo

			allocate (G(N,N),density(N,N)); G=mol%G; density=0

			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)

				if (mol%orb(k)%nels.EQ.2) then
					density(i,i)=1
					cycle
				endif

				density(i,i)=real(0.5,rglu); density(i,j)=real(0.5,rglu)
				density(j,i)=real(0.5,rglu); density(j,j)=real(0.5,rglu)
			enddo

			! correction to the calculation of two-electron integrals
			! in case of the atoms with unshared electron pairs.
			do k = 1,N
			do l = k,N
				if ((mol%atm(k)%nels.EQ.2).AND.(mol%atm(l)%nels.EQ.2)) then
					G(k,l)=4*G(k,l)
					G(l,k)=G(k,l)
					cycle
				endif
				
				if ((mol%atm(k)%nels.EQ.2).OR. (mol%atm(l)%nels.EQ.2)) then
					G(k,l)=2*G(k,l)
					G(l,k)=G(k,l)
				endif
			enddo
			enddo

			allocate (cueDistance(N,N)); cueDistance=0

			do i = 1,N
				do j = i,N
					cueDistance(i,j)=mol%cuedist(i,j)
					cueDistance(j,i)=mol%cuedist(j,i)
				enddo
			enddo

			allocate ( R(N,N,N,N),F(N,N) )

			!$omp parallel default(shared) private(i,j,k,l)
			!$omp do
			do i = 1,N
			do j = 1,N
			do k = 1,N
			do l = 1,N
				R(i,j,k,l)=spat_cue_int(i,j,k,l)/4
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			select case (uchGet(umethod))
				case ('cue-ccs')
					allocate (t1(Nocc,Nocc+1:N),&
							  d1(Nocc,Nocc+1:N) )

					t1=0; d1=0
				case ('cue-ccsd')
					allocate (t1(Nocc,Nocc+1:N),t2(Nocc,Nocc,Nocc+1:N,Nocc+1:N),&
							  d1(Nocc,Nocc+1:N),d2(Nocc,Nocc,Nocc+1:N,Nocc+1:N) )

					t1=0; t2=0; d1=0; d2=0

			end select

		case ('spin-cue-ccs','spin-cue-ccsd','spin-cue-ccsdt')
			No=2*N; Nocc=Nel/2; Ne=(N-Nocc)*Nocc
			allocate (cueIndex(2,No),V(0:N,No),iapairs(No))
			cueIndex=0; V=0; iapairs=0
			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)
				l=mol%orb(k)%ova

				if (mol%orb(k)%nels.EQ.2) then
					V(i,2*k-1)=int(1,kind=1); V(i,2*k)=int(1,kind=1)

					cueIndex(1,2*k-1)=i; cueIndex(1,2*k)=i
					cueIndex(2,2*k-1)=j; cueIndex(2,2*k)=j

					iapairs(2*k-1)=2*l-1; iapairs(2*k)=2*l
					cycle
				endif

				V(i,2*k-1)=1; V(i,2*l-1)=-1
				V(i,2*k  )=1; V(i,2*l  )=-1
				V(j,2*k-1)=1; V(j,2*l-1)= 1
				V(j,2*k  )=1; V(j,2*l  )= 1

				cueIndex(1,2*k-1)=i; cueIndex(1,2*l-1)=i
				cueIndex(1,2*k  )=i; cueIndex(1,2*l  )=i
				cueIndex(2,2*k-1)=j; cueIndex(2,2*l-1)=j
				cueIndex(2,2*k  )=j; cueIndex(2,2*l  )=j

				iapairs(2*k-1)=2*l-1; iapairs(2*k)=2*l
				iapairs(2*l-1)=2*k-1; iapairs(2*l)=2*k
			enddo

			allocate (G(N,N),density(N,N)); G=mol%G; density=0

			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)

				if (mol%orb(k)%nels.EQ.2) then
					density(i,i)=1
					cycle
				endif

				density(i,j)=real(0.5,rglu); density(j,i)=real(0.5,rglu)
				density(i,i)=real(0.5,rglu); density(j,j)=real(0.5,rglu)
			enddo

			do k = 1,N
			do l = k,N
				if ((mol%atm(k)%nels.GT.1).AND.(mol%atm(l)%nels.GT.1)) then
					G(k,l)=4*G(k,l)
					G(l,k)=G(k,l)
					cycle
				endif
				
				if ((mol%atm(k)%nels.GT.1).OR. (mol%atm(l)%nels.GT.1)) then
					G(k,l)=2*G(k,l)
					G(l,k)=G(k,l)
				endif
			enddo
			enddo

			allocate (cueDistance(No,No)); cueDistance=0

			do i = 1,N
				do j = i,N
					cueDistance(2*i-1,2*j-1)=mol%cuedist(i,j)
					cueDistance(2*i  ,2*j-1)=mol%cuedist(i,j)
					cueDistance(2*i-1,2*j  )=mol%cuedist(i,j)
					cueDistance(2*i  ,2*j  )=mol%cuedist(i,j)

					cueDistance(2*j-1,2*i-1)=mol%cuedist(i,j)
					cueDistance(2*j  ,2*i-1)=mol%cuedist(i,j)
					cueDistance(2*j-1,2*i  )=mol%cuedist(i,j)
					cueDistance(2*j  ,2*i  )=mol%cuedist(i,j)
				enddo
			enddo

			allocate ( R(No,No,No,No), F(No,No) )

			!$omp parallel default(shared) private(i,j,k,l)
			!$omp do
			do i = 1,No
			do j = 1,No
			do k = 1,No
			do l = 1,No
				R(i,j,k,l)=spin_cue_int(i,j,k,l)/4
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			select case (uchGet(umethod))
				case ('spin-cue-ccs')
					allocate (t1(Nel,Nel+1:No),&
							  d1(Nel,Nel+1:No) )

					t1=0; d1=0
				case ('spin-cue-ccsd')
					allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No),&
							  d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No) )

					t1=0; t2=0; d1=0; d2=0

				case ('spin-cue-ccsdt')
					allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No),t3(Nel,Nel,Nel,Nel+1:No,Nel+1:No,Nel+1:No),&
							  d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No),d3(Nel,Nel,Nel,Nel+1:No,Nel+1:No,Nel+1:No) )

					t1=0; t2=0; t3=0; d1=0; d2=0; d3=0

			end select

		case ('u-ccd','u-ccsd')
			No=N; Nocc=Nel/2; Ne=(N-Nocc)*Nocc

			allocate (hV(N,N))
			hV=0

			allocate ( excSet(Ne,2) )
			k=0
			do i = 1,Nocc
			do j = Nocc+1,N
				k=k+1; excSet(k,1)=i; excSet(k,2)=j
			enddo
			enddo

			! one- and two-electron integrals
			allocate (G(N,N),density(N,N)); G=mol%G; density=0
			
			allocate ( R(N,N,N,N),F(N,N) )

			call setSCFParameters
			call initSCF
			call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true)
			call getSCFResult(vectors=hV)

			!$omp parallel default(shared) private(i,j,k,l)
			!$omp do
			do i = 1,N
			do j = 1,N
			do k = 1,N
			do l = 1,N
				R(i,j,k,l)=spat_hf_int(i,j,k,l)
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			select case (uchGet(umethod))
				case ('u-ccd')
					allocate (t2(Nocc,Nocc,Nocc+1:N,Nocc+1:N),&
							  d2(Nocc,Nocc,Nocc+1:N,Nocc+1:N) )

					t2=0; d2=0
				case ('u-ccsd')
					allocate (t1(Nocc,Nocc+1:N),t2(Nocc,Nocc,Nocc+1:N,Nocc+1:N),&
							  d1(Nocc,Nocc+1:N),d2(Nocc,Nocc,Nocc+1:N,Nocc+1:N) )

					t1=0; t2=0; d1=0; d2=0

			end select

		case ('r-ccd','r-ccsd')
			No=N; Nocc=Nel/2; Ne=(N-Nocc)*Nocc

			allocate (hV(N,N))
			hV=0

			allocate ( excSet(Ne,2) )
			k=0
			do i = 1,Nocc
			do j = Nocc+1,N
				k=k+1; excSet(k,1)=i; excSet(k,2)=j
			enddo
			enddo

			! one- and two-electron integrals
			allocate (G(N,N),density(N,N)); G=mol%G; density=0
			
			allocate ( R(N,N,N,N),F(N,N) )
			call setSCFParameters

			select case (uchGet(umethod))
				case ('r-ccd')
					allocate (t2(Nocc,Nocc,Nocc+1:N,Nocc+1:N),&
							  d2(Nocc,Nocc,Nocc+1:N,Nocc+1:N) )

					t2=0; d2=0
				case ('r-ccsd')
					allocate (t1(Nocc,Nocc+1:N),t2(Nocc,Nocc,Nocc+1:N,Nocc+1:N),&
							  d1(Nocc,Nocc+1:N),d2(Nocc,Nocc,Nocc+1:N,Nocc+1:N) )

					t1=0; t2=0; d1=0; d2=0

			end select

		case ('spin-u-ccd','spin-u-ccsd','spin-u-ccsdt')
			No=2*N; Nocc=Nel/2; Ne=(N-Nocc)*Nocc
			allocate (hV(N,N),hVs(N,No))

			allocate (G(N,N),density(N,N)); G=mol%G; density=0
				
			allocate ( R(No,No,No,No), F(No,No) )

			call setSCFParameters
			call initSCF
			call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true)
			call getSCFResult(vectors=hV)

			do i = 1,N
			do j = 1,N
				hVs(i,2*j-1)=hV(i,j)
				hVs(i,2*j  )=hV(i,j)
			enddo
			enddo

			!$omp parallel default(shared) private(i,j,k,l)
			!$omp do
			do i = 1,No
			do j = 1,No
			do k = 1,No
			do l = 1,No
				R(i,j,k,l)=spin_hf_int(i,j,k,l)
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			select case (uchGet(umethod))
				case ('spin-u-ccd')
					allocate (t2(Nel,Nel,Nel+1:No,Nel+1:No),&
							  d2(Nel,Nel,Nel+1:No,Nel+1:No) )

					t2=0; d2=0
				case ('spin-u-ccsd')
					allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No),&
							  d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No) )

					t1=0; t2=0; d1=0; d2=0

				case ('spin-u-ccsdt')
					allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No),t3(Nel,Nel,Nel,Nel+1:No,Nel+1:No,Nel+1:No),&
							  d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No),d3(Nel,Nel,Nel,Nel+1:No,Nel+1:No,Nel+1:No) )

					t1=0; t2=0; t3=0; d1=0; d2=0; d3=0

			end select


		case ('spin-r-ccd','spin-r-ccsd','spin-r-ccsdt','spin-r-ccsd(t)')
			No=2*N; Nocc=Nel/2; Ne=(N-Nocc)*Nocc
			allocate (hVs(N,No),hV(N,N))

			allocate (G(N,N),density(N,N)); G=mol%G; density=0
				
			allocate ( R(No,No,No,No), F(No,No) )

			call setSCFParameters

			select case (uchGet(umethod))
				case ('spin-r-ccd')
					allocate (t2(Nel,Nel,Nel+1:No,Nel+1:No),&
							  d2(Nel,Nel,Nel+1:No,Nel+1:No) )

					t2=0; d2=0
				case ('spin-r-ccsd','spin-r-ccsd(t)')
					allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No),&
							  d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No) )

					t1=0; t2=0; d1=0; d2=0

				case ('spin-r-ccsdt')
					allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No),t3(Nel,Nel,Nel,Nel+1:No,Nel+1:No,Nel+1:No),&
							  d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No),d3(Nel,Nel,Nel,Nel+1:No,Nel+1:No,Nel+1:No) )

					t1=0; t2=0; t3=0; d1=0; d2=0; d3=0

			end select

	end select

	return
	end subroutine setCCParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine initCC
	implicit none

	integer (kind=iglu) :: i,j,a,b


	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			call prepareFockCUE('spin')

			NFnz=0
			do i = 1,No; do j = 1,No
				if (abs(F(i,j)).GT.gluCompare) NFnz=NFnz+1
			enddo; enddo

			if (allocated(Fnz)) deallocate (Fnz)
			allocate (Fnz(2,NFnz)); Fnz=0

			NFnz=0
			do i = 1,No; do j = 1,No
				if (abs(F(i,j)).GT.gluCompare) then
					NFnz=NFnz+1; Fnz(1,NFnz)=i; Fnz(2,NFnz)=j
				endif
			enddo
			enddo

			call initSpareCC

		case ('cue-ccs','cue-ccsd')
			call prepareFockCUE('spatial')

			NFnz=0
			do i = 1,Nocc; do j = Nocc+1,N
				if (abs(F(i,j)).GT.gluCompare) NFnz=NFnz+1
			enddo; enddo

			if (allocated(Fnz)) deallocate (Fnz)
			allocate (Fnz(2,NFnz)); Fnz=0

			NFnz=0
			do i = 1,Nocc; 	do j = Nocc+1,N
				if (abs(F(i,j)).GT.gluCompare) then
					NFnz=NFnz+1; Fnz(1,NFnz)=i; Fnz(2,NFnz)=j
				endif
			enddo; enddo

		case ('spin-cue-ccs','spin-cue-ccsd','spin-cue-ccsdt')
			call prepareFockCUE('spin')

			NFnz=0
			do i = 1,Nel; do j = Nel+1,No
				if (abs(F(i,j)).GT.gluCompare) NFnz=NFnz+1
			enddo; enddo

			if (allocated(Fnz)) deallocate (Fnz)
			allocate (Fnz(2,NFnz)); Fnz=0

			NFnz=0
			do i = 1,Nel; do j = Nel+1,No
				if (abs(F(i,j)).GT.gluCompare) then
					NFnz=NFnz+1; Fnz(1,NFnz)=i; Fnz(2,NFnz)=j
				endif
			enddo
			enddo

		case ('u-ccd','u-ccsd')
			call initSCF(hV)
			call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true)
			call getSCFResult(vectors=hV)
			call prepareDensity(hV)
			call prepareFock('spatial')
			!deallocate (G,hV,density)

		case ('r-ccd','r-ccsd')
			call initSCF
			call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true)
			call getSCFResult(vectors=hV)
			call prepareDensity(hV)
			call prepareFock('spatial')

			!$omp parallel default(shared) private(i,j,a,b)
			!$omp do
			do i = 1,N
			do j = 1,N
			do a = 1,N
			do b = 1,N
				R(i,j,a,b)=spat_hf_int(i,j,a,b)
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

		case ('spin-u-ccd','spin-u-ccsd','spin-u-ccsdt')
			call initSCF(hV)
			call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true)
			call getSCFResult(vectors=hV)
			call prepareDensity(hV)

			call prepareFock('spin')

		case ('spin-r-ccd','spin-r-ccsd','spin-r-ccsdt','spin-r-ccsd(t)')
			call initSCF
			call iterator(iterationSCF,energySCF,scfbd%maxiters,scfbd%accuracy,true)
			call getSCFResult(vectors=hV)
			call prepareDensity(hV)

			call prepareFock('spin')

			do i = 1,N
			do j = 1,N
				hVs(i,2*j-1)=hV(i,j)
				hVs(i,2*j  )=hV(i,j)
			enddo
			enddo

			!$omp parallel default(shared) private(i,j,a,b)
			!$omp do
			do i = 1,No
			do j = 1,No
			do a = 1,No
			do b = 1,No
				R(i,j,a,b)=spin_hf_int(i,j,a,b)
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

	end select

	if (diisbd%enabled) call prepareDIIS
	call guessCC

	return
	end subroutine initCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine guessCC

	use coupledClusterSparse, only: Indexs,cIndex,ftfm,spt1=>t1,spt2=>vt

	implicit none

	integer(kind=iglu) :: i,j,a,b,mm,nn


	select case (uchGet(umethod))
		case ('u-ccd','r-ccd')
			!$omp parallel default(shared) private(mm,nn,i,a,j,b)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					t2(i,j,a,b)=R(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
					t2(j,i,b,a)=t2(i,j,a,b)
				enddo
			enddo
			!$omp end parallel

		case ('u-ccsd','r-ccsd')
			!$omp parallel default(shared) private(mm,nn,i,a,j,b)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				t1(i,a)=F(i,a)/(F(i,i)-F(a,a))
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					t2(i,j,a,b)=R(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
					t2(j,i,b,a)=t2(i,j,a,b)
				enddo
			enddo
			!$omp end parallel

		case ('spin-u-ccd','spin-r-ccd')
			!$omp parallel default(shared) private(i,a,j,b)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				do j = 1,Nel
				do b = Nel+1,No
					t2(i,j,a,b)=R(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				enddo
				enddo
			enddo
			enddo
			!$omp end parallel

		case ('spin-u-ccsdt','spin-r-ccsdt','spin-u-ccsd','spin-r-ccsd','spin-r-ccsd(t)')
			!$omp parallel default(shared) private(i,a,j,b)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				t1(i,a)=F(i,a)/(F(i,i)-F(a,a))
				do j = 1,Nel
				do b = Nel+1,No
					t2(i,j,a,b)=R(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				enddo
				enddo
			enddo
			enddo
			!$omp end parallel

		case ('spare-cue-ccsd')
			!$omp parallel default(shared) private(i,a)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				if (notFitRadius(1,i,a)) cycle
				spt1(i,a)=F(i,a)/(F(i,i)-F(a,a))
			enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(mm,nn,i,a,j,b)
			!$omp do
			do mm = 1,Ne ! [ia||jb]
				i=Indexs(mm,1)
				a=Indexs(mm,2)

				j=iapairs(a); b=iapairs(i)

				nn=cIndex(j,b)
				spt2(ftfm(mm,nn))=( spin_cue_int(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b)) )/4

				if (b.EQ.a) then
					do j = 1,Nel
						b=iapairs(j)
						nn=cIndex(j,b)
						spt2(ftfm(mm,nn))=( spin_cue_int(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b)) )/4
					enddo
				endif
			enddo
			!$omp end parallel

		case ('cue-ccs')
			!$omp parallel default(shared) private(mm,i,a)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				if (notFitRadius(1,i,a)) cycle

				t1(i,a)=F(i,a)/(F(i,i)-F(a,a))
			enddo
			!$omp end parallel

		case ('cue-ccsd')
			!$omp parallel default(shared) private(mm,i,a)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				if (notFitRadius(1,i,a)) cycle

				t1(i,a)=F(i,a)/(F(i,i)-F(a,a))
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(mm,nn,i,a,j,b)
			!$omp do
			do mm = 1,Ne
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					if (notFitRadius(2,i,a,j,b)) cycle

					t2(i,j,a,b)=R(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
					t2(j,i,b,a)=t2(i,j,a,b)
				enddo
			enddo
			!$omp end parallel

		case ('spin-cue-ccs')
			!$omp parallel default(shared) private(i,a)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				if (notFitRadius(1,i,a)) cycle
				t1(i,a)=F(i,a)/(F(i,i)-F(a,a))
			enddo
			enddo
			!$omp end parallel

		case ('spin-cue-ccsd','spin-cue-ccsdt')
			!$omp parallel default(shared) private(i,a)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				if (notFitRadius(1,i,a)) cycle
				t1(i,a)=F(i,a)/(F(i,i)-F(a,a))
			enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(i,a,j,b)
			!$omp do
			do i = 1,Nel-1
			do a = Nel+1,No-1
				do j = i+1,Nel
				do b = a+1,No
					if (notFitRadius(2,i,a,j,b)) cycle
					t2(i,j,a,b)= R(i,a,j,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
					t2(j,i,a,b)=-t2(i,j,a,b)
					t2(i,j,b,a)=-t2(i,j,a,b)
					t2(j,i,b,a)= t2(i,j,a,b)
				enddo
				enddo
			enddo
			enddo
			!$omp end parallel

	end select

	return
	end subroutine guessCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine iterationCC(iteration,epsilon,saccuracy)
	implicit none

	integer(kind=iglu), intent(in)  :: iteration
	real   (kind=rglu), intent(in)  :: epsilon
	real   (kind=rglu), intent(out) :: saccuracy(5)


	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			call projection_ccsd_singles_spin_cue_spare
			call projection_ccsd_doubles_spin_cue_spare

		case ('cue-ccs')
			call projection_ccs_singles_spatial_cue

		case ('cue-ccsd')
			call projection_ccsd_singles_spatial_cue
			call projection_ccsd_doubles_spatial_cue

		case ('spin-cue-ccs')
			call projection_ccs_singles_spin_cue

		case ('spin-cue-ccsd')
			call projection_ccsd_singles_spin_cue
			call projection_ccsd_doubles_spin_cue

		case ('spin-cue-ccsdt')
			call projection_ccsdt_singles_spin_cue
			call projection_ccsdt_doubles_spin_cue
			call projection_ccsdt_triples_spin_cue

		case ('u-ccd','r-ccd')
			call projection_ccd_doubles_spatial_hf

		case ('u-ccsd','r-ccsd')
			call projection_ccsd_singles_spatial_hf
			call projection_ccsd_doubles_spatial_hf

		case ('spin-u-ccd','spin-r-ccd')
			call projection_ccd_doubles_spin_hf

		case ('spin-u-ccsd','spin-r-ccsd','spin-r-ccsd(t)')
			call projection_ccsd_singles_spin_hf
			call projection_ccsd_doubles_spin_hf

		case ('spin-u-ccsdt','spin-r-ccsdt')
			call projection_ccsdt_singles_spin_hf
			call projection_ccsdt_doubles_spin_hf
			call projection_ccsdt_triples_spin_hf

	end select

	call stepCC

	saccuracy=accuracy

	if (maxval(accuracy).LT.epsilon) return

	if (diisbd%enabled) then
		call pushCCVectors

		if (iteration.GT.Nd) then
			call newCCVectors(Nd)
		else
			call newCCVectors(iteration)
		endif
	endif

	return
	end subroutine iterationCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine stepCC

	use coupledClusterSparse, only: Indexs,numcol,erow,ftfm,cIndex
	use coupledClusterSparse, only: spt1=> t1, spd1=> d1, spt2=>vt, spd2=>vd

	implicit none

	integer(kind=iglu) :: i,j,k,a,b,c,mm,nn,vv,sta1,sto1,projtype
	real   (kind=rglu) :: rez,Ax


	select case (uchGet(ccbd%projType))
		case ('1'  ); projtype=1
		case ('2-1'); projtype=2
	end select

	accuracy=-1
	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			!$omp parallel default(shared) private(i,a)
			!$omp do
			do i = 1,Nel; do a = Nel+1,No
				if ( abs(spd1(i,a)).LT.1D-15) spd1(i,a)=0
				spt1(i,a)=spt1(i,a)+ccbd%iterStep(1)*spd1(i,a)/(F(i,i)-F(a,a))
			enddo; enddo
			!$omp end parallel

			spd2(0)=0; spt2(0)=0
			!$omp parallel default(shared) private(mm,sta1,sto1,i,a,vv,nn,j,b,Ax)
			!$omp do
			do mm = 1,Ne
				sta1=erow(mm)
				sto1=erow(mm+1)-1
				i=Indexs(mm,1)
				a=Indexs(mm,2)
				do vv = sta1,sto1
					nn=numcol(vv)

					j=Indexs(nn,1)
					b=Indexs(nn,2)

					if (abs(spd2(vv)).LT.1D-15) spd2(vv)=0

					spt2(vv)=spt2(vv)+ccbd%iterStep(2)*spd2(vv)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(i,j,a,b,Ax)
			!$omp do
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				Ax=spt2( ftfm( cIndex(i,a),cIndex(j,b) ) )
			
				spt2(ftfm(cIndex(i,a),cIndex(j,b)))= Ax
				spt2(ftfm(cIndex(j,a),cIndex(i,b)))=-Ax
				spt2(ftfm(cIndex(i,b),cIndex(j,a)))=-Ax
				spt2(ftfm(cIndex(j,b),cIndex(i,a)))= Ax
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(spd1))
			accuracy(2)=maxval(abs(spd2))

		case ('cue-ccsd','u-ccsd','r-ccsd')
			!$omp parallel default(shared) private(mm,i,a)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				t1(i,a)=t1(i,a)+ccbd%iterStep(1)*d1(i,a)/( F(i,i)-F(a,a) )
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(mm,i,a,nn,j,b)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					select case (projtype)
						case (1); t2(i,j,a,b)=t2(i,j,a,b)+ccbd%iterStep(2)*d2(i,j,a,b)/( F(i,i)-F(a,a)+F(j,j)-F(b,b) )
						case (2); t2(i,j,a,b)=t2(i,j,a,b)+ccbd%iterStep(2)*( 2*d2(i,j,a,b)+d2(i,j,b,a) )/( F(i,i)-F(a,a)+F(j,j)-F(b,b) )
					end select
					t2(j,i,b,a)=t2(i,j,a,b)
				enddo
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(d1))
			accuracy(2)=maxval(abs(d2))

		case ('spin-cue-ccsd','spin-u-ccsd','spin-r-ccsd','spin-r-ccsd(t)')
			!$omp parallel default(shared) private(i,a)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle

				t1(i,a)=t1(i,a)+ccbd%iterStep(1)*d1(i,a)/(F(i,i)-F(a,a))
			enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(i,j,a,b,rez)
			!$omp do
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle

				rez=t2(i,j,a,b)+ccbd%iterStep(2)*d2(i,j,a,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				t2(i,j,a,b)= rez; t2(j,i,b,a)= rez
				t2(i,j,b,a)=-rez; t2(j,i,a,b)=-rez
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(d1))
			accuracy(2)=maxval(abs(d2))

		case ('spin-cue-ccsdt','spin-r-ccsdt','spin-u-ccsdt')
			!$omp parallel default(shared) private(i,a)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle

				t1(i,a)=t1(i,a)+ccbd%iterStep(1)*d1(i,a)/(F(i,i)-F(a,a))
			enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(i,j,a,b,rez)
			!$omp do
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle

				rez=t2(i,j,a,b)+ccbd%iterStep(2)*d2(i,j,a,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				t2(i,j,a,b)= rez; t2(j,i,b,a)= rez
				t2(i,j,b,a)=-rez; t2(j,i,a,b)=-rez
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(i,j,k,a,b,c,rez)
			!$omp do
			do i = 1,Nel-2
			do j = i+1,Nel-1
			do k = j+1,Nel
			do a = Nel+1,No-2
			do b = a+1,No-1
			do c = b+1,No
				if (btest(i+j+k,0).NE.btest(a+b+c,0)) cycle

				rez=t3(i,j,k,a,b,c)+ccbd%iterStep(3)*d3(i,j,k,a,b,c)/(F(i,i)+F(j,j)+F(k,k)-F(a,a)-F(b,b)-F(c,c))
				t3(i,j,k,a,b,c)=+rez; t3(i,j,k,a,c,b)=-rez
				t3(i,j,k,b,a,c)=-rez; t3(i,j,k,b,c,a)=+rez
				t3(i,j,k,c,a,b)=+rez; t3(i,j,k,c,b,a)=-rez

				t3(i,k,j,a,b,c)=-rez; t3(i,k,j,a,c,b)=+rez
				t3(i,k,j,b,a,c)=+rez; t3(i,k,j,b,c,a)=-rez
				t3(i,k,j,c,a,b)=-rez; t3(i,k,j,c,b,a)=+rez

				t3(j,i,k,a,b,c)=-rez; t3(j,i,k,a,c,b)=+rez
				t3(j,i,k,b,a,c)=+rez; t3(j,i,k,b,c,a)=-rez
				t3(j,i,k,c,a,b)=-rez; t3(j,i,k,c,b,a)=+rez

				t3(j,k,i,a,b,c)=+rez; t3(j,k,i,a,c,b)=-rez
				t3(j,k,i,b,a,c)=-rez; t3(j,k,i,b,c,a)=+rez
				t3(j,k,i,c,a,b)=+rez; t3(j,k,i,c,b,a)=-rez

				t3(k,i,j,a,b,c)=+rez; t3(k,i,j,a,c,b)=-rez
				t3(k,i,j,b,a,c)=-rez; t3(k,i,j,b,c,a)=+rez
				t3(k,i,j,c,a,b)=+rez; t3(k,i,j,c,b,a)=-rez

				t3(k,j,i,a,b,c)=-rez; t3(k,j,i,a,c,b)=+rez
				t3(k,j,i,b,a,c)=+rez; t3(k,j,i,b,c,a)=-rez
				t3(k,j,i,c,a,b)=-rez; t3(k,j,i,c,b,a)=+rez
			enddo
			enddo
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(d1))
			accuracy(2)=maxval(abs(d2))
			accuracy(3)=maxval(abs(d3))

		case ('cue-ccs')
			!$omp parallel default(shared) private(mm,i,a)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				t1(i,a)=t1(i,a)+ccbd%iterStep(1)*d1(i,a)/( F(i,i)-F(a,a) )
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(d1))

		case ('spin-cue-ccs')
			!$omp parallel default(shared) private(i,a)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle

				t1(i,a)=t1(i,a)+ccbd%iterStep(1)*d1(i,a)/(F(i,i)-F(a,a))
			enddo
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(d1))

		case ('spin-r-ccd','spin-u-ccd')
			!$omp parallel default(shared) private(i,j,a,b,rez)
			!$omp do
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle

				rez=t2(i,j,a,b)+ccbd%iterStep(2)*d2(i,j,a,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				t2(i,j,a,b)= rez; t2(j,i,b,a)= rez
				t2(i,j,b,a)=-rez; t2(j,i,a,b)=-rez
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(d2))

		case ('r-ccd', 'u-ccd')
			!$omp parallel default(shared) private(mm,i,a,nn,j,b)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					select case (projtype)
						case (1); t2(i,j,a,b)=t2(i,j,a,b)+ccbd%iterStep(2)*d2(i,j,a,b)/( F(i,i)-F(a,a)+F(j,j)-F(b,b) )
						case (2); t2(i,j,a,b)=t2(i,j,a,b)+ccbd%iterStep(2)*( 2*d2(i,j,a,b)+d2(i,j,b,a) )/( F(i,i)-F(a,a)+F(j,j)-F(b,b) )
					end select
					t2(j,i,b,a)=t2(i,j,a,b)
				enddo
			enddo
			!$omp end parallel

			accuracy(1)=maxval(abs(d2))

	end select

	return
	end subroutine stepCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine pushCCVectors
	use coupledClusterSparse, only: vt,vd

	implicit none

	integer(kind=iglu) :: i,j,k,a,b,c,mm,nn,vv,projtype


	select case (uchGet(ccbd%projType))
		case ('1'  ); projtype=1
		case ('2-1'); projtype=2
	end select

	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			do k = Nd,2,-1
				do vv = 1,UBound(sd1,1)
					sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
				enddo

				do vv = 1,UBound(sd2,1)
					sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
				enddo
			enddo

			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				v=vv+1; sd1(vv,1)=d1(i,a); st1(vv,1)=t1(i,a)
			enddo
			enddo

			do vv = 1,UBound(vt,1)
				sd2(vv,1)=vd(vv); st2(vv,1)=vt(vv)
			enddo

		case ('cue-ccs')
			do k = Nd,2,-1
				do vv = 1,UBound(sd1,1)
					sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
				enddo
			enddo

			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				vv=vv+1; sd1(vv,1)=d1(i,a); st1(vv,1)=t1(i,a)
			enddo

		case ('spin-cue-ccs')
			do k = Nd,2,-1
				do vv = 1,UBound(sd1,1)
					sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
				enddo
			enddo

			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1; sd1(vv,1)=d1(i,a); st1(vv,1)=t1(i,a)
			enddo
			enddo

		case ('r-ccd','u-ccd')
			do k = Nd,2,-1
				do vv = 1,UBound(sd2,1)
					sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
				enddo
			enddo

			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					vv=vv+1
					st2(vv,1)=t2(i,j,a,b)
					select case (projtype)
						case (1); sd2(vv,1)=d2(i,j,a,b)
						case (2); sd2(vv,1)=2*d2(i,j,a,b)+d2(i,j,b,a)
					end select
				enddo
			enddo

		case ('spin-r-ccd','spin-u-ccd')
			do k = Nd,2,-1
				do vv = 1,UBound(sd2,1)
					sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
				enddo
			enddo

			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1; sd2(vv,1)=d2(i,j,a,b); st2(vv,1)=t2(i,j,a,b)
			enddo
			enddo
			enddo
			enddo

		case ('cue-ccsd','r-ccsd','u-ccsd')
			do k = Nd,2,-1
				do vv = 1,UBound(sd1,1)
					sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
				enddo

				do vv = 1,UBound(sd2,1)
					sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
				enddo
			enddo

			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				vv=vv+1; sd1(vv,1)=d1(i,a); st1(vv,1)=t1(i,a)
			enddo

			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					vv=vv+1
					st2(vv,1)=t2(i,j,a,b)
					select case (projtype)
						case (1); sd2(vv,1)=d2(i,j,a,b)
						case (2); sd2(vv,1)=2*d2(i,j,a,b)+d2(i,j,b,a)
					end select
				enddo
			enddo

		case ('spin-cue-ccsd','spin-r-ccsd','spin-u-ccsd','spin-r-ccsd(t)')
			do k = Nd,2,-1
				do vv = 1,UBound(sd1,1)
					sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
				enddo

				do vv = 1,UBound(sd2,1)
					sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
				enddo
			enddo

			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1; sd1(vv,1)=d1(i,a); st1(vv,1)=t1(i,a)
			enddo
			enddo

			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1; sd2(vv,1)=d2(i,j,a,b); st2(vv,1)=t2(i,j,a,b)
			enddo
			enddo
			enddo
			enddo

		case ('spin-cue-ccsdt','spin-r-ccsdt','spin-u-ccsdt')
			do k = Nd,2,-1
				do vv = 1,UBound(sd1,1)
					sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
				enddo

				do vv = 1,UBound(sd2,1)
					sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
				enddo

				do vv = 1,UBound(sd3,1)
					sd3(vv,k)=sd3(vv,k-1); st3(vv,k)=st3(vv,k-1)
				enddo
			enddo

			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1; sd1(vv,1)=d1(i,a); st1(vv,1)=t1(i,a)
			enddo
			enddo

			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1; sd2(vv,1)=d2(i,j,a,b); st2(vv,1)=t2(i,j,a,b)
			enddo
			enddo
			enddo
			enddo

			vv=0
			do i = 1,Nel-2
			do j = i+1,Nel-1
			do k = j+1,Nel
			do a = Nel+1,No-2
			do b = a+1,No-1
			do c = b+1,No
				if (btest(i+j+k,0).NE.btest(a+b+c,0)) cycle
				vv=vv+1; sd3(vv,1)=d3(i,j,k,a,b,c); st3(vv,1)=t3(i,j,k,a,b,c)
			enddo
			enddo
			enddo
			enddo
			enddo
			enddo

	end select

	return
	end subroutine pushCCVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine newCCVectors(iteration)

	use coupledClusterSparse, only: vt,vd
	use math, only: tred4

	implicit none

	integer(kind=iglu) :: iteration,k,l,vv,pp,i,j,a,b,c,mm,nn,ub
	real   (kind=rglu) :: sum


	do k = iteration,2,-1
	do l = iteration,2,-1
		diisMatrix(k,l)=diisMatrix(k-1,l-1)
	enddo
	enddo

	do k = 1,iteration
		call diisDotProduct(k,1)
	enddo

	if (iteration.EQ.1) return

	do k = 1,iteration+1
		diisMatrix(k,iteration+1)=1
		diisMatrix(iteration+1,k)=1
	enddo; diisMatrix(iteration+1,iteration+1)=0

	ub=iteration+1
	call tred4(diisMatrix(1:ub,1:ub),diisVectors(1:ub,1:ub),diisValues(1:ub),ub,real(1.d-100,rspu),real(1.d-300,rspu))

	do k = 1,iteration
		diisCoefficients(k)=0
		do l = 1,iteration+1
			diisCoefficients(k)=diisCoefficients(k)+diisVectors(k,l)*diisVectors(iteration+1,l)/diisValues(l)
		enddo
	enddo

	select case (uchGet(umethod))
		case ('spare-cue-ccsd')

			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st1(vv,pp)
				enddo
				t1(i,a)=sum
			enddo
			enddo

			do vv = 1,UBound(vt,1)
				sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st2(vv,pp)
				enddo
				vt(vv)=sum
			enddo

		case ('cue-ccs')
			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st1(vv,pp)
				enddo
				t1(i,a)=sum
			enddo

		case ('spin-cue-ccs')
			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st1(vv,pp)
				enddo
				t1(i,a)=sum
			enddo
			enddo

		case ('r-ccd','u-ccd')
			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					vv=vv+1; sum=0
					do pp = 1,Nd
						sum=sum+diisCoefficients(pp)*st2(vv,pp)
					enddo
					t2(i,j,a,b)=sum
				enddo
			enddo

		case ('spin-r-ccd','spin-u-ccd')
			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st2(vv,pp)
				enddo
				t2(i,j,a,b)=sum
			enddo
			enddo
			enddo
			enddo

		case ('cue-ccsd','r-ccsd','u-ccsd')
			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st1(vv,pp)
				enddo
				t1(i,a)=sum
			enddo

			vv=0
			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					vv=vv+1; sum=0
					do pp = 1,Nd
						sum=sum+diisCoefficients(pp)*st2(vv,pp)
					enddo
					t2(i,j,a,b)=sum
					t2(j,i,b,a)=sum
				enddo
			enddo

		case ('spin-cue-ccsd','spin-r-ccsd','spin-u-ccsd','spin-r-ccsd(t)')
			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st1(vv,pp)
				enddo
				t1(i,a)=sum
			enddo
			enddo

			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st2(vv,pp)
				enddo
				t2(i,j,a,b)= sum
				t2(j,i,a,b)=-sum
				t2(i,j,b,a)=-sum
				t2(j,i,b,a)= sum
			enddo
			enddo
			enddo
			enddo

		case ('spin-cue-ccsdt','spin-r-ccsdt','spin-u-ccsdt')
			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st1(vv,pp)
				enddo
				t1(i,a)=sum
			enddo
			enddo

			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st2(vv,pp)
				enddo
				t2(i,j,a,b)= sum
				t2(j,i,a,b)=-sum
				t2(i,j,b,a)=-sum
				t2(j,i,b,a)= sum
			enddo
			enddo
			enddo
			enddo

			vv=0
			do i = 1,Nel-2
			do j = i+1,Nel-1
			do k = j+1,Nel
			do a = Nel+1,No-2
			do b = a+1,No-1
			do c = b+1,No
				if (btest(i+j+k,0).NE.btest(a+b+c,0)) cycle

				vv=vv+1; sum=0
				do pp = 1,Nd
					sum=sum+diisCoefficients(pp)*st3(vv,pp)
				enddo
				t3(i,j,k,a,b,c)=+sum; t3(i,j,k,a,c,b)=-sum
				t3(i,j,k,b,a,c)=-sum; t3(i,j,k,b,c,a)=+sum
				t3(i,j,k,c,a,b)=+sum; t3(i,j,k,c,b,a)=-sum

				t3(i,k,j,a,b,c)=-sum; t3(i,k,j,a,c,b)=+sum
				t3(i,k,j,b,a,c)=+sum; t3(i,k,j,b,c,a)=-sum
				t3(i,k,j,c,a,b)=-sum; t3(i,k,j,c,b,a)=+sum

				t3(j,i,k,a,b,c)=-sum; t3(j,i,k,a,c,b)=+sum
				t3(j,i,k,b,a,c)=+sum; t3(j,i,k,b,c,a)=-sum
				t3(j,i,k,c,a,b)=-sum; t3(j,i,k,c,b,a)=+sum

				t3(j,k,i,a,b,c)=+sum; t3(j,k,i,a,c,b)=-sum
				t3(j,k,i,b,a,c)=-sum; t3(j,k,i,b,c,a)=+sum
				t3(j,k,i,c,a,b)=+sum; t3(j,k,i,c,b,a)=-sum

				t3(k,i,j,a,b,c)=+sum; t3(k,i,j,a,c,b)=-sum
				t3(k,i,j,b,a,c)=-sum; t3(k,i,j,b,c,a)=+sum
				t3(k,i,j,c,a,b)=+sum; t3(k,i,j,c,b,a)=-sum

				t3(k,j,i,a,b,c)=-sum; t3(k,j,i,a,c,b)=+sum
				t3(k,j,i,b,a,c)=+sum; t3(k,j,i,b,c,a)=-sum
				t3(k,j,i,c,a,b)=-sum; t3(k,j,i,c,b,a)=+sum
			enddo
			enddo
			enddo
			enddo
			enddo
			enddo

	end select

	return
	end subroutine newCCVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine energyCC(energy)

	use coupledClusterSparse, only: Indexs,numcol,erow
	use coupledClusterSparse, only: spt1=> t1, spd1=> d1, spt2=>vt, spd2=>vd

	implicit none

	real   (kind=rglu) :: energy(5)
	real   (kind=rglu) :: sum1,sum2,sum3,amp,c1,c2,c3,denom
	integer(kind=iglu) :: mm,nn,vv,i,j,k,l,a,b,c,d,sta1,sto1


	energy=0
	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			sum1=0
			!$omp parallel default(shared) private(i,a) reduction(+:sum1)
			!$omp do
			do i = 1,Nel
				do a = Nel+1,No
					sum1=sum1+F(i,a)*spt1(i,a)
				enddo
			enddo
			!$omp end parallel

			sum2=0
			!$omp parallel default(shared) private(mm,nn,vv,sta1,sto1,i,a,j,b) reduction(+:sum2)
			!$omp do
			do mm = 1,Ne
				i=indexs(mm,1); a=indexs(mm,2)
				sta1=erow(mm) ; sto1=erow(mm+1)-1		

				do vv = sta1,sto1
					nn=numcol(vv)
					j=indexs(nn,1)
					b=indexs(nn,2)

					sum2=sum2+spin_cue_int(a,i,b,j)*spt2(vv)
				enddo
			enddo
			!$omp end parallel

			sum3=0
			!$omp parallel default(shared) private(i,a,j,b) reduction(+:sum3)
			!$omp do	
			do i = 1,Nel
			do a = Nel+1,No
				do j = 1,Nel
				do b = Nel+1,No
					sum3=sum3+spin_cue_int(a,i,b,j)*(spt1(i,a)*spt1(j,b)-spt1(j,a)*spt1(i,b))
				enddo
				enddo
			enddo
			enddo
			!$omp end parallel
			energy(1)=refeEnergy+sum1+(sum2+sum3)/16

		case ('cue-ccs')
			sum1=0
			!$omp parallel default(shared) private(mm,i,a) reduction(+:sum1)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				sum1=sum1+F(i,a)*t1(i,a)
			enddo
			!$omp end parallel

			energy(1)=refeEnergy+2*sum1

		case ('spin-cue-ccs')
			sum1=0
			!$omp parallel default(shared) private(i,a) reduction(+:sum1)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				sum1=sum1+F(i,a)*t1(i,a)
			enddo
			enddo
			!$omp end parallel
			energy(1)=refeEnergy+sum1

		case ('cue-ccsd','u-ccsd','r-ccsd')
			sum1=0; sum2=0; sum3=0

			!$omp parallel default(shared) private(mm,i,a) reduction(+:sum1)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				sum1=sum1+F(i,a)*t1(i,a)
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(mm,i,a,nn,j,b) reduction(+:sum2)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				do nn = mm+1,Ne
					j=excSet(nn,1); b=excSet(nn,2)

					sum2=sum2+ R(a,i,b,j)*( t2(i,j,a,b) + t1(i,a)*t1(j,b) )
				enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(mm,i,j,a,b) reduction(+:sum3)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				j=excSet(mm,1); b=excSet(mm,2)

				sum3=sum3+ R(a,i,b,j)*( t2(i,j,a,b) + t1(i,a)*t1(j,b) )
			enddo
			!$omp end parallel
			energy(1)=refeEnergy+2*sum1+(2*sum2+sum3)

		case ('spin-cue-ccsd','spin-u-ccsd','spin-r-ccsd','spin-cue-ccsdt','spin-u-ccsdt','spin-r-ccsdt')
			sum1=0
			!$omp parallel default(shared) private(i,a) reduction(+:sum1)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				sum1=sum1+F(i,a)*t1(i,a)
			enddo
			enddo
			!$omp end parallel

			sum2=0
			!$omp parallel default(shared) private(i,j,a,b) reduction(+:sum2)
			!$omp do
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				sum2=sum2+R(i,a,j,b)*(t2(i,j,a,b)+t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a))
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel
			energy(1)=refeEnergy+sum1+sum2

		case ('u-ccd','r-ccd')
			sum1=0; sum2=0
			!$omp parallel default(shared) private(mm,i,a,nn,j,b) reduction(+:sum1)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				do nn = mm+1,Ne
					j=excSet(nn,1); b=excSet(nn,2)

					sum1=sum1+R(a,i,b,j)*t2(i,j,a,b)
				enddo
			enddo
			!$omp end parallel

			!$omp parallel default(shared) private(mm,i,j,a,b) reduction(+:sum2)
			!$omp do
			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				j=excSet(mm,1); b=excSet(mm,2)

				sum2=sum2+ R(a,i,b,j)*t2(i,j,a,b)
			enddo
			!$omp end parallel
			energy(1)=refeEnergy+2*sum1+sum2

		case ('spin-u-ccd','spin-r-ccd')
			sum1=0
			!$omp parallel default(shared) private(i,j,a,b) reduction(+:sum1)
			!$omp do
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				sum1=sum1+R(i,a,j,b)*(t2(i,j,a,b)+t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a))
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel
			energy(1)=refeEnergy+sum1

		case ('spin-r-ccsd(t)')
			sum1=0
			!$omp parallel default(shared) private(i,a) reduction(+:sum1)
			!$omp do
			do i = 1,Nel
			do a = Nel+1,No
				sum1=sum1+F(i,a)*t1(i,a)
			enddo
			enddo
			!$omp end parallel

			sum2=0
			!$omp parallel default(shared) private(i,j,a,b) reduction(+:sum2)
			!$omp do
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				sum2=sum2+R(i,a,j,b)*(t2(i,j,a,b)+t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a))
			enddo
			enddo
			enddo
			enddo
			!$omp end parallel
			energy(3)=refeEnergy+sum1+sum2

			c1=0; c2=0; c3=0
			!$omp parallel default(shared) private(i,j,k,a,b,c,l,d,sum1,sum2,denom,amp) reduction(+:c1,c2,c3)
			!$omp do
			do i = 1,Nel
			do j = 1,Nel
			do k = 1,Nel
			do a = Nel+1,No
			do b = Nel+1,No
			do c = Nel+1,No
					if (btest(i+j+k,0).NE.btest(a+b+c,0)) cycle

					sum1=0
					do l = 1,Nel
						sum1=sum1+R(j,c,k,l)*t2(i,l,a,b)& !<rst>
								 -R(i,c,k,l)*t2(j,l,a,b)&
								 -R(j,c,i,l)*t2(k,l,a,b)&
								 -R(j,a,k,l)*t2(i,l,c,b)&
								 +R(i,a,k,l)*t2(j,l,c,b)&
								 +R(j,a,i,l)*t2(k,l,c,b)&
								 -R(j,b,k,l)*t2(i,l,a,c)&
								 +R(i,b,k,l)*t2(j,l,a,c)&
								 +R(j,b,i,l)*t2(k,l,a,c)
					enddo

					sum2=0
					do d = Nel+1,No
						sum2=sum2+R(b,k,c,d)*t2(i,j,a,d)&
								 -R(b,i,c,d)*t2(k,j,a,d)&
								 -R(b,j,c,d)*t2(i,k,a,d)&
								 -R(a,k,c,d)*t2(i,j,b,d)&
								 +R(a,i,c,d)*t2(k,j,b,d)&
								 +R(a,j,c,d)*t2(i,k,b,d)&
								 -R(b,k,a,d)*t2(i,j,c,d)&
								 +R(b,i,a,d)*t2(k,j,c,d)&
								 +R(b,j,a,d)*t2(i,k,c,d)
					enddo

					denom=F(a,a)+F(b,b)+F(c,c)-F(i,i)-F(j,j)-F(k,k)

					amp=(sum1-sum2)/denom

					c1=c1+amp**2*denom
					c2=c2+amp*t1(i,a)*R(b,j,c,k)
					c3=c3+amp*t2(i,j,a,b)*F(k,c)
				enddo
				enddo
				enddo
			enddo
			enddo
			enddo
			!$omp end parallel
			c1=c1/36; c2=c2/4

			energy(2)=energy(3)+c1
			energy(1)=energy(3)+c1+c2+c3

	end select

	return
	end subroutine energyCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine finalizeCC
	implicit none


	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
		case ('cue-ccs')
		case ('cue-ccsd')
		case ('cue-ccsdt')

		case ('u-ccd')
		case ('u-ccsd')
		case ('u-ccsdt')

		case ('r-ccd')
		case ('r-ccsd')
		case ('r-ccsdt')
		case ('r-ccsd(t)')

	end select

	return
	end subroutine finalizeCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine prepareDIIS

	use coupledClusterSparse, only: vt,vd

	implicit none

	integer(kind=iglu) :: i,j,k,a,b,c,mm,nn,vv
	

	allocate (diisVectors(Nd+1,Nd+1),diisValues(Nd+1))
	allocate (diisMatrix (Nd+1,Nd+1),diisCoefficients(Nd))

	diisVectors=0; diisValues=0
	diisMatrix=0 ; diisCoefficients=0

	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1
			enddo
			enddo

			nn=UBound(sd2,1)

			allocate ( sd1(vv,Nd),st1(vv,Nd) ); sd1=0; st1=0	
			allocate ( sd2(nn,Nd),st2(nn,Nd) ); sd2=0; st2=0

		case ('cue-ccs')
			allocate ( sd1(Ne,Nd),st1(Ne,Nd) ); sd1=0; st1=0

		case ('spin-cue-ccs')
			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1
			enddo
			enddo

			allocate ( sd1(vv,Nd),st1(vv,Nd) ); sd1=0; st1=0

		case ('r-ccd','u-ccd')
			vv=0
			do mm = 1,Ne
			do nn = mm,Ne
				vv=vv+1
			enddo
			enddo

			allocate ( sd2(vv,Nd), st2(vv,Nd) ); sd2=0; st2=0

		case ('spin-r-ccd','spin-u-ccd')
			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1
			enddo
			enddo
			enddo
			enddo

			allocate ( sd2(vv,Nd), st2(vv,Nd) ); sd2=0; st2=0

		case ('cue-ccsd','r-ccsd','u-ccsd')
			vv=0
			do mm = 1,Ne
			do nn = mm,Ne
				vv=vv+1
			enddo
			enddo

			allocate ( sd1(Ne,Nd), st1(Ne,Nd) ); sd1=0; st1=0
			allocate ( sd2(vv,Nd), st2(vv,Nd) ); sd2=0; st2=0

		case ('spin-cue-ccsd','spin-r-ccsd','spin-u-ccsd','spin-r-ccsd(t)')
			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1
			enddo
			enddo
			allocate ( sd1(vv,Nd),st1(vv,Nd) ); sd1=0; st1=0

			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1
			enddo
			enddo
			enddo
			enddo
			allocate ( sd2(vv,Nd), st2(vv,Nd) ); sd2=0; st2=0

		case ('spin-cue-ccsdt','spin-r-ccsdt','spin-u-ccsdt')
			vv=0
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle
				vv=vv+1
			enddo
			enddo
			allocate ( sd1(vv,Nd),st1(vv,Nd) ); sd1=0; st1=0

			vv=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle
				vv=vv+1
			enddo
			enddo
			enddo
			enddo
			allocate ( sd2(vv,Nd), st2(vv,Nd) ); sd2=0; st2=0

			vv=0
			do i = 1,Nel-2
			do j = i+1,Nel-1
			do k = j+1,Nel
			do a = Nel+1,No-2
			do b = a+1,No-1
			do c = b+1,No
				if (btest(i+j+k,0).NE.btest(a+b+c,0)) cycle
				vv=vv+1
			enddo
			enddo
			enddo
			enddo
			enddo
			enddo

			allocate ( sd3(vv,Nd), st3(vv,Nd) ); sd3=0; st3=0

	end select

	return
	end subroutine prepareDIIS

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine diisDotProduct(k1,k2)

	implicit none

	integer(kind=iglu) :: k1,k2,vv
	real   (kind=rglu) :: sum1,sum2,sum3


	sum1=0 ! singles
	select case (uchGet(umethod))
		case ('cue-ccs','spin-cue-ccs',&
		      'cue-ccsd','r-ccsd','u-ccsd',&
		      'spin-cue-ccsd','spin-r-ccsd','spin-u-ccsd',&
		      'spare-cue-ccsd',&
		      'spin-r-ccsd(t)',&
		      'spin-cue-ccsdt','spin-r-ccsdt','spin-u-ccsdt')

			!$omp parallel default(shared) private(vv) reduction(+:sum1)
			!$omp do
			do vv = 1,UBound(sd1,1)
				sum1=sum1+sd1(vv,k1)*sd1(vv,k2)
			enddo
			!$omp end parallel

	end select

	sum2=0 ! doubles
	select case (uchGet(umethod))
		case ('r-ccd','u-ccd','spin-r-ccd','spin-u-ccd',&
		      'cue-ccsd','r-ccsd','u-ccsd',&
		      'spin-cue-ccsd','spin-r-ccsd','spin-u-ccsd',&
		      'spare-cue-ccsd',&
		      'spin-r-ccsd(t)',&
		      'spin-cue-ccsdt','spin-r-ccsdt','spin-u-ccsdt')

			!$omp parallel default(shared) private(vv) reduction(+:sum2)
			!$omp do
			do vv = 1,UBound(sd2,1)
				sum2=sum2+sd2(vv,k1)*sd2(vv,k2)
			enddo
			!$omp end parallel

	end select

	sum3=0 ! triples
	select case (uchGet(umethod))
		case('spin-cue-ccsdt','spin-r-ccsdt','spin-u-ccsdt')

			!$omp parallel default(shared) private(vv) reduction(+:sum3)
			!$omp do
			do vv = 1,UBound(sd3,1)
				sum3=sum3+sd3(vv,k1)*sd3(vv,k2)
			enddo
			!$omp end parallel

	end select

	diisMatrix(k1,k2)=sum1+sum2+sum3; diisMatrix(k2,k1)=diisMatrix(k1,k2)

	return
	end subroutine diisDotProduct

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	real(kind=rglu) function spat_cue_int(a,b,c,d) result(ret) !2*[ab|cd]-[ad|cb]
	implicit none

	integer(kind=iglu), intent(in) :: a,b,c,d
	real   (kind=rglu)             :: sum
	integer(kind=iglu)             :: k,l,mu,nu
	integer(kind=1)                :: irt,ir1,ir2
		

	sum=0
	do k = 1,2
		mu=cueIndex(k,a)
		ir1=2*V(mu,a)*V(mu,b)
		ir2=  V(mu,a)*V(mu,d)
		do l = 1,2
			nu=cueIndex(l,c)
			irt=( ir1*V(nu,d)-ir2*V(nu,b) )*V(nu,c)
			sum=sum+irt*G(mu,nu)
		enddo
	enddo

	ret=sum; return
	end function spat_cue_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	real(kind=rglu) function spin_cue_int(a,b,c,d) result(ret) ! [ab|cd]-[ad|cb]
	implicit none

	integer(kind=iglu), intent(in) :: a,b,c,d
	real   (kind=rglu)             :: sum1,sum2
	integer(kind=1)                :: irt,ir1
	integer(kind=iglu)             :: mu,nu

	sum1=0
	if ((.NOT.btest(a+b,0)).AND.(.NOT.btest(c+d,0))) then ! [ab|cd]
		mu=cueIndex(1,a); ir1=    V(mu,a)*V(mu,b)
		nu=cueIndex(1,c); irt=ir1*V(nu,d)*V(nu,c); sum1=     irt*G(mu,nu)
		nu=cueIndex(2,c); irt=ir1*V(nu,d)*V(nu,c); sum1=sum1+irt*G(mu,nu)

		mu=cueIndex(2,a); ir1=    V(mu,a)*V(mu,b)
		nu=cueIndex(1,c); irt=ir1*V(nu,d)*V(nu,c); sum1=sum1+irt*G(mu,nu)
		nu=cueIndex(2,c); irt=ir1*V(nu,d)*V(nu,c); sum1=sum1+irt*G(mu,nu)
	endif

	sum2=0
	if ((.NOT.btest(a+d,0)).AND.(.NOT.btest(c+b,0))) then ! [ad|cb]
		mu=cueIndex(1,a); ir1=    V(mu,a)*V(mu,d)
		nu=cueIndex(1,c); irt=ir1*V(nu,b)*V(nu,c); sum2=     irt*G(mu,nu)
		nu=cueIndex(2,c); irt=ir1*V(nu,b)*V(nu,c); sum2=sum2+irt*G(mu,nu)
				
		mu=cueIndex(2,a); ir1=    V(mu,a)*V(mu,d)
		nu=cueIndex(1,c); irt=ir1*V(nu,b)*V(nu,c); sum2=sum2+irt*G(mu,nu)
		nu=cueIndex(2,c); irt=ir1*V(nu,b)*V(nu,c); sum2=sum2+irt*G(mu,nu)
	endif

	ret=sum1-sum2; return
	end function spin_cue_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	real(kind=rglu) function spat_hf_int(a,b,c,d) result(ret) ! 2*[ab|cd]-[ad|cb]
	implicit none

	integer(kind=iglu), intent(in) :: a,b,c,d
	real   (kind=rglu)             :: sum,irt,ir1,ir2
	integer(kind=iglu)             :: mu,nu


	sum=0
	do mu = 1,N
		ir1=2*hV(mu,a)*hV(mu,b)
		ir2=  hV(mu,a)*hV(mu,d)
		do nu = 1,N
			irt=( ir1*hV(nu,d)-ir2*hV(nu,b) )*hV(nu,c)
			sum=sum+irt*G(mu,nu)
		enddo
	enddo

	ret=sum; return
	end function spat_hf_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	real(kind=rglu) function spin_hf_int(a,b,c,d) result(ret) ! [ab||cd]=[ab|cd]-[ad|bc]
	implicit none
	integer(kind=iglu), intent(in) :: a,b,c,d
	real   (kind=rglu)             :: sum1,sum2
	integer(kind=iglu)             :: mu,nu


	sum1=0
	if ((.NOT.btest(a+b,0)).AND.(.NOT.btest(c+d,0))) then ! [ab|cd]
		do mu = 1,N
		do nu = 1,N
			sum1=sum1+hVs(mu,a)*hVs(mu,b)*hVs(nu,c)*hVs(nu,d)*G(mu,nu)
		enddo
		enddo
	endif

	sum2=0
	if ((.NOT.btest(a+d,0)).AND.(.NOT.btest(b+c,0))) then ! [ad|bc]
		do mu = 1,N
		do nu = 1,N
			sum2=sum2+hVs(mu,a)*hVs(mu,d)*hVs(nu,b)*hVs(nu,c)*G(mu,nu)
		enddo
		enddo
	endif

	ret=sum1-sum2; return
	end function spin_hf_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine prepareDensity(V)

	implicit none

	real(kind=rglu)    :: V(N,N)
	real(kind=rglu)    :: sum
	integer(kind=iglu) :: i,mu,nu

	!$omp parallel default(shared) private(mu,nu,i,sum)
	!$omp do
	do mu = 1,N
		do nu = 1,N
			sum=0
			do i = 1,Nocc
				sum=sum+V(mu,i)*V(nu,i)
			enddo
			density(mu,nu)=sum
		enddo
	enddo
	!$omp end parallel

	end subroutine prepareDensity

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine prepareFock(orbitals)

	implicit none

	character (len=*) , intent(in)  :: orbitals
	real   (kind=rglu), allocatable :: X(:,:)
	real   (kind=rglu)              :: sum
	integer(kind=iglu)              :: mu,nu,i,j


	allocate (X(N,N)); X=0

	!$omp parallel default(shared) private(mu,nu,sum)
	!$omp do
	do mu = 1,N
		sum=0
		do nu = 1,N
			sum=sum+mol%G(mu,nu)*density(nu,nu)
		enddo 
		do nu = 1,N
			X(mu,nu)=mol%core(mu,nu)-mol%G(mu,nu)*density(mu,nu)
		enddo
		X(mu,mu)=X(mu,mu)+2*sum
	enddo
	!$omp end parallel

	call referenceEnergy(X)

	F=0
	!$omp parallel default(shared) private(i,j,mu,nu,sum)
	!$omp do
	do i = 1,N
		do j = i,N
			sum=0
			do mu = 1,N
			do nu = 1,N
				sum=sum+hV(mu,i)*hV(nu,j)*X(mu,nu)
			enddo
			enddo

			select case (orbitals)
				case ('spin')
					F(2*i-1,2*j-1)=sum; F(2*i,2*j)=sum
					F(2*j-1,2*i-1)=sum; F(2*j,2*i)=sum

				case ('spatial')
					F(i,j)=sum; F(j,i)=sum

			end select
		enddo
	enddo
	!$omp end parallel	
	
	deallocate (X)

	return
	end subroutine prepareFock

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine prepareFockCUE(orbitals)

	implicit none

	character (len=*) , intent(in)  :: orbitals
	real   (kind=rglu), allocatable :: X(:,:)
	real   (kind=rglu)              :: sum,c1,c2
	integer(kind=iglu)              :: mu,nu,i,j


	allocate (X(N,N)); X=0

	!$omp parallel default(shared) private(mu,nu,sum)
	!$omp do
	do mu = 1,N
		sum=0
		do nu = 1,N
			sum=sum+mol%G(mu,nu)*density(nu,nu)
		enddo 
		do nu = 1,N
			X(mu,nu)=mol%core(mu,nu)-mol%G(mu,nu)*density(mu,nu)
		enddo
		X(mu,mu)=X(mu,mu)+2*sum
	enddo
	!$omp end parallel

	call referenceEnergy(X)

	F=0
	!$omp parallel default(shared) private(i,j,c1,c2,mu,nu,sum)
	!$omp do
	do i = 1,N
		do j = i,N
			sum=0
			mu=mol%orb(i)%atoms(1); c1=mol%orb(i)%coef(1)
			nu=mol%orb(j)%atoms(1); c2=mol%orb(j)%coef(1)
			sum=sum+X(mu,nu)*c1*c2

			mu=mol%orb(i)%atoms(1); c1=mol%orb(i)%coef(1)
			nu=mol%orb(j)%atoms(2); c2=mol%orb(j)%coef(2)
			sum=sum+X(mu,nu)*c1*c2

			mu=mol%orb(i)%atoms(2); c1=mol%orb(i)%coef(2)
			nu=mol%orb(j)%atoms(1); c2=mol%orb(j)%coef(1)
			sum=sum+X(mu,nu)*c1*c2

			mu=mol%orb(i)%atoms(2); c1=mol%orb(i)%coef(2)
			nu=mol%orb(j)%atoms(2); c2=mol%orb(j)%coef(2)
			sum=sum+X(mu,nu)*c1*c2

			select case (orbitals)
				case ('spin')
					F(2*i-1,2*j-1)=sum; F(2*i,2*j)=sum
					F(2*j-1,2*i-1)=sum; F(2*j,2*i)=sum

				case ('spatial')
					F(i,j)=sum; F(j,i)=sum

			end select
		enddo
	enddo
	!$omp end parallel

	deallocate (X)

	return
	end subroutine prepareFockCUE

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine referenceEnergy(X)
	implicit none

	real   (kind=rglu) :: Eel,Enuc,X(:,:)
	integer(kind=iglu) :: mu,nu


	Eel=0
	!$omp parallel default(shared) private(mu,nu) reduction(+:Eel)
	!$omp do
	do mu = 1,N
		do nu = 1,N
			Eel=Eel+density(mu,nu)*( mol%core(mu,nu)+X(mu,nu) )
		enddo
	enddo
	!$omp end parallel

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
	refeEnergy=Eel+Enuc

	return
	end subroutine referenceEnergy

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	logical (kind=lglu) function notFitRadius(lvl,i,a,j,b,k,c) result(ret)
	implicit none

	integer(kind=iglu), intent(in)           :: lvl,i,a
	integer(kind=iglu), intent(in), optional :: j,b,k,c
	real   (kind=rglu)                       :: thCentroid(15)


	ret=false
	if (.NOT.cuebd%local(lvl)) return

	if (present(j).AND.present(k)) then
		thCentroid( 1)=cueDistance(i,a)
		thCentroid( 2)=cueDistance(i,b)
		thCentroid( 3)=cueDistance(i,c)
		thCentroid( 4)=cueDistance(i,j)
		thCentroid( 5)=cueDistance(i,k)
		thCentroid( 6)=cueDistance(j,a)
		thCentroid( 7)=cueDistance(j,b)
		thCentroid( 8)=cueDistance(j,c)
		thCentroid( 9)=cueDistance(j,k)
		thCentroid(10)=cueDistance(k,a)
		thCentroid(11)=cueDistance(k,b)
		thCentroid(12)=cueDistance(k,c)
		thCentroid(13)=cueDistance(a,b)
		thCentroid(14)=cueDistance(a,c)
		thCentroid(15)=cueDistance(b,c)
		ret=maxval(thCentroid(1:15)).GT.mol%cueLevel(lvl)
	elseif(present(j)) then
		thCentroid( 1)=cueDistance(i,a)
		thCentroid( 2)=cueDistance(i,b)
		thCentroid( 3)=cueDistance(i,j)
		thCentroid( 4)=cueDistance(j,a)
		thCentroid( 5)=cueDistance(j,b)
		thCentroid( 6)=cueDistance(a,b)
		ret=maxval(thCentroid(1:6)).GT.mol%cueLevel(lvl)
	else
		thCentroid( 1)=cueDistance(i,a)
		ret=maxval(thCentroid(1:1)).GT.mol%cueLevel(lvl)
	endif

	return
	end function notFitRadius

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module coupledCluster