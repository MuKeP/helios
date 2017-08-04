	module coupledCluster

	use hdb, only: uch,uchGet,uchSet,iglu,rglu,lglu,true,false
	use hdb, only: gluCompare
	use hdb, only: mol,ccbd,diisbd,cuebd,systembd,ou
	use hdb, only: operator(.in.),prMatrix

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	!use cc_sparse

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	! General CC arrays
	real   (kind=rglu), allocatable :: F(:,:),R(:,:,:,:)
	real   (kind=rglu), allocatable :: t1(:,:),t2(:,:,:,:),t3(:,:,:,:,:,:)
	real   (kind=rglu), allocatable :: d1(:,:),d2(:,:,:,:),d3(:,:,:,:,:,:)

	! Information on non-zero Fockian elements
	integer(kind=iglu), allocatable :: Fnz(:,:),excSet(:,:) ! (2,NFnz)

	real   (kind=rglu), allocatable :: density(:,:),cueDistance(:,:)

	! DIIS storage
!	real   (kind=rglu), allocatable :: st1(:,:),st2(:,:),st3(:,:)
!	real   (kind=rglu), allocatable :: sd1(:,:),sd2(:,:),sd3(:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type(uch)          :: umethod
	integer(kind=iglu) :: N,Ne,Nth,Nocc,Nel,No,NFnz
	logical(kind=lglu) :: dcue,dsparse
	real   (kind=rglu) :: accuracy(2)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	integer(kind=1)   , allocatable :: V(:,:)
	real   (kind=rglu), allocatable :: hV(:,:),G(:,:)
	integer(kind=iglu), allocatable :: cueIndex(:,:),iapairs(:) !,spinpairs(:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine setCCParameters(method)
	implicit none

	character (len=*), intent(in) :: method
	integer(kind=iglu)            :: i,j,k,l


	dcue=method .in. ['cue-ccs','cue-ccsd','cue-ccsdt']
	dsparse=cuebd%sparse .AND. (method.EQ.'cue-ccsd')
	Nth=systembd%nNodes

	if (dsparse) then
		umethod=uchSet('spare-'//method)
	else
		umethod=uchSet(method)
	endif

	N=mol%nAtoms; Nel=mol%nEls

	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			No=2*N; Nocc=mol%nEls/2
			allocate (cueIndex(2,No),V(0:N,No),iapairs(No)) !,spinpairs(No))
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

			!do k = 1,N
			!	spinpairs(2*k-1)=2*k; spinpairs(2*k)=2*k-1
			!enddo

			allocate (G(N,N),density(N,N)); G=mol%G; density=0

			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)

				if (mol%orb(k)%nels.EQ.2) then
					density(i,i)=1
					cycle
				endif

				density(i,j)=real(0.5,rglu)
				density(j,i)=real(0.5,rglu)
				density(i,i)=real(0.5,rglu)
				density(j,j)=real(0.5,rglu)
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
			No=N; Nocc=mol%nEls/2; Ne=(N-Nocc)*Nocc

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

			! one- and two-electron integrals
			allocate (G(N,N),density(N,N)); G=mol%G; density=0

			do k = 1,Nocc
				i=mol%orb(k)%atoms(1)
				j=mol%orb(k)%atoms(2)

				if (mol%orb(k)%nels.EQ.2) then
					density(i,i)=1
					cycle
				endif

				density(i,i)=real(0.5,rglu)
				density(i,j)=real(0.5,rglu)
				density(j,i)=real(0.5,rglu)
				density(j,j)=real(0.5,rglu)
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
			
			allocate ( R(N,N,N,N),F(N,N) )

			do i = 1,N
			do j = 1,N
			do k = 1,N
			do l = 1,N
				R(i,j,k,l)=spat_cue_int(i,j,k,l)*real(0.25,rglu)
			enddo
			enddo
			enddo
			enddo

			allocate (t1(Nocc,Nocc+1:N),t2(Nocc,Nocc,Nocc+1:N,Nocc+1:N),&
			          d1(Nocc,Nocc+1:N),d2(Nocc,Nocc,Nocc+1:N,Nocc+1:N) )

			t1=0; t2=0; d1=0; d2=0

		case ('cue-ccsdt')
			No=2*N; Nocc=mol%nEls/2; Ne=(N-Nocc)*Nocc
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

				density(i,j)=real(0.5,rglu)
				density(j,i)=real(0.5,rglu)
				density(i,i)=real(0.5,rglu)
				density(j,j)=real(0.5,rglu)
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
			
			allocate ( R(No,No,No,No), F(No,No) )

			do i = 1,No
			do j = 1,No
			do k = 1,No
			do l = 1,No
				R(i,j,k,l)=spin_cue_int(i,j,k,l)*real(0.25,rglu)
			enddo
			enddo
			enddo
			enddo

			allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No),&
			          d1(Nel,Nel+1:No),d2(Nel,Nel,Nel+1:No,Nel+1:No) )

			t1=0; t2=0; d1=0; d2=0

		case ('u-ccd','u-ccsd','u-ccsdt','r-ccd','r-ccsd','r-ccsdt','r-ccsd(t)','r-ccsd[t]')
			allocate (hV(N,No)); hV=0

	end select

	return
	end subroutine setCCParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine initCC
	implicit none

	integer (kind=iglu) :: i,j,a,b


	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			call prepareFockCUE('spin')

			!deallocate (G,V,density)

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
			
			!deallocate (G,V,density)

			! collecting non-zero Focking elements.
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

		case ('cue-ccsdt')
			call prepareFockCUE('spin')

			!deallocate (G,V,density)

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

		case ('u-ccd')
		case ('u-ccsd')
		case ('u-ccsdt')

		case ('r-ccd')
		case ('r-ccsd')
		case ('r-ccsdt')
		case ('r-ccsd(t)')
		case ('r-ccsd[t]')

		case default; stop
	end select

	return
	end subroutine initCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine iterationCC
	implicit none


	select case (uchGet(umethod))
		case ('spare-cue-ccsd')
			call projection_ccsd_singles_spin_cue_spare
			call projection_ccsd_doubles_spin_cue_spare
		case ('cue-ccs')
		case ('cue-ccsd')
			call projection_ccsd_singles_spatial_cue
			call projection_ccsd_doubles_spatial_cue

		case ('cue-ccsdt')
			call projection_ccsdt_singles_spin_cue
			call projection_ccsdt_doubles_spin_cue

		case ('u-ccd')
		case ('u-ccsd')
		case ('u-ccsdt')

		case ('r-ccd')
		case ('r-ccsd')
		case ('r-ccsdt')
		case ('r-ccsd(t)')
		case ('r-ccsd[t]')

		case default; stop
	end select

	!call stepCC
	!call pushCCVectors
	!call newCCVectors

	return
	end subroutine iterationCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine stepCC

	use coupledClusterSparse, only: Indexs,numcol,erow,ftfm,cIndex
	use coupledClusterSparse, only: spt1=> t1, spd1=> d1, spt2=>vt, spd2=>vd

	implicit none

	integer(kind=iglu) :: i,j,a,b,mm,nn,vv,sta1,sto1
	real   (kind=rglu) :: rez,Ax


	select case (uchGet(umethod))
		case ('spare-cue-ccsd')

			do i = 1,Nel; do a = Nel+1,No
				if ( abs(spd1(i,a)).LT.1D-15) spd1(i,a)=0
				spt1(i,a)=spt1(i,a)+spd1(i,a)/(F(i,i)-F(a,a))
			enddo; enddo

			spd2(0)=0
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

					spt2(vv)=spt2(vv)+spd2(vv)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				enddo
			enddo
			!$omp end parallel

			do i = 1,Nel
			do j = i+1,Nel
			do a = Nel+1,No
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

			accuracy(1)=maxval(abs(spd1))
			accuracy(2)=maxval(abs(spd2))
			write (*,'(2(2X,ES16.8)\)') accuracy; call energyCC

		case ('cue-ccs')
		case ('cue-ccsd')

			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)

				t1(i,a)=t1(i,a)+0.25d0*d1(i,a)/( F(i,i)-F(a,a) )
			enddo

			do mm = 1,Ne
				i=excSet(mm,1)
				a=excSet(mm,2)
				do nn = mm,Ne
					j=excSet(nn,1)
					b=excSet(nn,2)

					!t2(i,j,a,b)=t2(i,j,a,b)+0.25d0*( 2.0d0*d2(i,j,a,b)+d2(i,j,b,a) )/( F(i,i)-F(a,a)+F(j,j)-F(b,b) )
					t2(i,j,a,b)=t2(i,j,a,b)+0.25d0*d2(i,j,a,b)/( F(i,i)-F(a,a)+F(j,j)-F(b,b) )
					t2(j,i,b,a)=t2(i,j,a,b)
				enddo
			enddo
			!write (70,*) d1(1:Nocc,Nocc+1:N)
			!write (90,*) d2(1:Nocc,1:Nocc,Nocc+1:N,Nocc+1:N)
			!read (*,*)

			accuracy(1)=maxval(abs(d1))
			accuracy(2)=maxval(abs(d2))
			write (*,'(2(2X,ES16.8)\)') accuracy; call energyCC

		case ('cue-ccsdt')
			do i = 1,Nel
			do a = Nel+1,No
				if (btest(i,0).NE.btest(a,0)) cycle

				t1(i,a)=t1(i,a)+0.25d0*d1(i,a)/(F(i,i)-F(a,a))
			enddo
			enddo

			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				if (btest(i+j,0).NE.btest(a+b,0)) cycle

				rez=t2(i,j,a,b)+0.25d0*d2(i,j,a,b)/(F(i,i)+F(j,j)-F(a,a)-F(b,b))
				t2(i,j,a,b)= rez
				t2(j,i,b,a)= rez
				t2(i,j,b,a)=-rez
				t2(j,i,a,b)=-rez
			enddo
			enddo
			enddo
			enddo

			accuracy(1)=maxval(abs(d1))
			accuracy(2)=maxval(abs(d2))

			write (*,'(2(2X,ES16.8)\)') accuracy; call energyCC
		case ('u-ccd')
		case ('u-ccsd')
		case ('u-ccsdt')

		case ('r-ccd')
		case ('r-ccsd')
		case ('r-ccsdt')
		case ('r-ccsd(t)')
		case ('r-ccsd[t]')

		case default; stop
	end select





	return
	end subroutine stepCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine pushCCVectors
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
		case ('r-ccsd[t]')

		case default; stop
	end select

	return
	end subroutine pushCCVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine newCCVectors
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
		case ('r-ccsd[t]')

		case default; stop
	end select

	return
	end subroutine newCCVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine energyCC

	use coupledClusterSparse, only: Indexs,numcol,erow
	use coupledClusterSparse, only: spt1=> t1, spd1=> d1, spt2=>vt, spd2=>vd

	implicit none

	real(kind=rglu)    :: sum1,sum2,sum3
	integer(kind=iglu) :: mm,nn,vv,i,j,a,b,sta1,sto1


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
			write (*,'(F21.13\)') sum1+(sum2+sum3)/16
			!write (*,'(A,1X,F18.13\)') '  dE',-13.8722417885897d0-(sum1+(sum2+sum3)/16)

		case ('cue-ccs')
		case ('cue-ccsd')
			sum1=0; sum2=0; sum3=0

			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				sum1=sum1+F(i,a)*t1(i,a)
			enddo

			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)

				do nn = mm+1,Ne
					j=excSet(nn,1); b=excSet(nn,2)

					sum2=sum2+ R(a,i,b,j)*( t2(i,j,a,b) + t1(i,a)*t1(j,b) )
				enddo
			enddo

			do mm = 1,Ne
				i=excSet(mm,1); a=excSet(mm,2)
				j=excSet(mm,1); b=excSet(mm,2)

				sum3=sum3+ R(a,i,b,j)*( t2(i,j,a,b) + t1(i,a)*t1(j,b) )
			enddo

			write (*,'(F21.13\)') 2*sum1 + (2*sum2+sum3)

		case ('cue-ccsdt')
			sum1=0
			do i = 1,Nel
			do a = Nel+1,No
				sum1=sum1+F(i,a)*t1(i,a)
			enddo
			enddo

			sum2=0
			do i = 1,Nel-1
			do j = i+1,Nel
			do a = Nel+1,No-1
			do b = a+1,No
				sum2=sum2+R(i,a,j,b)*(t2(i,j,a,b)+t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a))
			enddo
			enddo
			enddo
			enddo
			write (*,'(F21.13\)') sum1+sum2
		case ('u-ccd')
		case ('u-ccsd')
		case ('u-ccsdt')

		case ('r-ccd')
		case ('r-ccsd')
		case ('r-ccsdt')
		case ('r-ccsd(t)')
		case ('r-ccsd[t]')

		case default; stop
	end select

	return
	end subroutine energyCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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
		case ('r-ccsd[t]')

		case default; stop
	end select

	return
	end subroutine finalizeCC

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu) function spat_cue_int(a,b,c,d) result(ret) !2*[ab|cd]-[ad|cb]
	implicit none

	integer(kind=1)   , parameter  :: btwo= int(2,kind=1)
	integer(kind=iglu), intent(in) :: a,b,c,d
	real   (kind=rglu)             :: sum
	integer(kind=iglu)             :: k,l,mu,nu
	integer(kind=1)                :: irt,ir1,ir2
		

	sum=0
	do k = 1,2
		mu=cueIndex(k,a)
		ir1=btwo*V(mu,a)*V(mu,b)
		ir2=     V(mu,a)*V(mu,d)
		do l = 1,2
			nu=cueIndex(l,c)
			irt=( ir1*V(nu,d)-ir2*V(nu,b) )*V(nu,c)
			sum=sum+irt*G(mu,nu)
		enddo
	enddo

	ret=sum; return
	end function spat_cue_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu) function spin_hf_int(a,b,c,d) result(ret) ! [ab||cd]=[ab|cd]-[ad|bc]
	implicit none
	integer(kind=iglu), intent(in) :: a,b,c,d
	real   (kind=rglu)             :: sum1,sum2
	integer(kind=iglu)             :: mu,nu


	sum1=0
	if ((.NOT.btest(a+b,0)).AND.(.NOT.btest(c+d,0))) then ! [ab|cd]
		do mu = 1,N
		do nu = 1,N
			sum1=sum1+hV(mu,a)*hV(mu,b)*hV(nu,c)*hV(nu,d)*G(mu,nu)
		enddo
		enddo
	endif

	sum2=0
	if ((.NOT.btest(a+d,0)).AND.(.NOT.btest(b+c,0))) then ! [ad|bc]
		do mu = 1,N
		do nu = 1,N
			sum2=sum2+hV(mu,a)*hV(mu,d)*hV(nu,b)*hV(nu,c)*G(mu,nu)
		enddo
		enddo
	endif

	ret=sum1-sum2; return
	end function spin_hf_int

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine prepareFock(orbitals)

	implicit none

	character (len=*) , intent(in)  :: orbitals
	real   (kind=rglu), allocatable :: X(:,:)
	real   (kind=rglu)              :: sum
	integer(kind=iglu)              :: mu,nu,i,j


	allocate (X(N,N)); X=0
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

	F=0
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
	enddo; deallocate (X)

	return
	end subroutine prepareFock

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	subroutine prepareFockCUE(orbitals)

	implicit none

	character (len=*) , intent(in)  :: orbitals
	real   (kind=rglu), allocatable :: X(:,:)
	real   (kind=rglu)              :: sum,c1,c2
	integer(kind=iglu)              :: mu,nu,i,j


	allocate (X(N,N)); X=0

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

	F=0
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
	enddo; deallocate (X)

	return
	end subroutine prepareFockCUE

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module coupledCluster