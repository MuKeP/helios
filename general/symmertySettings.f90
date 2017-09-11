	subroutine symmetrySettings

	use glob     , only: rglu,iglu,lglu,true,false,mid
	use glob     , only: uch,uchGet,uchSet
	use hdb      , only: mol,geometrybd,polarizbd,generalbd,ou,ouWidth
	use hdb      , only: GlEt,Et,MEt,MMEt
	use hdb      , only: pointAccordance,pointSet,pointToPut,pointToCalc
	use txtParser, only: tpAdjustc
	use printmod , only: prMatrix
	use math     , only: gltred4

	implicit none

	real   (kind=rglu) :: symmetryTolerance,gridRange,sum,ul,dl,diff
	integer(kind=iglu) :: N,M,Np,Ncue,dSize,i,j,k,l,a,b,c,sta,sto,cElements,ubndDistr,&
	                      ipos,jpos,zpos,uniqueAtoms,writeMypnts

	real   (kind=rglu), allocatable :: coords(:,:),distance(:,:),eigenVectors(:,:),eigenValues(:),&
	                                   d1gridEquivalence(:),&
									   d2gridEquivalence(:,:),evecDifference(:,:),&
									   d3gridEquivalence(:,:,:)

	integer(kind=iglu), allocatable :: d1grid(:),d2grid(:,:),d3grid(:,:,:),pntDistribution(:)
	integer(kind=iglu), allocatable :: atomEqu(:,:),bondEqu(:,:)



	write (ou,'(/A/)') tpAdjustc('Symmetry analysis',ouWidth,'=')

	symmetryTolerance=geometrybd%symmetryTolerance; N=mol%nAtoms; M=mol%nBonds; Ncue=mol%nEls/2

	write (ou,120) symmetryTolerance
	Np=polarizbd%nPoints; sta=-(Np-1)/2; sto=-sta
	select case( uchGet(generalbd%task) )
		case ('polarizability')
			dSize=N+Ncue+Np**3

		case ('density','coulson','hypercharges')
			dSize=N+Ncue+N*(N-1)/2

		case ('energy','wf-analize')
			return

	end select

	allocate (coords(3,dSize),distance(dSize,dSize),eigenVectors(dSize,dSize),eigenValues(dSize))

	cElements=0
	do i = 1,N
		cElements=cElements+1
		do k = 1,3
			coords(k,cElements)=mol%atm(i)%coords(k)
		enddo
	enddo

	do i = 1,Ncue
		cElements=cElements+1
		do k = 1,3
			coords(k,cElements)=mol%orb(i)%coords(k)
		enddo
	enddo

	select case( uchGet(generalbd%task) )
		case ('polarizability')
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

	call gltred4(distance,eigenVectors,eigenValues,dSize,real(1d-100,rglu),real(1d-100,rglu))

	select case( uchGet(generalbd%task) )
		case ('polarizability')
			allocate (d1gridEquivalence(sta:sto),&
			&         d2gridEquivalence(sta:sto,sta:sto),&
			&         d3gridEquivalence(sta:sto,sta:sto,sta:sto))

			allocate (GlEt(sta:sto,sta:sto,sta:sto))
			allocate (Et(sta:sto),MEt(sta:sto,sta:sto),MMEt(sta:sto,sta:sto,sta:sto))

			d3gridEquivalence=0
			do i = dSize-Np**3+1,dSize
				d3gridEquivalence( int(coords(1,i)/gridRange),int(coords(2,i)/gridRange),int(coords(3,i)/gridRange) )=eigenVectors(i,dSize)
			enddo

			ubndDistr=16; allocate (pntDistribution(0:ubndDistr))
			do l = 0,ubndDistr
				pntDistribution(l)=0; ul=real(10,rglu)**real(-l); dl=real(10,rglu)**(-l-1)

				if (l.EQ.ubndDistr) dl=-1

				do i = sta,sto
				do j = sta,sto
				do k = sta,sto
					if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle

					do a = sta,sto
					do b = sta,sto
					do c = sta,sto
						if ((a.NE.0).AND.(b.NE.0).AND.(c.NE.0)) cycle

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

			allocate (d1grid(sta:sto),d2grid(sta:sto,sta:sto),d3grid(sta:sto,sta:sto,sta:sto))

			pointAccordance=0; pointSet=0; pointToCalc=1
			select case (uchGet(polarizbd%scales))
				
				case ('x','y','z')

					select case (uchGet(polarizbd%scales))
						case ('x');	d1gridEquivalence=d3gridEquivalence(:,0,0); ipos=1
						case ('y'); d1gridEquivalence=d3gridEquivalence(0,:,0); ipos=2
						case ('z'); d1gridEquivalence=d3gridEquivalence(0,0,:); ipos=3
					end select

					d1grid=0; pointToPut=0
					do i = sto,sta,-1
						if (d1grid(i).NE.0) cycle

						do a = sto,sta,-1
							if (d1grid(a).NE.0) cycle

							diff=abs(d1gridEquivalence(i)-d1gridEquivalence(a))

							if (diff.LT.symmetryTolerance) then
								pointToPut=pointToPut+1
								pointAccordance(pointToPut,1,ipos)=i
								pointAccordance(pointToPut,2,ipos)=a
								d1grid(a)=-1
							endif
						enddo
						d1grid(i)=1
					enddo

					do i = sta,sto
						if (d1grid(i).EQ.1) then
							if (i.EQ.0) cycle
							pointToCalc=pointToCalc+1
							pointSet(ipos,pointToCalc)=i
						endif
					enddo

				case ('xy','yx','xz','zx','yz','zy')

					if (uchGet(polarizbd%scales).EQ.'yx') then
						polarizbd%scales=uchSet('xy')
					elseif(uchGet(polarizbd%scales).EQ.'zx') then
						polarizbd%scales=uchSet('xz')
					elseif(uchGet(polarizbd%scales).EQ.'zy') then
						polarizbd%scales=uchSet('yz')
					endif

					select case (uchGet(polarizbd%scales,1,1))
						case ('x'); ipos=1
						case ('y'); ipos=2
						case ('z'); ipos=3
					end select

					select case (uchGet(polarizbd%scales,2,2))
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

					select case (zpos)
						case (1); d2gridEquivalence=d3gridEquivalence(0,:,:)
						case (2); d2gridEquivalence=d3gridEquivalence(:,0,:)
						case (3); d2gridEquivalence=d3gridEquivalence(:,:,0)
					end select

					d2grid=0; pointToPut=0
					do i = sto,sta,-1
					do j = sto,sta,-1
						if (d2grid(i,j).NE.0) cycle

						do a = sto,sta,-1
						do b = sto,sta,-1
							if (d2grid(a,b).NE.0) cycle

							diff=abs(d2gridEquivalence(i,j)-d2gridEquivalence(a,b))

							if (diff.LT.symmetryTolerance) then
								pointToPut=pointToPut+1
								pointAccordance(pointToPut,1,ipos)=i; pointAccordance(pointToPut,1,jpos)=j
								pointAccordance(pointToPut,2,ipos)=a; pointAccordance(pointToPut,2,jpos)=b
								d2grid(a,b)=-1
							endif
						enddo
						enddo
						d2grid(i,j)=1
					enddo
					enddo

					do i = sta,sto
					do j = sta,sto
						if (d2grid(i,j).EQ.1) then
							if ((i.EQ.0).AND.(j.EQ.0)) cycle
							pointToCalc=pointToCalc+1
							pointSet(ipos,pointToCalc)=i
							pointSet(jpos,pointToCalc)=j
						endif
					enddo
					enddo

				case ('xyz')

					d3grid=0; pointToPut=0
					do i = sto,sta,-1
					do j = sto,sta,-1
					do k = sto,sta,-1
						if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle

						if (d3grid(i,j,k).NE.0) cycle

						do a = sto,sta,-1
						do b = sto,sta,-1
						do c = sto,sta,-1
							if (d3grid(a,b,c).NE.0) cycle

							diff=abs(d3gridEquivalence(i,j,k)-d3gridEquivalence(a,b,c))

							if (diff.LT.symmetryTolerance) then
								pointToPut=pointToPut+1
								pointAccordance(pointToPut,1,1)=i; pointAccordance(pointToPut,1,2)=j; pointAccordance(pointToPut,1,3)=k
								pointAccordance(pointToPut,2,1)=a; pointAccordance(pointToPut,2,2)=b; pointAccordance(pointToPut,2,3)=c
								d3grid(a,b,c)=-1
							endif
						enddo
						enddo
						enddo
						d3grid(i,j,k)=1
					enddo
					enddo
					enddo

					do i = sta,sto
					do j = sta,sto
					do k = sta,sto
						if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle

						if (d3grid(i,j,k).EQ.1) then
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

			deallocate ( d1grid,d2grid,d3grid )
			deallocate ( d1gridEquivalence,d2gridEquivalence,d3gridEquivalence)

			select case (len(uchGet(polarizbd%scales)))
				case (1); writeMypnts=Np
				case (2); writeMypnts=Np**2
				case (3); writeMypnts=3*(Np**2-2*Np+1)+3*(Np-1)+1 !unique+cross{x00,0y0,00z}+zero{000}
			end select

			write (ou,100)
			do i = 1,pointToPut
				write (ou,101) pointAccordance(i,1,:),pointAccordance(i,2,:)
			enddo
			write (ou,110) pointToCalc,writeMypnts
			write (ou,*)

			write (ou,130) (3*(Np**2-2*Np+1)+3*(Np-1)+1) * (3*(Np**2-2*Np+1)+3*(Np-1))
			do l = 0,ubndDistr
				ul=real(10,rglu)**(-l); dl=real(10,rglu)**(-l-1)
				if (l.EQ.ubndDistr) dl=0
				write (ou,131) ul,dl,pntDistribution(l)
			enddo
			write (ou,*)

			deallocate (pntDistribution)

		case ('density','coulson','hypercharges')
			allocate (atomEqu(N,N),bondEqu(M,M)); atomEqu=0; bondEqu=0
			allocate (evecDifference(dSize,dSize))

			do i = 1,dSize
			do j = i,dSize
				evecDifference(i,j)=abs(eigenVectors(i,dSize)-eigenVectors(j,dSize))
				evecDifference(j,i)=evecDifference(i,j)
			enddo
			enddo

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

			do l = 0,ubndDistr
				pntDistribution(l)=0; ul=real(10,rglu)**(-l); dl=real(10,rglu)**(-l-1)

				if (l.EQ.ubndDistr) dl=-1

				do i = 1,N
				do j = 1,N

					if ( (evecDifference(i,j).GT.dl) .AND. (evecDifference(i,j).LE.ul) ) then
						pntDistribution(l)=pntDistribution(l)+1
					endif
				enddo
				enddo
			enddo

			write (ou,132) N*N
			do l = 0,ubndDistr
				ul=real(10,rglu)**(-l); dl=real(10,rglu)**(-l-1)
				if (l.EQ.ubndDistr) dl=0
				write (ou,131) ul,dl,pntDistribution(l)
			enddo
			write (ou,*)
			call prMatrix(evecDifference,ou,'Atom equivalence','0.0000E00')
			deallocate (pntDistribution,evecDifference)

	end select

	write (ou,'(/A/)') tpAdjustc('The start of computation',ouWidth,'=')

100	format ( 4X,'Point accordance as a consequence of symmetry transformation.')
110 format ( 4X,'Points to be calculated: ',i<mid(pointToCalc)>,' from ',i<mid(writeMypnts)>)
101 format ( 4X,'(',i2,',',i2,',',i2,')--->(',i2,',',i2,',',i2,')')
120 format ( 4X,'Tollerance in grid symmetry definition:',1X,ES9.3)
121 format ( 4X,'Range for grid points:',1X,ES9.3)
130 format ( 4X,'Combinations in 3d grid',1X,i6,'. Distribution:')
131 format ( 4X,ES9.3,1X,'<-->',1X,ES9.3,1X,i5,1X,'combinations')
132 format ( 4X,'Atom-atom combinations ',1X,i6,'. Distribution:')
133 format (/4X,'Unique atom-atom pairs:',1X,i5)

	return

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine defAtomsEquivalence
		implicit none


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

				if (evecDifference(N+Ncue+k,N+Ncue+l).LE.symmetryTolerance) then
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

		return
		end subroutine defAtomsEquivalence

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine defBondsEquivalence
		implicit none

		logical(kind=lglu), allocatable :: d(:,:)
		integer(kind=iglu), allocatable :: bond(:,:)
		real   (kind=rglu), allocatable :: X(:),Y(:),Z(:)
		integer(kind=iglu)              :: i,j,a,b,k(4),l(4),aa,bb,cc,dd
		real   (kind=rglu)              :: centr(2)


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
		! omp private (j)
		do i = 1,N
			do j = 1,N
				d(i,j)=abs(atomEqu(i,i)).EQ.abs(atomEqu(j,j))
			enddo
		enddo

		do i = 1,M !!!
			write (ou,*) i,bond(1,i),bond(2,i)
		enddo

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

		deallocate (d,bond,X,Y,Z)

		return
		end subroutine defBondsEquivalence

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end subroutine symmetrySettings