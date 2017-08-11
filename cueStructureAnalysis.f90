	subroutine cueStructureAnalysis

	use hdb, only: iglu,rglu,lglu,true
	use hdb, only: mol,cuebd,ou,su,ouWidth,ouIndent
	use hdb, only: tpFill,tpShowRuler,tpAdjustc,prMatrix,mid

	implicit none

	integer(kind=iglu)              :: N,M,Nocc,NLayers,i,j,k,l,r,p
	logical(kind=lglu)              :: condit
	real   (kind=rglu)              :: sum

	integer(kind=iglu), allocatable :: idist(:,:),moDist(:,:)
	real   (kind=rglu), allocatable :: cueLayers(:),cueDistance(:,:)

	integer(kind=iglu), allocatable :: relLayerPos(:,:),relLayerPopulation(:,:),&
	&                                  absLayerPopulation(:,:),absMOEmployment(:)
	real   (kind=rglu), allocatable :: relMOImportance(:,:),absMOImportance(:)


	N=mol%nAtoms; Nocc=mol%nEls/2; M=mol%nBonds
	allocate (idist(N,N))       ; idist=N**2
	allocate (moDist(Nocc,Nocc)); modist=0
	allocate (cueDistance(N,N)) ; cueDistance=mol%cueDist

	write (ou,'(A/)') tpAdjustc('The start of CUE Analysis',ouWidth,'=')

	write (ou,'(A/)') tpAdjustc('CUE-MO Centroids',73)

	write (ou,'(5X,A,7X,A,3X,A,7X,A,11X,A,11X,A/)') 'CUE-MO','mu','nu','X','Y','Z'
	do k = 1,Nocc
		write (ou,'(2X,i4,",",1X,i4,4X,2(i4,1X),1X,3(F11.5,1X))') &
		k,mol%orb(k)%ova,mol%orb(k)%atoms(1),mol%orb(k)%atoms(2),mol%orb(k)%coords
	enddo

	call prMatrix(cueDistance(1:Nocc,1:Nocc),ou,'Distance matrix for MO','^.00000',maxwidth=ouWidth)

	do i = 1,N
		idist(i,i)=0
	enddo

	do k = 1,M
		i=mol%bnd(k)%atoms(1)
		j=mol%bnd(k)%atoms(2)

		idist(i,j)=1; idist(j,i)=1
	enddo

	! shortest path problem. Floyd-Warshall method.
	do k = 1,N
		do i = 1,N; do j = 1,N
			idist(i,j)=min(idist(i,j),idist(i,k)+idist(k,j))
		enddo; enddo
	enddo

	do p = 1,Nocc-1
		i=mol%orb(p)%atoms(1)
		j=mol%orb(p)%atoms(2)

		do r = p+1,Nocc
			k=mol%orb(r)%atoms(1)
			l=mol%orb(r)%atoms(2)

			moDist(p,r)=min(idist(i,k),idist(i,l),idist(j,k),idist(j,l))
			moDist(r,p)=moDist(p,r)
		enddo
	enddo
	NLayers=maxval(moDist)+1; allocate (cueLayers(0:NLayers))

	cueLayers=0
	do p = 1,Nocc-1
		do r = p+1,Nocc
			
			k=1
			if (moDist(p,r).EQ.k) then
				if (cueDistance(p,r).GT.cueLayers(k+1)) cueLayers(k+1)=cueDistance(p,r)
			endif
		enddo
	enddo

	l=2
	do k = 2,NLayers-1,2

		l=l+1
		do p = 1,Nocc-1
		do r = p+1,Nocc
			if ((moDist(p,r).EQ.k).OR.(moDist(p,r).EQ.k+1)) then
				if (cueDistance(p,r).GT.cueLayers(l)) then
					cueLayers(l)=cueDistance(p,r)
				endif
			endif
		enddo
		enddo

	enddo
	cueLayers=cueLayers+real(1e-5,rglu); cueLayers(0)=-real(1e-5,rglu)

	do k = 1,NLayers-1
		if (cueLayers(k).GT.cueLayers(k+1)) then
			NLayers=k
			exit
		endif
	enddo

	mol%nCUELayers=NLayers
	allocate (mol%cueLayers(0:NLayers)); mol%cueLayers=cueLayers(0:NLayers)

	do k = 0,3
		if ( (cuebd%radius(k).LE.0).OR.(cuebd%radius(k).GT.NLayers)) then
			cuebd%radius(k)=NLayers
		endif
	enddo

	do k = 1,mol%nCUELayers
		do l = 1,3
			if (cuebd%radius(l).EQ.k) then
				mol%cueLevel(l)=mol%cueLayers(k)
			endif
		enddo
	enddo
	mol%cueLevel(0)=mol%cueLayers( cuebd%radius(0) )

	do k = 1,3
		if (cuebd%radius(k).LT.mol%nCUELayers) then
			cuebd%local(k)=true
		endif
	enddo

	write (ou,'(/4X,A/)') 'Layer contents'
	do k = 1,NLayers

		write (ou,'(5X,A,1X,i<mid(NLayers)>,2X,A,F10.5,A,F10.5)') 'CUE Layer',k,'interval',cueLayers(k-1),'  <--->',cueLayers(k)
		l=0
		do i = 1,Nocc
		do j = i,Nocc
			condit=(cueDistance(i,j).GT.cueLayers(k-1)).AND.(cueDistance(i,j).LT.cueLayers(k))

			if (condit) then
				write (ou,'(4X,i4,2X,i4,2X,F10.5)') i,j,cueDistance(i,j)
				l=l+1
			endif
		enddo
		enddo
		write (ou,'(4X,A/22X,i4/)') tpFill(22,'_'),l

	enddo

	write (ou,'(/4X,A/)') 'The size of CUE Layers'
	do k = 1,NLayers
		write (ou,'(4X,i<mid(NLayers)>,1X,F10.5\)') k,cueLayers(k)
		do l = 1,3
			if (cuebd%radius(l).EQ.k) then
				write (ou,'(2X,"t",i1\)') l
			endif
		enddo
		write (ou,*)
	enddo

	allocate (relLayerPos(Nocc,Nocc)) !layer affiliation of "i-j" pairs
	do i = 1,Nocc
	do j = i,Nocc
		do k = 1,NLayers
			condit=(cueDistance(i,j).GT.cueLayers(k-1)).AND.(cueDistance(i,j).LT.cueLayers(k))

			if (condit) then
				relLayerPos(i,j)=k; relLayerPos(j,i)=k
			endif

		enddo
	enddo
	enddo

	allocate (relLayerPopulation(Nocc,NLayers),absLayerPopulation(NLayers,2))
	do k = 1,NLayers
		absLayerPopulation(k,1)=0
		do i = 1,Nocc
			l=0
			do j = 1,Nocc
				condit=(cueDistance(i,j).GT.cueLayers(k-1)).AND.(cueDistance(i,j).LT.cueLayers(k))
				if (condit) l=l+1
			enddo
			relLayerPopulation(i,k)=l                         !number of 'i'   inclusions in 'k' layer
			absLayerPopulation(k,1)=absLayerPopulation(k,1)+l !number of 'i-j' pairs inside  'k' layer.
		enddo
	enddo

	
	allocate (absMOEmployment(Nocc)) !number of layers, where does ethylene employ.
	do i = 1,Nocc
		absMOEmployment(i)=0
		do k = 1,NLayers
			do j = 1,Nocc
				condit=(cueDistance(i,j).GT.cueLayers(k-1)).AND.(cueDistance(i,j).LT.cueLayers(k))
				if (condit) then
					absMOEmployment(i)=absMOEmployment(i)+1; exit
				endif
			enddo
		enddo
	enddo

	!importance of current ethylene.
	allocate (relMOImportance(Nocc,NLayers),absMOImportance(Nocc))
	do i = 1,Nocc
		absMOImportance(i)=0
		do k = 1,NLayers

			sum=0
			do j = 1,Nocc
				if (i.EQ.j) cycle !remove infinity from all sums.
				condit=(cueDistance(i,j).GT.cueLayers(k-1)).AND.(cueDistance(i,j).LT.cueLayers(k))
				if (condit) sum=sum+1/cueDistance(i,j)
			enddo
		
			relMOImportance(i,k)=sum
			absMOImportance(i)=absMOImportance(i)+sum
		enddo
	enddo

	call prMatrix(relLayerPos       ,ou,'Relative MO-layer arrangement','^',maxwidth=ouWidth)
	call prMatrix(relLayerPopulation,ou,'Relative layer population'    ,'^',maxwidth=ouWidth,transpose=true)
	write (ou,'(/4X,A/)') 'Absolute layer population'
	do k = 1,NLayers
		write (ou,'(4X,A,i<mid(NLayers)>,A,1X,i<mid(Nocc**2)>)') 'Layer(',k,') ',absLayerPopulation(k,1)
	enddo

	call prMatrix(relMOImportance   ,ou,'Relative MO importance'       ,'^.00000',maxwidth=ouWidth,transpose=true)
	write (ou,'(/4X,A/)') 'Absolute MO employment and specific importance (by employment)'
	do i = 1,Nocc
		write (ou,'(4X,A,i<mid(Nocc)>,A,1X,i<mid(NLayers)>,1X,F10.5)') 'MO(',i,') ',absMOEmployment(i),absMOImportance(i)/absMOEmployment(i)
	enddo

	write (ou,'(/A/)') tpAdjustc('The end of CUE Analysis',ouWidth,'=')

	deallocate (idist,moDist,cueLayers,cueDistance)
	deallocate (relLayerPos,relLayerPopulation,absLayerPopulation)
	deallocate (absMOEmployment,relMOImportance,absMOImportance)

	return
	end subroutine cueStructureAnalysis