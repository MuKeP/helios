	subroutine projection_ccsd_singles_spatial_cue

	use glob          , only: rglu,iglu,lglu
	use coupledCluster, only: N,Nocc,Ne,NFnz
	use coupledCluster, only: F,R,Fnz
	use coupledCluster, only: t1,t2,d1,excSet,iapairs
	use coupledCluster, only: notFitRadius

	implicit none

	real   (kind=rglu) :: rez,sum
	integer(kind=iglu) :: u,i,a,j,b,k,c
	integer(kind=iglu) :: t,ab(2)


	do u = 1,Ne
		i=excSet(u,1)
		a=excSet(u,2)

		if (notFitRadius(1,i,a)) cycle

		rez=F(i,a)

		rez=rez+R(a,i,i,a)*t1(i,a)

		sum=0
		if (a.EQ.iapairs(i)) then
			!$omp parallel default(shared) private(j,b) reduction(+:sum)
			!$omp do
			do j = 1,Nocc
				b = iapairs(j)
				if ((i.EQ.j).AND.(a.EQ.b)) cycle
				sum=sum+R(a,i,j,b)*t1(j,b)
			enddo	
			!$omp end parallel
		endif
		rez=rez+sum

		sum=0
		!$omp parallel default(shared) private(j,b,c) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			b=a; c=iapairs(j)
			sum=sum+R(a,b,j,c)*(t2(i,j,b,c)+t1(i,b)*t1(j,c))

			if (a.NE.iapairs(j)) then
				b=iapairs(j); c=a
				sum=sum+R(a,b,j,c)*(t2(i,j,b,c)+t1(i,b)*t1(j,c))
			endif
		enddo	
		!$omp end parallel
		rez=rez+sum

		sum=0
		!$omp parallel default(shared) private(b,j,k) reduction(+:sum)
		!$omp do
		do b = Nocc+1,N
			k=i; j=iapairs(b)
			sum=sum+R(k,b,j,i)*(t2(j,k,a,b)+t1(j,a)*t1(k,b))

			if (i.NE.iapairs(b)) then
				k=iapairs(b); j=i
				sum=sum+R(k,b,j,i)*(t2(j,k,a,b)+t1(j,a)*t1(k,b))
			endif
		enddo	
		!$omp end parallel
		rez=rez-sum

		sum=0
		!$omp parallel default(shared) private(k,j,b,c,t,ab) reduction(+:sum)
		!$omp do
		do k = 1,Nocc
			do j = 1,Nocc
				ab(1)=iapairs(k); ab(2)=iapairs(j)
				do t = 1,2
					c=ab(t); b=ab(3-t)

					sum=sum+R(k,c,j,b)*(&
		                    t1(k,c)*( 2*t2(i,j,a,b)-t2(i,j,b,a)    )&
		                   -t1(k,a)*(   t2(i,j,c,b)+t1(j,b)*t1(i,c))&
		                   -t1(i,b)*    t2(j,k,a,c)&
		                                )

					if (k.EQ.j) exit
				enddo
			enddo
		enddo
		!$omp end parallel
		rez=rez+sum

		sum=0
		!$omp parallel default(shared) private(t,j,b) reduction(+:sum)
		!$omp do
		do t = 1,NFnz
			j=Fnz(1,t); b=Fnz(2,t)
			sum=sum+F(j,b)*( t1(i,b)*t1(j,a)-2*t2(i,j,a,b)+t2(i,j,b,a) )
		enddo
		!$omp end parallel
		rez=rez-sum

		sum=0
		!$omp parallel default(shared) private(b) reduction(+:sum)
		!$omp do
		do b = Nocc+1,N    
			sum=sum+F(a,b)*t1(i,b)
		enddo
		!$omp end parallel
		rez=rez+sum

		sum=0
		!$omp parallel default(shared) private(j) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			sum=sum+F(i,j)*t1(j,a)
		enddo
		!$omp end parallel
		rez=rez-sum

		d1(i,a)=rez
	enddo

	return
	end subroutine projection_ccsd_singles_spatial_cue