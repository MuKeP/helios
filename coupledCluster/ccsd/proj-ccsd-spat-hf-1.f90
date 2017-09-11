	subroutine projection_ccsd_singles_spatial_hf

	use glob          , only: rglu,iglu,lglu
	use coupledCluster, only: N,Nocc,Ne
	use coupledCluster, only: F,R
	use coupledCluster, only: t1,t2,d1,excSet

	implicit none

	real   (kind=rglu) :: sum,rez
	integer(kind=iglu) :: i,j,a,b, c,k,mm


	do mm = 1,Ne
		i=excSet(mm,1)
		a=excSet(mm,2)

		rez=0

		rez=F(i,a) ! -> 6

		sum=0 ! -> 8
		!$omp parallel default(shared) private(j,b) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			do b = Nocc+1,N
				sum=sum+R(a,i,j,b)*t1(j,b)
			enddo
		enddo
		!$omp end parallel
		rez=rez+sum

		sum=0 ! -> 7
		!$omp parallel default(shared) private(b) reduction(+:sum)
		!$omp do
		do b = Nocc+1,N
			sum=sum+F(a,b)*t1(i,b)
		enddo
		!$omp end parallel
		rez=rez+sum

		sum=0 ! -> 7
		!$omp parallel default(shared) private(j) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			sum=sum+F(i,j)*t1(j,a)
		enddo
		!$omp end parallel
		rez=rez-sum

		sum=0 ! -> 9
		!$omp parallel default(shared) private(j,b) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			do b = Nocc+1,N
				sum=sum+F(j,b)*(2*t2(i,j,a,b)-t2(j,i,a,b))
			enddo
		enddo
		!$omp end parallel
		rez=rez+sum

		sum=0 ! -> 10
		!$omp parallel default(shared) private(j,b,c) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			do b = Nocc+1,N
				do c = Nocc+1,N
					sum=sum+R(a,b,j,c)*(t2(i,j,b,c)+t1(i,b)*t1(j,c))
				enddo
			enddo
		enddo
		!$omp end parallel
		rez=rez+sum

		sum=0 ! -> 11
		!$omp parallel default(shared) private(j,k,b) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			do k = 1,Nocc
				do b = Nocc+1,N
					sum=sum+R(k,b,j,i)*(t2(j,k,a,b)+t1(j,a)*t1(k,b))
				enddo
			enddo
		enddo
		!$omp end parallel
		rez=rez-sum

		sum=0 ! -> 12
		!$omp parallel default(shared) private(j,b) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			do b = Nocc+1,N
				sum=sum+F(j,b)*t1(i,b)*t1(j,a)
			enddo
		enddo
		!$omp end parallel
		rez=rez-sum

		sum=0 ! -> 13
		!$omp parallel default(shared) private(j,k,b,c) reduction(+:sum)
		!$omp do
		do j = 1,Nocc
			do k = 1,Nocc
				do b = Nocc+1,N
					do c = Nocc+1,N
						sum=sum+R(k,c,j,b)*(&
											 +t1(k,c)*(2*t2(i,j,a,b)-t2(j,i,a,b))&
											 -t1(i,b)*t2(j,k,a,c)&
											 -t1(k,a)*(t2(j,i,b,c)+t1(j,b)*t1(i,c))&
						                   )
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel
		rez=rez+sum

		d1(i,a)=rez
	enddo

	return
	end subroutine projection_ccsd_singles_spatial_hf