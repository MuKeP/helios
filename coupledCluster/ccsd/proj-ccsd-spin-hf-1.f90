    subroutine projection_ccsd_singles_spin_hf

    use glob,           only: rglu,iglu,lglu
    use coupledCluster, only: No,Nel
    use coupledCluster, only: R,F
    use coupledCluster, only: t1,t2,d1

    implicit none

    real   (kind=rglu) :: rez,sum
    integer(kind=iglu) :: i,a,j,b,k,c


    do i = 1,Nel
    do a = Nel+1,No
        if (btest(i,0).NE.btest(a,0)) cycle

        rez=F(i,a)

        sum=0
        !$omp parallel default(shared) private(b) reduction(+:sum)
        !$omp do
        do b = Nel+1,No
            sum=sum+F(a,b)*t1(i,b)
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(j) reduction(+:sum)
        !$omp do
        do j = 1,Nel
            sum=sum+F(i,j)*t1(j,a)
        enddo
        !$omp end parallel
        rez=rez-sum

        sum=0
        !$omp parallel default(shared) private(j,b) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do b = Nel+1,No
            sum=sum+R(j,i,a,b)*t1(j,b)
        enddo
        enddo
        !$omp end parallel
        rez=rez-sum

        sum=0
        !$omp parallel default(shared) private(j,b) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do b = Nel+1,No
            sum=sum+F(j,b)*t2(i,j,a,b)
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(j,b,c) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do b = Nel+1,No
        do c = Nel+1,No
            sum=sum+R(a,b,j,c)*t2(i,j,b,c)
        enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum/2

        sum=0
        !$omp parallel default(shared) private(j,k,b) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do k = 1,Nel
        do b = Nel+1,No
            sum=sum+R(j,i,k,b)*t2(j,k,a,b)
        enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez-sum/2

        sum=0
        !$omp parallel default(shared) private(j,b) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do b = Nel+1,No
            sum=sum+F(j,b)*t1(i,b)*t1(j,a)
        enddo
        enddo
        !$omp end parallel
        rez=rez-sum

        sum=0
        !$omp parallel default(shared) private(j,b,c) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do b = Nel+1,No
        do c = Nel+1,No
            sum=sum+R(a,b,j,c)*t1(i,b)*t1(j,c)
        enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(j,k,b) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do k = 1,Nel
        do b = Nel+1,No
            sum=sum+R(j,i,k,b)*t1(j,a)*t1(k,b)
        enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez-sum

        sum=0
        !$omp parallel default(shared) private(j,k,b,c) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do k = 1,Nel
        do b = Nel+1,No
        do c = Nel+1,No
            sum=sum+ R(j,b,k,c)*(&

                  &+t1(k,c)*t2(i,j,a,b)&
                  &-t1(j,b)*t1(i,c)*t1(k,a)&
                 &-(t1(i,b)*t2(j,k,a,c)&
                   +t1(k,a)*t2(j,i,b,c))/2 &

            &)
        enddo
        enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        d1(i,a)=rez
    enddo
    enddo

    return
    end subroutine projection_ccsd_singles_spin_hf