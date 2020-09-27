    subroutine projection_ccd_doubles_spin_hf

    use glob,           only: rglu,iglu,lglu
    use coupledCluster, only: No,Nel
    use coupledCluster, only: R,F
    use coupledCluster, only: t2,d2

    implicit none

    real   (kind=rglu) :: rez,sum
    integer(kind=iglu) :: i,a,j,b,k,c,l,d


    do i = 1,Nel-1
    do a = Nel+1,No-1
    do j = i+1,Nel
    do b = a+1,No
        if (btest(i+j,0).NE.btest(a+b,0)) cycle

        rez=R(i,a,j,b)

        sum=0
        !$omp parallel default(shared) private(c) reduction(+:sum)
        !$omp do
        do c = Nel+1,No
            sum=sum+F(b,c)*t2(i,j,a,c)-F(a,c)*t2(i,j,b,c)
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(k) reduction(+:sum)
        !$omp do
        do k = 1,Nel
            sum=sum-F(k,j)*t2(i,k,a,b)+F(k,i)*t2(j,k,a,b)
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(k,c) reduction(+:sum)
        !$omp do
        do k = 1,Nel
        do c = Nel+1,No
            sum=sum-R(k,j,b,c)*t2(i,k,a,c)&
                   +R(k,j,a,c)*t2(i,k,b,c)&
                   +R(k,i,b,c)*t2(j,k,a,c)&
                   -R(k,i,a,c)*t2(j,k,b,c)

        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(k,l) reduction(+:sum)
        !$omp do
        do k = 1,Nel
        do l = 1,Nel
            sum =sum +R(k,i,l,j)*t2(k,l,a,b)/2
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(c,d) reduction(+:sum)
        !$omp do
        do c = Nel+1,No
        do d = Nel+1,No
            sum =sum +R(a,c,b,d)*t2(i,j,c,d)/2
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(k,l,c,d) reduction(+:sum)
        !$omp do
        do k = 1,Nel
        do l = 1,Nel
        do c = Nel+1,No
        do d = Nel+1,No
            sum=sum+R(k,c,l,d)*(&

             +t2(i,j,c,d)*t2(k,l,a,b)/4 &

                    +(&

                      -t2(k,i,c,d)*t2(l,j,a,b)&
                      +t2(k,j,c,d)*t2(l,i,a,b)&
                      -t2(i,j,a,c)*t2(k,l,b,d)&
                      +t2(i,j,b,c)*t2(k,l,a,d)&
                     )/2 &

             +t2(i,k,a,c)*t2(l,j,d,b)&
             -t2(j,k,a,c)*t2(l,i,d,b)&
             )
        enddo
        enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        d2(i,j,a,b)=rez
    enddo
    enddo
    enddo
    enddo

    return
    end subroutine projection_ccd_doubles_spin_hf