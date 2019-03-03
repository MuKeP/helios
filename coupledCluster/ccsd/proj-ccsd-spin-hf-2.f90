    subroutine projection_ccsd_doubles_spin_hf

    use glob,           only: rglu,iglu,lglu
    use coupledCluster, only: No,Nel
    use coupledCluster, only: R,F
    use coupledCluster, only: t1,t2,d2

    implicit none

    real   (kind=rglu) :: rez,sum,sum1,sum2
    integer(kind=iglu) :: i,a,j,b,k,c,l,d


    do i = 1,Nel-1
    do a = Nel+1,No-1
    do j = i+1,Nel
    do b = a+1,No
        if (btest(i+j,0).NE.btest(a+b,0)) cycle

        rez=R(i,a,j,b)

        sum=0; sum1=0
        !$omp parallel default(shared) private(c) reduction(+:sum,sum1)
        !$omp do
        do c = Nel+1,No
            sum=sum+F(b,c)*t2(i,j,a,c)-F(a,c)*t2(i,j,b,c)

            sum1=sum1+R(c,a,j,b)*t1(i,c)&
                     -R(c,a,i,b)*t1(j,c)
        enddo
        !$omp end parallel
        rez=rez+sum+sum1

        sum=0; sum1=0
        !$omp parallel default(shared) private(k) reduction(+:sum,sum1)
        !$omp do
        do k = 1,Nel
            sum=sum-F(k,j)*t2(i,k,a,b)+F(k,i)*t2(j,k,a,b)
            sum1=sum1+R(k,i,a,j)*t1(k,b)&
                     -R(k,i,b,j)*t1(k,a)
        enddo
        !$omp end parallel
        rez=rez+sum+sum1

        sum=0; sum1=0; sum2=0
        !$omp parallel default(shared) private(k,c) reduction(+:sum,sum1,sum2)
        !$omp do
        do k = 1,Nel
        do c = Nel+1,No
            sum=sum+t1(k,b)*(R(a,c,k,i)*t1(j,c)-R(a,c,k,j)*t1(i,c))&
                   +t1(k,a)*(R(b,c,k,j)*t1(i,c)-R(b,c,k,i)*t1(j,c))

            sum1=sum1+F(k,c)*(&

                                 -t2(k,j,a,b)*t1(i,c)&
                                 +t2(k,i,a,b)*t1(j,c)&
                                 -t2(i,j,a,c)*t1(k,b)&
                                 +t2(i,j,b,c)*t1(k,a)&
                                 )

            sum2=sum2-R(k,j,b,c)*t2(i,k,a,c)&
                     +R(k,j,a,c)*t2(i,k,b,c)&
                     +R(k,i,b,c)*t2(j,k,a,c)&
                     -R(k,i,a,c)*t2(j,k,b,c)

        enddo
        enddo
        !$omp end parallel
        rez=rez+sum+sum1+sum2

        sum=0
        !$omp parallel default(shared) private(k,l) reduction(+:sum)
        !$omp do
        do k = 1,Nel
        do l = 1,Nel
            sum =sum +R(k,i,l,j)*(t2(k,l,a,b)/2+t1(k,a)*t1(l,b))
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0
        !$omp parallel default(shared) private(c,d) reduction(+:sum)
        !$omp do
        do c = Nel+1,No
        do d = Nel+1,No
            sum =sum +R(a,c,b,d)*(t2(i,j,c,d)/2+t1(i,c)*t1(j,d))
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0; sum1=0
        !$omp parallel default(shared) private(k,c,d) reduction(+:sum,sum1)
        !$omp do
        do k = 1,Nel
        do c = Nel+1,No
        do d = Nel+1,No
            sum=sum+R(k,c,b,d)*(t1(j,d)*t2(i,k,a,c)-t1(i,d)*t2(j,k,a,c))&
                   +R(k,c,a,d)*(t1(i,d)*t2(j,k,b,c)-t1(j,d)*t2(i,k,b,c))


            sum1=sum1+R(b,c,k,d)*(t2(i,j,a,c)*t1(k,d)+t1(k,a)*t1(i,c)*t1(j,d)&
                                      +t1(k,a)*t2(i,j,c,d)/2)&

                     -R(a,c,k,d)*(t2(i,j,b,c)*t1(k,d)+t1(k,b)*t1(i,c)*t1(j,d)&
                                      +t1(k,b)*t2(i,j,c,d)/2)
        enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum+sum1

        sum=0
        !$omp parallel default(shared) private(k,l,c) reduction(+:sum)
        !$omp do
        do k = 1,Nel
        do l = 1,Nel
        do c = Nel+1,No
            sum=sum+R(k,i,l,c)*t1(l,c)*t2(j,k,a,b)&
                   -R(k,j,l,c)*t1(l,c)*t2(i,k,a,b)

            sum=sum+R(k,c,l,i)*(-t2(k,l,a,b)*t1(j,c)/2-t1(l,a)*t2(j,k,b,c)+t1(l,b)*(t2(j,k,a,c)-t1(j,c)*t1(k,a)))&
                   -R(k,c,l,j)*(-t2(k,l,a,b)*t1(i,c)/2-t1(l,a)*t2(i,k,b,c)+t1(l,b)*(t2(i,k,a,c)-t1(i,c)*t1(k,a)))

        enddo
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

             +t2(i,j,c,d)*t2(k,l,a,b)/4&

                    +(&

                      -t2(k,i,c,d)*t2(l,j,a,b)&
                      +t2(k,j,c,d)*t2(l,i,a,b)&
                      -t2(i,j,a,c)*t2(k,l,b,d)&
                      +t2(i,j,b,c)*t2(k,l,a,d)&
                      +t2(k,l,a,b)*t1(i,c)*t1(j,d)&
                      +t2(i,j,c,d)*t1(k,a)*t1(l,b)&
                     )/2&

             +t2(i,k,a,c)*t2(l,j,d,b)&
             -t2(j,k,a,c)*t2(l,i,d,b)&

             +t1(k,a)*(&
                         +t2(l,i,d,b)*t1(j,c)&
                         +t2(i,j,b,c)*t1(l,d)&
                         +t1(i,c)*(&
                                     +t1(j,d)*t1(l,b)&
                                     -t2(l,j,d,b)&
                                    )&
                        )&

             +t1(k,b)*(&
                         +t2(l,j,d,a)*t1(i,c)&
                         -t2(l,i,d,a)*t1(j,c)&
                         -t2(i,j,a,c)*t1(l,d)&
                        )&

             +t1(k,c)*(&
                         +t2(l,i,a,b)*t1(j,d)&
                         -t2(l,j,a,b)*t1(i,d)&
                        )&

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
    end subroutine projection_ccsd_doubles_spin_hf