    subroutine projection_ccsd_doubles_spatial_hf

    use glob,           only: rglu,iglu,lglu
    use coupledCluster, only: N,Nocc,Ne
    use coupledCluster, only: F,R
    use coupledCluster, only: t1,t2,d2,excSet

    implicit none

    real   (kind=rglu) :: sum,rez
    integer(kind=iglu) :: i,j,a,b, c,d,k,l,mm,nn


    do mm = 1,Ne
        i=excSet(mm,1)
        a=excSet(mm,2)

        do nn = mm,Ne

            j=excSet(nn,1)
            b=excSet(nn,2)

            rez=0

            rez=rez+R(a,i,b,j) ! -> 14

            sum=0 ! -> 25
            !$omp parallel default(shared) private(c) reduction(+:sum)
            !$omp do
            do c = Nocc+1,N
                sum=sum+R(b,j,a,c)*t1(i,c)&
                       +R(a,i,b,c)*t1(j,c)
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 26
            !$omp parallel default(shared) private(k) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                sum=sum+R(k,i,b,j)*t1(k,a)&
                       +R(k,j,a,i)*t1(k,b)
            enddo
            !$omp end parallel
            rez=rez-sum

            sum=0 ! -> 17
            !$omp parallel default(shared) private(k,l) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                do l = 1,Nocc
                    sum=sum+R(k,i,l,j)*(t2(k,l,a,b)+t1(k,a)*t1(l,b))
                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 18
            !$omp parallel default(shared) private(c,d) reduction(+:sum)
            !$omp do
            do c = Nocc+1,N
                do d = Nocc+1,N
                    sum=sum+R(b,d,a,c)*( t2(i,j,c,d)+t1(i,c)*t1(j,d) )
                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 19,20,21,22,27,28
            !$omp parallel default(shared) private(k,c) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                do c = Nocc+1,N

                    sum=sum&

                       +2*(&
                           +t2(i,k,a,c)*R(b,j,k,c)&
                             +t2(j,k,b,c)*R(a,i,k,c)&

                             +(&
                              -t2(i,k,c,a)*R(k,c,b,j)&
                                -t2(j,k,c,b)*R(k,c,a,i)&
                                -t2(j,k,a,c)*R(k,c,b,i)&
                                -t2(i,k,b,c)*R(k,c,a,j)&
                                -t2(k,j,a,c)*R(k,i,b,c)&
                                -t2(k,i,b,c)*R(k,j,a,c)&
                             )/2 &
                           )

                    sum=sum+&

                            (&
                           -R(a,c,k,j)*t1(i,c)*t1(k,b)&
                             -R(b,c,k,i)*t1(j,c)*t1(k,a)&
                             -R(k,c,b,j)*t1(i,c)*t1(k,a)&
                             -R(k,c,a,i)*t1(j,c)*t1(k,b)&
                           )

                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 33,36,38
            !$omp parallel default(shared) private(k,c,d) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                do c = Nocc+1,N
                    do d = Nocc+1,N

                    sum=sum+R(k,c,b,d)*&
                                         (t1(j,d)*( 2*t2(i,k,a,c)-t2(k,i,a,c) )-&
                                                      t2(j,k,a,c)*t1(i,d)-&
                                                      t2(k,j,a,d)*t1(i,c) )&

                           +R(k,c,a,d)*&
                                        (t1(i,d)*( 2*t2(j,k,b,c)-t2(k,j,b,c) )-&
                                                      t2(i,k,b,c)*t1(j,d)-&
                                                      t2(k,i,b,d)*t1(j,c) )

                    sum=sum+t1(k,d)*(R(k,d,b,c)*( 2*t2(i,j,a,c)-t2(i,j,c,a) )&
                           +         R(k,d,a,c)*( 2*t2(j,i,b,c)-t2(j,i,c,b) ) )

                    sum=sum-R(a,c,k,d)*t1(k,b)*( t2(i,j,c,d)+t1(i,c)*t1(j,d) )&
                            -R(b,c,k,d)*t1(k,a)*( t2(j,i,c,d)+t1(j,c)*t1(i,d) )

                    enddo
                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 34,35,37
            !$omp parallel default(shared) private(k,l,c) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                do c = Nocc+1,N
                    do l = 1,Nocc

                    sum=sum+R(k,c,l,j)*&
                                         (t1(l,b)*(-2*t2(i,k,a,c)+t2(k,i,a,c) )+&
                                             t2(i,k,b,c)*t1(l,a)+&
                                             t2(i,l,c,b)*t1(k,a) )&

                            +R(k,c,l,i)*&
                                         (t1(l,a)*(-2*t2(j,k,b,c)+t2(k,j,b,c) )+&
                                             t2(j,k,a,c)*t1(l,b)+&
                                             t2(j,l,c,a)*t1(k,b) )

                    sum=sum+t1(l,c)*(R(k,j,l,c)*( -2*t2(i,k,a,b)+t2(k,i,a,b) )&
                            +         R(k,i,l,c)*( -2*t2(j,k,b,a)+t2(k,j,b,a) ) )

                    sum=sum+R(k,c,l,j)*t1(i,c)*( t2(k,l,a,b)+t1(k,a)*t1(l,b) )&
                            +R(k,c,l,i)*t1(j,c)*( t2(k,l,b,a)+t1(k,b)*t1(l,a) )

                    enddo
                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 23,24,39
            !$omp parallel default(shared) private(k,l,c,d) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                do c = Nocc+1,N
                    do l = 1,Nocc
                        do d = Nocc+1,N

                        sum = sum + R(k,c,l,d)*(&

                           -2*(t2(l,j,d,b)*t1(i,c)*t1(k,a)+t2(l,i,d,a)*t1(j,c)*t1(k,b)&
                              +t2(l,j,a,b)*t1(k,c)*t1(i,d)+t2(l,i,b,a)*t1(k,c)*t1(j,d)&
                              +t2(i,j,a,c)*t1(k,b)*t1(l,d)+t2(j,i,b,c)*t1(k,a)*t1(l,d)&
                              +t2(k,i,c,d)*t2(l,j,a,b)+t2(k,j,c,d)*t2(l,i,b,a)&
                              +t2(i,j,a,c)*t2(k,l,b,d)+t2(j,i,b,c)*t2(k,l,a,d)&
                              +t2(i,k,c,a)*t2(l,j,d,b)+t2(j,k,c,b)*t2(l,i,d,a)&
                              -t2(i,k,a,c)*t2(l,j,d,b)-t2(j,k,b,c)*t2(l,i,d,a)&
                              )&

                              +t2(j,l,d,b)*t1(i,c)*t1(k,a)+t2(i,l,d,a)*t1(j,c)*t1(k,b)&
                              +t2(l,i,d,b)*t1(j,c)*t1(k,a)+t2(l,j,d,a)*t1(i,c)*t1(k,b)&
                              +t2(i,l,c,b)*t1(j,d)*t1(k,a)+t2(j,l,c,a)*t1(i,d)*t1(k,b)&
                              +t2(l,i,a,b)*t1(k,c)*t1(j,d)+t2(l,j,b,a)*t1(k,c)*t1(i,d)&
                               +t2(i,j,c,a)*t1(k,b)*t1(l,d)+t2(j,i,c,b)*t1(k,a)*t1(l,d)&

                              +t2(k,l,a,b)*t1(i,c)*t1(j,d)&
                              +t2(i,j,c,d)*(t2(k,l,a,b)+t1(k,a)*t1(l,b))&
                              +t1(i,c)*t1(k,a)*t1(j,d)*t1(l,b)&

                              +t2(k,j,c,d)*t2(l,i,a,b)+t2(k,i,c,d)*t2(l,j,b,a)&
                              +t2(i,j,b,c)*t2(k,l,a,d)+t2(j,i,a,c)*t2(k,l,b,d)&
                              -t2(l,i,d,b)*(t2(j,k,a,c)-t2(j,k,c,a))&
                              -t2(l,j,d,a)*(t2(i,k,b,c)-t2(i,k,c,b))&
                             +(t2(i,k,c,a)*t2(j,l,d,b)+t2(j,k,c,b)*t2(i,l,d,a)&
                             + t2(j,k,d,a)*t2(i,l,c,b)+t2(i,k,d,b)*t2(j,l,c,a))/2 &

                         )

                        enddo
                    enddo
                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 15
            !$omp parallel default(shared) private(c) reduction(+:sum)
            !$omp do
            do c = Nocc+1,N
                sum=sum+F(c,a)*(2*t2(i,j,c,b)-t2(j,i,c,b))&
                       +F(c,b)*(2*t2(j,i,c,a)-t2(i,j,c,a))
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0 ! -> 16
            !$omp parallel default(shared) private(k) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                sum=sum+F(k,i)*(2*t2(k,j,a,b)-t2(k,j,b,a))&
                       +F(k,j)*(2*t2(k,i,b,a)-t2(k,i,a,b))
            enddo
            !$omp end parallel
            rez=rez-sum

            sum=0 ! -> 29,31  30,32
            !$omp parallel default(shared) private(k,c) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                do c = Nocc+1,N
                    sum=sum+F(k,c)*(&

                                     -2*(&
                                         +t2(k,j,a,b)*t1(i,c)&
                                         +t2(k,i,b,a)*t1(j,c)&
                                         +t2(i,j,a,c)*t1(k,b)&
                                         +t2(j,i,b,c)*t1(k,a)&
                                         )&

                                         +t2(k,i,a,b)*t1(j,c)&
                                         +t2(k,j,b,a)*t1(i,c)&
                                         +t2(j,i,a,c)*t1(k,b)&
                                         +t2(i,j,b,c)*t1(k,a)&
                                    )

                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            d2(i,j,a,b)=rez; d2(j,i,b,a)=rez

        enddo
    enddo

    return
    end subroutine projection_ccsd_doubles_spatial_hf