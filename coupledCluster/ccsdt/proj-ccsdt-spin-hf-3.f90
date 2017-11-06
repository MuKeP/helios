    subroutine projection_ccsdt_triples_spin_hf

    use glob          , only: rglu,iglu,lglu
    use coupledCluster, only: No,Nel
    use coupledCluster, only: R,F
    use coupledCluster, only: t1,t2,t3,d3

    implicit none

    real   (kind=rglu) :: rez,sum,sum1,Ax,isum(12)
    integer(kind=iglu) :: i,a,j,b,k,c,l,d,m,e


    do i = 1,Nel-2
    do a = Nel+1,No-2
        do j = i+1,Nel-1
        do b = a+1,No-1
            do k = j+1,Nel
            do c = b+1,No
                if (btest(i+j+k,0).NE.btest(a+b+c,0)) cycle

                rez=0

                sum=0; sum1=0
                !$omp parallel default(shared) private(d) reduction(+:sum,sum1)
                !$omp do
                do d = Nel+1,No
                    sum=sum+ R(b,k,c,d)*t2(i,j,a,d)& !<rst>
                           - R(b,i,c,d)*t2(k,j,a,d)&
                           - R(b,j,c,d)*t2(i,k,a,d)&
                           - R(a,k,c,d)*t2(i,j,b,d)&
                           + R(a,i,c,d)*t2(k,j,b,d)&
                           + R(a,j,c,d)*t2(i,k,b,d)&
                           - R(b,k,a,d)*t2(i,j,c,d)&
                           + R(b,i,a,d)*t2(k,j,c,d)&
                           + R(b,j,a,d)*t2(i,k,c,d)

                    sum1=sum1+ F(c,d)*t3(i,j,k,a,b,d)&
                             - F(a,d)*t3(i,j,k,c,b,d)&
                             - F(b,d)*t3(i,j,k,a,c,d)
                enddo
                !$omp end parallel
                rez=rez-sum+sum1

                sum=0; sum1=0
                !$omp parallel default(shared) private(l) reduction(+:sum,sum1)
                !$omp do
                do l = 1,Nel
                    sum=sum+ R(j,c,k,l)*t2(i,l,a,b)& !<rst>
                           - R(i,c,k,l)*t2(j,l,a,b)&
                           - R(j,c,i,l)*t2(k,l,a,b)&
                           - R(j,a,k,l)*t2(i,l,c,b)&
                           + R(i,a,k,l)*t2(j,l,c,b)&
                           + R(j,a,i,l)*t2(k,l,c,b)&
                           - R(j,b,k,l)*t2(i,l,a,c)&
                           + R(i,b,k,l)*t2(j,l,a,c)&
                           + R(j,b,i,l)*t2(k,l,a,c)


                    sum1=sum1+ F(l,k)*t3(i,j,l,a,b,c)&
                             - F(l,i)*t3(k,j,l,a,b,c)&
                             - F(l,j)*t3(i,k,l,a,b,c)
                enddo
                !$omp end parallel
                rez=rez+sum-sum1

                sum=0; sum1=0
                !$omp parallel default(shared) private(l,d) reduction(+:sum,sum1)
                !$omp do
                do l = 1,Nel
                do d = Nel+1,No
                    sum=sum&

                           +R(i,a,d,l)*(& !<rs>
                                             +t3(k,j,l,c,b,d)&
                                             +t2(k,j,b,d)*t1(l,c)&
                                             -t2(k,j,c,d)*t1(l,b)&
                                             +t2(j,l,c,b)*t1(k,d)&
                                             -t2(k,l,c,b)*t1(j,d)&
                                            )&

                           +R(j,a,d,l)*(& !<rs>
                                             +t3(i,k,l,c,b,d)&
                                             -t2(i,l,c,b)*t1(k,d)&
                                             +t2(k,l,c,b)*t1(i,d)&
                                             +t2(i,k,b,d)*t1(l,c)&
                                             -t2(i,k,c,d)*t1(l,b)&
                                            )&

                           +R(k,a,d,l)*(& !<rs>
                                             -t3(i,j,l,c,b,d)&
                                             +t2(i,l,c,b)*t1(j,d)&
                                             -t2(j,l,c,b)*t1(i,d)&
                                             -t2(i,j,b,d)*t1(l,c)&
                                             +t2(i,j,c,d)*t1(l,b)&
                                            )&

                           +R(i,b,d,l)*(& !<rs>
                                             +t3(k,j,l,a,c,d)&
                                             +t2(j,l,a,c)*t1(k,d)&
                                             -t2(k,l,a,c)*t1(j,d)&
                                             -t2(k,j,a,d)*t1(l,c)&
                                             +t2(k,j,c,d)*t1(l,a)&
                                            )&

                           +R(j,b,d,l)*(& !<rs>
                                             +t3(i,k,l,a,c,d)&
                                             -t2(i,l,a,c)*t1(k,d)&
                                             +t2(k,l,a,c)*t1(i,d)&
                                             -t2(i,k,a,d)*t1(l,c)&
                                             +t2(i,k,c,d)*t1(l,a)&
                                            )&

                           +R(k,b,d,l)*(& !<rs>
                                             -t3(i,j,l,a,c,d)&
                                             +t2(i,j,a,d)*t1(l,c)&
                                             -t2(i,j,c,d)*t1(l,a)&
                                             +t2(i,l,a,c)*t1(j,d)&
                                             -t2(j,l,a,c)*t1(i,d)&
                                            )&

                           +R(i,c,d,l)*(& !<rs>
                                             -t3(k,j,l,a,b,d)&
                                             +t2(k,j,a,d)*t1(l,b)&
                                             -t2(k,j,b,d)*t1(l,a)&
                                             -t2(j,l,a,b)*t1(k,d)&
                                             +t2(k,l,a,b)*t1(j,d)&
                                            )&

                           +R(j,c,d,l)*(& !<rs>
                                             -t3(i,k,l,a,b,d)&
                                             +t2(i,k,a,d)*t1(l,b)&
                                             -t2(i,k,b,d)*t1(l,a)&
                                             +t2(i,l,a,b)*t1(k,d)&
                                             -t2(k,l,a,b)*t1(i,d)&
                                            )&

                           +R(k,c,d,l)*(& !<rs>
                                             +t3(i,j,l,a,b,d)&
                                             -t2(i,j,a,d)*t1(l,b)&
                                             +t2(i,j,b,d)*t1(l,a)&
                                             -t2(i,l,a,b)*t1(j,d)&
                                             +t2(j,l,a,b)*t1(i,d)&
                                            )

                    sum1=sum1&
                           +F(l,d)*(&
                                         +t3(i,j,k,a,b,d)*t1(l,c)&
                                         -t3(i,j,k,c,b,d)*t1(l,a)&
                                         -t3(i,j,k,a,c,d)*t1(l,b)&
                                         +t3(l,j,k,a,b,c)*t1(i,d)&
                                         -t3(l,i,k,a,b,c)*t1(j,d)&
                                         -t3(l,j,i,a,b,c)*t1(k,d)&
                                         +t2(i,j,a,d)*t2(l,k,b,c)&
                                         -t2(i,j,b,d)*t2(l,k,a,c)&
                                         -t2(i,j,c,d)*t2(l,k,b,a)&
                                         -t2(k,j,a,d)*t2(l,i,b,c)&
                                         +t2(k,j,b,d)*t2(l,i,a,c)&
                                         +t2(k,j,c,d)*t2(l,i,b,a)&
                                         -t2(i,k,a,d)*t2(l,j,b,c)&
                                         +t2(i,k,b,d)*t2(l,j,a,c)&
                                         +t2(i,k,c,d)*t2(l,j,b,a)&
                                       )

                enddo
                enddo
                !$omp end parallel
                rez=rez+sum-sum1

                sum=0
                !$omp parallel default(shared) private(l,m) reduction(+:sum)
                !$omp do
                do l = 1,Nel
                do m = 1,Nel
                    sum=sum+R(j,l,k,m)*(& !<rt>
                                             +t2(i,l,a,b)*t1(m,c)&
                                             -t2(i,l,c,b)*t1(m,a)&
                                             -t2(i,l,a,c)*t1(m,b)&
                                             +t3(i,l,m,a,b,c)/2 &
                                            )&

                           +R(i,l,k,m)*(& !<rt>
                                             -t2(j,l,a,b)*t1(m,c)&
                                             +t2(j,l,c,b)*t1(m,a)&
                                             +t2(j,l,a,c)*t1(m,b)&
                                             -t3(j,l,m,a,b,c)/2 &
                                            )&

                           +R(j,l,i,m)*(& !<rt>
                                             -t2(k,l,a,b)*t1(m,c)&
                                             +t2(k,l,c,b)*t1(m,a)&
                                             +t2(k,l,a,c)*t1(m,b)&
                                             -t3(k,l,m,a,b,c)/2 &
                                            )
                enddo
                enddo
                !$omp end parallel
                rez=rez+sum

                sum=0
                !$omp parallel default(shared) private(d,e) reduction(+:sum)
                !$omp do
                do d = Nel+1,No
                do e = Nel+1,No
                    sum=sum+R(b,d,c,e)*(& !<rt>
                                             +t2(i,j,a,d)*t1(k,e)&
                                             -t2(i,k,a,d)*t1(j,e)&
                                             -t2(k,j,a,d)*t1(i,e)&
                                             +t3(i,j,k,a,d,e)/2 &
                                            )&

                           +R(a,d,c,e)*(& !<rt>
                                             -t2(i,j,b,d)*t1(k,e)&
                                             +t2(k,j,b,d)*t1(i,e)&
                                             +t2(i,k,b,d)*t1(j,e)&
                                             -t3(i,j,k,b,d,e)/2 &
                                            )&

                           +R(b,d,a,e)*(& !<rt>
                                             -t2(i,j,c,d)*t1(k,e)&
                                             +t2(k,j,c,d)*t1(i,e)&
                                             +t2(i,k,c,d)*t1(j,e)&
                                             -t3(i,j,k,c,d,e)/2 &
                                            )
                enddo
                enddo
                !$omp end parallel
                rez=rez+sum

                sum=0
                !$omp parallel default(shared) private(l,d,e) reduction(+:sum)
                !$omp do
                do l = 1,Nel
                do d = Nel+1,No
                do e = Nel+1,No
                    sum=sum&

                           +R(a,d,l,e)*(& !<r>
                                             -t3(i,j,k,c,b,d)*t1(l,e)&
                                             -t2(i,j,b,d)*t2(l,k,e,c)&
                                             +t2(i,j,c,d)*t2(l,k,e,b)&
                                             +t2(k,j,b,d)*t2(l,i,e,c)&
                                             -t2(k,j,c,d)*t2(l,i,e,b)&
                                             +t2(i,k,b,d)*t2(l,j,e,c)&
                                             -t2(i,k,c,d)*t2(l,j,e,b)&

                                             +t1(j,e)*(&
                                                         +t2(i,k,c,d)*t1(l,b)&
                                                         -t2(i,k,b,d)*t1(l,c)&
                                                         -t2(l,k,b,c)*t1(i,d)&
                                                         +t2(l,i,b,c)*t1(k,d)&
                                                         -t3(i,k,l,c,b,d)&
                                                        )&

                                             +t1(k,e)*(&
                                                         +t2(i,j,b,d)*t1(l,c)&
                                                         -t2(i,j,c,d)*t1(l,b)&
                                                         +t2(l,j,b,c)*t1(i,d)&
                                                         +t3(i,j,l,c,b,d)&
                                                        )&

                                             +t1(i,e)*(&
                                                         +t2(k,j,c,d)*t1(l,b)&
                                                         -t2(k,j,b,d)*t1(l,c)&
                                                         -t3(k,j,l,c,b,d)&
                                                        )&

                                             +(&
                                                +t3(i,j,k,b,d,e)*t1(l,c)&
                                                -t3(i,j,k,c,d,e)*t1(l,b)&
                                                -t2(i,j,d,e)*t2(l,k,b,c)&
                                                +t2(k,j,d,e)*t2(l,i,b,c)&
                                                +t2(i,k,d,e)*t2(l,j,b,c)&
                                               )/2 &
                                            )&

                           +R(b,d,l,e)*(& !<r>

                                             -t3(i,j,k,a,c,d)*t1(l,e)&
                                             +t2(i,j,a,d)*t2(l,k,e,c)&
                                             -t2(i,j,c,d)*t2(l,k,e,a)&
                                             -t2(k,j,a,d)*t2(l,i,e,c)&
                                             +t2(k,j,c,d)*t2(l,i,e,a)&
                                             -t2(i,k,a,d)*t2(l,j,e,c)&
                                             +t2(i,k,c,d)*t2(l,j,e,a)&

                                             +t1(j,e)*(&
                                                         +t2(i,k,a,d)*t1(l,c)&
                                                         -t2(i,k,c,d)*t1(l,a)&
                                                         +t2(l,k,a,c)*t1(i,d)&
                                                         -t2(l,i,a,c)*t1(k,d)&
                                                         -t3(i,k,l,a,c,d)&
                                                        )&

                                             +t1(k,e)*(&
                                                         +t2(i,j,c,d)*t1(l,a)&
                                                         -t2(i,j,a,d)*t1(l,c)&
                                                         -t2(l,j,a,c)*t1(i,d)&
                                                         +t3(i,j,l,a,c,d)&
                                                        )&

                                             +t1(i,e)*(&
                                                         +t2(k,j,a,d)*t1(l,c)&
                                                         -t2(k,j,c,d)*t1(l,a)&
                                                         -t3(k,j,l,a,c,d)&
                                                        )&

                                             +(&
                                                +t2(i,j,d,e)*t2(l,k,a,c)&
                                                -t2(k,j,d,e)*t2(l,i,a,c)&
                                                -t2(i,k,d,e)*t2(l,j,a,c)&
                                                -t3(i,j,k,a,d,e)*t1(l,c)&
                                                +t3(i,j,k,c,d,e)*t1(l,a)&
                                                )/2 &
                                            )&

                           +R(c,d,l,e)*(& !<r>
                                             +t3(i,j,k,a,b,d)*t1(l,e)&
                                             -t2(i,j,a,d)*t2(l,k,e,b)&
                                             +t2(i,j,b,d)*t2(l,k,e,a)&
                                             +t2(k,j,a,d)*t2(l,i,e,b)&
                                             -t2(k,j,b,d)*t2(l,i,e,a)&
                                             +t2(i,k,a,d)*t2(l,j,e,b)&
                                             -t2(i,k,b,d)*t2(l,j,e,a)&

                                             +t1(j,e)*(&
                                                         +t2(i,k,b,d)*t1(l,a)&
                                                         -t2(i,k,a,d)*t1(l,b)&
                                                         +t2(l,k,b,a)*t1(i,d)&
                                                         -t2(l,i,b,a)*t1(k,d)&
                                                         +t3(i,k,l,a,b,d)&
                                                        )&

                                             +t1(k,e)*(&
                                                         +t2(i,j,a,d)*t1(l,b)&
                                                         -t2(i,j,b,d)*t1(l,a)&
                                                         -t2(l,j,b,a)*t1(i,d)&
                                                         -t3(i,j,l,a,b,d)&
                                                        )&

                                             +t1(i,e)*(&
                                                         +t2(k,j,b,d)*t1(l,a)&
                                                         -t2(k,j,a,d)*t1(l,b)&
                                                         +t3(k,j,l,a,b,d)&
                                                        )&

                                             +(&
                                                +t2(i,j,d,e)*t2(l,k,b,a)&
                                                -t2(k,j,d,e)*t2(l,i,b,a)&
                                                -t2(i,k,d,e)*t2(l,j,b,a)&
                                                +t3(i,j,k,a,d,e)*t1(l,b)&
                                                -t3(i,j,k,b,d,e)*t1(l,a)&
                                               )/2 &
                                            )

                enddo
                enddo
                enddo
                !$omp end parallel
                rez=rez+sum

                sum=0
                !$omp parallel default(shared) private(l,m,d) reduction(+:sum)
                !$omp do
                do l = 1,Nel
                do m = 1,Nel
                do d = Nel+1,No
                    sum=sum&
                           +R(k,l,d,m)*(& !<r>
                                             -t3(i,j,l,a,b,c)*t1(m,d)&
                                             +t2(i,l,a,b)*t2(m,j,d,c)&
                                             -t2(j,l,a,b)*t2(m,i,d,c)&
                                             -t2(i,l,c,b)*t2(m,j,d,a)&
                                             +t2(j,l,c,b)*t2(m,i,d,a)&
                                             -t2(i,l,a,c)*t2(m,j,d,b)&
                                             +t2(j,l,a,c)*t2(m,i,d,b)&

                                             +t1(m,c)*(&
                                                         -t2(i,l,a,b)*t1(j,d)&
                                                         +t2(j,l,a,b)*t1(i,d)&
                                                         -t2(i,j,a,d)*t1(l,b)&
                                                         +t2(i,j,b,d)*t1(l,a)&
                                                         +t3(i,j,l,a,b,d)&
                                                        )&

                                             +t1(m,a)*(&
                                                         +t2(i,l,c,b)*t1(j,d)&
                                                         -t2(j,l,c,b)*t1(i,d)&
                                                         +t2(i,j,c,d)*t1(l,b)&
                                                         -t3(i,j,l,c,b,d)&
                                                        )&

                                             +t1(m,b)*(&
                                                         +t2(i,l,a,c)*t1(j,d)&
                                                         -t2(j,l,a,c)*t1(i,d)&
                                                         -t3(i,j,l,a,c,d)&
                                                        )&

                                             +(&
                                                +t2(i,j,c,d)*t2(l,m,b,a)&
                                                +t2(i,j,b,d)*t2(l,m,a,c)&
                                                -t2(i,j,a,d)*t2(l,m,b,c)&
                                                -t3(l,m,i,a,b,c)*t1(j,d)&
                                                +t3(l,m,j,a,b,c)*t1(i,d)&
                                               )/2 &
                                            )&


                           +R(j,l,d,m)*(& !<r>
                                             +t3(i,k,l,a,b,c)*t1(m,d)&
                                             -t2(i,l,a,b)*t2(m,k,d,c)&
                                             +t2(k,l,a,b)*t2(m,i,d,c)&
                                             +t2(i,l,c,b)*t2(m,k,d,a)&
                                             -t2(k,l,c,b)*t2(m,i,d,a)&
                                             +t2(i,l,a,c)*t2(m,k,d,b)&
                                             -t2(k,l,a,c)*t2(m,i,d,b)&

                                             +t1(m,c)*(&
                                                         +t2(i,l,a,b)*t1(k,d)&
                                                         -t2(k,l,a,b)*t1(i,d)&
                                                         +t2(i,k,a,d)*t1(l,b)&
                                                         -t2(i,k,b,d)*t1(l,a)&
                                                         -t3(i,k,l,a,b,d)&
                                                        )&

                                             +t1(m,a)*(&
                                                         +t2(k,l,c,b)*t1(i,d)&
                                                         -t2(i,k,c,d)*t1(l,b)&
                                                         -t2(i,l,c,b)*t1(k,d)&
                                                         +t3(i,k,l,c,b,d)&
                                                        )&

                                             +t1(m,b)*(&
                                                         +t2(k,l,a,c)*t1(i,d)&
                                                         -t2(i,l,a,c)*t1(k,d)&
                                                         +t3(i,k,l,a,c,d)&
                                                        )&

                                             +(&
                                                +t3(l,m,i,a,b,c)*t1(k,d)&
                                                -t3(l,m,k,a,b,c)*t1(i,d)&
                                                +t2(i,k,a,d)*t2(l,m,b,c)&
                                                -t2(i,k,b,d)*t2(l,m,a,c)&
                                                -t2(i,k,c,d)*t2(l,m,b,a)&
                                               )/2 &
                                            )&


                           +R(i,l,d,m)*(& !<r>
                                             +t3(k,j,l,a,b,c)*t1(m,d)&
                                             +t2(j,l,a,b)*t2(m,k,d,c)&
                                             -t2(k,l,a,b)*t2(m,j,d,c)&
                                             -t2(j,l,c,b)*t2(m,k,d,a)&
                                             +t2(k,l,c,b)*t2(m,j,d,a)&
                                             -t2(j,l,a,c)*t2(m,k,d,b)&
                                             +t2(k,l,a,c)*t2(m,j,d,b)&

                                             +t1(m,c)*(&
                                                         +t2(k,l,a,b)*t1(j,d)&
                                                         +t2(k,j,a,d)*t1(l,b)&
                                                         -t2(k,j,b,d)*t1(l,a)&
                                                         -t2(j,l,a,b)*t1(k,d)&
                                                         -t3(k,j,l,a,b,d)&
                                                        )&

                                             +t1(m,a)*(&
                                                         +t2(j,l,c,b)*t1(k,d)&
                                                         -t2(k,l,c,b)*t1(j,d)&
                                                         -t2(k,j,c,d)*t1(l,b)&
                                                         +t3(k,j,l,c,b,d)&
                                                        )&

                                             +t1(m,b)*(&
                                                         +t2(j,l,a,c)*t1(k,d)&
                                                         -t2(k,l,a,c)*t1(j,d)&
                                                         +t3(k,j,l,a,c,d)&
                                                        )&

                                             +(&
                                                +t3(l,m,k,a,b,c)*t1(j,d)&
                                                -t3(l,m,j,a,b,c)*t1(k,d)&
                                                +t2(k,j,a,d)*t2(l,m,b,c)&
                                                -t2(k,j,b,d)*t2(l,m,a,c)&
                                                -t2(k,j,c,d)*t2(l,m,b,a)&
                                              )/2 &
                                            )
                enddo
                enddo
                enddo
                !$omp end parallel
                rez=rez+sum

                sum=0
                !$omp parallel default(shared) private(l,m,d,e) reduction(+:sum)
                !$omp do
                do l = 1,Nel
                do m = 1,Nel
                do d = Nel+1,No
                do e = Nel+1,No
                    Ax=R(l,d,m,e)

                    isum(1)=&
                         +        (&
                                   +t3(i,j,k,a,d,e)*t2(l,m,b,c)&
                                   -t3(i,j,k,b,d,e)*t2(l,m,a,c)&
                                   -t3(i,j,k,c,d,e)*t2(l,m,b,a)&
                                   +t3(l,m,k,a,b,c)*t2(i,j,d,e)&
                                   -t3(l,m,i,a,b,c)*t2(k,j,d,e)&
                                   -t3(l,m,j,a,b,c)*t2(i,k,d,e)&
                                  )/4

                    isum(2)=&
                         +       (&
                                  -t3(i,l,j,a,d,e)*t2(m,k,b,c)&
                                  +t3(k,l,j,a,d,e)*t2(m,i,b,c)&
                                  +t3(i,l,k,a,d,e)*t2(m,j,b,c)&
                                  +t3(i,l,j,b,d,e)*t2(m,k,a,c)&
                                  -t3(k,l,j,b,d,e)*t2(m,i,a,c)&
                                  -t3(i,l,k,b,d,e)*t2(m,j,a,c)&
                                  +t3(i,l,j,c,d,e)*t2(m,k,b,a)&
                                  -t3(k,l,j,c,d,e)*t2(m,i,b,a)&
                                  -t3(i,l,k,c,d,e)*t2(m,j,b,a)&
                                  -t3(l,m,k,b,e,c)*t2(i,j,a,d)&
                                  +t3(l,m,i,b,e,c)*t2(k,j,a,d)&
                                  +t3(l,m,j,b,e,c)*t2(i,k,a,d)&
                                  +t3(l,m,k,a,e,c)*t2(i,j,b,d)&
                                  -t3(l,m,i,a,e,c)*t2(k,j,b,d)&
                                  -t3(l,m,j,a,e,c)*t2(i,k,b,d)&
                                  +t3(l,m,k,b,e,a)*t2(i,j,c,d)&
                                  -t3(l,m,i,b,e,a)*t2(k,j,c,d)&
                                  -t3(l,m,j,b,e,a)*t2(i,k,c,d)&
                                  -t3(i,j,k,a,b,d)*t2(l,m,c,e)&
                                  +t3(i,j,k,c,b,d)*t2(l,m,a,e)&
                                  +t3(i,j,k,a,c,d)*t2(l,m,b,e)&
                                  -t3(m,j,k,a,b,c)*t2(l,i,d,e)&
                                  +t3(m,i,k,a,b,c)*t2(l,j,d,e)&
                                  +t3(m,j,i,a,b,c)*t2(l,k,d,e)&
                                  )/2

                    isum(3)=&
                         +t1(l,b)*(&
                                     +t2(k,j,d,e)*t2(m,i,a,c)&
                                     -t2(i,j,d,e)*t2(m,k,a,c)&
                                     +t2(i,k,d,e)*t2(m,j,a,c)&
                                     -t3(i,j,k,c,d,e)*t1(m,a)&
                                     +t3(i,j,k,a,d,e)*t1(m,c)&
                                    )/2

                    isum(4)=(&
                                     +t1(j,e)*(&
                                                +t2(i,k,c,d)*t2(l,m,b,a)&
                                                +t2(i,k,b,d)*t2(l,m,a,c)&
                                                -t2(i,k,a,d)*t2(l,m,b,c)&
                                                -t3(l,m,i,a,b,c)*t1(k,d)&
                                                +t3(l,m,k,a,b,c)*t1(i,d)&
                                                )&

                                     +t1(k,e)*(&
                                                +t2(i,j,a,d)*t2(l,m,b,c)&
                                                -t2(i,j,b,d)*t2(l,m,a,c)&
                                                -t2(i,j,c,d)*t2(l,m,b,a)&
                                                -t3(l,m,j,a,b,c)*t1(i,d)&
                                                )&

                                     +t1(i,e)*(&
                                                -t2(k,j,a,d)*t2(l,m,b,c)&
                                                +t2(k,j,b,d)*t2(l,m,a,c)&
                                                +t2(k,j,c,d)*t2(l,m,b,a)&
                                                )&

                                     +t1(l,a)*(&
                                                +t2(i,j,d,e)*t2(m,k,b,c)&
                                                -t2(k,j,d,e)*t2(m,i,b,c)&
                                                -t2(i,k,d,e)*t2(m,j,b,c)&
                                                -t3(i,j,k,b,d,e)*t1(m,c)&
                                                )&

                                     +t1(l,c)*(&
                                                -t2(i,j,d,e)*t2(m,k,b,a)&
                                                +t2(k,j,d,e)*t2(m,i,b,a)&
                                                +t2(i,k,d,e)*t2(m,j,b,a)&
                                                )&
                                 )/2

                    isum(5)=&
                        +t3(i,j,l,a,b,d)*t2(m,k,e,c)&
                        -t3(k,j,l,a,b,d)*t2(m,i,e,c)&
                        -t3(i,k,l,a,b,d)*t2(m,j,e,c)&
                        -t3(i,j,l,c,b,d)*t2(m,k,e,a)&
                        +t3(k,j,l,c,b,d)*t2(m,i,e,a)&
                        +t3(i,k,l,c,b,d)*t2(m,j,e,a)&
                        -t3(i,j,l,a,c,d)*t2(m,k,e,b)&
                        +t3(k,j,l,a,c,d)*t2(m,i,e,b)&
                        +t3(i,k,l,a,c,d)*t2(m,j,e,b)

                    isum(6)=&
                        +t1(l,a)*(&
                                    +t3(i,j,k,c,b,d)*t1(m,e)&
                                    +t2(i,j,b,d)*t2(m,k,e,c)&
                                    -t2(i,j,c,d)*t2(m,k,e,b)&
                                    -t2(k,j,b,d)*t2(m,i,e,c)&
                                    +t2(k,j,c,d)*t2(m,i,e,b)&
                                    -t2(i,k,b,d)*t2(m,j,e,c)&
                                    +t2(i,k,c,d)*t2(m,j,e,b)&
                                    +t1(j,e)*(&
                                                +t2(m,k,b,c)*t1(i,d)&
                                                -t2(m,i,b,c)*t1(k,d)&
                                                +t2(i,k,b,d)*t1(m,c)&
                                                )&

                                    +t1(m,c)*(&
                                                +t2(k,j,b,d)*t1(i,e)&
                                                -t2(i,j,b,d)*t1(k,e)&
                                                )&

                                    -t2(m,j,b,c)*t1(i,d)*t1(k,e)&
                                    )

                    isum(7)=&
                        +t1(l,b)*(&
                                    +t3(i,j,k,a,c,d)*t1(m,e)&
                                    -t2(i,j,a,d)*t2(m,k,e,c)&
                                    +t2(i,j,c,d)*t2(m,k,e,a)&
                                    +t2(k,j,a,d)*t2(m,i,e,c)&
                                    -t2(k,j,c,d)*t2(m,i,e,a)&
                                    +t2(i,k,a,d)*t2(m,j,e,c)&
                                    -t2(i,k,c,d)*t2(m,j,e,a)&
                                    +t1(j,e)*(&
                                                -t2(m,k,a,c)*t1(i,d)&
                                                +t2(m,i,a,c)*t1(k,d)&
                                                +t2(i,k,c,d)*t1(m,a)&
                                                -t2(i,k,a,d)*t1(m,c)&
                                                )&

                                    +t1(k,e)*(&
                                                +t2(m,j,a,c)*t1(i,d)&
                                                +t2(i,j,a,d)*t1(m,c)&
                                                -t2(i,j,c,d)*t1(m,a)&
                                                )&

                                    +t1(i,e)*(&
                                                -t2(k,j,a,d)*t1(m,c)&
                                                +t2(k,j,c,d)*t1(m,a)&
                                                )&

                                    )

                    isum(8)=&
                        +t1(l,c)*(&
                                    -t3(i,j,k,a,b,d)*t1(m,e)&
                                    -t2(i,j,b,d)*t2(m,k,e,a)&
                                    +t2(i,j,a,d)*t2(m,k,e,b)&
                                    +t2(k,j,b,d)*t2(m,i,e,a)&
                                    -t2(k,j,a,d)*t2(m,i,e,b)&
                                    +t2(i,k,b,d)*t2(m,j,e,a)&
                                    -t2(i,k,a,d)*t2(m,j,e,b)&
                                    +t1(j,e)*(&
                                                -t2(m,k,b,a)*t1(i,d)&
                                                +t2(m,i,b,a)*t1(k,d)&
                                                )&

                                    +t2(m,j,b,a)*t1(i,d)*t1(k,e)&
                                    )

                    isum(9)=&
                        +t1(l,d)*(&
                                    -t3(m,j,k,a,b,c)*t1(i,e)&
                                    +t3(m,i,k,a,b,c)*t1(j,e)&
                                    +t3(m,j,i,a,b,c)*t1(k,e)&
                                    +t2(k,j,a,e)*t2(m,i,b,c)&
                                    -t2(i,j,a,e)*t2(m,k,b,c)&
                                    +t2(i,k,a,e)*t2(m,j,b,c)&
                                    +t2(i,j,b,e)*t2(m,k,a,c)&
                                    -t2(k,j,b,e)*t2(m,i,a,c)&
                                    -t2(i,k,b,e)*t2(m,j,a,c)&
                                    +t2(i,j,c,e)*t2(m,k,b,a)&
                                    -t2(k,j,c,e)*t2(m,i,b,a)&
                                    -t2(i,k,c,e)*t2(m,j,b,a)&
                                    )

                    isum(10)=&
                        +t1(k,e)*(&
                                    +t3(i,j,l,c,b,d)*t1(m,a)&
                                    +t3(i,j,l,a,c,d)*t1(m,b)&
                                    -t3(i,j,l,a,b,d)*t1(m,c)&
                                    -t2(j,l,a,d)*t2(m,i,b,c)&
                                    +t2(j,l,c,d)*t2(m,i,b,a)&
                                    +t2(j,l,b,d)*t2(m,i,a,c)&
                                    +t2(i,l,a,d)*t2(m,j,b,c)&
                                    -t2(i,l,c,d)*t2(m,j,b,a)&
                                    -t2(i,l,b,d)*t2(m,j,a,c)&
                                    )

                    isum(11)=&
                        +t1(j,e)*(&
                                    -t3(i,k,l,c,b,d)*t1(m,a)&
                                    -t3(i,k,l,a,c,d)*t1(m,b)&
                                    +t3(i,k,l,a,b,d)*t1(m,c)&
                                    +t2(k,l,a,d)*t2(m,i,b,c)&
                                    -t2(k,l,c,d)*t2(m,i,b,a)&
                                    -t2(k,l,b,d)*t2(m,i,a,c)&
                                    -t2(i,l,a,d)*t2(m,k,b,c)&
                                    +t2(i,l,c,d)*t2(m,k,b,a)&
                                    +t2(i,l,b,d)*t2(m,k,a,c)&
                                    )

                    isum(12)=&
                        +t1(i,e)*(&
                                    -t3(k,j,l,c,b,d)*t1(m,a)&
                                    -t3(k,j,l,a,c,d)*t1(m,b)&
                                    +t3(k,j,l,a,b,d)*t1(m,c)&
                                    -t2(k,l,a,d)*t2(m,j,b,c)&
                                    +t2(k,l,b,d)*t2(m,j,a,c)&
                                    +t2(k,l,c,d)*t2(m,j,b,a)&
                                    +t2(j,l,a,d)*t2(m,k,b,c)&
                                    -t2(j,l,b,d)*t2(m,k,a,c)&
                                    -t2(j,l,c,d)*t2(m,k,b,a)&
                                    )
                    sum=sum+ Ax*(+isum( 1)+isum( 2)+isum( 3)+isum( 4)+isum( 5)+isum( 6)&
                                 +isum( 7)+isum( 8)+isum( 9)+isum(10)+isum(11)+isum(12))
                enddo
                enddo
                enddo
                enddo
                !$omp end parallel
                rez=rez+sum

                d3(i,j,k,a,b,c)=rez
            enddo
            enddo
        enddo
        enddo
    enddo
    enddo

    return
    end subroutine projection_ccsdt_triples_spin_hf