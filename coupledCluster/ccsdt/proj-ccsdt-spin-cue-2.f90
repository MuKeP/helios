    subroutine projection_ccsdt_doubles_spin_cue

    use glob          , only: rglu,iglu,lglu
    use coupledCluster, only: No,Nel,NFnz
    use coupledCluster, only: R,F,Fnz
    use coupledCluster, only: t1,t2,t3,d2,iapairs
    use coupledCluster, only: notFitRadius

    implicit none

    real   (kind=rglu) :: rez,sum
    integer(kind=iglu) :: ab(2),t,i,a,j,b,k,c,l,d


    do i = 1,Nel-1
    do a = Nel+1,No-1
    do j = i+1,Nel
    do b = a+1,No
        if (btest(i+j,0).NE.btest(a+b,0)) cycle

        if (notFitRadius(2,i,a,j,b)) cycle

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
        !$omp parallel default(shared) private(t,k,c) reduction(+:sum)
        !$omp do
        do t = 1,NFnz
            k=Fnz(1,t)
            c=Fnz(2,t)

            sum=sum+F(k,c)*(&
                                -t2(k,j,a,b)*t1(i,c)&
                                +t2(k,i,a,b)*t1(j,c)&
                                -t2(i,j,a,c)*t1(k,b)&
                                +t2(i,j,b,c)*t1(k,a)&
                                +t3(i,j,k,a,b,c)&
                                )
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0 ! opt
        sum=sum+R(b,j,a,a)*t1(i,a)
        sum=sum+R(b,j,a,b)*t1(i,b)
        sum=sum-R(b,i,a,a)*t1(j,a)
        sum=sum-R(b,i,a,b)*t1(j,b)
        sum=sum+R(j,a,i,i)*t1(i,b)
        sum=sum+R(j,a,i,j)*t1(j,b)
        sum=sum-R(j,b,i,i)*t1(i,a)
        sum=sum-R(j,b,i,j)*t1(j,a)
        rez=rez+sum

        sum=0 !opt
        !$omp parallel default(shared) private(k,c) reduction(+:sum)
        !$omp do
        do k = 1,Nel
            c=iapairs(k)
            sum=sum+R(i,a,c,k)*(-t1(k,b)*t1(j,c)+t2(j,k,b,c))
            sum=sum+R(j,a,c,k)*(+t1(k,b)*t1(i,c)-t2(i,k,b,c))
            sum=sum+R(j,b,c,k)*(-t1(k,a)*t1(i,c)+t2(i,k,a,c))
            sum=sum+R(i,b,c,k)*(+t1(k,a)*t1(j,c)-t2(j,k,a,c))
        enddo
        !$omp end parallel

        if (iapairs(i).NE.a) then
            c=a; k=i; sum=sum+R(i,a,c,k)*(-t1(k,b)*t1(j,c)+t2(j,k,b,c))
        endif

        if (iapairs(j).NE.a) then
            c=a; k=j; sum=sum+R(j,a,c,k)*(+t1(k,b)*t1(i,c)-t2(i,k,b,c))
        endif

        if (iapairs(j).NE.b) then
            c=b; k=j; sum=sum+R(j,b,c,k)*(-t1(k,a)*t1(i,c)+t2(i,k,a,c))
        endif

        if (iapairs(i).NE.b) then
            c=b; k=i; sum=sum+R(i,b,c,k)*(+t1(k,a)*t1(j,c)-t2(j,k,a,c))
        endif
        rez=rez+sum

        sum=0 !opt
        ab(1)=i; ab(2)=j
        do t = 1,2
            k=ab(t); l=ab(3-t)
            sum=sum+R(i,k,j,l)*(t2(k,l,a,b)/2+t1(k,a)*t1(l,b))
        enddo
        rez=rez+sum

        sum=0 !opt
        ab(1)=a; ab(2)=b
        do t = 1,2
            c=ab(t); d=ab(3-t)
            sum=sum+R(a,c,b,d)*(t2(i,j,c,d)/2+t1(i,c)*t1(j,d))
        enddo
        rez=rez+sum

        sum=0 !opt
        !$omp parallel default(shared) private(ab,t,k,c,d) reduction(+:sum)
        !$omp do
        do k = 1,Nel
            ab(1)=b; ab(2)=iapairs(k)
            do t = 1,2
                c=ab(t); d=ab(3-t)
                sum=sum+R(b,c,k,d)*(&
                                        -(t1(j,d)*t2(i,k,a,c)-t1(i,d)*t2(j,k,a,c))&
                                        +(t2(i,j,a,c)*t1(k,d)+t1(k,a)*t1(i,c)*t1(j,d)&
                                        +(t1(k,a)*t2(i,j,c,d)+t3(i,j,k,a,c,d))/2 )&
                                        )
            enddo

            ab(1)=a; ab(2)=iapairs(k)
            do t = 1,2
                c=ab(t); d=ab(3-t)
                sum=sum-R(a,c,k,d)*(&
                                        +(t1(i,d)*t2(j,k,b,c)-t1(j,d)*t2(i,k,b,c))&
                                        +(t2(i,j,b,c)*t1(k,d)+t1(k,b)*t1(i,c)*t1(j,d)&
                                        +(t1(k,b)*t2(i,j,c,d)+t3(i,j,k,b,c,d))/2 )&
                                        )
            enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0 !opt
        !$omp parallel default(shared) private(ab,t,c,k,l) reduction(+:sum)
        !$omp do
        do c = Nel+1,No
            ab(1)=i; ab(2)=iapairs(c)
            do t = 1,2
                k=ab(t); l=ab(3-t)
                sum=sum&
                       +R(i,k,c,l)*(&
                                         +      (&
                                                  +t3(j,k,l,a,b,c)&
                                                  +t2(k,l,a,b)*t1(j,c)&
                                                 )/2 &

                                         +t1(l,c)*t2(j,k,a,b)&
                                         +t1(l,a)*t2(j,k,b,c)&
                                        -t1(l,b)*(&
                                                     +t2(j,k,a,c)&
                                                     -t1(j,c)*t1(k,a)&
                                                    )&
                                        )
            enddo

            ab(1)=j; ab(2)=iapairs(c)
            do t = 1,2
                k=ab(t); l=ab(3-t)
                sum=sum&
                       -R(j,k,c,l)*(&
                                         +      (&
                                                  +t3(i,k,l,a,b,c)&
                                                  +t2(k,l,a,b)*t1(i,c)&
                                                 )/2 &

                                        +t1(l,c)*t2(i,k,a,b)&
                                        +t1(l,a)*t2(i,k,b,c)&
                                        -t1(l,b)*(&
                                                    +t2(i,k,a,c)&
                                                    -t1(i,c)*t1(k,a)&
                                                    )&
                                        )
            enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0 !opt
        !$omp parallel default(shared) private(ab,t,k,l,c,d) reduction(+:sum)
        !$omp do
        do k = 1,Nel
        do l = 1,Nel
            if (k.EQ.l) cycle

            ab(1)=iapairs(k); ab(2)=iapairs(l)
            do t = 1,2
                c=ab(t); d=ab(3-t)
                sum=sum+R(k,c,l,d)*(&

                +t2(i,j,c,d)*t2(k,l,a,b)/4&

                +       (&

                         -t2(k,i,c,d)*t2(l,j,a,b)&
                         +t2(k,j,c,d)*t2(l,i,a,b)&
                         -t2(i,j,a,c)*t2(k,l,b,d)&
                         +t2(i,j,b,c)*t2(k,l,a,d)&
                         +t2(k,l,a,b)*t1(i,c)*t1(j,d)&
                         +t2(i,j,c,d)*t1(k,a)*t1(l,b)&
                         -t1(l,b)*t3(i,k,j,a,c,d)&
                         +t1(l,a)*t3(i,k,j,b,c,d)&
                         -t1(i,c)*t3(k,l,j,a,d,b)&
                         +t1(j,c)*t3(k,l,i,a,d,b)&

                        )/2 &

                +t1(l,d)*t3(i,j,k,a,b,c)&

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
        !$omp end parallel
        rez=rez+sum

        d2(i,j,a,b)=rez
    enddo
    enddo
    enddo
    enddo

    return
    end subroutine projection_ccsdt_doubles_spin_cue