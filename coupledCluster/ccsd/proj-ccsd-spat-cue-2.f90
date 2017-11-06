    subroutine projection_ccsd_doubles_spatial_cue

    use glob          , only: rglu,iglu,lglu
    use coupledCluster, only: N,Nocc,Ne,NFnz
    use coupledCluster, only: F,R,Fnz !,mol
    use coupledCluster, only: t1,t2,d2,excSet,iapairs
    use coupledCluster, only: notFitRadius

    implicit none

    real   (kind=rglu) :: rez,sum !,oset(4)
    integer(kind=iglu) :: u,v,i,a,j,b,k,c,l,d
    integer(kind=iglu) :: t,ab(2)
    logical(kind=lglu) :: ai,aj,bi,bj,ba,ji,ki,kj,ka,kb


    d2=0.d0
    do u = 1,Ne
        i=excSet(u,1)
        a=excSet(u,2)
        do v = u,Ne
            j=excSet(v,1)
            b=excSet(v,2)

            if (notFitRadius(2,i,a,j,b)) cycle

            ai=a.EQ.iapairs(i); aj=a.EQ.iapairs(j)
            bi=b.EQ.iapairs(i); bj=b.EQ.iapairs(j)
            ba=a.EQ.b         ; ji=i.EQ.j

            rez=R(a,i,b,j)

            sum=0
            if ((ba.AND.aj).OR.bj) sum=sum+R(b,j,a,a)*t1(i,a)
            if ((ba.AND.bi).OR.ai) sum=sum+R(a,i,b,b)*t1(j,b)
            if ((ba.AND.bj).OR.aj) sum=sum+R(b,j,a,b)*t1(i,b)
            if ((ba.AND.ai).OR.bi) sum=sum+R(a,i,b,a)*t1(j,a)

            if ((ji.AND.bi).OR.bj) sum=sum-R(i,i,b,j)*t1(i,a)
            if ((ji.AND.aj).OR.ai) sum=sum-R(j,j,a,i)*t1(j,b)
            if ((ji.AND.bj).OR.bi) sum=sum-R(j,i,b,j)*t1(j,a)
            if ((ji.AND.ai).OR.aj) sum=sum-R(i,j,a,i)*t1(i,b)
            rez=rez+sum

            sum=0
            !$omp parallel default(shared) private(k,c,ki,kj,ka,kb) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                c=iapairs(k)
                ki = k.NE.i; kj = k.NE.j
                ka = c.EQ.a; kb = c.EQ.b

                if (.NOT.((.NOT.bj).AND.kj)) then
                    sum=sum+R(k,c,b,j)*(2*t2(i,k,a,c)-t2(i,k,c,a)-t1(i,c)*t1(k,a) )
                    if (.NOT.kb) sum=sum+R(k,b,b,j)*(2*t2(i,k,a,b)-t2(i,k,b,a)-t1(i,b)*t1(k,a) )
                endif

                if (.NOT.((.NOT.ai).AND.ki)) then
                    sum=sum+R(a,i,k,c)*(2*t2(j,k,b,c)-t2(j,k,c,b)-t1(j,c)*t1(k,b) )
                     if (.NOT.ka) sum=sum+R(a,i,k,a)*(2*t2(j,k,b,a)-t2(j,k,a,b)-t1(j,a)*t1(k,b) )
                endif

                if (.NOT.((.NOT.aj).AND.kj)) then
                    sum = sum-R(k,j,a,c)*(t1(i,c)*t1(k,b)+t2(i,k,c,b))-R(k,c,a,j)*t2(i,k,b,c)
                    if (.NOT.ka) sum=sum-R(k,j,a,a)*(t1(i,a)*t1(k,b)+t2(i,k,a,b))-R(k,a,a,j)*t2(i,k,b,a)
                endif

                if (.NOT.((.NOT.bi).AND.ki)) then
                    sum=sum-R(k,i,b,c)*(t1(j,c)*t1(k,a)+t2(j,k,c,a))-R(k,c,b,i)*t2(j,k,a,c)
                    if (.NOT.kb) sum = sum-R(k,i,b,b)*(t1(j,b)*t1(k,a)+t2(j,k,b,a))-R(k,b,b,i)*t2(j,k,a,b)
                endif
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0
            !$omp parallel default(shared) private(c) reduction(+:sum)
            !$omp do
            do c = Nocc+1,N
                sum=sum+F(c,a)*(2*t2(i,j,c,b)-t2(i,j,b,c))+&
                        F(c,b)*(2*t2(i,j,a,c)-t2(i,j,c,a))

            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0
            !$omp parallel default(shared) private(k) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                sum=sum+F(k,i)*(2*t2(k,j,a,b)-t2(k,j,b,a))+&
                        F(k,j)*(2*t2(k,i,b,a)-t2(k,i,a,b))
            enddo
            !$omp end parallel
            rez=rez-sum

            sum=0
            sum=sum+R(i,i,j,j)*(t2(i,j,a,b)+t1(i,a)*t1(j,b))
            if (.NOT.ji) sum=sum+R(j,i,i,j)*(t2(j,i,a,b)+t1(j,a)*t1(i,b))
            rez=rez+sum

            sum=0
            sum=sum+R(b,b,a,a)*( t2(i,j,a,b)+t1(i,a)*t1(j,b))
            if (.NOT.ba) sum=sum+R(b,a,a,b)*( t2(i,j,b,a)+t1(i,b)*t1(j,a))
            rez=rez+sum

            sum=0
            !$omp parallel default(shared) private(k,c,d,t,ab) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                ab(1)=b; ab(2)=iapairs(k)
                do t = 1,2
                    c=ab(t); d=ab(3-t)

                    sum=sum+R(k,c,b,d)*&
                                       ( t1(j,d)*( 2*t2(i,k,a,c)-t2(i,k,c,a) )&
                                                    -t2(j,k,a,c)*t1(i,d)&
                                                    -t2(j,k,d,a)*t1(i,c)&
                                        )&

                    +R(k,d,b,c)*t1(k,d)*(2*t2(i,j,a,c)-t2(i,j,c,a)     )&
                    -R(b,c,k,d)*t1(k,a)*(  t2(i,j,d,c)+t1(j,c)*t1(i,d) )
                enddo

                ab(1)=a; ab(2)=iapairs(k)
                do t = 1,2
                    c=ab(t); d=ab(3-t)

                    sum=sum+R(k,c,a,d)*&
                                       ( t1(i,d)*(2*t2(j,k,b,c)-t2(j,k,c,b) )&
                                                   -t2(i,k,b,c)*t1(j,d)&
                                                   -t2(i,k,d,b)*t1(j,c)&
                                        )&

                    +R(k,d,a,c)*t1(k,d)*(2*t2(j,i,b,c)-t2(j,i,c,b) )&
                    -R(a,c,k,d)*t1(k,b)*(  t2(i,j,c,d)+t1(i,c)*t1(j,d) )
                enddo

            enddo ! k
            !$omp end parallel
            rez=rez+sum

            sum=0
            !$omp parallel default(shared) private(k,l,c,t,ab) reduction(+:sum)
            !$omp do
            do c = Nocc+1,N
                ab(1)=j; ab(2)=iapairs(c)
                do t = 1,2
                    k=ab(t); l=ab(3-t)

                    sum=sum&
                            +R(k,c,l,j)*(&
                                          +t1(l,b)*(-2*t2(i,k,a,c)+t2(i,k,c,a) )&
                                          +t2(i,k,b,c)*t1(l,a)&
                                          +t2(i,l,c,b)*t1(k,a)&
                                          +t1(i,c)*(   t2(k,l,a,b)+t1(k,a)*t1(l,b) )&
                                         )&
                    +R(k,j,l,c)*t1(l,c)*(-2*t2(i,k,a,b)+t2(i,k,b,a) )
                enddo

                ab(1)=i; ab(2)=iapairs(c)
                do t = 1,2
                    k=ab(t); l=ab(3-t)

                    sum=sum&
                            +R(k,c,l,i)*(&
                                          +t1(l,a)*(-2*t2(j,k,b,c)+t2(j,k,c,b) )&
                                          +t2(j,k,a,c)*t1(l,b)&
                                          +t2(j,l,c,a)*t1(k,b)&
                                          +t1(j,c)*(   t2(k,l,b,a)+t1(k,b)*t1(l,a) )&
                                         )&
                     +R(k,i,l,c)* t1(l,c)*(-2*t2(j,k,b,a)+t2(j,k,a,b) )
                enddo

            enddo ! c
            !$omp end parallel
            rez=rez+sum

            sum=0
            !$omp parallel default(shared) private(k,c,l,d,t,ab) reduction(+:sum)
            !$omp do
            do k = 1,Nocc
                do l = 1,Nocc
                    ab(1)=iapairs(k); ab(2)=iapairs(l)
                    do t = 1,2
                        c=ab(t); d=ab(3-t)

                        sum=sum+R(k,c,l,d)*(&

                        -2*(&
                             +t2(j,l,b,d)*(t1(i,c)*t1(k,a)-t2(i,k,a,c)+t2(k,i,a,c) )&
                             +t2(i,l,a,d)*(t1(j,c)*t1(k,b)+t2(k,j,b,c)-t2(j,k,b,c) )&

                             +t1(k,c)*(t2(j,l,b,a)*t1(i,d)+t2(i,l,a,b)*t1(j,d))&
                             +t1(l,d)*(t2(i,j,a,c)*t1(k,b)+t2(i,j,c,b)*t1(k,a))&

                             +t2(i,k,d,c)*t2(j,l,b,a)+t2(j,k,d,c)*t2(i,l,a,b)&
                             +t2(i,j,a,c)*t2(l,k,d,b)+t2(i,j,c,b)*t2(k,l,a,d)&
                            )&

                        +t1(k,a)*(&
                                   +t2(l,j,b,d)*t1(i,c)&
                                   +t2(i,l,b,d)*t1(j,c)&
                                   +t2(l,i,b,c)*t1(j,d)&
                                   +t1(l,b)*( t2(i,j,c,d)+t1(i,c)*t1(j,d) )&
                                  )&

                        +t1(k,b)*(&
                                   +t2(l,j,a,c)*t1(i,d)&
                                   +t2(l,i,a,d)*t1(j,c)&
                                   +t2(j,l,a,d)*t1(i,c)&
                                  )&

                        +t1(k,c)*(t2(i,l,b,a)*t1(j,d)+t2(j,l,a,b)*t1(i,d))&
                        +t1(l,d)*(t2(i,j,c,a)*t1(k,b)+t2(i,j,b,c)*t1(k,a))&

                        +t2(j,k,d,c)*t2(i,l,b,a)+t2(i,k,d,c)*t2(j,l,a,b)&
                        +t2(i,j,b,c)*t2(k,l,a,d)+t2(i,j,c,a)*t2(l,k,d,b)&

                        -t2(i,l,b,d)*(t2(j,k,a,c)-t2(k,j,a,c))&
                        +t2(j,l,a,d)*(t2(k,i,b,c)-t2(i,k,b,c))&

                        +(+t2(k,i,a,c)*t2(l,j,b,d)+t2(k,j,b,c)*t2(l,i,a,d)&
                          +t2(j,k,d,a)*t2(l,i,b,c)+t2(i,k,d,b)*t2(l,j,a,c)&
                          )/2 &

                        +t2(k,l,a,b)*(t2(i,j,c,d)+t1(i,c)*t1(j,d))&

                                           &)
                        if (k.EQ.l) exit
                    enddo
                enddo
            enddo
            !$omp end parallel
            rez=rez+sum

            sum=0
            !$omp parallel default(shared) private(k,c,t) reduction(+:sum)
            !$omp do
            do t = 1,NFnz
                k=Fnz(1,t); c=Fnz(2,t)

                sum=sum+F(k,c)*(+t1(i,c)*(-2*t2(k,j,a,b)+t2(k,j,b,a))&
                                +t1(j,c)*(-2*t2(k,i,b,a)+t2(k,i,a,b))&
                                +t1(k,b)*(-2*t2(i,j,a,c)+t2(i,j,c,a))&
                                +t1(k,a)*(-2*t2(i,j,c,b)+t2(i,j,b,c))&
                                )


            enddo
            !$omp end parallel
            rez=rez+sum

            d2(i,j,a,b)=rez
            d2(j,i,b,a)=rez
        enddo
    enddo

    return
    end subroutine projection_ccsd_doubles_spatial_cue