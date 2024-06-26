    subroutine projection_ccsdt_singles_spin_cue

    use glob,           only: rglu,iglu,lglu
    use coupledCluster, only: No,Nel,NFnz
    use coupledCluster, only: R,F,Fnz
    use coupledCluster, only: t1,t2,t3,d1,iapairs
    use coupledCluster, only: notFitRadius

    implicit none

    real   (kind=rglu) :: sum,rez
    integer(kind=iglu) :: ab(2),t,i,a,j,b,k,c


    do i = 1,Nel
    do a = Nel+1,No
        if (btest(i,0).NE.btest(a,0)) cycle

        if (notFitRadius(1,i,a)) cycle

        rez=0

        rez=rez+F(i,a)

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
        !$omp parallel default(shared) private(t,j,b) reduction(+:sum)
        !$omp do
        do t = 1,NFnz
            j=Fnz(1,t)
            b=Fnz(2,t)

            sum=sum+F(j,b)*(t2(i,j,a,b)-t1(i,b)*t1(j,a))
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0 !opt
        !$omp parallel default(shared) private(j,b) reduction(+:sum)
        !$omp do
        do j = 1,Nel
            b=iapairs(j)
            sum=sum+R(i,a,b,j)*t1(j,b)
        enddo
        !$omp end parallel

        if (iapairs(i).NE.a) then
            j=i; b=a; sum=sum+R(i,a,b,j)*t1(j,b)
        endif
        rez=rez+sum

        sum=0 !opt
        !$omp parallel default(shared) private(ab,t,j,b,c) reduction(+:sum)
        !$omp do
        do j = 1,Nel
            ab(1)=a; ab(2)=iapairs(j)
            do t = 1,2
                b=ab(t); c=ab(3-t)

                sum=sum+R(a,b,j,c)*(t1(i,b)*t1(j,c)+t2(i,j,b,c)/2)
            enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        sum=0 !opt
        !$omp parallel default(shared) private(ab,t,b,j,k) reduction(+:sum)
        !$omp do
        do b = Nel+1,No
            ab(1)=i; ab(2)=iapairs(b)
            do t = 1,2
                j=ab(t); k=ab(3-t)
                sum=sum+R(i,j,b,k)*(t1(j,a)*t1(k,b)+t2(j,k,a,b)/2)
            enddo
        enddo
        !$omp end parallel
        rez=rez-sum

        sum=0 !opt
        !$omp parallel default(shared) private(ab,t,j,k,b,c) reduction(+:sum)
        !$omp do
        do j = 1,Nel
        do k = 1,Nel
            if (k.EQ.j) cycle
            ab(1)=iapairs(j); ab(2)=iapairs(k)
            do t = 1,2
                b=ab(t); c=ab(3-t)
                sum=sum+R(j,b,k,c)*(&
                                        +t3(i,j,k,a,b,c)/4 &

                                        -(&
                                           +t1(i,b)*t2(j,k,a,c)&
                                           +t1(k,a)*t2(j,i,b,c)&
                                          )/2 &

                                        +t1(k,c)*t2(i,j,a,b)&
                                        -t1(j,b)*t1(i,c)*t1(k,a)&
                                        )
            enddo
        enddo
        enddo
        !$omp end parallel
        rez=rez+sum

        d1(i,a)=rez
    enddo
    enddo

    return
    end subroutine projection_ccsdt_singles_spin_cue