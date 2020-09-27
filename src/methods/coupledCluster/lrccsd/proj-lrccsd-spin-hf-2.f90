    subroutine projection_lrccsd_doubles_spin_hf

    use glob,         only: iglu,rglu
    use lrccsdModule, only: Fab,Fij,ai1,ai2,ai3,ai4,ai5,at1,at2,at3,at4
    use lrccsdModule, only: Nel,No
    use lrccsdModule, only: r1,r2,d2,t2

    implicit none

    integer(kind=iglu) :: i,j,k,l,a,b,c,d
    real   (kind=rglu) :: rez,sum


    d2=0
    !$omp parallel default(shared) private (a,b,k,c,l,d,rez,sum)
    !$omp do
    do i=1, Nel-1
    do j=i+1, Nel
    do a=Nel+1,No-1
    do b=a+1, No
        if (btest(a+b,0).NE.btest(i+j,0)) cycle

        rez=0

        sum=0
        do c = Nel+1,No
            sum=sum+r1(i,c)*ai3(j,c,a,b)&
                   -r1(j,c)*ai3(i,c,a,b)&
                   +fab(b,c)*r2(i,j,a,c)&
                   -fab(a,c)*r2(i,j,b,c)
        enddo
        rez=rez+sum

        sum=0
        do k = 1,Nel
            sum=sum+r1(k,a)*ai2(i,j,k,b)&
                   -r1(k,b)*ai2(i,j,k,a)&
                   +fij(i,k)*r2(j,k,a,b)&
                   -fij(j,k)*r2(i,k,a,b)
        enddo
        rez=rez+sum

        sum=0
        do l = 1,Nel
            sum=sum+at1(i,l)*t2(j,l,a,b)&
                   -at1(j,l)*t2(i,l,a,b)&
                   +at3(i,l)*t2(j,l,a,b)&
                   -at3(j,l)*t2(i,l,a,b)
        enddo
        rez=rez+sum

        sum=0
        do d = Nel+1,No
            sum=sum+at2(a,d)*t2(i,j,b,d)&
                   -at2(b,d)*t2(i,j,a,d)&
                   +at4(a,d)*t2(i,j,b,d)&
                   -at4(b,d)*t2(i,j,a,d)
        enddo
        rez=rez+sum

        sum=0
        do k = 1,Nel
        do c = Nel+1,No
            sum=sum+r2(j,k,a,c)*ai1(i,c,k,b)&
                   -r2(i,k,a,c)*ai1(j,c,k,b)&
                   +r2(i,k,b,c)*ai1(j,c,k,a)&
                   -r2(j,k,b,c)*ai1(i,c,k,a)
        enddo
        enddo
        rez=rez+sum

        sum=0
        do k = 1,Nel
        do l = 1,Nel
            sum=sum+r2(k,l,a,b)*ai4(i,j,k,l)
        enddo
        enddo
        rez=rez+sum

        sum=0
        do c = Nel+1,No
        do d = Nel+1,No
            sum=sum+r2(i,j,c,d)*ai5(a,b,c,d)
        enddo
        enddo
        rez=rez+sum

        d2(i,j,a,b)=+rez
        d2(i,j,b,a)=-rez
        d2(j,i,a,b)=-rez
        d2(j,i,b,a)=+rez
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel

    return
    end subroutine projection_lrccsd_doubles_spin_hf