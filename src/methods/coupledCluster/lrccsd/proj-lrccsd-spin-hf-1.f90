    subroutine projection_lrccsd_singles_spin_hf

    use glob,         only: iglu,rglu
    use lrccsdModule, only: Fab,Fij,ai1,Fia,ai6,ai7
    use lrccsdModule, only: Nel,No
    use lrccsdModule, only: r1,r2,d1

    implicit none

    integer(kind=iglu) :: i,j,k,a,b,c
    real   (kind=rglu) :: rez,sum


    d1=0
    !$omp parallel default(shared) private(i,a,j,b,k,c,rez,sum)
    !$omp do
    do i = 1,Nel
    do a = Nel+1,No
        if (btest(a,0).NE.btest(i,0)) cycle

        rez=0

        sum=0
        do b = Nel+1,No
            sum=sum+Fab(a,b)*r1(i,b)
        enddo
        rez=rez+sum

        sum=0
        do j = 1,Nel
            sum=sum-Fij(i,j)*r1(j,a)
        enddo
        rez=rez+sum

        sum=0
        do j = 1,Nel
        do b = Nel+1,No
            sum=sum-ai1(i,b,j,a)*r1(j,b)
        enddo
        enddo
        rez=rez+sum

        sum=0
        do j = 1,Nel
        do b = Nel+1,No
            sum=sum+Fia(j,b)*r2(i,j,a,b)
        enddo
        enddo
        rez=rez+sum

        sum=0
        do j = 1,Nel
        do k = 1,Nel
        do b = Nel+1,No
            sum=sum+ai6(j,k,i,b)*r2(j,k,a,b)
        enddo
        enddo
        enddo
        rez=rez-sum/2

        sum=0
        do j = 1,Nel
        do b = Nel+1,No
        do c = Nel+1,No
            sum=sum+ai7(j,a,b,c)*r2(i,j,b,c)
        enddo
        enddo
        enddo
        rez=rez-sum/2

        d1(i,a)=rez
    enddo
    enddo
    !$omp end parallel

    return
    end subroutine projection_lrccsd_singles_spin_hf