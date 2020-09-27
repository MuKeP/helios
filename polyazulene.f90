    subroutine analize_azulenes

    use glob,                  only: rglu,iglu,lglu,void,gluCompare,mid,true,false
    use hdb,                   only: mol,ou,ccbd,ouWidth
    use printmod,              only: prMatrix,prVector
    use txtParser,             only: tpAdjustc
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu)              :: i,j,a,b,k,l,c,d,m,uc(4),holder(4),count,mm,nn,pp,vv
    real   (kind=rglu), allocatable :: holdv(:),holdr(:,:)
    real   (kind=rglu)              :: sum3,rest,sumw,sss1,sss2
    real   (kind=rglu)              :: Ax,Bx,Cx,Dx,Ex,Fx,Gx,Hx,Ix,Jx


    umethod='cue-ccsd'
    void=prepareForAnalize()

    m=N/10
    allocate(holdv(m),holdr(m,m))

    ! do i = 1,2*N
    !     write (ou,*) i,ucell(i)
    ! enddo

    ! return

    holdv=0; holdr=0; sumw=0
    write (ou,*) Nse
    do mm = 1,Nse
        write (ou, '(1X,i4\)') mm
        i=spinInds(mm,1)
        a=spinInds(mm,2)

        Ax=cc_c1_t1([i,a])**2
        sumw=sumw+Ax
        call put_matrix(Ax, [i,a,0,0,0,0,0,0])


        !!$omp parallel default(shared) private(nn,pp,vv,j,b,k,c,l,d,Bx,Cx,Dx,Ex,Fx,Gx,Hx,Ix,Jx) reduction(+:hold,sumw)
        !!$omp do
        do nn = 1,Nse
            j=spinInds(nn,1)
            b=spinInds(nn,2)

            if ((j.LE.i).OR.(b.LE.a)) cycle

            Bx=cc_c2_t2([i,j,a,b])**2
            Cx=cc_c2_t1t1([i,j,a,b])**2
            sumw=sumw+Bx+Cx
            call put_matrix(Bx+Cx, [i,a,j,b,0,0,0,0])

            do pp = 1,Nse
                k=spinInds(pp,1)
                c=spinInds(pp,2)

                if ((k.LE.j).OR.(c.LE.b)) cycle

                Dx = cc_c3_t2t1([i,j,k,a,b,c])**2
                Ex = cc_c3_t1t1t1([i,j,k,a,b,c])**2
                Fx = cc_c3_t3([i,j,k,a,b,c])**2
                sumw=sumw+Dx+Ex+Fx
                call put_matrix(Dx+Ex, [i,a,j,b,k,c,0,0])

                do vv = 1,Nse
                    l=spinInds(vv,1)
                    d=spinInds(vv,2)

                    if ((l.LE.k).OR.(d.LE.c)) cycle

                    Gx=cc_c4_t2t2([i,j,k,l,a,b,c,d])**2
                    Hx=cc_c4_t2t1t1([i,j,k,l,a,b,c,d])**2
                    Ix=cc_c4_t1t1t1t1([i,j,k,l,a,b,c,d])**2
                    Jx=cc_c4_t3t1([i,j,k,l,a,b,c,d])**2
                    sumw=sumw+Gx+Hx+Ix+Jx
                    call put_matrix(Gx+Hx+Ix+Jx, [i,a,j,b,k,c,l,d])

                enddo !vv
            enddo !pp
        enddo !nn
        !!$omp end parallel

    enddo !mm

    sss1=0; sss2=0
    do i = 1,m
        sss2=sss2+holdv(i)
        do j = 1,m
            sss1=sss1+holdr(i,j)
        enddo
    enddo

    holdr=100*holdr/sss1
    holdv=100*holdv/sss2
    rest=100*rest/sumw
    call prVector(holdv,ou,'Unitcell employment','^.00','horizontal',ouWidth)
    call prMatrix(holdr,ou,'Inter-uc transition','^.00',ouWidth)
    write (ou,'(2X,A,2X,F6.2)') 'The rest.',rest
    write (ou,'(2X,A,2X,F6.2)') 'W-F norm',sumw

!    stop
    return

    contains

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    integer(kind=iglu) function ucell(orbital) result(ret)
    implicit none

    integer(kind=iglu) :: orbital,uorbital

    uorbital=orbital

    if (uorbital.EQ.0) then
        ret=0
        return
    endif
    if (uorbital.GT.N) then
        uorbital=uorbital-N
    endif

    ret=int((uorbital-1)/10)+1

    return
    end function ucell

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    integer(kind=iglu) function set_len(mons,set,vals) result(ret)
    implicit none

    integer(kind=iglu) :: i,k,mons(8),set(8),vals
    logical(kind=lglu) :: pres

    set=0; ret=0; vals=0
    do i = 1,8
        if (mons(i).EQ.0) exit
        vals=vals+1
        pres=false
        do k = 1,ret
            if (set(k).EQ.mons(i)) then
                pres=true
                exit
            endif
        enddo
        if (.NOT.pres) then
            ret=ret+1
            set(ret)=mons(i)
        endif
    enddo

    return
    end function set_len

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    subroutine put_vector(value,arr)
    implicit none

    integer(kind=iglu), intent(in) :: arr(8)
    real   (kind=rglu), intent(in) :: value
    integer(kind=iglu)             :: i,k


    do i = 1,8
        k=ucell(arr(i))
        holdv(k)=holdv(k)+value
    enddo


    return
    end subroutine put_vector

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    subroutine put_matrix(value,arr)
    implicit none

    integer(kind=iglu), intent(in) :: arr(8)
    real   (kind=rglu), intent(in) :: value
    integer(kind=iglu)             :: i,k,uc(8),holder(8),vals


    call put_vector(value,arr)

    do i = 1,8
        uc(i)=ucell(arr(i))
    enddo

    count=set_len(uc,holder,vals)

    select case(count)
        case(1)
            k=holder(1); l=k
            holdr(k,l)=holdr(k,l)+value

        case(2)
            k=holder(1); l=holder(2)
            holdr(k,l)=holdr(k,l)+value
            holdr(l,k)=holdr(l,k)+value

        case default
            rest=rest+value

    end select

    write (100+vals/2,100) uc(1:vals),value

100 format(<vals>(1X,i2),1X,F20.12)

    return
    end subroutine put_matrix

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    end subroutine analize_azulenes