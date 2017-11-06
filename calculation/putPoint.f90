    subroutine putPoint

    use glob, only: assignment(=)
    use glob, only: iglu,uch
    use hdb , only: pointToPut,pointAccordance,GlEt,Et,MEt,MMEt
    use hdb , only: polarizbd,geometrybd

    implicit none

    integer(kind=iglu) :: ll,i,j,k,i1,j1,k1, ipos,jpos,zpos,sta,sto


    sta=-(polarizbd%nPoints-1)/2; sto=-sta
    !write (70,*) GlEt
    !stop

    do ll = 1,pointToPut
        i =pointAccordance(ll,1,1); j =pointAccordance(ll,1,2); k =pointAccordance(ll,1,3)
        i1=pointAccordance(ll,2,1); j1=pointAccordance(ll,2,2); k1=pointAccordance(ll,2,3)

        !write (*,'(3(1X,i2)," ->",3(1X,i2))') i,j,k, i1,j1,k1

        GlEt(i1,j1,k1)=GlEt(i,j,k)
    enddo

    !write (*,*) GlEt
    !stop

    if (geometrybd%searchLinear(2)) then
        do i = sta,sto
            GlEt(i,:,:)=GlEt(i,0,0)
        enddo
        polarizbd%scales='xyz'
    endif

    if (geometrybd%searchPlanar(2)) then
        do i = sta,sto
        do j = sta,sto
            GlEt(i,j,:)=GlEt(i,j,0)
        enddo
        enddo
        polarizbd%scales='xyz'
    endif

    MMEt=GlEt
    select case (polarizbd%scales%ln)

        case (1)
            select case (polarizbd%scales%get())
                case ('x'); Et=GlEt(:,0,0)
                case ('y'); Et=GlEt(0,:,0)
                case ('z'); Et=GlEt(0,0,:)
            end select

        case (2)
            select case (polarizbd%scales%get(1,1))
                case ('x'); ipos=1
                case ('y'); ipos=2
                case ('z'); ipos=3
            end select

            select case (polarizbd%scales%get(2,2))
                case ('x'); jpos=1
                case ('y'); jpos=2
                case ('z'); jpos=3
            end select

            if (ipos.GT.jpos) then
                i=jpos; jpos=ipos; ipos=i
            endif

            if (ipos.EQ.1) then
                if (jpos.EQ.2) zpos=3
                if (jpos.EQ.3) zpos=2
            else;              zpos=1
            endif

            select case (zpos)
                case (1); MEt=GlEt(0,:,:)
                case (2); MEt=GlEt(:,0,:)
                case (3); MEt=GlEt(:,:,0)
            end select

        case (3)
            MMEt=GlEt

    end select

    return
    end subroutine putPoint