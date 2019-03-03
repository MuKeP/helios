    subroutine putPoint

    use glob, only: assignment(=)
    use glob, only: iglu,uch
    use hdb , only: pointToPut,pointAccordance,GlEt,Et,MEt,MMEt
    use hdb , only: polarizbd,geometrybd

    implicit none

    integer(kind=iglu) :: ll,i,j,k,i1,j1,k1, ipos,jpos,zpos,sta,sto


    sta=-(polarizbd%nPoints-1)/2; sto=-sta

    do ll = 1,pointToPut
        i =pointAccordance(ll,1,1); j =pointAccordance(ll,1,2); k =pointAccordance(ll,1,3)
        i1=pointAccordance(ll,2,1); j1=pointAccordance(ll,2,2); k1=pointAccordance(ll,2,3)

        !write (*,'(3(1X,i2)," ->",3(1X,i2))') i,j,k, i1,j1,k1

        GlEt(i1,j1,k1)=GlEt(i,j,k)
    enddo

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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine putAtom(array)

    use glob, only: iglu,rglu
    use hdb,  only: atomEqu,mol

    implicit none

    real   (kind=rglu) :: array(1:mol%nAtoms,1:mol%nAtoms)
    integer(kind=iglu) :: N,i,j,k,l


    N=mol%nAtoms

    do i = 1,N
        j=abs(atomEqu(i,i))

        array(i,i)=array(j,j)
    enddo

    do i = 1,N-1
    do j = i+1,N

        k=int( abs(atomEqu(i,j))/N )
        l=abs(atomEqu(i,j))-k*N

        if (l.EQ.0) then; k=k-1; l=N; endif

        array(i,j)=array(k,l); array(j,i)=array(k,l)
    enddo
    enddo

    return
    end subroutine putAtom

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine putAtomDiagonal(array, transformation)

    use glob, only: iglu,rglu,lglu,true,false
    use hdb,  only: atomEqu,mol

    implicit none

    real   (kind=rglu) :: array(1:mol%nAtoms)
    character  (len=*) :: transformation
    integer(kind=iglu) :: N,i,j,k,l
    logical(kind=lglu) :: tcase


    N=mol%nAtoms

    do i = 1,N
        j=abs(atomEqu(i,i))

        select case (transformation)

            case('i')
                tcase=(mol%atm(i)%coords(1)*mol%atm(j)%coords(1).LT.0).AND.&
                      (mol%atm(i)%coords(2)*mol%atm(j)%coords(2).LT.0).AND.&
                      (mol%atm(i)%coords(3)*mol%atm(j)%coords(3).LT.0)

            case('x')
                tcase=mol%atm(i)%coords(1)*mol%atm(j)%coords(1).LT.0

            case('y')
                tcase=mol%atm(i)%coords(2)*mol%atm(j)%coords(2).LT.0

            case('z')
                tcase=mol%atm(i)%coords(3)*mol%atm(j)%coords(3).LT.0

            case default
                tcase=false

        endselect

        if (tcase) then
            array(i)=-array(j)
        else
            array(i)=+array(j)
        endif
    enddo

    return
    end subroutine putAtomDiagonal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine putBond(array)

    use glob, only: iglu,rglu
    use hdb,  only: bondEqu,mol

    implicit none

    real   (kind=rglu) :: array(1:mol%nBonds,1:mol%nBonds)
    integer(kind=iglu) :: M,i,j,k,l


    M=UBound(array,1)

    do i = 1,M
        j=abs(bondEqu(i,i))

        array(i,i)=array(j,j)
    enddo

    do i = 1,M-1
    do j = i+1,M

        k=int( abs(bondEqu(i,j))/M )
        l=abs(bondEqu(i,j))-k*M

        if (l.EQ.0) then; k=k-1; l=M; endif

        array(i,j)=array(k,l); array(j,i)=array(k,l)
    enddo
    enddo

    return
    end subroutine putBond