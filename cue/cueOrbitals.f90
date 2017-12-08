    subroutine cueOrbitals

    use glob    , only: iglu,rglu,void,i8kind,glControlMemory
    use hdb     , only: eu,ou,cueConstant1
    use hdb     , only: mol,cuebd
    use printmod, only: prMatrix

    implicit none

    integer(kind=iglu)  :: N,M,Nel,Nocc,i,j,k,l,v,p
    real(kind=rglu)     :: coords(3),val
    integer(kind=iglu), allocatable :: connectivity(:)



    N=mol%nAtoms; M=mol%nBonds
    Nel=mol%nEls; Nocc=mol%nEls/2

    void=glControlMemory(int(iglu*N,kind=i8kind),'tmp. Connectivity check')
    allocate (connectivity(N)); connectivity=0

    do k = 1,N
        if (mol%atm(k)%nels.EQ.2) then
            connectivity(k)=connectivity(k)+1
        endif
    enddo

    do k = 1,M
        if (mol%bnd(k)%kind.EQ.'double') then
            i=mol%bnd(k)%atoms(1); j=mol%bnd(k)%atoms(2)
            connectivity(i)=connectivity(i)+1
            connectivity(j)=connectivity(j)+1
        endif
    enddo

    if ( .NOT.( (maxval(connectivity).EQ.1).AND.(minval(connectivity).EQ.1) ) ) then
        write (eu,*) 'Error. Proposed CUE basis is not correct.'; stop
    endif

    deallocate (connectivity)
    void=glControlMemory(int(sizeof(connectivity),kind=i8kind),'tmp. Connectivity check','free')

    l=0
    do k = 1,M
        if (mol%bnd(k)%kind.EQ.'double') then
            l=l+1; v=l+Nocc
            i=mol%bnd(k)%atoms(1); j=mol%bnd(k)%atoms(2)

            mol%orb(l)%ova=v;      mol%orb(v)%ova=l

            mol%orb(l)%atoms(1)=i; mol%orb(v)%atoms(1)=i
            mol%orb(l)%atoms(2)=j; mol%orb(v)%atoms(2)=j
            mol%orb(l)%nels=1;     mol%orb(v)%nels=1

            mol%orb(l)%coef(1)=cueConstant1;  mol%orb(v)%coef(1)=-cueConstant1
            mol%orb(l)%coef(2)=cueConstant1;  mol%orb(v)%coef(2)= cueConstant1

            do p = 1,3
                coords(p)=( mol%atm(i)%coords(p)+mol%atm(j)%coords(p) )/2
            enddo
            mol%orb(l)%coords=coords; mol%orb(v)%coords=coords
        endif
    enddo

    do j = 1,N
        if (mol%atm(j)%nels.EQ.2) then
            l=l+1

            mol%orb(l)%ova=Nocc+1 !any orthogonal vacant orbital

            mol%orb(l)%atoms(1)=j
            if (j.EQ.N) then !any non-equal to j
                mol%orb(l)%atoms(2)=1
            else
                mol%orb(l)%atoms(2)=N
            endif

            k=mol%orb(l)%ova

            !write (*,*) mol%orb(k)%atoms(1),mol%orb(k)%atoms(2),mol%orb(l)%atoms(2)

            mol%orb(l)%nels=2

            mol%orb(l)%coef(1)=1
            mol%orb(l)%coef(2)=0

            mol%orb(l)%coords=mol%atm(j)%coords
        endif
    enddo

    !$omp parallel default(shared) private(i,j,k,val)
    !$omp do
    do i = 1,N-1
        do j = i+1,N
            val=0
            do k = 1,3
                val=val+( mol%orb(i)%coords(k)-mol%orb(j)%coords(k) )**2
            enddo
            mol%cueDist(i,j)=sqrt(val)
            mol%cueDist(j,i)=sqrt(val)
        enddo
    enddo
    !$omp end parallel

    call cueStructureAnalysis

    return
    end subroutine cueOrbitals