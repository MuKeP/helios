    module coupledClusterAnalize

    use glob,           only: rglu,iglu,lglu,i8kind,uch,true,false,void,glControlMemory
    use hdb,            only: mol
    use coupledCluster, only: hV,t3
    use txtParser,      only: operator(.in.)

    integer(kind=iglu), parameter :: chapsLen=32
    real   (kind=rglu), parameter :: prBar=0.05_rglu,ouThreshold=0.001_rglu

    real   (kind=rglu), allocatable                :: V(:,:),G(:,:),distMO(:,:),F(:,:),cMO(:,:)
    real   (kind=rglu), allocatable, dimension (:) :: sumt1,sumt2,sumt3,maxt1,maxt2,maxt3,mean,minv,&
                                                      maxv,inteContr,inteAmp,efft1,efft2,efft3,sd
    real   (kind=rglu), allocatable                :: t1(:,:),t2(:,:,:,:),pvec(:)
    integer(kind=iglu), allocatable                :: spinInds(:,:),spatInds(:,:),countt(:,:)
    real   (kind=rglu), allocatable                :: holdContributions(:,:,:),holdValues(:,:)

    logical(kind=lglu) :: pattern(3),dCUE
    integer(kind=iglu) :: N,Nel,No,Nocc,maxNei,Ne,Nse
    real(kind=rglu)    :: multiplier

    type(uch)          :: umethod

    private :: rglu,iglu,lglu,true,false,void,mol,hV,operator(.in.)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function prepareForAnalize() result(ret)
    use coupledCluster, only: cct1=>t1,cct2=>t2,Fcc=>F
    use lrccsdModule,   only: r1,r2
    use coupledCluster, only: convertSpatialToSpinCC
    implicit none

    integer(kind=iglu) :: k,l,i,j,a,mm
    real(kind=rglu)    :: Ax


    N=mol%nAtoms; Nel=mol%nEls; No=2*Nel; Nocc=Nel/2
    Nse=Nel*(No-Nel); Ne=Nocc*(N-Nocc)

    void=glControlMemory(int( rglu*(Nel*Nel*Nel*Nel+Nel*Nel), kind=i8kind),'Coupled cluster wave-function analize')

    allocate (t1(Nel,Nel+1:No),t2(Nel,Nel,Nel+1:No,Nel+1:No)); t1=0; t2=0
    if ('lr' .in. umethod%get()) then
        t1=r1
        t2=r2
    elseif(umethod%get().EQ.'cue-ccs') then
        void=convertSpatialToSpinCC(Nocc, Nel, cct1,cct2,t1,t2,'s')
    elseif((umethod%get().EQ.'r-ccd').OR.(umethod%get().EQ.'u-ccd')) then
        void=convertSpatialToSpinCC(Nocc, Nel, cct1,cct2,t1,t2,'d')
    else
        void=convertSpatialToSpinCC(Nocc, Nel, cct1,cct2,t1,t2,'sd')
    endif

    pattern=false
    if ('ccs' .in. umethod%get()) then
        pattern(1)=true
    endif

    if (('ccd' .in. umethod%get()).OR.('ccsd' .in. umethod%get())) then
        pattern(2)=true
    endif

    if (('ccsdt' .in. umethod%get())) then
        pattern(3)=true
    endif

    maxNei=mol%nCUELayers

    void=glControlMemory(int( rglu*(N*No+Nel/2) ,kind=i8kind),'Coupled cluster wave-function analize')
    allocate(V(N,No),pvec(Nel/2)); V=0; pvec=0

    if ('cue' .in. umethod%get()) then
        dCUE=true
        multiplier=1._rglu/4._rglu
        do k = 1,Nocc
            i=mol%orb(k)%atoms(1)
            j=mol%orb(k)%atoms(2)
            l=mol%orb(k)%ova

            if (mol%orb(k)%nels.EQ.2) then
                V(i,2*k-1)=1._rglu
                V(i,2*k  )=1._rglu
                cycle
            endif

            Ax=1._rglu/sqrt(2._rglu)

            V(i,2*k-1)=Ax; V(i,2*l-1)=-Ax
            V(i,2*k  )=Ax; V(i,2*l  )=-Ax
            V(j,2*k-1)=Ax; V(j,2*l-1)= Ax
            V(j,2*k  )=Ax; V(j,2*l  )= Ax
        enddo
    else
        dCUE=false
        multiplier=1
        if ('spin' .in. umethod%get()) then
            V=hV
        else
            do i = 1,N
            do j = 1,N
                V(i,2*j-1)=hV(i,j)
                V(i,2*j  )=hV(i,j)
            enddo
            enddo
        endif
    endif

    void=glControlMemory(int( rglu*(N*N+No*No+3*No+No*No) ,kind=i8kind),'Coupled cluster wave-function analize')

    allocate(G(N,N)); G=0
    G=mol%G

    allocate(distMO(No,No)); distMO=0
    do i = 1,N
        do j = i,N
            distMO(2*i-1,2*j-1)=mol%cuedist(i,j)
            distMO(2*i  ,2*j-1)=mol%cuedist(i,j)
            distMO(2*i-1,2*j  )=mol%cuedist(i,j)
            distMO(2*i  ,2*j  )=mol%cuedist(i,j)

            distMO(2*j-1,2*i-1)=mol%cuedist(i,j)
            distMO(2*j  ,2*i-1)=mol%cuedist(i,j)
            distMO(2*j-1,2*i  )=mol%cuedist(i,j)
            distMO(2*j  ,2*i  )=mol%cuedist(i,j)
        enddo
    enddo

    allocate(cMO(No,3),F(No,No)); cMO=0; F=0
    do i = 1,N
        cMO(2*i-1,:) = mol%orb(i)%coords
        cMO(2*i  ,:) = mol%orb(i)%coords
    enddo
    if ('spin' .in. umethod%get()) then
        F=Fcc
    else
        do i = 1,N
        do j = 1,N
            F(2*i-1,2*j-1)=Fcc(i,j); F(2*i,2*j)=Fcc(i,j)
            F(2*j-1,2*i-1)=Fcc(i,j); F(2*j,2*i)=Fcc(i,j)
        enddo
        enddo
    endif

    void=glControlMemory(int( iglu*(2*Nse+2*Ne) ,kind=i8kind),'Coupled cluster wave-function analize')

    allocate ( spinInds(Nse,2),spatInds(Ne,2) )
    allocate ( holdValues(maxNei+1,5) ); holdValues=0

    allocate ( sumt1(maxNei),sumt2(maxNei),sumt3(maxNei),maxt1(maxNei),maxt2(maxNei),maxt3(maxNei),&
    &          inteContr(maxNei),inteAmp(maxNei),efft1(maxNei),efft2(maxNei),efft3(maxNei) )
    allocate ( countt(maxNei,3),mean(maxNei),sd(maxNei),minv(maxNei),maxv(maxNei) )

    void=glControlMemory(int( rglu*(Nel*maxNei*2) ,kind=i8kind),'Coupled cluster wave-function analize')

    allocate ( holdContributions(Nel,maxNei,2) ); holdContributions=0

    mm=0
    do i = 1,Nel
    do a = Nel+1,No
        mm=mm+1; spinInds(mm,1)=i; spinInds(mm,2)=a
    enddo
    enddo

    mm=0
    do i = 1,Nocc
    do a = Nocc+1,N
        mm=mm+1; spatInds(mm,1)=i; spatInds(mm,2)=a
    enddo
    enddo

    ret=0; return
    end function prepareForAnalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function finishAnalize() result(ret)
    implicit none

    deallocate (sumt1,sumt2,sumt3,maxt1,maxt2,maxt3,inteContr,inteAmp,efft1,efft2,efft3,countt)
    deallocate (spinInds,spatInds,holdValues,holdContributions,mean,sd,minv,maxv)
    deallocate (t1,t2,F,cMO,distMO,G,V,pvec)

    void=glControlMemory(int( sizeof(t1)+sizeof(t2)+&
                              sizeof(F)+sizeof(cMO)+sizeof(pvec)+&
                              sizeof(distMO)+sizeof(G)+sizeof(V),&
                              kind=i8kind),'Coupled cluster wave-function analize', 'free')

    ret=0
    return
    end function finishAnalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical(kind=lglu) function inLayer(dist,layer) result(ret)
    implicit none

    real(kind=rglu)    :: dist
    integer(kind=iglu) :: layer


    ret=(dist.GT.mol%cueLayers(layer-1)).AND.(dist.LE.mol%cueLayers(layer))
    return
    end function inLayer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function twoeintegr(a,b,c,d) result(ret)
    implicit none
    integer(kind=iglu), intent(in) :: a,b,c,d
    real   (kind=rglu)             :: sum1,sum2
    integer(kind=iglu)             :: mu,nu


    sum1=0
    if ((.NOT.btest(a+b,0)).AND.(.NOT.btest(c+d,0))) then ! [ab|cd]
        do mu = 1,N
        do nu = 1,N
            sum1=sum1+V(mu,a)*V(mu,b)*V(nu,c)*V(nu,d)*G(mu,nu)
        enddo
        enddo
    endif

    sum2=0
    if ((.NOT.btest(a+d,0)).AND.(.NOT.btest(b+c,0))) then ! [ad|bc]
        do mu = 1,N
        do nu = 1,N
            sum2=sum2+V(mu,a)*V(mu,d)*V(nu,b)*V(nu,c)*G(mu,nu)
        enddo
        enddo
    endif

    ret=sum1-sum2; return
    end function twoeintegr

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c1_t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1


    ret=0
    if (.NOT.pattern(1)) return

    i1=arr(1)
    a1=arr(2)


    ret=ret+t1(i1,a1)

    return
    end function c1_t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c2_t1t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2


    ret=0
    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2)
    a1=arr(3); a2=arr(4)

    ret=0

    ret=ret+t1(i1,a1)*t1(i2,a2)-t1(i2,a1)*t1(i1,a2)

    return
    end function c2_t1t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c2_t2(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2


    ret=0
    if (.NOT.pattern(2)) return

    i1=arr(1); i2=arr(2)
    a1=arr(3); a2=arr(4)

    ret=ret+t2(i1,i2,a1,a2)

    return
    end function c2_t2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c3_t1t1t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2,i3,a3


    ret=0
    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2); i3=arr(3)
    a1=arr(4); a2=arr(5); a3=arr(6)

    ret=ret+t1(i1,a1)*t1(i2,a2)*t1(i3,a3)&
           -t1(i2,a1)*t1(i1,a2)*t1(i3,a3)&
           -t1(i1,a1)*t1(i3,a2)*t1(i2,a3)&
           +t1(i3,a1)*t1(i1,a2)*t1(i2,a3)&
           -t1(i3,a1)*t1(i2,a2)*t1(i1,a3)&
           +t1(i2,a1)*t1(i3,a2)*t1(i1,a3)

    return
    end function c3_t1t1t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c3_t1t2(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2,i3,a3


    ret=0
    if (.NOT.(pattern(1).AND.pattern(2))) return

    i1=arr(1); i2=arr(2); i3=arr(3)
    a1=arr(4); a2=arr(5); a3=arr(6)

    ret=ret+t1(i1,a1)*t2(i2,i3,a2,a3)&
           -t1(i2,a1)*t2(i1,i3,a2,a3)&
           +t1(i3,a1)*t2(i1,i2,a2,a3)&
           -t1(i1,a2)*t2(i2,i3,a1,a3)&
           +t1(i2,a2)*t2(i1,i3,a1,a3)&
           -t1(i3,a2)*t2(i1,i2,a1,a3)&
           +t1(i1,a3)*t2(i2,i3,a1,a2)&
           -t1(i2,a3)*t2(i1,i3,a1,a2)&
           +t1(i3,a3)*t2(i1,i2,a1,a2)

    return
    end function c3_t1t2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c3_t3(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2,i3,a3


    ret=0
    if (.NOT.pattern(3)) return

    i1=arr(1); i2=arr(2); i3=arr(3)
    a1=arr(4); a2=arr(5); a3=arr(6)

    ret=ret+t3(i1,i2,i3,a1,a2,a3)

    return
    end function c3_t3

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c4_t1t1t1t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2,i3,a3,i4,a4


    ret=0
    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+t1(i1,a1)*t1(i2,a2)*t1(i3,a3)*t1(i4,a4)&
           -t1(i2,a1)*t1(i1,a2)*t1(i3,a3)*t1(i4,a4)&
           -t1(i1,a1)*t1(i3,a2)*t1(i2,a3)*t1(i4,a4)&
           +t1(i3,a1)*t1(i1,a2)*t1(i2,a3)*t1(i4,a4)&
           -t1(i3,a1)*t1(i2,a2)*t1(i1,a3)*t1(i4,a4)&
           +t1(i2,a1)*t1(i3,a2)*t1(i1,a3)*t1(i4,a4)&
           -t1(i1,a1)*t1(i2,a2)*t1(i4,a3)*t1(i3,a4)&
           +t1(i2,a1)*t1(i1,a2)*t1(i4,a3)*t1(i3,a4)&
           +t1(i1,a1)*t1(i4,a2)*t1(i2,a3)*t1(i3,a4)&
           -t1(i4,a1)*t1(i1,a2)*t1(i2,a3)*t1(i3,a4)&
           +t1(i4,a1)*t1(i2,a2)*t1(i1,a3)*t1(i3,a4)&
           -t1(i2,a1)*t1(i4,a2)*t1(i1,a3)*t1(i3,a4)&
           -t1(i1,a1)*t1(i4,a2)*t1(i3,a3)*t1(i2,a4)&
           +t1(i4,a1)*t1(i1,a2)*t1(i3,a3)*t1(i2,a4)&
           +t1(i1,a1)*t1(i3,a2)*t1(i4,a3)*t1(i2,a4)&
           -t1(i3,a1)*t1(i1,a2)*t1(i4,a3)*t1(i2,a4)&
           +t1(i3,a1)*t1(i4,a2)*t1(i1,a3)*t1(i2,a4)&
           -t1(i4,a1)*t1(i3,a2)*t1(i1,a3)*t1(i2,a4)&
           -t1(i4,a1)*t1(i2,a2)*t1(i3,a3)*t1(i1,a4)&
           +t1(i2,a1)*t1(i4,a2)*t1(i3,a3)*t1(i1,a4)&
           +t1(i4,a1)*t1(i3,a2)*t1(i2,a3)*t1(i1,a4)&
           -t1(i3,a1)*t1(i4,a2)*t1(i2,a3)*t1(i1,a4)&
           +t1(i3,a1)*t1(i2,a2)*t1(i4,a3)*t1(i1,a4)&
           -t1(i2,a1)*t1(i3,a2)*t1(i4,a3)*t1(i1,a4)

    return
    end function c4_t1t1t1t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c4_t1t1t2(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2,i3,a3,i4,a4


    ret=0
    if (.NOT.(pattern(1).AND.pattern(2))) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+t1(i1,a1)*t1(i2,a2)*t2(i3,i4,a3,a4)&
           -t1(i2,a1)*t1(i1,a2)*t2(i3,i4,a3,a4)&
           -t1(i1,a1)*t1(i3,a2)*t2(i2,i4,a3,a4)&
           +t1(i3,a1)*t1(i1,a2)*t2(i2,i4,a3,a4)&
           -t1(i3,a1)*t1(i2,a2)*t2(i1,i4,a3,a4)&
           +t1(i2,a1)*t1(i3,a2)*t2(i1,i4,a3,a4)&
           +t1(i1,a1)*t1(i4,a2)*t2(i2,i3,a3,a4)&
           -t1(i4,a1)*t1(i1,a2)*t2(i2,i3,a3,a4)&
           +t1(i4,a1)*t1(i2,a2)*t2(i1,i3,a3,a4)&
           -t1(i2,a1)*t1(i4,a2)*t2(i1,i3,a3,a4)&
           +t1(i3,a1)*t1(i4,a2)*t2(i1,i2,a3,a4)&
           -t1(i4,a1)*t1(i3,a2)*t2(i1,i2,a3,a4)&
           -t1(i1,a1)*t1(i2,a3)*t2(i3,i4,a2,a4)&
           +t1(i2,a1)*t1(i1,a3)*t2(i3,i4,a2,a4)&
           +t1(i1,a1)*t1(i3,a3)*t2(i2,i4,a2,a4)&
           -t1(i3,a1)*t1(i1,a3)*t2(i2,i4,a2,a4)&
           +t1(i3,a1)*t1(i2,a3)*t2(i1,i4,a2,a4)&
           -t1(i2,a1)*t1(i3,a3)*t2(i1,i4,a2,a4)&
           -t1(i1,a1)*t1(i4,a3)*t2(i2,i3,a2,a4)&
           +t1(i4,a1)*t1(i1,a3)*t2(i2,i3,a2,a4)&
           -t1(i4,a1)*t1(i2,a3)*t2(i1,i3,a2,a4)&
           +t1(i2,a1)*t1(i4,a3)*t2(i1,i3,a2,a4)&
           -t1(i3,a1)*t1(i4,a3)*t2(i1,i2,a2,a4)&
           +t1(i4,a1)*t1(i3,a3)*t2(i1,i2,a2,a4)&
           -t1(i1,a3)*t1(i2,a2)*t2(i3,i4,a1,a4)&
           +t1(i2,a3)*t1(i1,a2)*t2(i3,i4,a1,a4)&
           +t1(i1,a3)*t1(i3,a2)*t2(i2,i4,a1,a4)&
           -t1(i3,a3)*t1(i1,a2)*t2(i2,i4,a1,a4)&
           +t1(i3,a3)*t1(i2,a2)*t2(i1,i4,a1,a4)&
           -t1(i2,a3)*t1(i3,a2)*t2(i1,i4,a1,a4)&
           -t1(i1,a3)*t1(i4,a2)*t2(i2,i3,a1,a4)&
           +t1(i4,a3)*t1(i1,a2)*t2(i2,i3,a1,a4)&
           -t1(i4,a3)*t1(i2,a2)*t2(i1,i3,a1,a4)&
           +t1(i2,a3)*t1(i4,a2)*t2(i1,i3,a1,a4)&
           -t1(i3,a3)*t1(i4,a2)*t2(i1,i2,a1,a4)&
           +t1(i4,a3)*t1(i3,a2)*t2(i1,i2,a1,a4)&
           +t1(i1,a1)*t1(i2,a4)*t2(i3,i4,a2,a3)&
           -t1(i2,a1)*t1(i1,a4)*t2(i3,i4,a2,a3)&
           -t1(i1,a1)*t1(i3,a4)*t2(i2,i4,a2,a3)&
           +t1(i3,a1)*t1(i1,a4)*t2(i2,i4,a2,a3)&
           -t1(i3,a1)*t1(i2,a4)*t2(i1,i4,a2,a3)&
           +t1(i2,a1)*t1(i3,a4)*t2(i1,i4,a2,a3)&
           +t1(i1,a1)*t1(i4,a4)*t2(i2,i3,a2,a3)&
           -t1(i4,a1)*t1(i1,a4)*t2(i2,i3,a2,a3)&
           +t1(i4,a1)*t1(i2,a4)*t2(i1,i3,a2,a3)&
           -t1(i2,a1)*t1(i4,a4)*t2(i1,i3,a2,a3)&
           +t1(i3,a1)*t1(i4,a4)*t2(i1,i2,a2,a3)&
           -t1(i4,a1)*t1(i3,a4)*t2(i1,i2,a2,a3)&
           +t1(i1,a4)*t1(i2,a2)*t2(i3,i4,a1,a3)&
           -t1(i2,a4)*t1(i1,a2)*t2(i3,i4,a1,a3)&
           -t1(i1,a4)*t1(i3,a2)*t2(i2,i4,a1,a3)&
           +t1(i3,a4)*t1(i1,a2)*t2(i2,i4,a1,a3)&
           -t1(i3,a4)*t1(i2,a2)*t2(i1,i4,a1,a3)&
           +t1(i2,a4)*t1(i3,a2)*t2(i1,i4,a1,a3)&
           +t1(i1,a4)*t1(i4,a2)*t2(i2,i3,a1,a3)&
           -t1(i4,a4)*t1(i1,a2)*t2(i2,i3,a1,a3)&
           +t1(i4,a4)*t1(i2,a2)*t2(i1,i3,a1,a3)&
           -t1(i2,a4)*t1(i4,a2)*t2(i1,i3,a1,a3)&
           +t1(i3,a4)*t1(i4,a2)*t2(i1,i2,a1,a3)&
           -t1(i4,a4)*t1(i3,a2)*t2(i1,i2,a1,a3)&
           +t1(i1,a3)*t1(i2,a4)*t2(i3,i4,a1,a2)&
           -t1(i2,a3)*t1(i1,a4)*t2(i3,i4,a1,a2)&
           -t1(i1,a3)*t1(i3,a4)*t2(i2,i4,a1,a2)&
           +t1(i3,a3)*t1(i1,a4)*t2(i2,i4,a1,a2)&
           -t1(i3,a3)*t1(i2,a4)*t2(i1,i4,a1,a2)&
           +t1(i2,a3)*t1(i3,a4)*t2(i1,i4,a1,a2)&
           +t1(i1,a3)*t1(i4,a4)*t2(i2,i3,a1,a2)&
           -t1(i4,a3)*t1(i1,a4)*t2(i2,i3,a1,a2)&
           +t1(i4,a3)*t1(i2,a4)*t2(i1,i3,a1,a2)&
           -t1(i2,a3)*t1(i4,a4)*t2(i1,i3,a1,a2)&
           +t1(i3,a3)*t1(i4,a4)*t2(i1,i2,a1,a2)&
           -t1(i4,a3)*t1(i3,a4)*t2(i1,i2,a1,a2)

    return
    end function c4_t1t1t2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c4_t2t2(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2,i3,a3,i4,a4


    ret=0
    if (.NOT.pattern(2)) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+t2(i1,i2,a1,a2)*t2(i3,i4,a3,a4)&
           -t2(i1,i3,a1,a2)*t2(i2,i4,a3,a4)&
           -t2(i3,i2,a1,a2)*t2(i1,i4,a3,a4)&
           +t2(i1,i4,a1,a2)*t2(i2,i3,a3,a4)&
           +t2(i4,i2,a1,a2)*t2(i1,i3,a3,a4)&
           +t2(i3,i4,a1,a2)*t2(i1,i2,a3,a4)&
           -t2(i1,i2,a1,a3)*t2(i3,i4,a2,a4)&
           +t2(i1,i3,a1,a3)*t2(i2,i4,a2,a4)&
           +t2(i3,i2,a1,a3)*t2(i1,i4,a2,a4)&
           -t2(i1,i4,a1,a3)*t2(i2,i3,a2,a4)&
           -t2(i4,i2,a1,a3)*t2(i1,i3,a2,a4)&
           -t2(i3,i4,a1,a3)*t2(i1,i2,a2,a4)&
           -t2(i1,i2,a3,a2)*t2(i3,i4,a1,a4)&
           +t2(i1,i3,a3,a2)*t2(i2,i4,a1,a4)&
           +t2(i3,i2,a3,a2)*t2(i1,i4,a1,a4)&
           -t2(i1,i4,a3,a2)*t2(i2,i3,a1,a4)&
           -t2(i4,i2,a3,a2)*t2(i1,i3,a1,a4)&
           -t2(i3,i4,a3,a2)*t2(i1,i2,a1,a4)

    return
    end function c4_t2t2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function c4_t1t3(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)

    integer(kind=iglu)             :: i1,a1,i2,a2,i3,a3,i4,a4


    ret=0
    if (.NOT.(pattern(1).AND.pattern(3))) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+t1(i1,a1)*t3(i2,i3,i4,a2,a3,a4)&
           -t1(i2,a1)*t3(i1,i3,i4,a2,a3,a4)&
           +t1(i3,a1)*t3(i1,i2,i4,a2,a3,a4)&
           -t1(i4,a1)*t3(i1,i2,i3,a2,a3,a4)&
           -t1(i1,a2)*t3(i2,i3,i4,a1,a3,a4)&
           +t1(i2,a2)*t3(i1,i3,i4,a1,a3,a4)&
           -t1(i3,a2)*t3(i1,i2,i4,a1,a3,a4)&
           +t1(i4,a2)*t3(i1,i2,i3,a1,a3,a4)&
           +t1(i1,a3)*t3(i2,i3,i4,a1,a2,a4)&
           -t1(i2,a3)*t3(i1,i3,i4,a1,a2,a4)&
           +t1(i3,a3)*t3(i1,i2,i4,a1,a2,a4)&
           -t1(i4,a3)*t3(i1,i2,i3,a1,a2,a4)&
           -t1(i1,a4)*t3(i2,i3,i4,a1,a2,a3)&
           +t1(i2,a4)*t3(i1,i3,i4,a1,a2,a3)&
           -t1(i3,a4)*t3(i1,i2,i4,a1,a2,a3)&
           +t1(i4,a4)*t3(i1,i2,i3,a1,a2,a3)

    return
    end function c4_t1t3

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module coupledClusterAnalize