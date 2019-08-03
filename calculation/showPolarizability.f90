    subroutine showPolarizability(method,state,ou,width)

    use glob     , only: rspu,r16kind,r8kind,i4kind,iglu,rglu,lglu,void,true,false,uch,purifyValues
    use txtParser, only: tpLowerCase,tpUpperCase,tpFill,tpAdjustc
    use printmod , only: prStrByVal
    use harvester, only: hpolariz,init_harvester,store_harvester
    use hdb      , only: HartreeEnergy,BohrRadius,dipoleToDeby,gammaToesu,GlEt
    use hdb      , only: Et,MEt,MMEt,mol,polarizbd,heVersion,generalbd,systembd
    use derivat  , only: deShareParams,deLagDeriv,deLSMDeriv,deFinalize,crt

    implicit none

    character (len=*) , intent(in) :: method
    integer(kind=iglu), intent(in) :: ou,width,state
    type(hpolariz)                 :: harv

    real   (kind=rglu), parameter  :: sigmaBondAlpha=3.48_rglu,CCsigmaBondGamma=220,CHsigmaBondGamma=240
    real   (kind=rglu), parameter  :: derivThreshold(4)=[0.001_rglu,0.01_rglu,0.1_rglu,1._rglu]

    character (len=1)  ::    sX,   sY
    character (len=2)  ::   sXX,  sYY
    character (len=3)  ::  sXXX, sYYY, sXYY, sYXX
    character (len=4)  :: sXXXX,sYYYY,sXXYY,sYYXX
    character (len=3)  :: scale

    integer(kind=iglu) :: i,j,k,l,Nel,Nat,Np,pp,sta,sto
    logical(kind=lglu) :: doharvest

    real   (kind=rglu) :: bstDeriv(4,10)
    real   (kind=rglu) :: dx,dy,dz,dmod
    real   (kind=rglu) :: axx,ayy,azz,avA
    real   (kind=rglu) :: bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod
    real   (kind=rglu) :: gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averG
    real   (kind=rglu) :: sigmaContributionAlpha,sigmaContributionGamma
    real   (kind=rglu) :: esu,PDip,PPol,PHpl1,PHpl2,transition(0:8)


    doharvest = systembd%harvest

    Nel=mol%nEls; Nat=mol%uniqueAtoms; scale=polarizbd%scales%get(); Np=polarizbd%nPoints

    sigmaContributionAlpha=  sigmaBondAlpha*(mol%uniqueBonds+mol%chBonds)
    sigmaContributionGamma=CCsigmaBondGamma*mol%uniqueBonds+CHsigmaBondGamma*mol%chBonds

    sta=-(Np-1)/2; sto= (Np-1)/2

    transition(0)=1; transition(1)=-dipoleToDeby
    do i = 2,UBound(transition,1)
        transition(i)=-HartreeEnergy**(i-1)/BohrRadius**(i)
    enddo
    PDip=transition(1); PPol=transition(2); PHpl1=transition(3); PHpl2=transition(4)
    esu=gammaToesu*1D-36

    void=deShareParams(Np,polarizbd%derivStep,MMEt)

    if (state.NE.0) then
        write (ou,'(/A//)') tpAdjustc(' Polarizability results for '//method//' State #'//prStrByVal(state)//' ', width, '#')
    else
        write (ou,'(/A//)') tpAdjustc(' Polarizability results for '//method//' ', width, '#')
    endif

    if (doharvest) then
        call init_harvester(harv,heVersion,&
                            method,int(state,i4kind),&
                            int(mol%nAtoms,i4kind),int(mol%nEls,i4kind),&
                            int(polarizbd%nPoints,i4kind),real(polarizbd%derivStep,r8kind),&
                            polarizbd%scales%get())

        harv%Energies=GlEt

        harv%sigmaContributionAlpha=sigmaContributionAlpha
        harv%sigmaContributionGamma=sigmaContributionGamma
    endif


    bstDeriv=0

    select case ( len_trim(polarizbd%scales%get()) )

        case (1)
            do
                sX=scale(1:1); sXX=sX//sX; sXXX=sX//sXX; sXXXX=sXX//sXX
                do i = sta,sto
                    write (ou,101) i,Et(i)
                enddo

                write (ou,200) sX
                write (ou,201) (l, l=3,Np,2)
                write (ou,211) [(deLagDeriv(j,tpLowerCase(sX)   ),j=3,Np,2)]*PDip
                write (ou,212) [(deLagDeriv(j,tpLowerCase(sXX)  ),j=3,Np,2)]*PPol
                write (ou,213) [(deLagDeriv(j,tpLowerCase(sXXX) ),j=5,Np,2)]*PHpl1
                write (ou,214) [(deLagDeriv(j,tpLowerCase(sXXXX)),j=5,Np,2)]*PHpl2,&
                               [(deLagDeriv(j,tpLowerCase(sXXXX)),j=5,Np,2)]*PHpl2*esu

                write (ou,281) deLagDeriv(Np,tpLowerCase(sXX))*PPol/3,&
                               deLagDeriv(Np,tpLowerCase(sXX))*PPol/(3*Nel)

                if (doharvest) then
                    select case(polarizbd%scales%get())
                        case('x')
                            harv%m(1)=deLagDeriv(Np,tpLowerCase(sX))*PDip
                            harv%a(1)=deLagDeriv(Np,tpLowerCase(sXX))*PPol
                        case('y')
                            harv%m(2)=deLagDeriv(Np,tpLowerCase(sX))*PDip
                            harv%a(2)=deLagDeriv(Np,tpLowerCase(sXX))*PPol
                        case('z')
                            harv%m(3)=deLagDeriv(Np,tpLowerCase(sX))*PDip
                            harv%a(3)=deLagDeriv(Np,tpLowerCase(sXX))*PPol
                    end select
                endif

                if (Np.LE.3) exit

                if (Nat.NE.Nel) write (ou,282) deLagDeriv(Np,tpLowerCase(sXX))*PPol/(3*Nat)
                write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
                write (ou,283) deLagDeriv(Np,tpLowerCase(sXXXX))*PHpl2/5,&
                               deLagDeriv(Np,tpLowerCase(sXXXX))*PHpl2/(5*Nel)

                if (Nat.NE.Nel) write (ou,284) deLagDeriv(Np,tpLowerCase(sXXXX))*PHpl2/(5*Nat)
                write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma

                if (doharvest) then
                    select case(polarizbd%scales%get())
                        case('x')
                            harv%b(1)=deLagDeriv(Np,tpLowerCase(sXXX))*PHpl1
                            harv%g(1)=deLagDeriv(Np,tpLowerCase(sXXXX))*PHpl2
                        case('y')
                            harv%b(2)=deLagDeriv(Np,tpLowerCase(sXXX))*PHpl1
                            harv%g(2)=deLagDeriv(Np,tpLowerCase(sXXXX))*PHpl2
                        case('z')
                            harv%b(3)=deLagDeriv(Np,tpLowerCase(sXXX))*PHpl1
                            harv%g(3)=deLagDeriv(Np,tpLowerCase(sXXXX))*PHpl2
                    end select
                endif

                exit
            enddo

        case (2)
            scale=tpUpperCase(scale); sX=scale(1:1); sY=scale(2:2)
            sXX=sX//sX; sYY=sY//sY
            sXXX=sX//sXX; sYYY=sY//sYY
            sXYY=sX//sYY; sYXX=sY//sXX
            sXXXX=sXX//sXX; sYYYY=sYY//sYY
            sXXYY=sXX//sYY; sYYXX=sYY//sXX
            scale=tpLowerCase(scale)

            do i = sta,sto
                do j = sta,sto
                    write (ou,102) i,j,MEt(i,j)
                enddo
                write (ou,*)
            enddo

            do ! ~~~~~~~~~ Dipole moment ~~~~~~~~~ !
                if (Np.LT.3) exit
                write (ou,301)

                do k = 3,Np,2
                    dx=deLagDeriv(k,tpLowerCase(sX))*PDip
                    dy=deLagDeriv(k,tpLowerCase(sY))*PDip
                    dmod=sqrt(dx**2+dy**2)

                    void=purifyValues(derivThreshold(1),dx,dy,dmod)

                    write (ou,311) k,sX,dx,sY,dy,dmod

                    bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,4)=dmod
                enddo
                exit
            enddo

            do ! ~~~~~~~~~ Polarizability ~~~~~~~~~ !
                if (Np.LT.3) exit
                write (ou,321)

                do k = 3,Np,2
                    axx=deLagDeriv(k,tpLowerCase(sXX))*PPol
                    ayy=deLagDeriv(k,tpLowerCase(sYY))*PPol
                    avA=(axx+ayy)/3

                    void=purifyValues(derivThreshold(2),axx,ayy,ava)

                    write (ou,331) k,sXX,axx,sYY,ayy,avA

                    bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,4)=avA
                enddo
                if (Nel.NE.Nat) write (ou,391) avA/Nat

                write (ou,392) avA/Nel,myCont([axx,ayy],3)
                write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
                exit
            enddo

            do ! ~~~~~~~~~ 1st Hyperpolarizability ~~~~~~~~~ !
                if (Np.LT.5) exit
                write (ou,341)

                do k = 5,Np,2
                    bxxx=deLagDeriv(k,tpLowerCase(sXXX))*PHpl1
                    byyy=deLagDeriv(k,tpLowerCase(sYYY))*PHpl1
                    bxyy=deLagDeriv(k,tpLowerCase(sXYY))*PHpl1
                    byxx=deLagDeriv(k,tpLowerCase(sYXX))*PHpl1
                    bmod=sqrt((bxxx+bxyy)**2+(byyy+byxx)**2)

                    void=purifyValues(derivThreshold(3),bxxx,byyy,bxyy,byxx,bmod)

                    write (ou,351) k,sXXX,bxxx,sYYY,byyy,sXYY,bxyy,sYXX,byxx,bmod

                    bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,4)=bxyy
                    bstDeriv(3,5)=byxx; bstDeriv(3,10)=bmod
                enddo
                exit
            enddo

            do ! ~~~~~~~~~ 2nd Hyperpolarizability ~~~~~~~~~ !
                if (Np.LT.5) exit
                write (ou,361)

                do k = 5,Np,2
                    gxxxx=deLagDeriv(k,tpLowerCase(sXXXX))*PHpl2
                    gyyyy=deLagDeriv(k,tpLowerCase(sYYYY))*PHpl2
                    gxxyy=deLagDeriv(k,tpLowerCase(sXXYY))*PHpl2
                    gyyxx=deLagDeriv(k,tpLowerCase(sYYXX))*PHpl2
                    averg=(gxxxx+gyyyy+gxxyy+gyyxx)/5

                    void=purifyValues(derivThreshold(4),gxxxx,gyyyy,gxxyy,gyyxx,averg)

                    bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,4)=gxxyy
                    bstDeriv(4,5)=gyyxx; bstDeriv(4,10)=averg

                    write (ou,371) k,sXXXX,gxxxx,sYYYY,gyyyy,sXXYY,gxxyy,sYYXX,gyyxx,averg,averg*esu
                enddo

                if (Nel.NE.Nat) write (ou,393) averg/Nat
                write (ou,394) averg/Nel,myCont([gxxxx,gyyyy,gxxyy,gyyxx],5)
                write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
                exit
            enddo

            if (doharvest) then
                do
                    dx=deLagDeriv(Np,tpLowerCase(sX))*PDip
                    dy=deLagDeriv(Np,tpLowerCase(sY))*PDip
                    dmod=sqrt(dx**2+dy**2)
                    void=purifyValues(derivThreshold(1),dx,dy,dmod)

                    axx=deLagDeriv(Np,tpLowerCase(sXX))*PPol
                    ayy=deLagDeriv(Np,tpLowerCase(sYY))*PPol
                    avA=(axx+ayy)/3
                    void=purifyValues(derivThreshold(2),axx,ayy,ava)

                    select case(polarizbd%scales%get())
                        case('xy')
                            harv%m(1)=dx;  harv%m(2)=dy;  harv%m(4)=dmod
                            harv%a(1)=axx; harv%a(2)=ayy; harv%a(4)=ava

                        case('xz')
                            harv%m(1)=dx;  harv%m(3)=dy;  harv%m(4)=dmod
                            harv%a(1)=axx; harv%a(3)=ayy; harv%a(4)=ava

                        case('yz')
                            harv%m(2)=dx;  harv%m(3)=dy;  harv%m(4)=dmod
                            harv%a(2)=axx; harv%a(3)=ayy; harv%a(4)=ava

                    end select

                    if (Np.LT.5) exit

                    bxxx=deLagDeriv(Np,tpLowerCase(sXXX))*PHpl1
                    byyy=deLagDeriv(Np,tpLowerCase(sYYY))*PHpl1
                    bxyy=deLagDeriv(Np,tpLowerCase(sXYY))*PHpl1
                    byxx=deLagDeriv(Np,tpLowerCase(sYXX))*PHpl1
                    bmod=sqrt((bxxx+bxyy)**2+(byyy+byxx)**2)

                    void=purifyValues(derivThreshold(3),bxxx,byyy,bxyy,byxx,bmod)

                    gxxxx=deLagDeriv(Np,tpLowerCase(sXXXX))*PHpl2
                    gyyyy=deLagDeriv(Np,tpLowerCase(sYYYY))*PHpl2
                    gxxyy=deLagDeriv(Np,tpLowerCase(sXXYY))*PHpl2
                    gyyxx=deLagDeriv(Np,tpLowerCase(sYYXX))*PHpl2
                    averg=(gxxxx+gyyyy+gxxyy+gyyxx)/5

                    void=purifyValues(derivThreshold(4),gxxxx,gyyyy,gxxyy,gyyxx,averg)

                    select case(polarizbd%scales%get())
                        case('xy')
                            harv%b(1)=bxxx;  harv%b(2)=byyy;  harv%b(4)=bxyy;  harv%b(5)=byxx;  harv%b(10)=bmod
                            harv%g(1)=gxxxx; harv%g(2)=gyyyy; harv%g(4)=gxxyy; harv%g(5)=gyyxx; harv%g(10)=averg

                        case('xz')
                            harv%b(1)=bxxx;  harv%b(3)=byyy;  harv%b(6)=bxyy;  harv%b(7)=byxx;  harv%b(10)=bmod
                            harv%g(1)=gxxxx; harv%g(3)=gyyyy; harv%g(6)=gxxyy; harv%g(7)=gyyxx; harv%g(10)=averg

                        case('yz')
                            harv%b(2)=bxxx;  harv%b(3)=byyy;  harv%b(8)=bxyy;  harv%b(9)=byxx;  harv%b(10)=bmod
                            harv%g(2)=gxxxx; harv%g(3)=gyyyy; harv%g(8)=gxxyy; harv%g(9)=gyyxx; harv%g(10)=averg

                    end select

                    exit
                enddo
            endif

        case (3)
            do k = sta,sto
            do i = sta,sto
            do j = sta,sto
                if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle
                write (ou,104) i,j,k,MMEt(i,j,k)
            enddo
            enddo
            enddo
            write (ou,105)

            do ! ~~~~~~~~~ Dipole moment ~~~~~~~~~ !
                if (Np.LT.3) exit
                write (ou,401)

                do k = 3,Np,2
                    dx=deLagDeriv(k,'x')*PDip
                    dy=deLagDeriv(k,'y')*PDip
                    dz=deLagDeriv(k,'z')*PDip
                    dmod=sqrt(dx**2+dy**2+dz**2)

                    void=purifyValues(derivThreshold(1),dx,dy,dz,dmod)

                    write (ou,411) k,dx,dy,dz,dmod
                    bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,3)=dz; bstDeriv(1,4)=dmod
                enddo
                exit
            enddo

            do ! ~~~~~~~~~ Polarizability ~~~~~~~~~ !
                if (Np.LT.3) exit
                write (ou,421)

                do k = 3,Np,2
                    axx=deLagDeriv(k,'xx')*PPol
                    ayy=deLagDeriv(k,'yy')*PPol
                    azz=deLagDeriv(k,'zz')*PPol
                    avA=(axx+ayy+azz)/3

                    void=purifyValues(derivThreshold(2),axx,ayy,azz,ava)

                    write (ou,431) k,axx,ayy,azz,avA

                    bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,3)=azz; bstDeriv(2,4)=avA
                enddo

                if (Nel.NE.Nat) then
                    write (ou,491) avA/Nat
                endif

                write (ou,492) avA/Nel,myCont([axx,ayy,azz],3)
                write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
                exit
            enddo

            do ! ~~~~~~~~~ 1st Hyperpolarizability ~~~~~~~~~ !
                if (Np.LT.5) exit
                write (ou,441)

                do k = 5,Np,2
                    bxxx = deLagDeriv(k,'xxx')*PHpl1
                    byyy = deLagDeriv(k,'yyy')*PHpl1
                    bzzz = deLagDeriv(k,'zzz')*PHpl1
                    bxyy = deLagDeriv(k,'xyy')*PHpl1
                    byxx = deLagDeriv(k,'yxx')*PHpl1
                    bxzz = deLagDeriv(k,'xzz')*PHpl1
                    bzxx = deLagDeriv(k,'zxx')*PHpl1
                    byzz = deLagDeriv(k,'yzz')*PHpl1
                    bzyy = deLagDeriv(k,'zyy')*PHpl1
                    bmod=sqrt((bxxx+bxyy+bxzz)**2+(byyy+byxx+byzz)**2+(bzzz+bzxx+bzyy)**2)

                    void=purifyValues(derivThreshold(3),bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod)

                    write (ou,451) k,bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod

                    bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,3)=bzzz; bstDeriv(3,4)=bxyy
                    bstDeriv(3,5)=byxx; bstDeriv(3,6)=bxzz; bstDeriv(3,7)=bzxx; bstDeriv(3,8)=byzz
                    bstDeriv(3,9)=bzyy; bstDeriv(3,10)=bmod
                enddo
                exit
            enddo

            do ! ~~~~~~~~~ 2nd Hyperpolarizability ~~~~~~~~~ !
                if (Np.LT.5) exit
                write (ou,461)

                do k = 5,Np,2
                    gxxxx=deLagDeriv(k,'xxxx')*PHpl2
                    gyyyy=deLagDeriv(k,'yyyy')*PHpl2
                    gzzzz=deLagDeriv(k,'zzzz')*PHpl2
                    gxxyy=deLagDeriv(k,'xxyy')*PHpl2
                    gyyxx=deLagDeriv(k,'yyxx')*PHpl2
                    gxxzz=deLagDeriv(k,'xxzz')*PHpl2
                    gzzxx=deLagDeriv(k,'zzxx')*PHpl2
                    gyyzz=deLagDeriv(k,'yyzz')*PHpl2
                    gzzyy=deLagDeriv(k,'zzyy')*PHpl2
                    averg=(gxxxx+gyyyy+gzzzz+gxxyy+gyyxx+gxxzz+gzzxx+gyyzz+gzzyy)/5

                    void=purifyValues(derivThreshold(4),gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averg)

                    write (ou,471) k,gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averg,averg*esu

                    bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,3)=gzzzz; bstDeriv(4,4)=gxxyy
                    bstDeriv(4,5)=gyyxx; bstDeriv(4,6)=gxxzz; bstDeriv(4,7)=gzzxx; bstDeriv(4,8)=gyyzz
                    bstDeriv(4,9)=gzzyy; bstDeriv(4,10)=averg
                enddo

                if (Nel.NE.Nat) then
                    write (ou,493) averg/Nat
                endif
                write (ou,494) averg/Nel,myCont([gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy],5)
                write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
                exit
            enddo

            if (doharvest) then
                do
                    dx=deLagDeriv(Np,'x')*PDip
                    dy=deLagDeriv(Np,'y')*PDip
                    dz=deLagDeriv(Np,'z')*PDip
                    dmod=sqrt(dx**2+dy**2+dz**2)

                    void=purifyValues(derivThreshold(1),dx,dy,dz,dmod)

                    axx=deLagDeriv(Np,'xx')*PPol
                    ayy=deLagDeriv(Np,'yy')*PPol
                    azz=deLagDeriv(Np,'zz')*PPol
                    avA=(axx+ayy+azz)/3

                    void=purifyValues(derivThreshold(2),axx,ayy,azz,ava)

                    harv%m(1)=dx;  harv%m(2)=dy;  harv%m(3)=dz;  harv%m(4)=dmod
                    harv%a(1)=axx; harv%a(2)=ayy; harv%a(3)=azz; harv%a(4)=ava

                    if (Np.LT.5) exit

                    bxxx = deLagDeriv(Np,'xxx')*PHpl1
                    byyy = deLagDeriv(Np,'yyy')*PHpl1
                    bzzz = deLagDeriv(Np,'zzz')*PHpl1
                    bxyy = deLagDeriv(Np,'xyy')*PHpl1
                    byxx = deLagDeriv(Np,'yxx')*PHpl1
                    bxzz = deLagDeriv(Np,'xzz')*PHpl1
                    bzxx = deLagDeriv(Np,'zxx')*PHpl1
                    byzz = deLagDeriv(Np,'yzz')*PHpl1
                    bzyy = deLagDeriv(Np,'zyy')*PHpl1
                    bmod=sqrt((bxxx+bxyy+bxzz)**2+(byyy+byxx+byzz)**2+(bzzz+bzxx+bzyy)**2)

                    void=purifyValues(derivThreshold(3),bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod)

                    gxxxx=deLagDeriv(Np,'xxxx')*PHpl2
                    gyyyy=deLagDeriv(Np,'yyyy')*PHpl2
                    gzzzz=deLagDeriv(Np,'zzzz')*PHpl2
                    gxxyy=deLagDeriv(Np,'xxyy')*PHpl2
                    gyyxx=deLagDeriv(Np,'yyxx')*PHpl2
                    gxxzz=deLagDeriv(Np,'xxzz')*PHpl2
                    gzzxx=deLagDeriv(Np,'zzxx')*PHpl2
                    gyyzz=deLagDeriv(Np,'yyzz')*PHpl2
                    gzzyy=deLagDeriv(Np,'zzyy')*PHpl2
                    averg=(gxxxx+gyyyy+gzzzz+gxxyy+gyyxx+gxxzz+gzzxx+gyyzz+gzzyy)/5

                    void=purifyValues(derivThreshold(4),gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averg)

                    harv%b( 1)=bxxx; harv%b( 2)=byyy; harv%b( 3)=bzzz; harv%b( 4)=bxyy; harv%b( 5)=byxx
                    harv%b( 6)=bxzz; harv%b( 7)=bzxx; harv%b( 8)=byzz; harv%b( 9)=bzyy; harv%b(10)=bmod

                    harv%g( 1)=gxxxx; harv%g( 2)=gyyyy; harv%g( 3)=gzzzz; harv%g( 4)=gxxyy; harv%g( 5)=gyyxx
                    harv%g( 6)=gxxzz; harv%g( 7)=gzzxx; harv%g( 8)=gyyzz; harv%g( 9)=gzzyy; harv%g(10)=averg
                    exit
                enddo
            endif

    end select

    if (Np.GE.5) void=lsmInterpolation( len_trim(polarizbd%scales%get()) )
    if (Np.GE.3) write (ou,502)
    if (Np.GE.5) write (ou,503)

    call deFinalize

    write (ou,'(/A)') tpFill(width, '#')

    if (doharvest) then
        call store_harvester(harv,generalbd%harvestfile%get())
    endif
    return

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

101 format ( 4X,'Field',1X,i2,1X            ,F25.16)
102 format ( 4X,'Field',1X,i2,1X,i2,1X      ,F25.16)
104 format ( 4X,'Field',1X,i2,1X,i2,1X,i2,1X,F25.16)
105 format (/4X,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'/)
111 format (//4X,A/4X,A,1X,'polarizability',1X,A,1X,A/4X,A//)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

200 format (///4X,'~~~~~~~~~~ Axis ',A1,' ~~~~~~~~~~'/)

201 format (2X,24X,9(6X,i2,' points'))
211 format (2X,'Dipole moment           deby',1X,ES14.7,1X,ES14.7,1X,ES14.7,1X,ES14.7)
212 format (2X,'Polarizability          au  ',1X,ES14.7,1X,ES14.7,1X,ES14.7,1X,ES14.7)
213 format (2X,'1st hyperpolarizability au  ',1X,15X,      ES14.7,1X,ES14.7)
214 format (2X,'2nd hyperpolarizability au  ',1X,15X,      ES14.7,1X,ES14.7/&
            2X,'2nd hyperpolarizability esu ',1X,15X,      ES14.7,1X,ES14.7)

281 format (/'         Average polarizability                =',1X,F13.4,1X,'au'/&
             'Specific average polarizability (per electron) =',1X,F13.4,1X,'au' )
282 format ( 'Specific average polarizability (per atom)     =',1X,F13.4,1X,'au'/)

283 format (/'         Average 2nd hyperpolarizability                =',1X,ES14.7,1X,'au'/&
             'Specific average 2nd hyperpolarizability (per electron) =',1X,ES14.7,1X,'au')
284 format ( 'Specific average 2nd hyperpolarizability (per atom)     =',1X,ES14.7,1X,'au')

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

301 format (/ 7X,'Dipole moment (deby)'/)
311 format (  2X,i2,' points',3X,A1,' =',1X,F7.3,4X,A1,' =',1X,F7.3,4X,'|D| =',1X,F7.3)

321 format (/ 7X,'Polarizability (au)'/)
331 format (  2X,i2,' points',2X,A2,' =',1X,ES14.7,2X,A2,' =',1X,ES14.7,2X,'<A> =',1X,ES14.7)

391 format (/2X,'Specific average polarizability (per atom)     =',1X,ES14.7,1X,'au')
392 format (/2X,'Specific average polarizability (per electron) =',1X,ES14.7,1X,'au'/&
             2X,'Contribution of general component =',3X,F6.2,'%')

341 format (//7X,'1st hyperpolarizability (au)')
351 format ( /2X,i2,' points',2X,A3,' =',1X,ES14.7,3X,A3,' =',1X,ES14.7/&
             13X,                A3,' =',1X,ES14.7,3X,A3,' =',1X,ES14.7,3X,'|B| =',1X,ES14.7)

361 format (//7X,'2nd hyperpolarizability (au)')
371 format ( /2X,i2,' points',2X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
                             13X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
                             14X,   '<G> =',2X,ES14.7,10X,'=',2X,ES14.7,' esu')

393 format (/2X,'Specific average 2nd hyperpolarizability (per atom)     =',9X,ES14.7,1X,'au')
394 format (/2X,'Specific average 2nd hyperpolarizability (per electron) =',9X,ES14.7,1X,'au'/&
             2X,'Contribution of general component =',3X,F6.2,'%')

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

401 format (/ 7X,'Dipole moment (deby)'/)
411 format (  2X,i2,' points',3X,'X = ',F7.3,4X,'Y = ',F7.3,4X,'Z = ',F7.3,4X,'|D| = ',F7.3)

421 format (/ 7X,   'Polarizability (au)'/)
431 format (  2X,i2,' points  XX = ',ES14.7,'  YY = ',ES14.7,'  ZZ = ',ES14.7,'  <A> = ',ES14.7)

491 format (/2X,'Specific Average Polarizability (per atom)     =',1X,ES14.7,1X,'au')
492 format (/2X,'Specific Average Polarizability (per electron) =',1X,ES14.7,1X,'au'/&
             2X,'Contribution of general component =',1X,F6.2,'%')

441 format (//7X,   '1st hyperpolarizability (au)')
451 format ( /2X,i2,' points',3X,'XXX =',1X,ES14.7,3X,'YYY =',1X,ES14.7,3X,'ZZZ =',1X,ES14.7/&
              2X,            12X,'XYY =',1X,ES14.7,3X,'YXX =',1X,ES14.7/&
              2X,            12X,'XZZ =',1X,ES14.7,3X,'ZXX =',1X,ES14.7/&
              2X,            12X,'YZZ =',1X,ES14.7,3X,'ZYY =',1X,ES14.7,3X,'|B| =',1X,ES14.7)

461 format (//7X   ,'2nd hyperpolarizability (au)')
471 format ( /2X,i2,' points',2X,'XXXX =',2X,ES14.7,5X,'YYYY =',2X,ES14.7,5X,'ZZZZ =',2X,ES14.7/&
                             13X,'XXYY =',2X,ES14.7,5X,'YYXX =',2X,ES14.7/&
                             13X,'XXZZ =',2X,ES14.7,5X,'ZZXX =',2X,ES14.7/&
                             13X,'YYZZ =',2X,ES14.7,5X,'ZZYY =',2X,ES14.7/&
                             13X,' <G> =',2X,ES14.7,10X,    '=',2X,ES14.7,' esu')

493 format (/2X,'Specific average 2nd hyperpolarizability (per atom)     =',1X,ES14.7,1X,'au')
494 format (/2X,'Specific average 2nd hyperpolarizability (per electron) =',1X,ES14.7,1X,'au'/&
             2X,'Contribution of general component =',1X,F6.2,'%')

500 format (/2X,'Sigma contribution to average polarizability =',1X,F11.2/&
             2X,F5.3,' au per C-C and C-H bonds [1]')
501 format (/2X,'Sigma contribution to average 2nd hyperpolarizability =',1X,ES10.3/&
             2X,F5.0,' au per C-C bond [2]'/&
             2X,F5.0,' au per C-H bond [2]')
502 format (/2X,'________________________________________________________'/&
             2X,'[1] Yu.B. Visotski and B.C. Briantsev, Quantum Chemistry'/&
             2X,'    of Radicals and Ion-Radicals with Conjugated Bonds'/&
             2X,'    (Donetsk, DonGuet, 2004), p. 207, in Russian.')
503 format ( 2X,'[2] B.M. Pierce, J. Chem. Phys. 91, 791 (1989).')

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function lsmInterpolation(scl) result(ret)

    implicit none

    integer(kind=iglu), intent(in) :: scl
    character                      :: line*510,value*100,pFormat*9,tmplt*4

    integer(kind=iglu)             :: j,k,l,m,sc1,sc2,sc3,pAccuracy,overAllLen,cLen,scSet(3,3)
    real   (kind=rglu)             :: cofs(0:4),avar


    ret=0; if (rspu.NE.r16kind) return

    write (ou,'(/A)') tpAdjustc('LSM INSET TO CONTROL NUMERICAL STEADINESS',width,'=')

    cofs(0)=1; cofs(1)=PDip; cofs(2)=PPol; cofs(3)=PHpl1; cofs(4)=PHpl2

    pAccuracy=5; overAllLen=4+4+pAccuracy
    pFormat='(ES'//prStrByVal(overAllLen)//'.'//prStrByVal(pAccuracy)//')'
    overAllLen=overAllLen+1

    select case (scl)

        case(1)
            sc1=index( crt,scale(1:1) )
            line=tpFill(line)
            do k = 0,4
                line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
            enddo
            cLen=5+5*overAllLen
            write (ou,'(A)') line(1:cLen)

            do k = 0,4
                avar=transition(k)*deLSMDeriv(repeat(crt(sc1:sc1),k))
                if (abs(avar).LT.1D-5) avar=0
                value=tpFill(value); write (value,fmt=pFormat) avar
                line(5+1+k*overAllLen:5+(k+1)*overAllLen)=value(1:overAllLen)
            enddo
            write (ou,'(A)') line(1:cLen)

        case(2)
            sc1=index( crt,scale(1:1) )
            sc2=index( crt,scale(2:2) )

            select case (scale(1:2))
                case ('xy'); sc3=3
                case ('xz'); sc3=2
                case ('yz'); sc3=1
            end select

            line=tpFill(line)

            do k = 0,4
                line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
            enddo
            cLen=5+5*overAllLen

            write (ou,'(A)') line(1:cLen)

            do l = 0,4
                line=tpFill(line); tmplt=tpFill(tmplt)
                tmplt=tpFill(l,crt(sc2:sc2))

                line(2:5)=adjustr( tmplt )
                do k = 0,4
                    if (k+l.LE.4) then
                        avar=transition(k+l)*deLSMDeriv(repeat(scale(1:1),k)//repeat(scale(2:2),l))
                        if (abs(avar).LT.1D-5) avar=0
                        value=tpFill(value); write (value,fmt=pFormat) avar
                        line(5+1+k*overAllLen:5+(k+1)*overAllLen)=value(1:overAllLen)
                    else
                        line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpFill(overAllLen)
                    endif
                enddo

                write (ou,'(A)') line(1:cLen)
            enddo

        case(3)
            scSet(1,1)=1; scSet(2,1)=2; scSet(3,1)=3 !xy
            scSet(1,2)=1; scSet(2,2)=3; scSet(3,2)=2 !xz
            scSet(1,3)=2; scSet(2,3)=3; scSet(3,3)=1 !yz
            do j = 1,3
                sc1=scSet(1,j); sc2=scSet(2,j); sc3=scSet(3,j); line=tpFill(line)

                do k = 0,4
                    line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
                enddo
                cLen=5+5*overAllLen

                write (ou,'(A)') line(1:cLen)

                m=0
                do l = 0,4
                    line=tpFill(line); tmplt=tpFill(tmplt)
                    tmplt=tpFill(l,crt(sc2:sc2))

                    line(2:5)=adjustr( tmplt )
                    do k = 0,4
                        if (k+l.LE.4) then
                            avar=transition(k+l)*deLSMDeriv(tpFill(k,crt(sc1:sc1))//tpFill(l,crt(sc2:sc2))//tpFill(m,crt(sc3:sc3)))
                            if (abs(avar).LT.1D-5) avar=0
                            value=tpFill(value); write (value,fmt=pFormat) avar
                            line(5+1+k*overAllLen:5+(k+1)*overAllLen)=value(1:overAllLen)
                        else
                            line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpFill(overAllLen)
                        endif
                    enddo

                    write (ou,'(A)') line(1:cLen)
                enddo
            enddo

    end select

    write (ou,'(A)') tpFill(width,'=')

    return
    end function lsmInterpolation

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function myCont(arr,divisor) result(ret)
    implicit none

    real   (kind=rglu), intent(in) :: arr(:)
    integer(kind=iglu), intent(in) :: divisor


    ret=100*maxval(abs(arr))/ sum(abs(arr))

    return
    end function myCont

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end subroutine showPolarizability