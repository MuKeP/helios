    subroutine wfAnalize(method, tdlist, prHeader)

    use glob,                  only: rglu,iglu,lglu,void,gluCompare,mid,true
    use hdb,                   only: mol,ou,setfLineLen,ccbd
    use printmod,              only: prMatrix,prVector
    use txtParser,             only: tpAdjustc
    use coupledClusterAnalize

    implicit none

    character(len=*)  , intent(in) :: method
    integer(kind=iglu), intent(in) :: tdlist
    logical(kind=lglu)             :: chaps(chapsLen)

    logical(kind=lglu), optional   :: prHeader
    logical(kind=lglu)             :: dprHeader

    character (len=1)  :: sh(6)

    integer(kind=iglu) :: ilen,ulen,zlen,t1plen,t2plen,dlayer,cChap
    integer(kind=iglu) :: i,j,a,b,k,l,c,d,e,ui,ua,uj,ub,uk,uc,jb(2,2),kc(2,2),mm,nn,pp,vv,gg
    integer(kind=iglu) :: count1,mcountt1,mcountt2,mcountt3,currl1,currl2,currl3

    real(kind=rglu)    :: Ax,Bx,Cx,Dx,Ex,Fx,Gx,Hx,Ix,Jx,cvalue,value,&
                          holdDist(28),hold(10),sumw(4),sss(4)
    real(kind=rglu)    :: fullt(3),energyContrib(2),cval(2),cor(9),cmv(4),birad(4)
    real(kind=rglu)    :: lav,rav,birav,birav2,dist12,dist23,dist13,divisor,gsum

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ pars todolist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! 2#01111111
    !    7654321

    ! 1 | print t1/t2/t3 by layers    | 1
    ! 2 | contribution to energy      | 2
    ! 3 | structure analize by layers | 4
    ! 4 | C1,C2,C3,C4 estimation      | 8
    ! 5 | biradical character         | 16
    ! 6 | statistics on t1/t2/t3      | 32
    ! 7 | configurational compound    | 64
    !     1&2&3&4&5&6&7 = 127 = 2#01111111

    dprHeader=true; if (present(prHeader)) dprHeader=prHeader

    do i = 1,chapsLen
        chaps(i)=btest(tdlist,i-1)
    enddo

    if (tdlist.NE.0) then
        if (dprHeader) write (ou,'(//A/)') tpAdjustc('Wave-function analysis',setfLineLen,'=')
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    umethod=method
    void=prepareForAnalize()

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=1
    if (chaps(cChap).AND.pattern(1)) then

    if (dCUE) then
        write (ou,100) Nse ! print t1 by layers
        do k = 1,maxNei
            count1=0
            write (ou,107) k

            do i = 1,Nel
            do a = Nel+1,No
                if (inLayer(distMO(i,a),k)) then
                    if (.NOT.btest(i+a,0)) then
                        if (btest(i,0)) then; sh(1)='a'; else; sh(1)='b'; endif
                        ui=i/2+mod(i,2); ua=(a-Nel)/2+mod(a-Nel,2)

                        if (abs(t1(i,a)).LT.gluCompare) then
                              write (ou,105) sh(1),ui,ua,     0.
                        else; write (ou,105) sh(1),ui,ua,t1(i,a)
                        endif
                        count1=count1+1
                    endif
                end if
            enddo
            enddo

            write (ou,108) count1
        enddo
    else
        write (ou,190) Nse,ccbd%printThreshold
        write (ou,'(/)')
        do i = 1,Nel
        do a = Nel+1,No
            if (.NOT.btest(i+a,0)) then
                if (btest(i,0)) then; sh(1)='a'; else; sh(1)='b'; endif
                ui=i/2+mod(i,2); ua=(a-Nel)/2+mod(a-Nel,2)

                if (abs(t1(i,a)).GT.ccbd%printThreshold) then
                    write (ou,105) sh(1),ui,ua,t1(i,a)
                endif
                count1=count1+1
            endif
        enddo
        enddo
    endif

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=1
    if (chaps(cChap).AND.pattern(2)) then

    if (dCUE) then
        write (ou,101) (Nse*(Nse-1))/2 ! print t2 by layers
        do l = 1,maxNei
            count1=0

            write (ou,107) l

            do mm = 1,Nse
                i=spinInds(mm,1)
                a=spinInds(mm,2)
                do nn = 1,Nse
                    j=spinInds(nn,1)
                    b=spinInds(nn,2)

                    holdDist(1)=distMO(i,j); holdDist(2)=distMO(i,a)
                    holdDist(3)=distMO(i,b); holdDist(4)=distMO(j,a)
                    holdDist(5)=distMO(j,b); holdDist(6)=distMO(a,b)
                    cmv(1)=distMO(i,a); cmv(2)=distMO(j,b); cmv(3)=maxval(holdDist(1:6))

                    do k = 1,maxNei
                        if (cmv(1).LT.mol%cueLayers(k)) then
                            currl1=k; exit
                        endif
                    enddo
                    do k = 1,maxNei
                        if (cmv(2).LT.mol%cueLayers(k)) then
                            currl2=k; exit
                        endif
                    enddo

                    cor(1)=(cMO(i,1)+cMO(a,1))/2
                    cor(2)=(cMO(i,2)+cMO(a,2))/2
                    cor(3)=(cMO(i,3)+cMO(a,3))/2

                    cor(4)=(cMO(j,1)+cMO(b,1))/2
                    cor(5)=(cMO(j,2)+cMO(b,2))/2
                    cor(6)=(cMO(j,3)+cMO(b,3))/2

                    dist12=sqrt( (cor(1)-cor(4))**2 + (cor(2)-cor(5))**2 + (cor(3)-cor(6))**2 )

                    do k = 1,maxNei
                        if (dist12.LT.mol%cueLayers(k)) then
                            d=k; exit
                        endif
                    enddo

                    if (inLayer(cmv(3),l)) then
                        if (.NOT.btest(i+a+j+b,0)) then
                            if (btest(i,0)) then; sh(1)='a'; else; sh(1)='b'; endif
                            if (btest(a,0)) then; sh(2)='a'; else; sh(2)='b'; endif
                            if (btest(j,0)) then; sh(3)='a'; else; sh(3)='b'; endif
                            if (btest(b,0)) then; sh(4)='a'; else; sh(4)='b'; endif

                            ui=i/2+mod(i,2); ua=(a-Nel)/2+mod(a-Nel,2)
                            uj=j/2+mod(j,2); ub=(b-Nel)/2+mod(b-Nel,2)

                            cval(1)=t2(i,j,a,b)
                            cval(2)=twoeintegr(i,j,a,b)*multiplier
                            if (abs(cval(1)*cval(2)).GT.gluCompare) then
                                if (abs(cval(1)).LT.gluCompare) then
                                      write (ou,106) sh(1),sh(2),sh(3),sh(4),ui,ua,uj,ub,   0   ,cval(2),           0   ,currl1,currl2,d
                                else; write (ou,106) sh(1),sh(2),sh(3),sh(4),ui,ua,uj,ub,cval(1),cval(2),cval(1)*cval(2),currl1,currl2,d
                                endif
                                count1=count1+1
                            endif
                        endif
                    end if
                enddo
            enddo

            write (ou,109) count1
        enddo
    else
        write (ou,191) (Nse*(Nse-1))/2,ccbd%printThreshold
        do mm = 1,Nse
            i=spinInds(mm,1)
            a=spinInds(mm,2)
            do nn = 1,Nse
                j=spinInds(nn,1)
                b=spinInds(nn,2)

                if (.NOT.btest(i+a+j+b,0)) then
                    if (btest(i,0)) then; sh(1)='a'; else; sh(1)='b'; endif
                    if (btest(a,0)) then; sh(2)='a'; else; sh(2)='b'; endif
                    if (btest(j,0)) then; sh(3)='a'; else; sh(3)='b'; endif
                    if (btest(b,0)) then; sh(4)='a'; else; sh(4)='b'; endif

                    ui=i/2+mod(i,2); ua=(a-Nel)/2+mod(a-Nel,2)
                    uj=j/2+mod(j,2); ub=(b-Nel)/2+mod(b-Nel,2)

                    cval(1)=t2(i,j,a,b)
                    cval(2)=twoeintegr(i,j,a,b)*multiplier
                    if ((abs(cval(1)).GT.ccbd%printThreshold).AND.(abs(cval(2)).GT.ccbd%printThreshold)) then
                        write (ou,196) sh(1),sh(2),sh(3),sh(4),ui,ua,uj,ub,cval(1),cval(2),cval(1)*cval(2)
                    endif
                end if
            enddo
        enddo
    endif

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=1
    if (chaps(cChap).AND.pattern(3)) then

    write (ou,111) (Nse*(Nse-1)*(Nse-1))/2,ouThreshold ! print t3 by layers
    do l = 2,maxNei
        count1=0

        write (ou,117) l

        do mm = 1,Nse
            i=spinInds(mm,1)
            a=spinInds(mm,2)
            do nn = 1,Nse
                j=spinInds(nn,1)
                b=spinInds(nn,2)
                do pp = 1,Nse
                    k=spinInds(pp,1)
                    c=spinInds(pp,2)

                    holdDist( 1)=distMO(i,j); holdDist( 2)=distMO(i,a)
                    holdDist( 3)=distMO(i,b); holdDist( 4)=distMO(j,a)
                    holdDist( 5)=distMO(j,b); holdDist( 6)=distMO(a,b)
                    holdDist( 7)=distMO(i,k); holdDist( 8)=distMO(i,c)
                    holdDist( 9)=distMO(j,k); holdDist(10)=distMO(j,c)
                    holdDist(11)=distMO(a,k); holdDist(12)=distMO(a,c)
                    holdDist(13)=distMO(b,k); holdDist(14)=distMO(b,c)
                    holdDist(15)=distMO(k,c)

                    cmv(1)=distMO(i,a); cmv(2)=distMO(j,b); cmv(3)=distMO(k,c); cmv(4)=maxval(holdDist(1:15))

                    do e = 1,maxNei
                        if (cmv(1).LT.mol%cueLayers(e)) then
                            currl1=e; exit
                        endif
                    enddo
                    do e = 1,maxNei
                        if (cmv(2).LT.mol%cueLayers(e)) then
                            currl2=e; exit
                        endif
                    enddo
                    do e = 1,maxNei
                        if (cmv(3).LT.mol%cueLayers(e)) then
                            currl3=e; exit
                        endif
                    enddo

                    if (inLayer(cmv(4),l)) then
                        if (.NOT.btest(i+a+j+b+k+c,0)) then
                            if (btest(i,0)) then; sh(1)='a'; else; sh(1)='b'; endif
                            if (btest(a,0)) then; sh(2)='a'; else; sh(2)='b'; endif
                            if (btest(j,0)) then; sh(3)='a'; else; sh(3)='b'; endif
                            if (btest(b,0)) then; sh(4)='a'; else; sh(4)='b'; endif
                            if (btest(k,0)) then; sh(5)='a'; else; sh(5)='b'; endif
                            if (btest(c,0)) then; sh(6)='a'; else; sh(6)='b'; endif

                            ui=i/2+mod(i,2); ua=(a-Nel)/2+mod(a-Nel,2)
                            uj=j/2+mod(j,2); ub=(b-Nel)/2+mod(b-Nel,2)
                            uk=k/2+mod(k,2); uc=(c-Nel)/2+mod(c-Nel,2)

                            value=t3(i,j,k,a,b,c)
                            if (abs(value).GT.ouThreshold) then
                                write (ou,116) sh(1),sh(2),sh(3),sh(4),sh(5),sh(6),ui,ua,uj,ub,uk,uc,value,currl1,currl2,currl3
                                count1=count1+1
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo
        write (ou,119) count1
    enddo

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=2 ! contribution to the correlation energy by layers
    if (chaps(cChap)) then

    energyContrib=0
    write (ou,103)
    do k = 1,maxNei

        cval=0
        do mm = 1,Nse
            i=spinInds(mm,1)
            a=spinInds(mm,2)

            if (pattern(1)) then
                if (inLayer(distMO(i,a),k)) then
                    cval(1)=cval(1)+F(i,a)*t1(i,a)
                endif
            endif
            do nn = 1,Nse
                j=spinInds(nn,1)
                b=spinInds(nn,2)

                holdDist(1)=distMO(i,j); holdDist(2)=distMO(i,a)
                holdDist(3)=distMO(i,b); holdDist(4)=distMO(j,a)
                holdDist(5)=distMO(j,b); holdDist(6)=distMO(a,b)

                if (pattern(2)) then
                    if (inLayer(maxval(holdDist(1:6)),k)) then
                        cval(2)=cval(2)+twoeintegr(a,i,b,j)*(t2(i,j,a,b)+t1(i,a)*t1(j,b)-t1(j,a)*t1(i,b))
                    endif
                endif
            enddo
        enddo

        cval(2)=0.25_rglu*cval(2)*multiplier

        if (abs(cval(1)).LT.1D-14) cval(1)=0
        if (abs(cval(2)).LT.1D-14) cval(2)=0
        write (ou,104) k,cval(1),cval(2)
        energyContrib(1)=energyContrib(1)+cval(1)
        energyContrib(2)=energyContrib(2)+cval(2)
    enddo
    write (ou,114) energyContrib

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=3 ! t1 structure analize by layers
    if (chaps(cChap).AND.pattern(1)) then

    fullt(1)=0; sumt1=0; maxt1=0; countt=0; mcountt1=0
    do mm = 1,Ne
        ui=spatInds(mm,1)
        ua=spatInds(mm,2)

        i=ui*2-1; a=ua*2-1

        cmv(1)=distMO(i,a)
        do k = 1,maxNei
            if (cmv(1).LT.mol%cueLayers(k)) then
                if (abs( t1(i,a) ).GT.abs( maxt1(k) )) maxt1(k)=t1(i,a)
                countt(k,1)=countt(k,1)+1
                sumt1(k)=sumt1(k)+t1(i,a)**2
                exit
            end if
        enddo
        fullt(1)=fullt(1)+t1(i,a)**2
        mcountt1=mcountt1+1
    enddo

    countt(:,1)=countt(:,1); mcountt1=mcountt1
    sumt1      =sumt1      ; fullt(1)=fullt(1)

    sss=0
    do k = 1,maxNei
        if (abs(maxt1(k)).LT.gluCompare) maxt1(k)=0
        sss(1)=sss(1)+100.*sumt1(k)   /fullt(1)
        sss(2)=sss(2)+100.*countt(k,1)/mcountt1
        inteContr(k)=sss(1); inteAmp(k)=sss(2)
    enddo

    do k = 1,maxNei
        efft1(k)=inteContr(k)/inteAmp(k)
    enddo

    write (ou,198)
    write (ou,200)
    do k = 1,maxNei
        write (ou,201) k,sumt1(k),100.*sumt1(k)/fullt(1),maxt1(k),countt(k,1),efft1(k)
    enddo
    write (ou,202)
    write (ou,203) fullt(1),100.,maxval(abs(maxt1)),mcountt1,1.

    write (ou,204) 't1',sqrt(fullt(1))/sqrt(real(N, rglu))

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=3 ! t2 structure analize by layers
    if (chaps(cChap).AND.pattern(2)) then

    fullt(2)=0; sumt2=0; maxt2=0; countt=0; mcountt2=0
    do mm = 1,Ne
        ui=spatInds(mm,1)
        ua=spatInds(mm,2)
        do nn = mm,Ne
            uj=spatInds(nn,1)
            ub=spatInds(nn,2)

            i=ui*2-1; a=ua*2-1
            j=uj*2  ; b=ub*2

            holdDist(1)=distMO(i,j); holdDist(2)=distMO(i,a)
            holdDist(3)=distMO(i,b); holdDist(4)=distMO(j,a)
            holdDist(5)=distMO(j,b); holdDist(6)=distMO(a,b)

            do k = 1,maxNei
                if (maxval(holdDist(1:6)).LT.mol%cueLayers(k)) then
                    if ( abs( t2(i,j,a,b) ).GT.abs( maxt2(k) )) maxt2(k)=t2(i,j,a,b)
                    countt(k,2)=countt(k,2)+1
                    sumt2(k)=sumt2(k)+t2(i,j,a,b)**2
                    exit
                end if
            enddo
            fullt(2)=fullt(2)+t2(i,j,a,b)**2
            mcountt2=mcountt2+1
        enddo
    enddo

    sss=0
    do k = 1,maxNei
        if (abs(maxt2(k)).LT.gluCompare) maxt2(k)=0
        sss(1)=sss(1)+100.*sumt2(k)   /fullt(2)
        sss(2)=sss(2)+100.*countt(k,2)/mcountt2
        inteContr(k)=sss(1); inteAmp(k)=sss(2)
    enddo

    do k = 1,maxNei
        efft2(k)=inteContr(k)/inteAmp(k)
    enddo

    write (ou,199)
    write (ou,200)
    do k = 1,maxNei
        write (ou,201) k,sumt2(k),100.*sumt2(k)/fullt(2),maxt2(k),countt(k,2),efft2(k)
    enddo
    write (ou,202)
    write (ou,203) fullt(2),100.,maxval(abs(maxt2)),mcountt2,1.

    write (ou,204) 't2',sqrt(fullt(2))/sqrt(real(N, rglu))

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=3 ! t3 structure analize by layers
    if (chaps(cChap).AND.pattern(3)) then

    fullt(3)=0; sumt3=0; maxt3=0; countt=0; mcountt3=0
    do mm = 1,Ne
        ui=spatInds(mm,1)
        ua=spatInds(mm,2)
        do nn = mm,Ne
            uj=spatInds(nn,1)
            ub=spatInds(nn,2)
            do pp = nn+1,Ne
                uk=spatInds(pp,1)
                uc=spatInds(pp,2)

                i=ui*2-1; a=ua*2-1

                jb(1,1)=uj*2; jb(1,2)=uj*2-1
                jb(2,1)=ub*2; jb(2,2)=ub*2-1

                kc(1,1)=uk*2; kc(1,2)=uk*2-1
                kc(2,1)=uc*2; kc(2,2)=uc*2-1

                do vv = 1,2
                    do gg = 1,2
                        j=jb(1,vv); b=jb(2,vv)
                        k=kc(1,gg); c=kc(2,gg)

                        holdDist( 1)=distMO(i,j); holdDist( 2)=distMO(i,a)
                        holdDist( 3)=distMO(i,b); holdDist( 4)=distMO(j,a)
                        holdDist( 5)=distMO(j,b); holdDist( 6)=distMO(a,b)
                        holdDist( 7)=distMO(i,k); holdDist( 8)=distMO(i,c)
                        holdDist( 9)=distMO(j,k); holdDist(10)=distMO(j,c)
                        holdDist(11)=distMO(a,k); holdDist(12)=distMO(a,c)
                        holdDist(13)=distMO(b,k); holdDist(14)=distMO(b,c)
                        holdDist(15)=distMO(k,c)

                        do l = 2,maxNei
                            if (maxval(holdDist(1:15)).LT.mol%cueLayers(l)) then
                                if ( abs( t3(i,j,k,a,b,c) ).GT.abs( maxt3(l) )) maxt3(l)=t3(i,j,k,a,b,c)
                                countt(l,3)=countt(l,3)+1
                                sumt3(l)=sumt3(l)+t3(i,j,k,a,b,c)**2
                                exit
                            end if
                        enddo
                        fullt(3)=fullt(3)+t3(i,j,k,a,b,c)**2
                        mcountt3=mcountt3+1
                    enddo
                enddo
            enddo
        enddo
    enddo

    sss=0
    do k = 2,maxNei
        if (abs(maxt3(k)).LT.gluCompare) maxt3(k)=0
        sss(1)=sss(1)+100.*sumt3(k)   /fullt(3)
        sss(2)=sss(2)+100.*countt(k,3)/mcountt3
        inteContr(k)=sss(1); inteAmp(k)=sss(2)
    enddo

    do k = 2,maxNei
        efft3(k)=inteContr(k)/inteAmp(k)
    enddo

    write (ou,399)
    write (ou,200)
    do k = 2,maxNei
        write (ou,201) k,sumt3(k),100.*sumt3(k)/fullt(3),maxt3(k),countt(k,3),efft3(k)
    enddo
    write (ou,202)
    write (ou,203) fullt(3),100.,maxval(abs(maxt3)),mcountt3,1.

    write (ou,204) 't3',sqrt(fullt(3))/sqrt(real(N, rglu))

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=4 ! C1, C2, C3, C4 estimation
    if (chaps(cChap)) then

    write (ou,'(//)')

    sumw=0; hold=0
    do dlayer = 1,maxNei

        if ((ou.NE.0).AND.(ou.EQ.6)) then
            backspace(ou)
        endif
        write (ou,'(A,1X,i2,A,i2,1X\)') 'Layer',dlayer,'/',maxNei

        sss=0; uj=5
        do mm = 1,Nse
            i=spinInds(mm,1)
            a=spinInds(mm,2)
            ui=int( 100.*float(mm)/float(Nse) )
            if (ui.EQ.uj) then
                uj=uj+5; write (ou,'(1X,i<mid(ui)>\)') ui
            endif

            if (distMO(i,a).GT.mol%cueLayers(dlayer)) cycle

            if (inLayer(distMO(i,a),dlayer)) then
                Ax=c1_t1([i,a])
                hold(1)=hold(1)+Ax*Ax
                cvalue=Ax; value=cvalue*cvalue
                sss(1)=sss(1)+value

                holdContributions(i    ,dlayer,1)=holdContributions(i    ,dlayer,1)+value
                holdContributions(a-Nel,dlayer,2)=holdContributions(a-Nel,dlayer,2)+value
            endif

            !$omp parallel default(shared) private(holdDist,nn,pp,vv,j,b,k,c,l,d,Bx,Cx,Dx,Ex,Fx,Gx,Hx,Ix,Jx,cvalue,value) reduction(+:sss,hold,holdContributions)
            !$omp do
            do nn = 1,Nse
                j=spinInds(nn,1)
                b=spinInds(nn,2)

                if ((j.LE.i).OR.(b.LE.a)) cycle

                holdDist(1)=distMO(i,a); holdDist(2)=distMO(i,j)
                holdDist(3)=distMO(i,b); holdDist(4)=distMO(a,j)
                holdDist(5)=distMO(a,b); holdDist(6)=distMO(j,b)

                if (maxval(holdDist(1:6)).GT.mol%cueLayers(dlayer)) cycle

                if (inLayer(maxval(holdDist(1:6)),dlayer)) then
                    Bx=c2_t2([i,j,a,b])
                    Cx=c2_t1t1([i,j,a,b])
                    hold(2)=hold(2)+Bx*Bx
                    hold(3)=hold(3)+Cx*Cx
                    cvalue=Bx+Cx; value=cvalue*cvalue
                    sss(2)=sss(2)+value

                    holdContributions(i    ,dlayer,1)=holdContributions(i    ,dlayer,1)+value
                    holdContributions(j    ,dlayer,1)=holdContributions(j    ,dlayer,1)+value
                    holdContributions(a-Nel,dlayer,2)=holdContributions(a-Nel,dlayer,2)+value
                    holdContributions(b-Nel,dlayer,2)=holdContributions(b-Nel,dlayer,2)+value
                endif

                do pp = 1,Nse
                    k=spinInds(pp,1)
                    c=spinInds(pp,2)

                    if ((k.LE.j).OR.(c.LE.b)) cycle
                    if (distMO(k,c).GT.mol%cueLayers(dlayer)) cycle

                    holdDist( 7)=distMO(i,k); holdDist( 8)=distMO(i,c)
                    holdDist( 9)=distMO(j,k); holdDist(10)=distMO(j,c)
                    holdDist(11)=distMO(a,k); holdDist(12)=distMO(a,c)
                    holdDist(13)=distMO(b,k); holdDist(14)=distMO(b,c)
                    holdDist(15)=distMO(k,c)

                    if (inLayer(maxval(holdDist(1:15)),dlayer)) then
                        Dx = c3_t1t2([i,j,k,a,b,c])
                        Ex = c3_t1t1t1([i,j,k,a,b,c])
                        Fx = c3_t3([i,j,k,a,b,c])

                        hold(4)=hold(4)+Dx*Dx
                        hold(5)=hold(5)+Ex*Ex
                        hold(9)=hold(9)+Fx*Fx
                        cvalue=Dx+Ex+Fx; value=cvalue*cvalue

                        sss(3)=sss(3)+value
                        holdContributions(i    ,dlayer,1)=holdContributions(i    ,dlayer,1)+value
                        holdContributions(j    ,dlayer,1)=holdContributions(j    ,dlayer,1)+value
                        holdContributions(k    ,dlayer,1)=holdContributions(k    ,dlayer,1)+value
                        holdContributions(a-Nel,dlayer,2)=holdContributions(a-Nel,dlayer,2)+value
                        holdContributions(b-Nel,dlayer,2)=holdContributions(b-Nel,dlayer,2)+value
                        holdContributions(c-Nel,dlayer,2)=holdContributions(c-Nel,dlayer,2)+value
                    endif

                    do vv = 1,Nse
                        l=spinInds(vv,1)
                        d=spinInds(vv,2)

                        if ((l.LE.k).OR.(d.LE.c)) cycle
                        if (distMO(l,d).GT.mol%cueLayers(dlayer)) cycle

                        holdDist(16)=distMO(i,l); holdDist(17)=distMO(i,d)
                        holdDist(18)=distMO(j,l); holdDist(19)=distMO(j,d)
                        holdDist(20)=distMO(k,l); holdDist(21)=distMO(k,d)
                        holdDist(22)=distMO(a,l); holdDist(23)=distMO(a,d)
                        holdDist(24)=distMO(b,l); holdDist(25)=distMO(b,d)
                        holdDist(26)=distMO(c,l); holdDist(27)=distMO(c,d)
                        holdDist(28)=distMO(l,d)

                        if (inLayer(maxval(holdDist(1:28)),dlayer)) then
                            Gx=c4_t2t2([i,j,k,l,a,b,c,d])
                            Hx=c4_t1t1t2([i,j,k,l,a,b,c,d])
                            Ix=c4_t1t1t1t1([i,j,k,l,a,b,c,d])
                            Jx=c4_t1t3([i,j,k,l,a,b,c,d])

                            hold(6) =hold(6) +Gx*Gx
                            hold(7) =hold(7) +Hx*Hx
                            hold(8) =hold(8) +Ix*Ix
                            hold(10)=hold(10)+Jx*Jx
                            cvalue=Gx+Hx+Ix+Jx; value=cvalue*cvalue

                            sss(4)=sss(4)+value
                            holdContributions(i    ,dlayer,1)=holdContributions(i    ,dlayer,1)+value
                            holdContributions(j    ,dlayer,1)=holdContributions(j    ,dlayer,1)+value
                            holdContributions(k    ,dlayer,1)=holdContributions(k    ,dlayer,1)+value
                            holdContributions(l    ,dlayer,1)=holdContributions(l    ,dlayer,1)+value
                            holdContributions(a-Nel,dlayer,2)=holdContributions(a-Nel,dlayer,2)+value
                            holdContributions(b-Nel,dlayer,2)=holdContributions(b-Nel,dlayer,2)+value
                            holdContributions(c-Nel,dlayer,2)=holdContributions(c-Nel,dlayer,2)+value
                            holdContributions(d-Nel,dlayer,2)=holdContributions(d-Nel,dlayer,2)+value
                        endif

                    enddo !vv
                enddo !pp
            enddo !nn
            !$omp end parallel
        enddo !mm

        holdValues(dlayer,1)=sss(1); holdValues(dlayer,2)=sss(2)
        holdValues(dlayer,3)=sss(3); holdValues(dlayer,4)=sss(4)
        holdValues(dlayer,5)=sum(sss)

        sumw(1)=sumw(1)+sss(1); sumw(2)=sumw(2)+sss(2)
        sumw(3)=sumw(3)+sss(3); sumw(4)=sumw(4)+sss(4)
        write (ou,*)
    enddo !dlayer
    write (ou,*)

    holdValues(maxNei+1,1)=sumw(1); holdValues(maxNei+1,2)=sumw(2)
    holdValues(maxNei+1,3)=sumw(3); holdValues(maxNei+1,4)=sumw(4)
    holdValues(maxNei+1,5)=sum(sumw)

    birav =(holdValues(maxNei+1,5)-holdValues(maxNei+1,1))/(1.+holdValues(maxNei+1,5))
    birav2= holdValues(maxNei+1,2)/(1+holdValues(maxNei+1,1)+holdValues(maxNei+1,2))

    sss(1)=holdValues(maxNei+1,5)
    holdValues=holdValues/sss(1)
    holdValues(maxNei+1,5)=sss(1)

    lav=0
    do dlayer = 1,maxNei
        lav=lav+dlayer*holdValues(dlayer,5)
    enddo
    write (ou,'(4X,"Average L:",2X,F7.4)') lav
    write (ou,'(4X,"Nodd (WF):",2X,F7.4)') birav
    write (ou,'(4X,"Nodd2(WF):",2X,F7.4)') birav2

    write (ou,300)
    do dlayer = 1,maxNei
        write (ou,301) dlayer,holdValues(dlayer,:)
    enddo
    write (ou,302) holdValues(maxNei+1,:)
    write (ou,303) (1-holdValues(1,5))/holdValues(1,5)

    divisor=1
    divisor=sum(hold)
    !divisor=sum(sumw)

    write (ou,304) 'Contribution by type:'
    if (pattern(1)) write (ou,305) '                  t1:',100.*hold(1) /divisor
    if (pattern(2)) write (ou,305) '                  t2:',100.*hold(2) /divisor
    if (pattern(1)) write (ou,305) '               t1*t1:',100.*hold(3) /divisor
    if (pattern(1)) write (ou,305) '               t1*t2:',100.*hold(4) /divisor
    if (pattern(1)) write (ou,305) '            t1*t1*t1:',100.*hold(5) /divisor
    if (pattern(3)) write (ou,305) '                  t3:',100.*hold(9) /divisor
    if (pattern(2)) write (ou,305) '               t2*t2:',100.*hold(6) /divisor
    if (pattern(1)) write (ou,305) '            t2*t1*t1:',100.*hold(7) /divisor
    if (pattern(1)) write (ou,305) '         t1*t1*t1*t1:',100.*hold(8) /divisor
    if (pattern(3)) write (ou,305) '               t1*t3:',100.*hold(10)/divisor
    write (ou,304) repeat('_',28)
    if (pattern(1)) write (ou,305) '             pure t1:',100.*(hold(1)+hold(3)+hold(5)+hold(8))/divisor
    if (pattern(2)) write (ou,305) '             pure t2:',100.*(hold(2)+hold(6))/divisor
    if (pattern(3)) write (ou,305) '             pure t3:',100.*(hold(9))/divisor
    write (ou,305) '               mixed:',100.*(hold(4)+hold(7)+hold(10))/divisor
    write (ou,304) repeat('_',28)
    write (ou,305) ' overall uncertainty:',100.*abs(sum(sumw)-sum(hold))                     /sum(sumw)
    write (ou,305) '   uncertainty in C1:',100.*abs(sumw(1)-hold(1))                         /sumw(1)
    write (ou,305) '   uncertainty in C2:',100.*abs(sumw(2)-hold(2)-hold(3))                 /sumw(2)
    write (ou,305) '   uncertainty in C3:',100.*abs(sumw(3)-hold(4)-hold(5)-hold(9))         /sumw(3)
    write (ou,305) '   uncertainty in C4:',100.*abs(sumw(4)-hold(6)-hold(7)-hold(8)-hold(10))/sumw(4)

    gsum=0
    do dLayer = 1,maxNei
        gsum=gsum+sum(holdContributions(:,dLayer,1))
        gsum=gsum+sum(holdContributions(:,dLayer,2))
    enddo

    do i = 2,Nel,2
        holdContributions(i/2,:,:)=holdContributions(i,:,:)+holdContributions(i-1,:,:)
    enddo

    holdContributions=100.*holdContributions/gsum

    call prMatrix(holdContributions(1:Nel/2,:,1),ou,'Donate','00.00',transpose=true,maxwidth=setfLineLen)
    call prMatrix(holdContributions(1:Nel/2,:,2),ou,'Accept','00.00',transpose=true,maxwidth=setfLineLen)
    call prMatrix(holdContributions(1:Nel/2,:,1)+holdContributions(:,:,2),ou,'Total','00.00',transpose=true,maxwidth=setfLineLen)

    pvec=0
    do i = 1,Nel/2
        pvec(i)=sum(holdContributions(i,:,1))+sum(holdContributions(i,:,2))
    enddo

    call prVector(pvec,ou,'Summated','00.00','horizontal',setfLineLen)

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=5 ! biradical character
    if (chaps(cChap)) then

    birad=0; rav=0
    do mm = 1,Nse
        i=spinInds(mm,1)
        a=spinInds(mm,2)

        do nn = 1,Nse
            j=spinInds(nn,1)
            b=spinInds(nn,2)

            if ((j.LE.i).OR.(b.LE.a)) cycle

            Bx=t2(i,j,a,b); Cx=t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a)
            cvalue=(Bx+Cx)*(Bx+Cx)

            rav=rav+abs(distMO(i,j)-distMO(a,b))*cvalue

            ui=i/2+mod(i,2); ua=(a-Nel)/2+mod(a-Nel,2)
            uj=j/2+mod(j,2); ub=(b-Nel)/2+mod(b-Nel,2)

            if ((ui.EQ.uj).AND.(ua.NE.ub)) then
                birad(1)=birad(1)+cvalue
            elseif ((ui.NE.uj).AND.(ua.EQ.ub)) then
                birad(2)=birad(2)+cvalue
            elseif ((ui.EQ.uj).AND.(ua.EQ.ub)) then
                birad(3)=birad(3)+cvalue
            elseif ((ui.NE.uj).AND.(ua.NE.ub)) then
                birad(4)=birad(4)+cvalue
            endif

        enddo
    enddo

    cvalue=sum(birad)

    rav=rav/cvalue
    write (ou,500) ( birad(i),100.*birad(i)/cvalue, i=1,4 ),cvalue,rav
500 format (//'C(iiab):  ',F8.4,1X,F6.2,'%'&
             /'C(ijaa):  ',F8.4,1X,F6.2,'%'&
             /'C(iiaa):  ',F8.4,1X,F6.2,'%'&
             /'C(ijab):  ',F8.4,1X,F6.2,'%'&
             /'sum:      ',F8.4            &
             /'Average R:',F8.4/)

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=6 ! statistics for t1
    if (chaps(cChap).AND.pattern(1)) then

    do k = 1,maxNei

        mean(k)=0; count1=0; minv(k)=1; maxv(k)=0
        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)

            i=ui*2-1; a=ua*2-1

            cmv(1)=distMO(i,a)
            if (inLayer(distMO(i,a),k)) then
                if (abs(t1(i,a)).GT.maxv(k)) maxv(k)=abs(t1(i,a))
                if (abs(t1(i,a)).LT.minv(k)) minv(k)=abs(t1(i,a))
                mean(k)=mean(k)+abs(t1(i,a))
                count1=count1+1
            end if
        enddo

        if (count1.EQ.0) count1=1

        mean(k)=mean(k)/dble(count1)
    enddo

    write (ou,'(/4X,A)'          ) ' Statistics for t1'
    write (ou,'(4X,A,4(3X,A,3X),1X,A)') ' k',' maxval','minval',' mean ','  sd  ','  number  '
    do k = 1,maxNei

        sd(k)=0; count1=0

        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)

            i=ui*2-1; a=ua*2-1
            cmv(1)=distMO(i,a)
            if (inLayer(distMO(i,a),k)) then
                sd(k)=sd(k)+(abs(t1(i,a))-mean(k))**2
                count1=count1+1
            end if
        enddo

        if (count1.EQ.0) count1=1

        sd(k)=sqrt(sd(k)/dble(count1))
        write (ou,'(4X,i2,4(1X,ES11.4),1X,i10)') k,maxv(k),minv(k),mean(k),sd(k),count1
    enddo

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=6 ! statistics for t2
    if (chaps(cChap).AND.pattern(2)) then

    do k = 1,maxNei

        mean(k)=0; minv(k)=1; maxv(k)=0; count1=0
        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)
            do nn = mm,Ne
                uj=spatInds(nn,1)
                ub=spatInds(nn,2)

                i=ui*2-1; a=ua*2-1
                j=uj*2  ; b=ub*2

                holdDist(1)=distMO(i,j)
                holdDist(2)=distMO(i,a)
                holdDist(3)=distMO(i,b)
                holdDist(4)=distMO(j,a)
                holdDist(5)=distMO(j,b)
                holdDist(6)=distMO(a,b)

                if (inLayer(maxval(holdDist(1:6)),k)) then
                    if (abs(t2(i,j,a,b)).GT.maxv(k)) maxv(k)=abs(t2(i,j,a,b))
                    if (abs(t2(i,j,a,b)).LT.minv(k)) minv(k)=abs(t2(i,j,a,b))
                    mean(k)=mean(k)+abs(t2(i,j,a,b))
                    count1=count1+1
                end if
            enddo
        enddo

        if (count1.EQ.0) count1=1

        mean(k)=mean(k)/dble(count1)
    enddo

    write (ou,'(/4X,A)'          ) ' Statistics for t2'
    write (ou,'(4X,A,4(3X,A,3X),1X,A)') ' k',' maxval','minval',' mean ','  sd  ','  number  '
    do k = 1,maxNei

        sd(k)=0; count1=0
        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)
            do nn = mm,Ne
                uj=spatInds(nn,1)
                ub=spatInds(nn,2)

                i=ui*2-1; a=ua*2-1
                j=uj*2  ; b=ub*2

                holdDist(1)=distMO(i,j)
                holdDist(2)=distMO(i,a)
                holdDist(3)=distMO(i,b)
                holdDist(4)=distMO(j,a)
                holdDist(5)=distMO(j,b)
                holdDist(6)=distMO(a,b)

                if (inLayer(maxval(holdDist(1:6)),k)) then
                    sd(k)=sd(k)+(abs(t2(i,j,a,b))-mean(k))**2
                    count1=count1+1
                end if
            enddo
        enddo

        if (count1.EQ.0) count1=1

        sd(k)=sqrt(sd(k)/dble(count1))
        write (ou,'(4X,i2,4(1X,ES11.4),1X,i10)') k,maxv(k),minv(k),mean(k),sd(k),count1
    enddo

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=6 ! statistics for t3
    if (chaps(cChap).AND.pattern(3)) then

    do l = 2,maxNei

        mean(l)=0; minv(l)=1; maxv(l)=0; count1=0
        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)
            do nn = mm,Ne
                uj=spatInds(nn,1)
                ub=spatInds(nn,2)
                do pp = nn+1,Ne
                    uk=spatInds(pp,1)
                    uc=spatInds(pp,2)

                    i=ui*2-1; a=ua*2-1

                    jb(1,1)=uj*2; jb(1,2)=uj*2-1
                    jb(2,1)=ub*2; jb(2,2)=ub*2-1

                    kc(1,1)=uk*2; kc(1,2)=uk*2-1
                    kc(2,1)=uc*2; kc(2,2)=uc*2-1

                    do vv = 1,2
                        do gg = 1,2
                            j=jb(1,vv); b=jb(2,vv)
                            k=kc(1,gg); c=kc(2,gg)

                            holdDist( 1)=distMO(i,j); holdDist( 2)=distMO(i,a)
                            holdDist( 3)=distMO(i,b); holdDist( 4)=distMO(j,a)
                            holdDist( 5)=distMO(j,b); holdDist( 6)=distMO(a,b)
                            holdDist( 7)=distMO(i,k); holdDist( 8)=distMO(i,c)
                            holdDist( 9)=distMO(j,k); holdDist(10)=distMO(j,c)
                            holdDist(11)=distMO(a,k); holdDist(12)=distMO(a,c)
                            holdDist(13)=distMO(b,k); holdDist(14)=distMO(b,c)
                            holdDist(15)=distMO(k,c)

                            if (inLayer(maxval(holdDist(1:15)),l)) then
                                if (abs(t3(i,j,k,a,b,c)).GT.maxv(l)) maxv(l)=abs(t3(i,j,k,a,b,c))
                                if (abs(t3(i,j,k,a,b,c)).LT.minv(l)) minv(l)=abs(t3(i,j,k,a,b,c))
                                mean(l)=mean(l)+abs(t3(i,j,k,a,b,c))
                                count1=count1+1
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

        if (count1.EQ.0) count1=1

        mean(l)=mean(l)/real(count1, rglu)
    enddo

    write (ou,'(/4X,A)'          ) ' Statistics for t3'
    write (ou,'(4X,A,4(3X,A,3X),1X,A)') ' k',' maxval','minval',' mean ','  sd  ','  number  '
    do l = 2,maxNei

        sd(l)=0; count1=0
        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)
            do nn = mm,Ne
                uj=spatInds(nn,1)
                ub=spatInds(nn,2)
                do pp = nn+1,Ne
                    uk=spatInds(pp,1)
                    uc=spatInds(pp,2)

                    i=ui*2-1; a=ua*2-1
                    jb(1,1)=uj*2; jb(1,2)=uj*2-1
                    jb(2,1)=ub*2; jb(2,2)=ub*2-1

                    kc(1,1)=uk*2; kc(1,2)=uk*2-1
                    kc(2,1)=uc*2; kc(2,2)=uc*2-1

                    do vv = 1,2
                        do gg = 1,2
                            j=jb(1,vv); b=jb(2,vv)
                            k=kc(1,gg); c=kc(2,gg)

                            holdDist( 1)=distMO(i,j); holdDist( 2)=distMO(i,a)
                            holdDist( 3)=distMO(i,b); holdDist( 4)=distMO(j,a)
                            holdDist( 5)=distMO(j,b); holdDist( 6)=distMO(a,b)
                            holdDist( 7)=distMO(i,k); holdDist( 8)=distMO(i,c)
                            holdDist( 9)=distMO(j,k); holdDist(10)=distMO(j,c)
                            holdDist(11)=distMO(a,k); holdDist(12)=distMO(a,c)
                            holdDist(13)=distMO(b,k); holdDist(14)=distMO(b,c)
                            holdDist(15)=distMO(k,c)

                            if (inLayer(maxval(holdDist(1:15)),l)) then
                                sd(l)=sd(l)+(abs(t3(i,j,k,a,b,c))-mean(l))**2
                                count1=count1+1
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

        if (count1.EQ.0) count1=1

        sd(l)=sqrt(sd(l)/real(count1, rglu))
        write (ou,'(4X,i2,4(1X,ES11.4),1X,i10)') l,maxv(l),minv(l),mean(l),sd(l),count1
    enddo

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    cChap=7 ! configurational compound
    if (chaps(cChap)) then

    ilen=mid(Nel/2)
    zlen=7
    t1plen=15+2*ilen; t2plen=19+4*ilen

    write (ou,497)
    if (dCUE) then
        write (ou,498)
    else
        write (ou,499)
    endif

    write (ou,400); ulen=zlen

    do i = 1,Nel,2
    do a = 1,Nel,2
        if (abs(t1(i,a+Nel)).GT.prBar) then
            if (dCUE) then
                ui=(i+1)/2
                ua=(a+1)/2
            else
                ui=Nel/2-(i+1)/2+1
                ua=(a+1)/2
            endif

            if (ulen+t1plen.GT.setfLineLen) then
                write (ou,401) t1(i,a+Nel),ui,ua; ulen=zlen+t1plen
            else
                write (ou,402) t1(i,a+Nel),ui,ua; ulen=ulen+t1plen
            endif
        endif
    enddo
    enddo

    write (ou,*); ulen=setfLineLen

    do i = 1,Nel,2
    do a = 1,Nel,2
    do j = i,Nel
    do b = a,Nel

        if (abs(t2(i,j,a+Nel,b+Nel)).GT.prBar) then
            if (dCUE) then
                ui=(i+1)/2
                ua=(a+1)/2
                uj=(j+1)/2
                ub=(b+1)/2
            else
                ui=Nel/2-(i+1)/2+1
                ua=(a+1)/2
                uj=Nel/2-(j+1)/2+1
                ub=(b+1)/2
            endif

            if (ulen+t2plen.GT.setfLineLen) then
                write (ou,411) t2(i,j,a+Nel,b+Nel),ui,ua,uj,ub; ulen=zlen+t2plen
            else
                write (ou,412) t2(i,j,a+Nel,b+Nel),ui,ua,uj,ub; ulen=ulen+t2plen
            endif
        endif

    enddo
    enddo
    enddo
    enddo

    write (ou,*)

    endif !cChap

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    void=finishAnalize()

    if (tdlist.NE.0) then
        if (dprHeader) write (ou,'(/A/)') tpAdjustc('End of analysis',setfLineLen,'=')
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

100 format (/4X,'Number of unique one-electron excitations:',1X,1i12)
190 format (/4X,'Number of unique one-electron excitations:',1X,1i12,1X,'Print threshold:',ES8.1)


105 format ( 2X,A1,2(1X,i4),2X,ES16.9)
108 format ( 31('_')/22X,i9)

107 format (/4X,'Excitations in l=',1X,i3,1X,'ethylenes:'/)

101 format (/4X,'Number of unique two-electron excitations:',1X,1i12/&
             4X,'Only non-zero contribution to correlation energy is shown.'//&
             4X,'In case of CUE, indexes correspond to ethylenes.'//&
             7X,'i',4X,'a',4X,'j',4X,'b',9X,'t2',14X,'integral',7X,'contribution',5X,'ia',2X,'jb',1X,'exc')

191 format (/4X,'Number of unique two-electron excitations:',1X,1i12,1X,'Print threshold:',ES8.1/&
             7X,'i',4X,'a',4X,'j',4X,'b',9X,'t2',14X,'integral',7X,'contribution')

106 format (4A1,4(i4,1X),1X,ES16.9,2X,ES16.9,2X,ES16.9,1X,3(1X,i3))
196 format (4A1,4(i4,1X),1X,ES16.9,2X,ES16.9,2X,ES16.9)

109 format ( 90('_')/81X,i9)

117 format (/4X,'Excitations in l=',1X,i3,1X,'ethylenes:'/)
111 format (/4X,'Number of unique three-electron excitations:',1X,1i14/&
             4X,'Print threshold',1X,ES9.1//&
             9X,'i',4X,'a',4X,'j',4X,'b',4X,'k',4X,'c',10X,'t3',9X,'ia',2X,'jb',2X,'kc')
116 format ( 6A1,6(i4,1X),1X,ES16.9,1X,3(1X,i3))
119 format ( 66('_')/57X,i9)

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
103 format (/4X,'Contribution to correlation energy from C1 and C2 by layers:'/)
104 format ( 8X,i<mid(maxNei)>,2X,ES16.9,2X,ES16.9)
114 format ( 8X,<mid(maxNei)+2*18>('_')/&
             8X,<mid(maxNei)>X,2X,ES16.9,2X,ES16.9)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
198 format (/4X,13X,'-------------t1 NormAnalize-------------')
199 format (/4X,13X,'-------------t2 NormAnalize-------------')
399 format (/4X,13X,'-------------t3 NormAnalize-------------')
200 format ( 2X,14X,'|',4X,'||T||',4X,'|',1X,'Contrib. %',1X,'|',3X,'max(T)',3X,'|',3X,'N(T)',3X,'|',2X,'Efficacy')
201 format ( 3X,'||cue(',i<mid(maxNei)>,')||',<4-mid(maxNei)>X,'|', 2X,F10.5,1X, '|',   1X,F10.5,1X,    '|', 1X,ES10.3,1X, '|',  1X,i8,1X,  '|',1X,F9.3)
202 format ( 2X,76('_'))
203 format ( 2X,14X,'|', 2X,F10.5,1X, '|',   1X,F10.5,1X,    '|', 1X,ES10.3,1X, '|',  1X,i8,1X,  '|',1X,F9.3)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
204 format (/4X,A,'-diagnostics:',1X,F9.6)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
300 format (/2X,'All values (except the last) are normalized on {W[1]+W[2]+W[3]+W[4]}.'/&
             2X,'L',5X,'Vl[1]',7X,'Vl[2]',7X,'Vl[3]',7X,'Vl[4]',9X,'Vl')
301 format ( 1X,i2,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7)
302 format ( 2X,61('-')/&
             2X,'W',1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7/)

303 format (/4X,'Theta:',1X,F8.3/)
304 format ( 4X,A)
305 format ( 4X,A,1X,F6.3)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
497 format (//4X,'Configurational compound of wave-function')
498 format (  4X,"Orbitals correspond to ethylenes (i): i = occupied, i' = vacant.")
499 format (  4X,"Occupied orbitals: 1=HOMO, 2=HOMO-1, etc. Vacant orbitals: 1'=LUMO, 2'=LUMO+1, etc.")

400 format ('|PSI> ='\)

401 format (/<zlen>X,1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"'>"\)
402 format (         1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"'>"\)

411 format (/<zlen>X,1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"',",i<ilen>,'->',i<ilen>,"'>"\)
412 format (         1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"',",i<ilen>,'->',i<ilen>,"'>"\)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    return
    end subroutine wfAnalize