
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    logical(kind=lglu) function inLayer(dist,layer) result(ret)
    use glob, only: rglu,iglu,lglu
    use hdb,  only: mol
    implicit none

    real(kind=rglu)    :: dist
    integer(kind=iglu) :: layer


    ret=(dist.GT.mol%cueLayers(layer-1)).AND.(dist.LE.mol%cueLayers(layer))
    return
    end function inLayer

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ AMPLITUDES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine print_t1_by_layer_nonlocal

    use glob,     only: rglu,iglu,lglu,true,false
    use hdb,      only: ccbd,ou
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    character (len=1)  :: sh
    integer(kind=iglu) :: i, a, ui, ua, count1


    write (ou,100) prStrByVal(Nel*(No-Nel)), ccbd%printThreshold
    count1=0
    do i = 1,Nel
    do a = Nel+1,No
        if (.NOT.btest(i+a,0)) then
            if (btest(i,0)) then; sh='a'; else; sh='b'; endif
            ui=i/2+mod(i,2); ua=a/2+mod(a,2)

            if (abs(t1(i,a)).GE.ccbd%printThreshold) then
                write (ou,101) sh,ui,ua,t1(i,a),f(i,a),t1(i,a)*f(i,a)
                count1=count1+1
            endif
        endif
    enddo
    enddo
    write (ou,102) count1

100 format (/2X,'Number of unique 1-electron excitations (',A,') [output |t1| >',ES8.1,']:'//&
             2X,'Shell',2X,'i',4X,'a',8X,'t(i,a)',11X,'f(i,a)',9X,'Contribution'/&
             2X, 67('-'))
101 format ( 4X,A1,2(1X,i4),3(2X,ES16.9))
102 format ( 4X,65('_')/60X,i9)

    return
    end subroutine print_t1_by_layer_nonlocal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine print_t1_by_layer_local

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd, mol, ou
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    character (len=1)  :: sh
    integer(kind=iglu) :: i, a, ui, ua, k, count


    write (ou,100) prStrByVal(Nel*(No-Nel)), ccbd%printThreshold
    count=0
    do k = 1,maxNei
        do i = 1,Nel
        do a = Nel+1,No
            if (inLayer(distMO(i,a),k)) then
                if (.NOT.btest(i+a,0)) then
                    if (btest(i,0)) then; sh='a'; else; sh='b'; endif
                    ui=i/2+mod(i,2); ua=a/2+mod(a,2)

                    if (abs(t1(i,a)).GE.ccbd%printThreshold) then
                        write (ou,101) sh,ui,ua,t1(i,a),f(i,a),t1(i,a)*f(i,a),k
                        count=count+1
                    endif
                endif
            endif
        enddo
        enddo
    enddo
    write (ou,102) count

100 format (/2X,'Number of unique 1-electron excitations (',A,') [output |t1| >',ES8.1,']:'//&
             2X,'Shell',2X,'i',4X,'a',8X,'t(i,a)',11X,'f(i,a)',9X,'Contribution',5X,'Layer'/&
             2X,75('-'))
101 format ( 4X,A1,2(1X,i4),3(2X,ES16.9),4X,i2)
102 format ( 4X,71('_')/66X,i9)

    return
    end subroutine print_t1_by_layer_local

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine print_t2_by_layer_nonlocal

    use glob,     only: rglu,iglu,lglu,true,false
    use hdb,      only: ccbd,ou
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu)             :: count,nn,mm,i,a,j,b,ui,uj,ua,ub
    real(kind=rglu)                :: val1,val2
    character (len=1)              :: sh(4)


    write (ou,100) prStrByVal((Nse*(Nse-1))/2), ccbd%printThreshold
    count = 0
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

                ui=i/2+mod(i,2); ua=a/2+mod(a,2)
                uj=j/2+mod(j,2); ub=b/2+mod(b,2)

                val1=t2(i,j,a,b)
                val2=twoeintegr(i,j,a,b)
                if (abs(val1).GT.ccbd%printThreshold) then
                    write (ou,101) sh,ui,ua,uj,ub,val1,val2,val1*val2
                    count = count + 1
                endif
            endif
        enddo
    enddo
    write (ou,102) count

100 format (/2X,'Number of unique 2-electron excitations (',A,') [output |t2| >',ES8.1,']:'//&
             2X,'Shell',4X,'i',4X,'a',4X,'j',4X,'b',5X,'t(i,a,j,b)',9X,'[ia||jb]',8X,'Contribution'/&
             2X,79('-'))
101 format ( 3X,4A1,4(1X,i4),3(2X,ES16.9))
102 format ( 3X,78('_')/72X,i9)

    return
    end subroutine print_t2_by_layer_nonlocal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine print_t2_by_layer_local

    use glob,     only: rglu,iglu,lglu,true,false
    use hdb,      only: ccbd,ou,mol
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu)             :: count,nn,mm,i,a,j,b,ui,uj,ua,ub,k,l,currl1,currl2,d
    real(kind=rglu)                :: val1,val2,holdDist(6),cmv(3),cor(6),dist12
    character (len=1)              :: sh(4)


    write (ou,100) prStrByVal((Nse*(Nse-1))/2), ccbd%printThreshold
    count=0
    do l = 1,maxNei
        do mm = 1,Nse
            i=spinInds(mm,1)
            a=spinInds(mm,2)
            do nn = 1,Nse
                j=spinInds(nn,1)
                b=spinInds(nn,2)

                holdDist(1)=distMO(i,j); holdDist(2)=distMO(i,a)
                holdDist(3)=distMO(i,b); holdDist(4)=distMO(j,a)
                holdDist(5)=distMO(j,b); holdDist(6)=distMO(a,b)
                cmv(1)=distMO(i,a); cmv(2)=distMO(j,b); cmv(3)=maxval(holdDist)

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

                dist12=sqrt((cor(1)-cor(4))**2+(cor(2)-cor(5))**2+(cor(3)-cor(6))**2)

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

                        ui=i/2+mod(i,2); ua=a/2+mod(a,2)
                        uj=j/2+mod(j,2); ub=b/2+mod(b,2)

                        val1=t2(i,j,a,b)
                        val2=twoeintegr(i,j,a,b)
                        if (abs(val1).GT.ccbd%printThreshold) then
                            write (ou,101) sh(1:4),ui,ua,uj,ub,val1,val2,val1*val2,currl1,currl2,d
                            count=count+1
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
    write (ou,102) count

100 format (/2X,'Number of unique 2-electron excitations (',A,') [output |t2| >',ES8.1,']:'//&
             2X,'Shell',4X,'i',4X,'a',4X,'j',4X,'b',5X,'t(i,a,j,b)',9X,'[ia||jb]',8X,'Contribution',&
             3X,'ia',1X,'jb',1X,'exc'/&
             2X,89('-'))
101 format ( 3X,4A1,4(1X,i4),3(2X,ES16.9),3(1X,i2))
102 format ( 3X,87('_')/81X,i9)

    return
    end subroutine print_t2_by_layer_local

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine print_t3_by_layer_nonlocal

    use glob,     only: rglu,iglu,lglu,true,false
    use hdb,      only: ccbd,ou
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu)             :: count,nn,mm,pp,i,a,j,b,k,c,ui,uj,uk,ua,ub,uc
    real(kind=rglu)                :: value
    character (len=1)              :: sh(6)


    write (ou,100) prStrByVal((Nse*(Nse-1)*(Nse-2))/2),ccbd%printThreshold
    count=0
    do mm = 1,Nse
        i=spinInds(mm,1)
        a=spinInds(mm,2)
        do nn = 1,Nse
            j=spinInds(nn,1)
            b=spinInds(nn,2)
            do pp = 1,Nse
                k=spinInds(pp,1)
                c=spinInds(pp,2)

                if (.NOT.btest(i+a+j+b+k+c,0)) then
                    if (btest(i,0)) then; sh(1)='a'; else; sh(1)='b'; endif
                    if (btest(a,0)) then; sh(2)='a'; else; sh(2)='b'; endif
                    if (btest(j,0)) then; sh(3)='a'; else; sh(3)='b'; endif
                    if (btest(b,0)) then; sh(4)='a'; else; sh(4)='b'; endif
                    if (btest(k,0)) then; sh(5)='a'; else; sh(5)='b'; endif
                    if (btest(c,0)) then; sh(6)='a'; else; sh(6)='b'; endif

                    ui=i/2+mod(i,2); ua=a/2+mod(a,2)
                    uj=j/2+mod(j,2); ub=b/2+mod(b,2)
                    uk=k/2+mod(k,2); uc=c/2+mod(c,2)

                    value=t3(i,j,k,a,b,c)
                    if (abs(value).GT.ccbd%printThreshold) then
                        write (ou,101) sh(1:6),ui,ua,uj,ub,uk,uc,value
                        count=count+1
                    endif
                endif
            enddo
        enddo
    enddo
    write (ou,102) count

100 format (/2X,'Number of unique 3-electron excitations (',A,') [output |t3| >',ES8.1,']:'//&
             3X,'Shell',5X,'i',4X,'a',4X,'j',4X,'b',4X,'k',4X,'c',3X,'t(i,a,j,b,k,c)'/&
             3X,54('-'))
101 format ( 3X,6A1,6(1X,i4),2X,ES16.9)
102 format ( 3X,54('_')/48X,i9)

    return
    end subroutine print_t3_by_layer_nonlocal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine print_t3_by_layer_local

    use glob,     only: rglu,iglu,lglu,true,false
    use hdb,      only: ccbd,ou,mol
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu)             :: count,nn,mm,pp,i,a,j,b,k,c,e,ui,uj,uk,ua,ub,uc,currl1,currl2,currl3,l,d
    real(kind=rglu)                :: value,holdDist(15),cmv(4)
    character (len=1)              :: sh(6)


    write (ou,100) prStrByVal((Nse*(Nse-1)*(Nse-2))/2),ccbd%printThreshold
    count=0
    do l = 2,maxNei
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

                    cmv(1)=distMO(i,a)
                    cmv(2)=distMO(j,b)
                    cmv(3)=distMO(k,c)
                    cmv(4)=maxval(holdDist)

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
                            if (abs(value).GT.ccbd%printThreshold) then
                                write (ou,101) sh(1:6),ui,ua,uj,ub,uk,uc,value,currl1,currl2,currl3
                                count=count+1
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo
    write (ou,102) count


100 format (/2X,'Number of unique 3-electron excitations (',A,') [output |t3| >',ES8.1,']:'//&
             3X,'Shell',5X,'i',4X,'a',4X,'j',4X,'b',4X,'k',4X,'c',3X,'t(i,a,j,b,k,c)',3X,'ia',2X,'jb',2X,'kc'/&
             3X,66('-'))
101 format ( 3X,6A1,6(1X,i4),2X,ES16.9,3(2X,i2))
102 format ( 3X,66('_')/60X,i9)

    return
    end subroutine print_t3_by_layer_local

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CORRELATION ENERGY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine correlation_energy_contribution_local

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu)           :: mm,nn,k,i,j,a,b
    real(kind=rglu)              :: holdDist(6),overall,pval1,pval2,sum1,sum2
    real(kind=rglu), allocatable :: contrib(:,:)


    allocate(contrib(maxNei, 2))

    write (ou,100)

    contrib=0
    do k = 1,maxNei

        do i = 1,Nel
        do a = Nel+1,No
            if (inLayer(distMO(i,a),k)) then
                contrib(k,1)=contrib(k,1)+f(i,a)*t1(i,a)
            endif
        enddo
        enddo

        do i = 1,Nel-1
        do a = Nel+1,No-1
            do j = i+1,Nel
            do b = a+1,No
                holdDist(1)=distMO(i,j); holdDist(2)=distMO(i,a)
                holdDist(3)=distMO(i,b); holdDist(4)=distMO(j,a)
                holdDist(5)=distMO(j,b); holdDist(6)=distMO(a,b)

                if (inLayer(maxval(holdDist),k)) then
                    contrib(k,2)=contrib(k,2)+twoeintegr(i,a,j,b)*&
                                 (t2(i,j,a,b)+t1(i,a)*t1(j,b)-t1(j,a)*t1(i,b))
                endif

            enddo
            enddo
        enddo
        enddo

    enddo

    sum1=0; sum2=0
    do k = 1,maxNei
        sum1=sum1+contrib(k,1)
        sum2=sum2+contrib(k,2)
    enddo
    overall=sum1+sum2

    do k = 1,maxNei
        pval1=0
        if (abs(contrib(k,1)).GT.gluCompare) then
            pval1=contrib(k,1)
        endif

        pval2=0
        if (abs(contrib(k,2)).GT.gluCompare) then
            pval2=contrib(k,2)
        endif
        write(ou,101) k,pval1,100*pval1/overall,pval2,100*pval2/overall,pval1+pval2,100*(pval1+pval2)/overall
    enddo
    write(ou,102) sum1,100*sum1/overall,sum2,100*sum2/overall,overall

    deallocate(contrib)

100 format (/4X,'Contribution to CC correction from C1 and C2 by layers:'//&
             2X,'Layer',11X,'C1[L]',22X,'C2[L]',18X,'C1[L]+C2[L]')
101 format (5X,i2,3(1X,ES16.9,1X,'(',F6.2,'%)'))
102 format (5X,83('-')/&
            7X,   2(1X,ES16.9,1X,'(',F6.2,'%)'),1X,ES16.9)

    return
    end subroutine correlation_energy_contribution_local

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine correlation_energy_contribution_nonlocal

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu)           :: k,i,j,a,b
    real(kind=rglu)              :: overall,pval1,pval2,sum1,sum2
    real(kind=rglu), allocatable :: contrib(:,:)


    allocate(contrib(0:No/4-1, 2))

    write (ou,100)

    contrib=0
    do k = 0,No/4-1
        do i = Nel-2*k-1, Nel
        do a = Nel+1, Nel+2*k+2
            contrib(k,1)=contrib(k,1)+f(i,a)*t1(i,a)
        enddo
        enddo

        do i = Nel-2*k-1, Nel-1
        do a = Nel+1, Nel+2*k+2-1
            do j = i+1,Nel
            do b = a+1,Nel+2*k+2
                contrib(k,2)=contrib(k,2)+twoeintegr(i,a,j,b)*(t2(i,j,a,b)+t1(i,a)*t1(j,b)-&
                                                                           t1(j,a)*t1(i,b))
            enddo
            enddo
        enddo
        enddo
    enddo

    sum1=contrib(No/4-1,1); sum2=contrib(No/4-1,2)
    overall=sum1+sum2
    do k = No/4-1,1,-1
        contrib(k,1)=contrib(k,1)-contrib(k-1,1)
        contrib(k,2)=contrib(k,2)-contrib(k-1,2)
    enddo

    do k = 0,No/4-1
        pval1=0
        if (abs(contrib(k,1)).GT.1D-10) then
            pval1=contrib(k,1)
        endif

        pval2=0
        if (abs(contrib(k,2)).GT.1D-10) then
            pval2=contrib(k,2)
        endif
        write(ou,101) k,pval1,100*pval1/overall,pval2,100*pval2/overall,pval1+pval2,100*(pval1+pval2)/overall
    enddo
    if (abs(sum1).LT.1D-10) sum1=0
    if (abs(sum2).LT.1D-10) sum2=0

    write(ou,102) sum1,100*sum1/overall,sum2,100*sum2/overall,overall

    deallocate(contrib)

100 format (/11X,'Contribution to CC correction from C1 and C2. Differential contribution'/&
             10X,'from every next pair of MOs accounted (HOMO-x ... HOMO, LUMO ... LUMO+x):'//&
              6X,'x',13X,'C1[x]',22X,'C2[x]',18X,'C1[x]+C2[x]')
101 format (  5X,i2,3(1X,ES16.9,1X,'(',F6.2,'%)'))
102 format (  5X,83('-')/,7X,2(1X,ES16.9,1X,'(',F6.2,'%)'),1X,ES16.9)

    return
    end subroutine correlation_energy_contribution_nonlocal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine t1_structure_analysis_by_layers

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare,mid
    use hdb,      only: ccbd,ou,mol
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu) :: mm,k,i,a,ui,ua,mcountt1
    real(kind=rglu)    :: cmv,fullt,sss(2)


    fullt=0; sumt1=0; maxt1=0; countt=0; mcountt1=0
    do mm = 1,Ne
        ui=spatInds(mm,1)
        ua=spatInds(mm,2)

        i=ui*2-1; a=ua*2-1

        cmv=distMO(i,a)
        do k = 1,maxNei
            if (cmv.LT.mol%cueLayers(k)) then
                if (abs(t1(i,a)).GT.abs(maxt1(k))) maxt1(k)=t1(i,a)
                countt(k,1)=countt(k,1)+1
                sumt1(k)=sumt1(k)+t1(i,a)**2
                exit
            end if
        enddo
        fullt=fullt+t1(i,a)**2
        mcountt1=mcountt1+1
    enddo

    sss=0
    do k = 1,maxNei
        if (abs(maxt1(k)).LT.gluCompare) maxt1(k)=0
        sss(1)=sss(1)+100*sumt1(k)   /fullt
        sss(2)=sss(2)+100*countt(k,1)/mcountt1
        inteContr(k)=sss(1); inteAmp(k)=sss(2)
    enddo

    do k = 1,maxNei
        efft1(k)=inteContr(k)/inteAmp(k)
    enddo

    write (ou,100)
    write (ou,101)
    do k = 1,maxNei
        write (ou,102) k,sumt1(k),100.*sumt1(k)/fullt,maxt1(k),countt(k,1),efft1(k)
    enddo
    write (ou,103)
    write (ou,104) fullt,100.,maxval(abs(maxt1)),mcountt1,1

    write (ou,105) 't1',sqrt(fullt/N)

100 format (/20X,'------------- T1 Norm Analysis -------------'/)
101 format (  2X,14X,'|',4X,'||T||',4X,'|',1X,'Contrib. %',1X,'|',3X,'max(T)',3X,'|',3X,'N(T)',3X,'|',2X,'Efficacy')
102 format (  3X,'||cue(',i<mid(maxNei)>,')||',<4-mid(maxNei)>X,'|',2X,F10.5,1X,'|',1X,F10.5,1X,'|',1X,ES10.3,1X,'|',1X,i8,1X,'|',1X,F9.3)
103 format (  2X,76('_'))
104 format (  2X,14X,'|', 2X,F10.5,1X, '|',1X,F10.5,1X,'|', 1X,ES10.3,1X,'|',1X,i8,1X,'|',1X,F9.3)
105 format (/ 4X,A,'-diagnostics:',1X,F9.6)

    return
    end subroutine t1_structure_analysis_by_layers

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine t2_structure_analysis_by_layers

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare,mid
    use hdb,      only: ccbd,ou,mol
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu) :: mm,nn,k,i,a,j,b,ui,ua,uj,ub,mcountt2
    real(kind=rglu)    :: cmv,fullt,sss(2),holdDist(6)


    fullt=0; sumt2=0; maxt2=0; countt=0; mcountt2=0
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
                    if (abs(t2(i,j,a,b)).GT.abs(maxt2(k))) maxt2(k)=t2(i,j,a,b)
                    countt(k,2)=countt(k,2)+1
                    sumt2(k)=sumt2(k)+t2(i,j,a,b)**2
                    exit
                end if
            enddo
            fullt=fullt+t2(i,j,a,b)**2
            mcountt2=mcountt2+1
        enddo
    enddo

    sss=0
    do k = 1,maxNei
        if (abs(maxt2(k)).LT.gluCompare) maxt2(k)=0
        sss(1)=sss(1)+100.*sumt2(k)   /fullt
        sss(2)=sss(2)+100.*countt(k,2)/mcountt2
        inteContr(k)=sss(1); inteAmp(k)=sss(2)
    enddo

    do k = 1,maxNei
        efft2(k)=inteContr(k)/inteAmp(k)
    enddo

    write (ou,100)
    write (ou,101)
    do k = 1,maxNei
        write (ou,102) k,sumt2(k),100.*sumt2(k)/fullt,maxt2(k),countt(k,2),efft2(k)
    enddo
    write (ou,103)
    write (ou,104) fullt,100,maxval(abs(maxt2)),mcountt2,1

    write (ou,105) 't2',sqrt(fullt/N)

100 format (/20X,'------------- T2 Norm Analysis -------------'/)
101 format (  2X,14X,'|',4X,'||T||',4X,'|',1X,'Contrib. %',1X,'|',3X,'max(T)',3X,'|',3X,'N(T)',3X,'|',2X,'Efficacy')
102 format (  3X,'||cue(',i<mid(maxNei)>,')||',<4-mid(maxNei)>X,'|',2X,F10.5,1X,'|',1X,F10.5,1X,'|',1X,ES10.3,1X,'|',1X,i8,1X,'|',1X,F9.3)
103 format (  2X,76('_'))
104 format (  2X,14X,'|', 2X,F10.5,1X, '|',1X,F10.5,1X,'|', 1X,ES10.3,1X,'|',1X,i8,1X,'|',1X,F9.3)
105 format (/ 4X,A,'-diagnostics:',1X,F9.6)

    return
    end subroutine t2_structure_analysis_by_layers

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine t3_structure_analysis_by_layers

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare,mid
    use hdb,      only: ccbd,ou,mol
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu) :: mm,nn,pp,vv,gg,i,a,j,b,k,c,l,ui,ua,uj,ub,uk,uc,jb(2,2),kc(2,2),mcountt3
    real(kind=rglu)    :: cmv,fullt,sss(2),holdDist(15)


    fullt=0; sumt3=0; maxt3=0; countt=0; mcountt3=0
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
                                if (abs(t3(i,j,k,a,b,c)).GT.abs(maxt3(l))) maxt3(l)=t3(i,j,k,a,b,c)
                                countt(l,3)=countt(l,3)+1
                                sumt3(l)=sumt3(l)+t3(i,j,k,a,b,c)**2
                                exit
                            end if
                        enddo
                        fullt=fullt+t3(i,j,k,a,b,c)**2
                        mcountt3=mcountt3+1
                    enddo
                enddo
            enddo
        enddo
    enddo

    sss=0
    do k = 2,maxNei
        if (abs(maxt3(k)).LT.gluCompare) maxt3(k)=0
        sss(1)=sss(1)+100*sumt3(k)   /fullt
        sss(2)=sss(2)+100*countt(k,3)/mcountt3
        inteContr(k)=sss(1); inteAmp(k)=sss(2)
    enddo

    do k = 2,maxNei
        efft3(k)=inteContr(k)/inteAmp(k)
    enddo

    write (ou,100)
    write (ou,101)
    do k = 2,maxNei
        write (ou,102) k,sumt3(k),100*sumt3(k)/fullt,maxt3(k),countt(k,3),efft3(k)
    enddo
    write (ou,103)
    write (ou,104) fullt,100,maxval(abs(maxt3)),mcountt3,1

    write (ou,105) 't3',sqrt(fullt/N)

100 format (/20X,'------------- T3 Norm Analysis -------------'/)
101 format (  2X,14X,'|',4X,'||T||',4X,'|',1X,'Contrib. %',1X,'|',3X,'max(T)',3X,'|',3X,'N(T)',3X,'|',2X,'Efficacy')
102 format (  3X,'||cue(',i<mid(maxNei)>,')||',<4-mid(maxNei)>X,'|',2X,F10.5,1X,'|',1X,F10.5,1X,'|',1X,ES10.3,1X,'|',1X,i8,1X,'|',1X,F9.3)
103 format (  2X,76('_'))
104 format (  2X,14X,'|', 2X,F10.5,1X, '|',1X,F10.5,1X,'|', 1X,ES10.3,1X,'|',1X,i8,1X,'|',1X,F9.3)
105 format (/ 4X,A,'-diagnostics:',1X,F9.6)

    return
    end subroutine t3_structure_analysis_by_layers

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine compound_by_excitation_level

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol,setfLineLen
    use printmod, only: prStrByVal,prMatrix,prVector
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu) :: dLayer, i,a,j,b,k,l,c,d,mm,nn,pp,vv,cper,lper
    real(kind=rglu)    :: VL(10),holdDist(28),hold(10),cvalue,value,sumw(4),sss(4),divisor,gsum,&
                          birav,birav2,lav


    write (ou,'(/)')

    sumw=0; hold=0
    do dlayer = 1,maxNei

        if ((ou.NE.0).AND.(ou.EQ.6)) then
            backspace(ou)
        endif
        write (ou,'(A,1X,i2,A,i2,1X,"|"\)') 'Layer',dlayer,'/',maxNei

        sss=0; cper=5
        do mm = 1,Nse
            i=spinInds(mm,1)
            a=spinInds(mm,2)
            lper=int(100.*float(mm)/float(Nse))
            if (lper.EQ.cper) then
                cper=cper+5; write (ou,'(1X,A\)') prStrByVal(lper)
            endif

            if (distMO(i,a).GT.mol%cueLayers(dlayer)) cycle

            if (inLayer(distMO(i,a),dlayer)) then
                VL(1)=cc_c1_t1([i,a])
                hold(1)=hold(1)+VL(1)*VL(1)
                cvalue=VL(1); value=cvalue*cvalue
                sss(1)=sss(1)+value

                holdContributions(i    ,dlayer,1)=holdContributions(i    ,dlayer,1)+value
                holdContributions(a-Nel,dlayer,2)=holdContributions(a-Nel,dlayer,2)+value
            endif

            !$omp parallel default(shared) private(holdDist,nn,pp,vv,j,b,k,c,l,d,VL,cvalue,value) reduction(+:sss,hold,holdContributions)
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
                    VL(2)=cc_c2_t2([i,j,a,b])
                    VL(3)=cc_c2_t1t1([i,j,a,b])
                    hold(2)=hold(2)+VL(2)*VL(2)
                    hold(3)=hold(3)+VL(3)*VL(3)
                    cvalue=VL(2)+VL(3); value=cvalue*cvalue
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
                        VL(4) = cc_c3_t2t1([i,j,k,a,b,c])
                        VL(5) = cc_c3_t1t1t1([i,j,k,a,b,c])
                        VL(6) = cc_c3_t3([i,j,k,a,b,c])

                        hold(4)=hold(4)+VL(4)*VL(4)
                        hold(5)=hold(5)+VL(5)*VL(5)
                        hold(9)=hold(9)+VL(6)*VL(6)
                        cvalue=VL(4)+VL(5)+VL(6); value=cvalue*cvalue

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
                            VL(7) =cc_c4_t2t2([i,j,k,l,a,b,c,d])
                            VL(8) =cc_c4_t2t1t1([i,j,k,l,a,b,c,d])
                            VL(9) =cc_c4_t1t1t1t1([i,j,k,l,a,b,c,d])
                            VL(10)=cc_c4_t3t1([i,j,k,l,a,b,c,d])

                            hold(6) =hold(6) +VL(7)*VL(7)
                            hold(7) =hold(7) +VL(8)*VL(8)
                            hold(8) =hold(8) +VL(9)*VL(9)
                            hold(10)=hold(10)+VL(10)*VL(10)
                            cvalue=VL(7)+VL(8)+VL(9)+VL(10); value=cvalue*cvalue

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

    write (ou,100)
    do dlayer = 1,maxNei
        write (ou,101) dlayer,holdValues(dlayer,:)
    enddo
    write (ou,102) holdValues(maxNei+1,:)
    write (ou,103) (1-holdValues(1,5))/holdValues(1,5)

    divisor=1
    divisor=sum(hold)
    !divisor=sum(sumw)

    write (ou,104) 'Contribution by type:'
    if (pattern(1)) write (ou,105) '                  t1:',100*hold(1) /divisor
    if (pattern(2)) write (ou,105) '                  t2:',100*hold(2) /divisor
    if (pattern(1)) write (ou,105) '               t1*t1:',100*hold(3) /divisor
    if (pattern(1)) write (ou,105) '               t1*t2:',100*hold(4) /divisor
    if (pattern(1)) write (ou,105) '            t1*t1*t1:',100*hold(5) /divisor
    if (pattern(3)) write (ou,105) '                  t3:',100*hold(9) /divisor
    if (pattern(2)) write (ou,105) '               t2*t2:',100*hold(6) /divisor
    if (pattern(1)) write (ou,105) '            t2*t1*t1:',100*hold(7) /divisor
    if (pattern(1)) write (ou,105) '         t1*t1*t1*t1:',100*hold(8) /divisor
    if (pattern(3)) write (ou,105) '               t1*t3:',100*hold(10)/divisor
    write (ou,104) repeat('_',28)
    if (pattern(1)) write (ou,105) '             pure t1:',100*(hold(1)+hold(3)+hold(5)+hold(8))/divisor
    if (pattern(2)) write (ou,105) '             pure t2:',100*(hold(2)+hold(6))/divisor
    if (pattern(3)) write (ou,105) '             pure t3:',100*(hold(9))/divisor
    write (ou,105)                 '               mixed:',100*(hold(4)+hold(7)+hold(10))/divisor
    write (ou,104) repeat('_',28)
    write (ou,105) ' overall uncertainty:',100*abs(sum(sumw)-sum(hold))                     /sum(sumw)
    write (ou,105) '   uncertainty in C1:',100*abs(sumw(1)-hold(1))                         /sumw(1)
    write (ou,105) '   uncertainty in C2:',100*abs(sumw(2)-hold(2)-hold(3))                 /sumw(2)
    write (ou,105) '   uncertainty in C3:',100*abs(sumw(3)-hold(4)-hold(5)-hold(9))         /sumw(3)
    write (ou,105) '   uncertainty in C4:',100*abs(sumw(4)-hold(6)-hold(7)-hold(8)-hold(10))/sumw(4)

    gsum=0
    do dLayer = 1,maxNei
        gsum=gsum+sum(holdContributions(:,dLayer,1))
        gsum=gsum+sum(holdContributions(:,dLayer,2))
    enddo

    do i = 2,Nel,2
        holdContributions(i/2,:,:)=holdContributions(i,:,:)+holdContributions(i-1,:,:)
    enddo

    holdContributions=100*holdContributions/gsum

    call prMatrix(holdContributions(1:Nel/2,:,1),ou,'Donate','00.00',transpose=true,maxwidth=setfLineLen)
    call prMatrix(holdContributions(1:Nel/2,:,2),ou,'Accept','00.00',transpose=true,maxwidth=setfLineLen)
    call prMatrix(holdContributions(1:Nel/2,:,1)+holdContributions(:,:,2),ou,'Total','00.00',transpose=true,maxwidth=setfLineLen)

    pvec=0
    do i = 1,Nel/2
        pvec(i)=sum(holdContributions(i,:,1))+sum(holdContributions(i,:,2))
    enddo

    call prVector(pvec,ou,'Summated','00.00','horizontal',setfLineLen)


100 format (/2X,'All values (except the last) are normalized on {W[1]+W[2]+W[3]+W[4]}.'/&
             2X,'L',5X,'Vl[1]',7X,'Vl[2]',7X,'Vl[3]',7X,'Vl[4]',9X,'Vl')
101 format ( 1X,i2,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7)
102 format ( 2X,61('-')/&
             2X,'W',1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7/)

103 format (/4X,'Theta:',1X,F8.3/)
104 format ( 4X,A)
105 format ( 4X,A,1X,F6.3)

    return
    end subroutine compound_by_excitation_level

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine biradical_character

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol,setfLineLen
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu) :: mm,nn,i,a,j,b,ui,ua,uj,ub
    real(kind=rglu)    :: birad(4),rav,value,Ax,Bx


    birad=0; rav=0
    do mm = 1,Nse
        i=spinInds(mm,1)
        a=spinInds(mm,2)

        do nn = 1,Nse
            j=spinInds(nn,1)
            b=spinInds(nn,2)

            if ((j.LE.i).OR.(b.LE.a)) cycle

            Ax=t2(i,j,a,b); Bx=t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a)
            value=(Ax+Bx)*(Ax+Bx)

            rav=rav+abs(distMO(i,j)-distMO(a,b))*value

            ui=i/2+mod(i,2); ua=(a-Nel)/2+mod(a-Nel,2)
            uj=j/2+mod(j,2); ub=(b-Nel)/2+mod(b-Nel,2)

            if ((ui.EQ.uj).AND.(ua.NE.ub)) then
                birad(1)=birad(1)+value
            elseif ((ui.NE.uj).AND.(ua.EQ.ub)) then
                birad(2)=birad(2)+value
            elseif ((ui.EQ.uj).AND.(ua.EQ.ub)) then
                birad(3)=birad(3)+value
            elseif ((ui.NE.uj).AND.(ua.NE.ub)) then
                birad(4)=birad(4)+value
            endif

        enddo
    enddo

    value=sum(birad)

    rav=rav/value
    write (ou,100) ( birad(i),100.*birad(i)/value, i=1,4 ),value,rav
100 format (//'C(iiab):  ',F8.4,1X,F6.2,'%'&
             /'C(ijaa):  ',F8.4,1X,F6.2,'%'&
             /'C(iiaa):  ',F8.4,1X,F6.2,'%'&
             /'C(ijab):  ',F8.4,1X,F6.2,'%'&
             /'sum:      ',F8.4            &
             /'Average R:',F8.4/)

    return
    end subroutine biradical_character

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine statistics_on_t1

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol,setfLineLen
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu) :: k,mm,ui,ua,i,a,count



    do k = 1,maxNei
        mean(k)=0; count=0; minv(k)=1; maxv(k)=0
        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)

            i=ui*2-1; a=ua*2-1
            if (inLayer(distMO(i,a),k)) then
                if (abs(t1(i,a)).GT.maxv(k)) maxv(k)=abs(t1(i,a))
                if (abs(t1(i,a)).LT.minv(k)) minv(k)=abs(t1(i,a))
                mean(k)=mean(k)+abs(t1(i,a))
                count=count+1
            end if
        enddo
        if (count.EQ.0) count=1
        mean(k)=mean(k)/count
    enddo

    write (ou,100) 'Statistics on t1'
    write (ou,101) 'k',' maxval','minval',' mean ','  sd  ','  number  '
    do k = 1,maxNei
        sd(k)=0; count=0
        do mm = 1,Ne
            ui=spatInds(mm,1)
            ua=spatInds(mm,2)

            i=ui*2-1; a=ua*2-1
            if (inLayer(distMO(i,a),k)) then
                sd(k)=sd(k)+(abs(t1(i,a))-mean(k))**2
                count=count+1
            end if
        enddo
        if (count.EQ.0) count=1

        sd(k)=sqrt(sd(k)/count)
        write (ou,102) k,maxv(k),minv(k),mean(k),sd(k),count
    enddo

100 format(/5X,A)
101 format( 5X,A,4(3X,A,3X),1X,A)
102 format( 4X,i2,4(1X,ES11.4),1X,i10)

    return
    end subroutine statistics_on_t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine statistics_on_t2

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol,setfLineLen
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu) :: k,mm,nn,ui,ua,uj,ub,i,a,j,b,count
    real(kind=rglu)    :: holdDist(6)


    do k = 1,maxNei
        mean(k)=0; minv(k)=1; maxv(k)=0; count=0
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
                    count=count+1
                end if
            enddo
        enddo

        if (count.EQ.0) count=1

        mean(k)=mean(k)/count
    enddo

    write (ou,100) 'Statistics for t2'
    write (ou,101) 'k',' maxval','minval',' mean ','  sd  ','  number  '
    do k = 1,maxNei
        sd(k)=0; count=0
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
                    count=count+1
                end if
            enddo
        enddo

        if (count.EQ.0) count=1

        sd(k)=sqrt(sd(k)/count)
        write (ou,102) k,maxv(k),minv(k),mean(k),sd(k),count
    enddo

100 format(/5X,A)
101 format( 5X,A,4(3X,A,3X),1X,A)
102 format( 4X,i2,4(1X,ES11.4),1X,i10)

    return
    end subroutine statistics_on_t2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine statistics_on_t3

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol,setfLineLen
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu) :: l,mm,nn,pp,gg,vv,ui,ua,uj,ub,uk,uc,i,a,j,b,k,c,count,jb(2,2),kc(2,2)
    real(kind=rglu)    :: holdDist(15)


    do l = 2,maxNei
        mean(l)=0; minv(l)=1; maxv(l)=0; count=0
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
                                count=count+1
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

        if (count.EQ.0) count=1

        mean(l)=mean(l)/count
    enddo

    write (ou,100) 'Statistics for t3'
    write (ou,101) 'k',' maxval','minval',' mean ','  sd  ','  number  '
    do l = 2,maxNei

        sd(l)=0; count=0
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
                                count=count+1
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

        if (count.EQ.0) count=1

        sd(l)=sqrt(sd(l)/count)
        write (ou,102) l,maxv(l),minv(l),mean(l),sd(l),count
    enddo

100 format(/5X,A)
101 format( 5X,A,4(3X,A,3X),1X,A)
102 format( 4X,i2,4(1X,ES11.4),1X,i10)

    return
    end subroutine statistics_on_t3

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine configuration_compound

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare,mid
    use hdb,      only: ccbd,ou,mol,setfLineLen
    use printmod, only: prStrByVal
    use coupledClusterAnalize

    implicit none

    integer(kind=iglu) :: i,a,j,b,ui,ua,uj,ub,zlen,ulen,t1plen,t2plen,ilen


    ilen=mid(Nel/2)
    zlen=7
    t1plen=15+2*ilen; t2plen=19+4*ilen

    write (ou,100)
    if (dCUE) then
        write (ou,101)
    else
        write (ou,102)
    endif

    write (ou,103); ulen=zlen

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
                write (ou,104) t1(i,a+Nel),ui,ua; ulen=zlen+t1plen
            else
                write (ou,105) t1(i,a+Nel),ui,ua; ulen=ulen+t1plen
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
                write (ou,106) t2(i,j,a+Nel,b+Nel),ui,ua,uj,ub; ulen=zlen+t2plen
            else
                write (ou,107) t2(i,j,a+Nel,b+Nel),ui,ua,uj,ub; ulen=ulen+t2plen
            endif
        endif

    enddo
    enddo
    enddo
    enddo

    write (ou,*)

100 format (//4X,'Configurational compound of wave-function')
101 format (  4X,"Orbitals correspond to ethylenes (i): i = occupied, i' = vacant.")
102 format (  4X,"Occupied orbitals: 1=HOMO, 2=HOMO-1, etc. Vacant orbitals: 1'=LUMO, 2'=LUMO+1, etc.")

103 format ('|PSI> ='\)

104 format (/<zlen>X,1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"'>"\)
105 format (         1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"'>"\)

106 format (/<zlen>X,1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"',",i<ilen>,'->',i<ilen>,"'>"\)
107 format (         1X,sp,ES9.2,ss,"|",i<ilen>,'->',i<ilen>,"',",i<ilen>,'->',i<ilen>,"'>"\)

    return
    end subroutine configuration_compound

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine compound_by_excitation_level_lr

    use glob,     only: rglu,iglu,lglu,true,false,gluCompare
    use hdb,      only: ccbd,ou,mol,setfLineLen
    use printmod, only: prStrByVal,prMatrix,prVector
    use coupledClusterAnalize

    implicit none

    logical(kind=lglu), external :: inLayer

    integer(kind=iglu) :: dLayer, i,a,j,b,k,l,c,d,mm,nn,pp,vv,cper,lper
    real(kind=rglu)    :: VL(4),holdDist(28),hold(10),cvalue,value,sumw(4),sss(4),divisor,gsum,&
                          birav,birav2,lav


    write (ou,'(//)')

    sumw=0; hold=0
    do dlayer = 1,maxNei

        if ((ou.NE.0).AND.(ou.EQ.6)) then
            backspace(ou)
        endif
        write (ou,'(A,1X,i2,A,i2,1X\)') 'Layer',dlayer,'/',maxNei

        sss=0; cper=5
        do mm = 1,Nse
            i=spinInds(mm,1)
            a=spinInds(mm,2)
            lper=int( 100.*float(mm)/float(Nse) )
            if (lper.EQ.cper) then
                cper=cper+5; write (ou,'(1X,A\)') prStrByVal(lper)
            endif

            if (distMO(i,a).GT.mol%cueLayers(dlayer)) cycle

            if (inLayer(distMO(i,a),dlayer)) then

                ! r0*T1 + R1
                VL(1)=r0*cc_c1_t1([i,a])+lr_c1_r1([i,a])
                sss(1)=sss(1)+VL(1)*VL(1)
            endif

            !$omp parallel default(shared) private(holdDist,nn,pp,vv,j,b,k,c,l,d,VL) reduction(+:sss)
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

                    ! r0*T2 + r0*T1T1 + R1T1 + R2
                    VL(2)=r0*(cc_c2_t2([i,j,a,b])+cc_c2_t1t1([i,j,a,b]))+lr_c2_r1t1([i,j,a,b])+&
                          lr_c2_r2([i,j,a,b])
                    sss(2)=sss(2)+VL(2)*VL(2)
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

                        ! r0*T1T1T1 + r0*T1T2 + R1T1T1 + R1T2 + R2T1
                        VL(3) = r0*(cc_c3_t1t1t1([i,j,k,a,b,c])+cc_c3_t2t1([i,j,k,a,b,c]))+&
                                lr_c3_r1t1t1([i,j,k,a,b,c])+lr_c3_r1t2([i,j,k,a,b,c])+&
                                lr_c3_r2t1([i,j,k,a,b,c])

                        sss(3)=sss(3)+VL(3)*VL(3)
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
                            ! r0*T1T1T1T1 + r0*T1T1T2 + r0*T2T2 + R1T1T1T1 + R1T1T2 + R2T1T1 + R2T2
                            VL(4)=r0*(cc_c4_t1t1t1t1([i,j,k,l,a,b,c,d])+cc_c4_t2t1t1([i,j,k,l,a,b,c,d])+&
                                  cc_c4_t2t2([i,j,k,l,a,b,c,d]))+lr_c4_r1t1t1t1([i,j,k,l,a,b,c,d])+&
                                  lr_c4_r1t2t1([i,j,k,l,a,b,c,d])+lr_c4_r2t1t1([i,j,k,l,a,b,c,d])+&
                                  lr_c4_r2t2([i,j,k,l,a,b,c,d])

                            sss(4)=sss(4)+VL(4)*VL(4)
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

    write (ou,100)
    do dlayer = 1,maxNei
        write (ou,101) dlayer,holdValues(dlayer,:)
    enddo
    write (ou,102) holdValues(maxNei+1,:)
    write (ou,103) (1-holdValues(1,5))/holdValues(1,5)


100 format (/2X,'All values (except the last) are normalized on {W[1]+W[2]+W[3]+W[4]}.'/&
             2X,'L',5X,'Vl[1]',7X,'Vl[2]',7X,'Vl[3]',7X,'Vl[4]',9X,'Vl')
101 format ( 1X,i2,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7)
102 format ( 2X,61('-')/&
             2X,'W',1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7,1X,F11.7/)

103 format (/4X,'Theta:',1X,F8.3/)

    return
    end subroutine compound_by_excitation_level_lr

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !
