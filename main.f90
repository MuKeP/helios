    program HELIOS

    use glob    , only: timecontrol,glMemoryLeft,timeControlCheckpoint,timestamp
    use glob    , only: void,signal,false,true,rglu,iglu,nullSub,pi,uch
    use hdb     , only: onLoad,trapSignals,generalbd,setParams,onTrap,ccbd,ou,finalizeHelios
    use hdb     , only: sighup,sigabrt,sigint,sigterm,sigcont,sigstop,sigusr1
    use property, only: getPolarizability,getEnergy,getRDM

    use sorts

    use hdb, only: mol,polarizbd,perturbate
    use printmod, only: prMatrix,prStrByVal

    use coupledCluster
    use lrccsdmodule

    implicit none

    real(kind=rglu) :: Ax,Bx,Eref,field

    integer(kind=iglu) :: mu,pX,pY,pZ,k
    real(kind=rglu)    :: sta,sto


    call onLoad; call trapSignals
    void=signal(SIGHUP , ontrap, -1)
    void=signal(SIGABRT, ontrap, -1)
    void=signal(SIGINT , ontrap, -1)
    void=signal(SIGTERM, ontrap, -1)
    void=signal(SIGCONT, ontrap, -1)
    void=signal(SIGSTOP, ontrap, -1)
    void=signal(SIGUSR1, ontrap, -1)

    call parseInput
    call setParams
    call readMoleculeInformation

!    ou=6

!    call setCCParameters('spin-cue-ccsd')
!    call initCC
!    call iterator(iterationCC,energyCC,ccbd%maxiters,ccbd%accuracy,false)
!
!    call setLRParameters('spin-cue-ccsd')
!    call initLR
!
!    stop

    field=0.01_rglu
    mol%perturbation=0
    pX=1; pY=2; pZ=3
    do mu = 1,mol%nAtoms
        mol%perturbation(mu,mu)=+pX*field*mol%atm(mu)%coords(1)&
                                +pY*field*mol%atm(mu)%coords(2)&
                                +pZ*field*mol%atm(mu)%coords(3)
    enddo
    !write (*,*) mol%perturbation
    !call perturbate


    !call prMatrix(mol%core,6,'Perturbated core (outside)','^.00000',maxwidth=79)

    void=timeControlCheckpoint('Started at',drop=true)
    void=timeControlCheckpoint('',raw=true)
    void=timeControlCheckpoint('Huckel    '//prStrByVal(getEnergy('huckel'   ),4,12),raw=true)
    void=timeControlCheckpoint('cue-ccs   '//prStrByVal(getEnergy('cue-ccs'  ),4,12),raw=true)
    void=timeControlCheckpoint('hf        '//prStrByVal(getEnergy('hf'       ),4,12),raw=true)
    void=timeControlCheckpoint('mp2       '//prStrByVal(getEnergy('mp2'      ),4,12),raw=true)
    void=timeControlCheckpoint('mp3       '//prStrByVal(getEnergy('mp3'      ),4,12),raw=true)
    void=timeControlCheckpoint('r-ccd     '//prStrByVal(getEnergy('r-ccd'    ),4,12),raw=true)
    void=timeControlCheckpoint('u-ccd     '//prStrByVal(getEnergy('u-ccd'    ),4,12),raw=true)
    void=timeControlCheckpoint('cue-ccsd  '//prStrByVal(getEnergy('cue-ccsd' ),4,12),raw=true)
    void=timeControlCheckpoint('r-ccsd    '//prStrByVal(getEnergy('r-ccsd'   ),4,12),raw=true)
    void=timeControlCheckpoint('u-ccsd    '//prStrByVal(getEnergy('u-ccsd'   ),4,12),raw=true)
    void=timeControlCheckpoint('r-ccsd(t) '//prStrByVal(getEnergy('r-ccsd(t)'),4,12),raw=true)
    !void=timeControlCheckpoint('Starting CCSDT part'//' (Dropped '//timestamp()//')',drop=true,raw=true)
    void=timeControlCheckpoint('cue-ccsdt '//prStrByVal(getEnergy('cue-ccsdt'),4,12),raw=true)
    void=timeControlCheckpoint('u-ccsdt   '//prStrByVal(getEnergy('u-ccsdt'  ),4,12),raw=true)
    void=timeControlCheckpoint('r-ccsdt   '//prStrByVal(getEnergy('r-ccsdt'  ),4,12),raw=true)
    void=timeControlCheckpoint('fci       '//prStrByVal(getEnergy('fci'      ),4,12),raw=true)
    void=timeControlCheckpoint('',raw=true)
    void=timeControlCheckpoint('Finished at',drop=true)
    stop



    !Ax=getEnergy('huckel')   ; write (*,*) 'huckel   ',Ax
    !Ax=getEnergy('cue-ccs')  ; write (*,*) 'cue-ccs  ',Ax
    !Ax=getEnergy('hf')       ; write (*,*) 'hf       ',Ax
    !Ax=getEnergy('mp2')      ; write (*,*) 'mp2      ',Ax
    !Ax=getEnergy('mp3')      ; write (*,*) 'mp3      ',Ax
    !Ax=getEnergy('r-ccd')    ; write (*,*) 'r-ccd    ',Ax
    !Ax=getEnergy('u-ccd')    ; write (*,*) 'u-ccd    ',Ax
    !Ax=getEnergy('cue-ccsd') ; write (*,*) 'cue-ccsd ',Ax
    !Ax=getEnergy('r-ccsd')   ; write (*,*) 'r-ccsd   ',Ax
    !Ax=getEnergy('u-ccsd')   ; write (*,*) 'u-ccsd   ',Ax
    !Ax=getEnergy('r-ccsd(t)'); write (*,*) 'r-ccsd(t)',Ax
    !Ax=getEnergy('cue-ccsdt'); write (*,*) 'cue-ccsdt',Ax
    !Ax=getEnergy('u-ccsdt')  ; write (*,*) 'u-ccsdt  ',Ax
    !Ax=getEnergy('r-ccsdt')  ; write (*,*) 'r-ccsdt  ',Ax
    !Ax=getEnergy('fci')      ; write (*,*) 'fci      ',Ax
    !stop

    sta=timecontrol()
    Ax=getEnergy('hf')         ; write (*,*) '0 hf       ',Ax; Bx=Ax
    Ax=getEnergy('hf',1)       ; write (*,*) '1 hf       ',Ax !+Bx
    Ax=getEnergy('hf',2)       ; write (*,*) '2 hf       ',Ax !+Bx
    !stop
    Ax=getEnergy('cue-ccsd')   ; write (*,*) '0 cue-ccsd ',Ax; Bx=Ax
    Ax=getEnergy('cue-ccsd',1) ; write (*,*) '1 cue-ccsd ',Ax !+Bx
    Ax=getEnergy('cue-ccsd',2) ; write (*,*) '2 cue-ccsd ',Ax !+Bx
    !stop
    Ax=getEnergy('r-ccsd')     ; write (*,*) '0 r-ccsd   ',Ax
    Ax=getEnergy('r-ccsd',1)   ; write (*,*) '1 r-ccsd   ',Ax
    Ax=getEnergy('r-ccsd',2)   ; write (*,*) '2 r-ccsd   ',Ax
    !stop
    Ax=getEnergy('u-ccsd')     ; write (*,*) '0 u-ccsd   ',Ax
    Ax=getEnergy('u-ccsd',1)   ; write (*,*) '1 u-ccsd   ',Ax
    Ax=getEnergy('u-ccsd',2)   ; write (*,*) '2 u-ccsd   ',Ax
    !stop
    Ax=getEnergy('fci')        ; write (*,*) '0 fci      ',Ax
    Ax=getEnergy('fci',1)      ; write (*,*) '1 fci      ',Ax
    Ax=getEnergy('fci',2)      ; write (*,*) '2 fci      ',Ax

    stop
    select case( generalbd%task%get() )
        case ('energy')        ; call nullSub
        case ('wf-analize')    ; call nullSub
        case ('density')       ; call getRDM
        case ('coulson')       ; call nullSub
        case ('polarizability'); call getPolarizability
        case ('hypercharges')  ; call nullSub
    end select

    call primaryInformation('end')

    call finalizeHelios

    stop
    end program HELIOS
