!###################################################################################################
! This is the part of HELIOS
!
! Copyright (c) 2015-2020 V.N. Karazin Kharkiv National University, Anton B. Zakharov
!
! Author:
!   Anton Zakharov <abzakharov(at)karazin(dot)ua>
!
! Contribution:
!   Vladimir Ivanov <vivanov(at)karazin(dot)ua>
!###################################################################################################

    program HELIOS

    use glob,     only: void,signal,nullSub

    use hdb,      only: onLoad,trapSignals,generalbd,setParams,onTrap,ou
    use hdb,      only: sighup,sigabrt,sigint,sigterm,sigcont,sigstop,sigusr1

    use property, only: getPolarizability,getEnergy,getRDM,getCoulson
    use property, only: getWaveFunctionAnalize,getHypercharges,getOnlyEnergies

    implicit none


    call onLoad

    call trapSignals
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

    select case(generalbd%task%get())
        case ('energy')        ; call getOnlyEnergies
        case ('wf-analysis')   ; call getWaveFunctionAnalize
        case ('polarizability'); call getPolarizability
        case ('density')       ; call getRDM
        case ('coulson')       ; call getCoulson
        case ('hypercharges')  ; call getHypercharges
    end select

    call primaryInformation('end')

    end program HELIOS
