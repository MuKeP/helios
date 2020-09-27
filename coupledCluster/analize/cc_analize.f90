    subroutine wfAnalize(method, tdlist, prHeader)

    use glob,                  only: iglu,lglu,void,true,false
    use hdb,                   only: ou,setfLineLen
    use txtParser,             only: tpAdjustc
    use coupledClusterAnalize, only: chapsLen,umethod,dCUE,prepareForAnalize,finishAnalize,pattern,&
                                     mode

    implicit none

    character(len=*)  , intent(in)           :: method
    integer(kind=iglu), intent(in)           :: tdlist
    logical(kind=lglu), intent(in), optional :: prHeader

    logical(kind=lglu) :: chaps(chapsLen), dprHeader
    integer(kind=iglu) :: i

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ parse todolist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! 2#01111111
    !    7654321

    ! 1 | print t1/t2/t3 by layers    | 1
    ! 2 | contribution to energy      | 2
    ! 3 | structure analize by layers | 4
    ! 4 | C1,C2,C3,C4 estimation      | 8
    ! 5 | biradical character         | 16
    ! 6 | statistics on t1/t2/t3      | 32
    ! 7 | configurational compound    | 64
    !   1 | 2 | 3 | 4| 5 | 6 | 7 = 127 = 2#01111111

    dprHeader=true; if (present(prHeader)) dprHeader=prHeader

    do i = 1,chapsLen
        chaps(i)=btest(tdlist,i-1)
    enddo

    if (tdlist.NE.0) then
        if (dprHeader) write (ou,'(//A/)') tpAdjustc(' Wave-function analysis ',setfLineLen,'=')
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    umethod=method
    void=prepareForAnalize()

    if (mode.EQ.'lr') then
        chaps=false
        call compound_by_excitation_level_lr
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! print t1 amplitudes greater than threshold
    if (chaps(1).AND.pattern(1)) then
        if (dCUE) then
            call print_t1_by_layer_local
        else
            call print_t1_by_layer_nonlocal
        endif
    endif

    ! print t2 amplitudes greater than threshold
    if (chaps(1).AND.pattern(2)) then
        if (dCUE) then
            call print_t2_by_layer_local
        else
            call print_t2_by_layer_nonlocal
        endif
    endif

    ! print t3 amplitudes greater than threshold
    if (chaps(1).AND.pattern(3)) then
        if (dCUE) then
            call print_t3_by_layer_local
        else
            call print_t3_by_layer_nonlocal
        endif
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! contribution to the correlation energy by layers
    if (chaps(2)) then
        if (dCUE) then
            call correlation_energy_contribution_local
        else
            call correlation_energy_contribution_nonlocal
        endif
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! t1 structure analize by layers
    if (chaps(3).AND.pattern(1)) then
        if (dCUE) then
            call t1_structure_analysis_by_layers
        endif
    endif

    ! t2 structure analize by layers
    if (chaps(3).AND.pattern(2)) then
        if (dCUE) then
            call t2_structure_analysis_by_layers
        endif
    endif

    ! t3 structure analize by layers
    if (chaps(3).AND.pattern(3)) then
        if (dCUE) then
            call t3_structure_analysis_by_layers
        endif
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! C1, C2, C3, C4 estimation
    if (chaps(4)) then
        call compound_by_excitation_level
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! biradical character
    if (chaps(5)) then
        call biradical_character
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! statistics for t1
    if (chaps(6).AND.pattern(1)) then
        if (dCUE) then
            call statistics_on_t1
        endif
    endif

    ! statistics for t2
    if (chaps(6).AND.pattern(2)) then
        if (dCUE) then
            call statistics_on_t2
        endif
    endif

    ! statistics for t3
    if (chaps(6).AND.pattern(3)) then
        if (dCUE) then
            call statistics_on_t3
        endif
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! configurational compound
    if (chaps(7)) then
        call configuration_compound
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    void=finishAnalize()
    if (tdlist.NE.0) then
        if (dprHeader) write (ou,'(/A/)') tpAdjustc(' End of analysis ',setfLineLen,'=')
    endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    return
    end subroutine wfAnalize
