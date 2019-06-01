    subroutine checkInput

    use glob,       only: assignment(=),i8kind
    use glob,       only: true,false,void,iglu,glControlMemory
    use hdb,        only: generalbd,geometrybd,systembd,densitybd,mol,methodNames,hyperchargesbd
    use hdb,        only: coulsonbd,bdNames,eu,ou,init,ouWidth,controlFileName,ierror
    use printmod,   only: prJoin
    use txtParser,  only: tpIndex,tpCount,tpIntByStr,tpSplit,tpSplitLen,tpSplitHold,operator(.in.)
    use txtParser,  only: tpIntArrayByStr,tpReplace,tpDigits,tpSetAccordance

    implicit none

    integer(kind=iglu)              :: i,j,k,l,nobjs,count,pair(2)
    integer(kind=iglu), allocatable :: array(:)


    ! if ((densitybd%dtype%get() .in. ['all','charges','orders', 'scharges', 'sorders']).AND.(generalbd%task%get().NE.'density')) then
    !     write (eu,'(2X,A)') 'rdm:elements = '//densitybd%dtype%get()//' supposes general:jobtype = density.'; stop
    ! endif

    if (generalbd%task%get().EQ.'density') then

        if (densitybd%dtype%get().EQ.'none') then
            ierror='For general:jobtype = density choose any of rdm:elements = [all, charges, orders, scharges, sorders]'
            call primaryInformation('error')
        endif

        if (densitybd%NOAnalize .AND. (densitybd%dtype%get().NE.'all')) then
            ierror='For rdm:natural-orbitals = True, required rdm:elements = all'
            call primaryInformation('error')
        endif

        if (densitybd%dtype%get().EQ.'scharges') then
            if (densitybd%gcharges%get().NE.'none') then
                ierror='If rdm:elements = scharges, required rdm:collect-charges = none'
                call primaryInformation('error')
            endif

            if (densitybd%scharges%get().EQ.'none') then
                ierror='For rdm:elements = scharges, required rdm:selected-charges != none'
                call primaryInformation('error')
            endif
            void=tpSplit(densitybd%scharges%get(),';'); nobjs=tpSplitLen

            if (nobjs.EQ.0) then
                ierror='For rdm:selected-charges must contain at least one value'
                call primaryInformation('error')
            endif

            do i = 1,nobjs
                k=tpIntByStr(tpSplitHold(i)%get())
                if ((k.LE.0).OR.(k.GT.mol%nAtoms)) then
                    ierror='For rdm:selected-charges values must be: N >= value > 0'
                    call primaryInformation('error')
                endif
            enddo
        endif

        if (densitybd%dtype%get().EQ.'sorders') then
            if (densitybd%gcharges%get().NE.'none') then
                ierror='If rdm:elements = sorders, required rdm:collect-charges = none'
                call primaryInformation('error')
            endif

            void=tpSplit(densitybd%sorders%get(),';'); nobjs=tpSplitLen
            do k = 1,nobjs
                pair=tpIntArrayByStr(tpSplitHold(k)%get(),2)
                i=pair(1); j=pair(2)

                if (((i.LE.0).OR.(i.GT.mol%nAtoms)).OR.((j.LE.0).OR.(j.GT.mol%nAtoms))) then
                    ierror='For rdm:selected-orders values must be: N >= value > 0'
                    call primaryInformation('error')
                endif
            enddo
        endif

        if (densitybd%gcharges%get().NE.'none') then
            void=tpSplit(densitybd%gcharges%get(),';'); nobjs=tpSplitLen
            do k = 1,nobjs
                count=tpCount(tpSplitHold(k)%get(),',')+1
                allocate (array(count)); array=tpIntArrayByStr(tpSplitHold(k)%get(),count)
                do l = 1,count
                    if ((array(l).GT.mol%nAtoms).OR.(array(l).LE.0)) then
                        ierror='For rdm:collect-charges values must be: N >= value > 0'
                        call primaryInformation('error')
                    endif
                enddo
                deallocate (array)
            enddo
        endif

    endif

    if (generalbd%task%get().EQ.'hypercharges') then

        if (hyperchargesbd%dtype%get().NE.'scharges') then

            if (hyperchargesbd%gcharges%get().NE.'none') then
                ierror='If hypercharges:elements = scharges, required rdm:collect-charges = none'
                call primaryInformation('error')
            endif

            if (hyperchargesbd%scharges%get().EQ.'none') then
                ierror='For hypercharges:elements = scharges, required rdm:selected-charges != none'
                call primaryInformation('error')
            endif

            void=tpSplit(hyperchargesbd%scharges%get(),';'); nobjs=tpSplitLen

            if (nobjs.EQ.0) then
                ierror='For hypercharges:selected-charges must contain at least one value or be equal "none"'
                call primaryInformation('error')
            endif

            do i = 1,nobjs
                k=tpIntByStr(tpSplitHold(i)%get())
                if ((k.LE.0).OR.(k.GT.mol%nAtoms)) then
                    ierror='For hypercharges:selected-charges values must be: N >= value > 0'
                    call primaryInformation('error')
                endif
            enddo
        endif

        if (hyperchargesbd%gcharges%get().NE.'none') then
            void=tpSplit(hyperchargesbd%gcharges%get(),';'); nobjs=tpSplitLen
            do k = 1,nobjs
                count=tpCount(tpSplitHold(k)%get(),',')+1
                allocate (array(count)); array=tpIntArrayByStr(tpSplitHold(k)%get(),count)
                do l = 1,count
                    if ((array(l).GT.mol%nAtoms).OR.(array(l).LE.0)) then
                        ierror='For hypercharges:collect-charges values must be: N >= value > 0'
                        call primaryInformation('error')
                    endif
                enddo
                deallocate (array)
            enddo
        endif

    endif

    if (generalbd%task%get().EQ.'coulson') then

        if (coulsonbd%ctype%get().EQ.'atom-bond') then
            if ((coulsonbd%selected%get().EQ.'diagonals').OR.(coulsonbd%selected%get().EQ.'offdiagonals')) then
                ierror='For coulson:elements = atom-bond, any coulson:selected != none is not acceptable'
                call primaryInformation('error')
            endif
        endif

        if (coulsonbd%ctype%get().EQ.'none') then
            ierror='For general:jobtype = coulson choose any of coulson:elements = [atom-atom, atom-bond, bond-bond]'
            call primaryInformation('error')
        endif

        if (coulsonbd%selected%get().EQ.'diagonals') then
            if (coulsonbd%soffdiagonals%get().NE.'none') then
                ierror='If coulson:elements = diagonals, required coulson:selected-offdiagonals = none'
                call primaryInformation('error')
            endif

            if (coulsonbd%sdiagonals%get().EQ.'none') then
                ierror='For coulson:elements = diagonals, required coulson:selected-diagonals != none'
                call primaryInformation('error')
            endif

            if (.NOT.tpSetAccordance(coulsonbd%sdiagonals%get(),tpDigits()//';')) then
                ierror='For coulson:elements = diagonals, incorrect syntax'
                call primaryInformation('error')
            endif
            void=tpSplit(coulsonbd%sdiagonals%get(),';'); nobjs=tpSplitLen

            if (nobjs.EQ.0) then
                ierror='For coulson:selected-diagonals must contain at least one value'
                call primaryInformation('error')
            endif

            do i = 1,nobjs
                k=tpIntByStr(tpSplitHold(i)%get())
                if (coulsonbd%ctype%get().EQ.'atom-atom') then
                    if ((k.LE.0).OR.(k.GT.mol%nAtoms)) then
                        ierror='For coulson:selected-diagonals values must be: Natoms >= value > 0'
                        call primaryInformation('error')
                    endif
                elseif (coulsonbd%ctype%get().EQ.'bond-bond') then
                    if ((k.LE.0).OR.(k.GT.mol%nBonds)) then
                        ierror='For coulson:selected-diagonals values must be: Nbonds >= value > 0'
                        call primaryInformation('error')
                    endif
                endif
            enddo
        endif

        if (coulsonbd%selected%get().EQ.'offdiagonals') then
            if (coulsonbd%sdiagonals%get().NE.'none') then
                ierror='If coulson:elements = offdiagonals, required coulson:collect-diagonals = none'
                call primaryInformation('error')
            endif

            void=tpSplit(coulsonbd%soffdiagonals%get(),';'); nobjs=tpSplitLen
            do k = 1,nobjs
                pair=tpIntArrayByStr(tpSplitHold(k)%get(),2)
                i=pair(1); j=pair(2)

                if (coulsonbd%ctype%get().EQ.'atom-atom') then
                    if (((i.LE.0).OR.(i.GT.mol%nAtoms)).OR.((j.LE.0).OR.(j.GT.mol%nAtoms))) then
                        ierror='For coulson:selected-offdiagonals values must be: Natoms >= value > 0'
                        call primaryInformation('error')
                    endif

                elseif (coulsonbd%ctype%get().EQ.'bond-bond') then
                    if (((i.LE.0).OR.(i.GT.mol%nBonds)).OR.((j.LE.0).OR.(j.GT.mol%nBonds))) then
                        ierror='For coulson:selected-offdiagonals values must be: Nbonds >= value > 0'
                        call primaryInformation('error')
                    endif
                endif

            enddo
        endif
    endif

    if (generalbd%methods%get().EQ.'all') then
        generalbd%methods=prJoin(methodNames,'+')
    endif

    ! if ('%methods%' .in. systembd%throughHeader%get()) then
    !     systembd%throughHeader=tpReplace(systembd%throughHeader%get(), '%methods%', generalbd%methods%get())
    ! endif

    generalbd%bondsAlternated=abs(generalbd%alternation).GT.0
    geometrybd%searchLinear(2)=false
    geometrybd%searchPlanar(2)=false

    return
    end subroutine checkInput