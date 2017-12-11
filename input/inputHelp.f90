    subroutine inputHelp

    use glob,      only: iglu,void
    use hdb,       only: su,showInputHelp,bdNames
    use datablock, only: bdPrintHelp
    use txtParser, only: operator(.in.)

    implicit none

    integer(kind=iglu) :: k


    if (showInputHelp%get().EQ.'all') then
        write (su,'(/2X,A/)') 'Following blocks are available:'
        do k = 1,UBound(bdNames,1)
            write (su,'(2X,i2,")",1X,A)') k,bdNames(k)
        enddo
        write (su,'(/2X,A)') 'To show help on every block type "-i %blockname"'
    elseif (.NOT. (showInputHelp%get() .in. bdNames)) then
        write (su,'(/2X,A)') 'Undefined block name '//showInputHelp%get()//' found. type "-i all" to see available.'
    else
        void=bdPrintHelp(showInputHelp%get(),su,78)
    endif

    stop
    end subroutine inputHelp