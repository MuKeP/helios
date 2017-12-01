    module sorts

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob, only: true,false,void,uch,lglu,rglu,iglu,compareStrings
    use glob, only: r16kind, r8kind, r4kind
    use glob, only: i8kind , i4kind, i2kind, i1kind

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    ! todo: more sorts.
    character (len=*), parameter :: srVersion='1.000'
    character (len=*), parameter :: srDate   ='2017.12.01'
    character (len=*), parameter :: srAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    interface qsort
        module procedure qsort_wrapper_i8 ,qsort_wrapper_i4,qsort_wrapper_i2,qsort_wrapper_i1,&
                         qsort_wrapper_r16,qsort_wrapper_r8,qsort_wrapper_r4,&
                         qsort_wrapper_ch,qsort_wrapper_uc
    end interface qsort

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    ! module attributes
    public :: srVersion,srDate,srAuthor
    public :: qsort

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ QUICK SORT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_i8(array,rev) result(ret)
    implicit none
    integer(kind=i8kind)         :: array(:),swap
    logical(kind=lglu), optional :: rev
    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    asize=size(array)

    call qsort_i8(array,1,asize)
    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_i8(array,first,last)
    implicit none
    integer(kind=i8kind) :: array(:),middle,swap
    integer(kind=iglu)   :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i).LT.middle)
            i=i+1
        enddo
        do while(middle.LT.array(j))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_i8(array,first,i-1)
    if (j+1.LT.last)  call qsort_i8(array,j+1,last)

    return
    end subroutine qsort_i8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_i4(array,rev) result(ret)
    implicit none
    integer(kind=i4kind)         :: array(:),swap
    logical(kind=lglu), optional :: rev
    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    asize=size(array)

    call qsort_i4(array,1,asize)
    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_i4(array,first,last)
    implicit none
    integer(kind=i4kind) :: array(:),middle,swap
    integer(kind=iglu)   :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i).LT.middle)
            i=i+1
        enddo
        do while(middle.LT.array(j))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_i4(array,first,i-1)
    if (j+1.LT.last)  call qsort_i4(array,j+1,last)

    return
    end subroutine qsort_i4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_i2(array,rev) result(ret)
    implicit none
    integer(kind=i2kind)         :: array(:),swap
    logical(kind=lglu), optional :: rev
    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    asize=size(array)

    call qsort_i2(array,1,asize)
    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_i2(array,first,last)
    implicit none
    integer(kind=i2kind) :: array(:),middle,swap
    integer(kind=iglu)   :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i).LT.middle)
            i=i+1
        enddo
        do while(middle.LT.array(j))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_i2(array,first,i-1)
    if (j+1.LT.last)  call qsort_i2(array,j+1,last)

    return
    end subroutine qsort_i2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_i1(array,rev) result(ret)
    implicit none
    integer(kind=i1kind)         :: array(:),swap
    logical(kind=lglu), optional :: rev
    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    asize=size(array)

    call qsort_i1(array,1,asize)
    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_i1(array,first,last)
    implicit none
    integer(kind=i1kind) :: array(:),middle,swap
    integer(kind=iglu)   :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i).LT.middle)
            i=i+1
        enddo
        do while(middle.LT.array(j))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_i1(array,first,i-1)
    if (j+1.LT.last)  call qsort_i1(array,j+1,last)

    return
    end subroutine qsort_i1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_r16(array,rev) result(ret)
    implicit none
    real(kind=r16kind)           :: array(:),swap
    logical(kind=lglu), optional :: rev
    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    asize=size(array)

    call qsort_r16(array,1,asize)
    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_r16(array,first,last)
    implicit none
    real(kind=r16kind) :: array(:),middle,swap
    integer(kind=iglu) :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i).LT.middle)
            i=i+1
        enddo
        do while(middle.LT.array(j))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_r16(array,first,i-1)
    if (j+1.LT.last)  call qsort_r16(array,j+1,last)

    return
    end subroutine qsort_r16

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_r8(array,rev) result(ret)
    implicit none
    real(kind=r8kind)            :: array(:),swap
    logical(kind=lglu), optional :: rev
    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    asize=size(array)

    call qsort_r8(array,1,asize)
    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_r8(array,first,last)
    implicit none
    real(kind=r8kind)  :: array(:),middle,swap
    integer(kind=iglu) :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i).LT.middle)
            i=i+1
        enddo
        do while(middle.LT.array(j))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_r8(array,first,i-1)
    if (j+1.LT.last)  call qsort_r8(array,j+1,last)

    return
    end subroutine qsort_r8

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_r4(array,rev) result(ret)
    implicit none
    real(kind=r4kind)            :: array(:),swap
    logical(kind=lglu), optional :: rev
    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    asize=size(array)

    call qsort_r4(array,1,asize)
    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_r4(array,first,last)
    implicit none
    real(kind=r4kind)  :: array(:),middle,swap
    integer(kind=iglu) :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i).LT.middle)
            i=i+1
        enddo
        do while(middle.LT.array(j))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_r4(array,first,i-1)
    if (j+1.LT.last)  call qsort_r4(array,j+1,last)

    return
    end subroutine qsort_r4

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_ch(array,rev,attribute) result(ret)
    implicit none
    character (len=*)              :: array(:)
    character (len=:), allocatable :: swap

    logical(kind=lglu), optional :: rev
    character (len=*) , optional :: attribute
    character (len=8)            :: uattr

    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    uattr='alphabet'; if (present(attribute)) uattr=attribute
    asize=size(array)

    select case (trim(uattr))
        case ('alphabet'); call qsort_ch_alpha(array,1,asize)
        case ('length')  ; call qsort_ch_len(array,1,asize)
    end select

    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_ch

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_ch_len(array,first,last)
    implicit none
    character (len=*)              :: array(:)
    character (len=:), allocatable :: swap,middle

    integer(kind=iglu) :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(len_trim(array(i)).LT.len_trim(middle))
            i=i+1
        enddo
        do while(len_trim(middle).LT.len_trim(array(j)))
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_ch_len(array,first,i-1)
    if (j+1.LT.last)  call qsort_ch_len(array,j+1,last)

    return
    end subroutine qsort_ch_len

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_ch_alpha(array,first,last)
    implicit none
    character (len=*)              :: array(:)
    character (len=:), allocatable :: swap,middle

    integer(kind=iglu) :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(compareStrings(array(i),middle).LT.0)
            i=i+1
        enddo
        do while(compareStrings(middle,array(j)).LT.0)
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_ch_alpha(array,first,i-1)
    if (j+1.LT.last)  call qsort_ch_alpha(array,j+1,last)

    return
    end subroutine qsort_ch_alpha

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) function qsort_wrapper_uc(array,rev,attribute) result(ret)
    implicit none
    type(uch)                    :: array(:),swap
    logical(kind=lglu), optional :: rev
    character (len=*) , optional :: attribute
    character (len=8)            :: uattr

    logical(kind=lglu)           :: urev
    integer(kind=iglu)           :: k,asize


    urev=false; if (present(rev)) urev=rev
    uattr='alphabet'; if (present(attribute)) uattr=attribute
    asize=size(array)

    select case (trim(uattr))
        case ('alphabet'); call qsort_uc_alpha(array,1,asize)
        case ('length')  ; call qsort_uc_len(array,1,asize)
    end select

    if (urev) then
        do k = 1,asize/2
            swap=array(k)
            array(k)=array(asize-k+1)
            array(asize-k+1)=swap
        enddo
    endif

    ret=0; return
    end function qsort_wrapper_uc

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_uc_len(array,first,last)
    implicit none
    type(uch)          :: array(:),swap,middle
    integer(kind=iglu) :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(array(i)%ln.LT.middle%ln)
            i=i+1
        enddo
        do while(middle%ln.LT.array(j)%ln)
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_uc_len(array,first,i-1)
    if (j+1.LT.last)  call qsort_uc_len(array,j+1,last)

    return
    end subroutine qsort_uc_len

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    recursive subroutine qsort_uc_alpha(array,first,last)
    implicit none
    type(uch)          :: array(:),swap,middle
    integer(kind=iglu) :: first,last,i,j


    middle=array((first+last)/2)
    i=first; j=last
    do
        do while(compareStrings(array(i)%get(),middle%get()).LT.0)
            i=i+1
        enddo
        do while(compareStrings(middle%get(),array(j)%get()).LT.0)
            j=j-1
        enddo
        if (i.GE.j) exit
        swap=array(i); array(i)=array(j); array(j)=swap
        i=i+1
        j=j-1
    enddo
    if (first.LT.i-1) call qsort_uc_alpha(array,first,i-1)
    if (j+1.LT.last)  call qsort_uc_alpha(array,j+1,last)

    return
    end subroutine qsort_uc_alpha

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

end module sorts