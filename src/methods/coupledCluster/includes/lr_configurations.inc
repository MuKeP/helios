
    real(kind=rglu) function lr_c1_r1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,a1


    ret=0

    i1=arr(1)
    a1=arr(2)

    ret=ret+r1(i1,a1)

    return
    end function lr_c1_r1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c2_r2(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,a1,a2


    ret=0

    i1=arr(1); i2=arr(2)
    a1=arr(3); a2=arr(4)

    ret=ret+r2(i1,i2,a1,a2)

    return
    end function lr_c2_r2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c2_r1t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,a1,a2


    ret=0

    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2)
    a1=arr(3); a2=arr(4)

    ret=ret+r1(i1,a1)*t1(i2,a2)&
           -r1(i1,a2)*t1(i2,a1)&
           -r1(i2,a1)*t1(i1,a2)&
           +r1(i2,a2)*t1(i1,a1)

    return
    end function lr_c2_r1t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c3_r2t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,i3,a1,a2,a3


    ret=0

    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2); i3=arr(3)
    a1=arr(4); a2=arr(5); a3=arr(6)

    ret=ret+r2(i1,i2,a1,a2)*t1(i3,a3)&
           -r2(i1,i2,a1,a3)*t1(i3,a2)&
           +r2(i1,i2,a2,a3)*t1(i3,a1)&
           -r2(i1,i3,a1,a2)*t1(i2,a3)&
           +r2(i1,i3,a1,a3)*t1(i2,a2)&
           -r2(i1,i3,a2,a3)*t1(i2,a1)&
           +r2(i2,i3,a1,a2)*t1(i1,a3)&
           -r2(i2,i3,a1,a3)*t1(i1,a2)&
           +r2(i2,i3,a2,a3)*t1(i1,a1)

    return
    end function lr_c3_r2t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c3_r1t2(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,i3,a1,a2,a3


    ret=0

    if (.NOT.pattern(2)) return

    i1=arr(1); i2=arr(2); i3=arr(3)
    a1=arr(4); a2=arr(5); a3=arr(6)

    ret=ret+r1(i1,a1)*t2(i2,i3,a2,a3)&
           -r1(i1,a2)*t2(i2,i3,a1,a3)&
           +r1(i1,a3)*t2(i2,i3,a1,a2)&
           -r1(i2,a1)*t2(i1,i3,a2,a3)&
           +r1(i2,a2)*t2(i1,i3,a1,a3)&
           -r1(i2,a3)*t2(i1,i3,a1,a2)&
           +r1(i3,a1)*t2(i1,i2,a2,a3)&
           -r1(i3,a2)*t2(i1,i2,a1,a3)&
           +r1(i3,a3)*t2(i1,i2,a1,a2)

    return
    end function lr_c3_r1t2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c3_r1t1t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,i3,a1,a2,a3


    ret=0

    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2); i3=arr(3)
    a1=arr(4); a2=arr(5); a3=arr(6)

    ret=ret+r1(i1,a1)*t1(i2,a2)*t1(i3,a3)&
           -r1(i1,a1)*t1(i2,a3)*t1(i3,a2)&
           -r1(i1,a2)*t1(i2,a1)*t1(i3,a3)&
           +r1(i1,a2)*t1(i2,a3)*t1(i3,a1)&
           +r1(i1,a3)*t1(i2,a1)*t1(i3,a2)&
           -r1(i1,a3)*t1(i2,a2)*t1(i3,a1)&
           -r1(i2,a1)*t1(i1,a2)*t1(i3,a3)&
           +r1(i2,a1)*t1(i1,a3)*t1(i3,a2)&
           +r1(i2,a2)*t1(i1,a1)*t1(i3,a3)&
           -r1(i2,a2)*t1(i1,a3)*t1(i3,a1)&
           -r1(i2,a3)*t1(i1,a1)*t1(i3,a2)&
           +r1(i2,a3)*t1(i1,a2)*t1(i3,a1)&
           +r1(i3,a1)*t1(i1,a2)*t1(i2,a3)&
           -r1(i3,a1)*t1(i1,a3)*t1(i2,a2)&
           -r1(i3,a2)*t1(i1,a1)*t1(i2,a3)&
           +r1(i3,a2)*t1(i1,a3)*t1(i2,a1)&
           +r1(i3,a3)*t1(i1,a1)*t1(i2,a2)&
           -r1(i3,a3)*t1(i1,a2)*t1(i2,a1)

    return
    end function lr_c3_r1t1t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c4_r2t2(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,i3,i4,a1,a2,a3,a4


    ret=0

    if (.NOT.pattern(2)) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+r2(i1,i2,a1,a2)*t2(i3,i4,a3,a4)&
           -r2(i1,i2,a1,a3)*t2(i3,i4,a2,a4)&
           +r2(i1,i2,a1,a4)*t2(i3,i4,a2,a3)&
           +r2(i1,i2,a2,a3)*t2(i3,i4,a1,a4)&
           -r2(i1,i2,a2,a4)*t2(i3,i4,a1,a3)&
           +r2(i1,i2,a3,a4)*t2(i3,i4,a1,a2)&
           -r2(i1,i3,a1,a2)*t2(i2,i4,a3,a4)&
           +r2(i1,i3,a1,a3)*t2(i2,i4,a2,a4)&
           -r2(i1,i3,a1,a4)*t2(i2,i4,a2,a3)&
           -r2(i1,i3,a2,a3)*t2(i2,i4,a1,a4)&
           +r2(i1,i3,a2,a4)*t2(i2,i4,a1,a3)&
           -r2(i1,i3,a3,a4)*t2(i2,i4,a1,a2)&
           +r2(i1,i4,a1,a2)*t2(i2,i3,a3,a4)&
           -r2(i1,i4,a1,a3)*t2(i2,i3,a2,a4)&
           +r2(i1,i4,a1,a4)*t2(i2,i3,a2,a3)&
           +r2(i1,i4,a2,a3)*t2(i2,i3,a1,a4)&
           -r2(i1,i4,a2,a4)*t2(i2,i3,a1,a3)&
           +r2(i1,i4,a3,a4)*t2(i2,i3,a1,a2)&
           +r2(i2,i3,a1,a2)*t2(i1,i4,a3,a4)&
           -r2(i2,i3,a1,a3)*t2(i1,i4,a2,a4)&
           +r2(i2,i3,a1,a4)*t2(i1,i4,a2,a3)&
           +r2(i2,i3,a2,a3)*t2(i1,i4,a1,a4)&
           -r2(i2,i3,a2,a4)*t2(i1,i4,a1,a3)&
           +r2(i2,i3,a3,a4)*t2(i1,i4,a1,a2)&
           -r2(i2,i4,a1,a2)*t2(i1,i3,a3,a4)&
           +r2(i2,i4,a1,a3)*t2(i1,i3,a2,a4)&
           -r2(i2,i4,a1,a4)*t2(i1,i3,a2,a3)&
           -r2(i2,i4,a2,a3)*t2(i1,i3,a1,a4)&
           +r2(i2,i4,a2,a4)*t2(i1,i3,a1,a3)&
           -r2(i2,i4,a3,a4)*t2(i1,i3,a1,a2)&
           +r2(i3,i4,a1,a2)*t2(i1,i2,a3,a4)&
           -r2(i3,i4,a1,a3)*t2(i1,i2,a2,a4)&
           +r2(i3,i4,a1,a4)*t2(i1,i2,a2,a3)&
           +r2(i3,i4,a2,a3)*t2(i1,i2,a1,a4)&
           -r2(i3,i4,a2,a4)*t2(i1,i2,a1,a3)&
           +r2(i3,i4,a3,a4)*t2(i1,i2,a1,a2)

    return
    end function lr_c4_r2t2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c4_r2t1t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,i3,i4,a1,a2,a3,a4


    ret=0

    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+r2(i1,i2,a1,a2)*t1(i3,a3)*t1(i4,a4)&
           -r2(i1,i2,a1,a2)*t1(i3,a4)*t1(i4,a3)&
           -r2(i1,i2,a1,a3)*t1(i3,a2)*t1(i4,a4)&
           +r2(i1,i2,a1,a3)*t1(i3,a4)*t1(i4,a2)&
           +r2(i1,i2,a1,a4)*t1(i3,a2)*t1(i4,a3)&
           -r2(i1,i2,a1,a4)*t1(i3,a3)*t1(i4,a2)&
           +r2(i1,i2,a2,a3)*t1(i3,a1)*t1(i4,a4)&
           -r2(i1,i2,a2,a3)*t1(i3,a4)*t1(i4,a1)&
           -r2(i1,i2,a2,a4)*t1(i3,a1)*t1(i4,a3)&
           +r2(i1,i2,a2,a4)*t1(i3,a3)*t1(i4,a1)&
           +r2(i1,i2,a3,a4)*t1(i3,a1)*t1(i4,a2)&
           -r2(i1,i2,a3,a4)*t1(i3,a2)*t1(i4,a1)&
           -r2(i1,i3,a1,a2)*t1(i2,a3)*t1(i4,a4)&
           +r2(i1,i3,a1,a2)*t1(i2,a4)*t1(i4,a3)&
           +r2(i1,i3,a1,a3)*t1(i2,a2)*t1(i4,a4)&
           -r2(i1,i3,a1,a3)*t1(i2,a4)*t1(i4,a2)&
           -r2(i1,i3,a1,a4)*t1(i2,a2)*t1(i4,a3)&
           +r2(i1,i3,a1,a4)*t1(i2,a3)*t1(i4,a2)&
           -r2(i1,i3,a2,a3)*t1(i2,a1)*t1(i4,a4)&
           +r2(i1,i3,a2,a3)*t1(i2,a4)*t1(i4,a1)&
           +r2(i1,i3,a2,a4)*t1(i2,a1)*t1(i4,a3)&
           -r2(i1,i3,a2,a4)*t1(i2,a3)*t1(i4,a1)&
           -r2(i1,i3,a3,a4)*t1(i2,a1)*t1(i4,a2)&
           +r2(i1,i3,a3,a4)*t1(i2,a2)*t1(i4,a1)&
           +r2(i1,i4,a1,a2)*t1(i2,a3)*t1(i3,a4)&
           -r2(i1,i4,a1,a2)*t1(i2,a4)*t1(i3,a3)&
           -r2(i1,i4,a1,a3)*t1(i2,a2)*t1(i3,a4)&
           +r2(i1,i4,a1,a3)*t1(i2,a4)*t1(i3,a2)&
           +r2(i1,i4,a1,a4)*t1(i2,a2)*t1(i3,a3)&
           -r2(i1,i4,a1,a4)*t1(i2,a3)*t1(i3,a2)&
           +r2(i1,i4,a2,a3)*t1(i2,a1)*t1(i3,a4)&
           -r2(i1,i4,a2,a3)*t1(i2,a4)*t1(i3,a1)&
           -r2(i1,i4,a2,a4)*t1(i2,a1)*t1(i3,a3)&
           +r2(i1,i4,a2,a4)*t1(i2,a3)*t1(i3,a1)&
           +r2(i1,i4,a3,a4)*t1(i2,a1)*t1(i3,a2)&
           -r2(i1,i4,a3,a4)*t1(i2,a2)*t1(i3,a1)&
           +r2(i2,i3,a1,a2)*t1(i1,a3)*t1(i4,a4)&
           -r2(i2,i3,a1,a2)*t1(i1,a4)*t1(i4,a3)&
           -r2(i2,i3,a1,a3)*t1(i1,a2)*t1(i4,a4)&
           +r2(i2,i3,a1,a3)*t1(i1,a4)*t1(i4,a2)&
           +r2(i2,i3,a1,a4)*t1(i1,a2)*t1(i4,a3)&
           -r2(i2,i3,a1,a4)*t1(i1,a3)*t1(i4,a2)&
           +r2(i2,i3,a2,a3)*t1(i1,a1)*t1(i4,a4)&
           -r2(i2,i3,a2,a3)*t1(i1,a4)*t1(i4,a1)&
           -r2(i2,i3,a2,a4)*t1(i1,a1)*t1(i4,a3)&
           +r2(i2,i3,a2,a4)*t1(i1,a3)*t1(i4,a1)&
           +r2(i2,i3,a3,a4)*t1(i1,a1)*t1(i4,a2)&
           -r2(i2,i3,a3,a4)*t1(i1,a2)*t1(i4,a1)&
           -r2(i2,i4,a1,a2)*t1(i1,a3)*t1(i3,a4)&
           +r2(i2,i4,a1,a2)*t1(i1,a4)*t1(i3,a3)&
           +r2(i2,i4,a1,a3)*t1(i1,a2)*t1(i3,a4)&
           -r2(i2,i4,a1,a3)*t1(i1,a4)*t1(i3,a2)&
           -r2(i2,i4,a1,a4)*t1(i1,a2)*t1(i3,a3)&
           +r2(i2,i4,a1,a4)*t1(i1,a3)*t1(i3,a2)&
           -r2(i2,i4,a2,a3)*t1(i1,a1)*t1(i3,a4)&
           +r2(i2,i4,a2,a3)*t1(i1,a4)*t1(i3,a1)&
           +r2(i2,i4,a2,a4)*t1(i1,a1)*t1(i3,a3)&
           -r2(i2,i4,a2,a4)*t1(i1,a3)*t1(i3,a1)&
           -r2(i2,i4,a3,a4)*t1(i1,a1)*t1(i3,a2)&
           +r2(i2,i4,a3,a4)*t1(i1,a2)*t1(i3,a1)&
           +r2(i3,i4,a1,a2)*t1(i1,a3)*t1(i2,a4)&
           -r2(i3,i4,a1,a2)*t1(i1,a4)*t1(i2,a3)&
           -r2(i3,i4,a1,a3)*t1(i1,a2)*t1(i2,a4)&
           +r2(i3,i4,a1,a3)*t1(i1,a4)*t1(i2,a2)&
           +r2(i3,i4,a1,a4)*t1(i1,a2)*t1(i2,a3)&
           -r2(i3,i4,a1,a4)*t1(i1,a3)*t1(i2,a2)&
           +r2(i3,i4,a2,a3)*t1(i1,a1)*t1(i2,a4)&
           -r2(i3,i4,a2,a3)*t1(i1,a4)*t1(i2,a1)&
           -r2(i3,i4,a2,a4)*t1(i1,a1)*t1(i2,a3)&
           +r2(i3,i4,a2,a4)*t1(i1,a3)*t1(i2,a1)&
           +r2(i3,i4,a3,a4)*t1(i1,a1)*t1(i2,a2)&
           -r2(i3,i4,a3,a4)*t1(i1,a2)*t1(i2,a1)

    return
    end function lr_c4_r2t1t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c4_r1t2t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,i3,i4,a1,a2,a3,a4


    ret=0

    if (.NOT.pattern(1)) return
    if (.NOT.pattern(2)) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+r1(i1,a1)*t2(i2,i3,a2,a3)*t1(i4,a4)&
           -r1(i1,a1)*t2(i2,i3,a2,a4)*t1(i4,a3)&
           +r1(i1,a1)*t2(i2,i3,a3,a4)*t1(i4,a2)&
           -r1(i1,a1)*t2(i2,i4,a2,a3)*t1(i3,a4)&
           +r1(i1,a1)*t2(i2,i4,a2,a4)*t1(i3,a3)&
           -r1(i1,a1)*t2(i2,i4,a3,a4)*t1(i3,a2)&
           +r1(i1,a1)*t2(i3,i4,a2,a3)*t1(i2,a4)&
           -r1(i1,a1)*t2(i3,i4,a2,a4)*t1(i2,a3)&
           +r1(i1,a1)*t2(i3,i4,a3,a4)*t1(i2,a2)&
           -r1(i1,a2)*t2(i2,i3,a1,a3)*t1(i4,a4)&
           +r1(i1,a2)*t2(i2,i3,a1,a4)*t1(i4,a3)&
           -r1(i1,a2)*t2(i2,i3,a3,a4)*t1(i4,a1)&
           +r1(i1,a2)*t2(i2,i4,a1,a3)*t1(i3,a4)&
           -r1(i1,a2)*t2(i2,i4,a1,a4)*t1(i3,a3)&
           +r1(i1,a2)*t2(i2,i4,a3,a4)*t1(i3,a1)&
           -r1(i1,a2)*t2(i3,i4,a1,a3)*t1(i2,a4)&
           +r1(i1,a2)*t2(i3,i4,a1,a4)*t1(i2,a3)&
           -r1(i1,a2)*t2(i3,i4,a3,a4)*t1(i2,a1)&
           +r1(i1,a3)*t2(i2,i3,a1,a2)*t1(i4,a4)&
           -r1(i1,a3)*t2(i2,i3,a1,a4)*t1(i4,a2)&
           +r1(i1,a3)*t2(i2,i3,a2,a4)*t1(i4,a1)&
           -r1(i1,a3)*t2(i2,i4,a1,a2)*t1(i3,a4)&
           +r1(i1,a3)*t2(i2,i4,a1,a4)*t1(i3,a2)&
           -r1(i1,a3)*t2(i2,i4,a2,a4)*t1(i3,a1)&
           +r1(i1,a3)*t2(i3,i4,a1,a2)*t1(i2,a4)&
           -r1(i1,a3)*t2(i3,i4,a1,a4)*t1(i2,a2)&
           +r1(i1,a3)*t2(i3,i4,a2,a4)*t1(i2,a1)&
           -r1(i1,a4)*t2(i2,i3,a1,a2)*t1(i4,a3)&
           +r1(i1,a4)*t2(i2,i3,a1,a3)*t1(i4,a2)&
           -r1(i1,a4)*t2(i2,i3,a2,a3)*t1(i4,a1)&
           +r1(i1,a4)*t2(i2,i4,a1,a2)*t1(i3,a3)&
           -r1(i1,a4)*t2(i2,i4,a1,a3)*t1(i3,a2)&
           +r1(i1,a4)*t2(i2,i4,a2,a3)*t1(i3,a1)&
           -r1(i1,a4)*t2(i3,i4,a1,a2)*t1(i2,a3)&
           +r1(i1,a4)*t2(i3,i4,a1,a3)*t1(i2,a2)&
           -r1(i1,a4)*t2(i3,i4,a2,a3)*t1(i2,a1)&
           -r1(i2,a1)*t2(i1,i3,a2,a3)*t1(i4,a4)&
           +r1(i2,a1)*t2(i1,i3,a2,a4)*t1(i4,a3)&
           -r1(i2,a1)*t2(i1,i3,a3,a4)*t1(i4,a2)&
           +r1(i2,a1)*t2(i1,i4,a2,a3)*t1(i3,a4)&
           -r1(i2,a1)*t2(i1,i4,a2,a4)*t1(i3,a3)&
           +r1(i2,a1)*t2(i1,i4,a3,a4)*t1(i3,a2)&
           -r1(i2,a1)*t2(i3,i4,a2,a3)*t1(i1,a4)&
           +r1(i2,a1)*t2(i3,i4,a2,a4)*t1(i1,a3)&
           -r1(i2,a1)*t2(i3,i4,a3,a4)*t1(i1,a2)&
           +r1(i2,a2)*t2(i1,i3,a1,a3)*t1(i4,a4)&
           -r1(i2,a2)*t2(i1,i3,a1,a4)*t1(i4,a3)&
           +r1(i2,a2)*t2(i1,i3,a3,a4)*t1(i4,a1)&
           -r1(i2,a2)*t2(i1,i4,a1,a3)*t1(i3,a4)&
           +r1(i2,a2)*t2(i1,i4,a1,a4)*t1(i3,a3)&
           -r1(i2,a2)*t2(i1,i4,a3,a4)*t1(i3,a1)&
           +r1(i2,a2)*t2(i3,i4,a1,a3)*t1(i1,a4)&
           -r1(i2,a2)*t2(i3,i4,a1,a4)*t1(i1,a3)&
           +r1(i2,a2)*t2(i3,i4,a3,a4)*t1(i1,a1)&
           -r1(i2,a3)*t2(i1,i3,a1,a2)*t1(i4,a4)&
           +r1(i2,a3)*t2(i1,i3,a1,a4)*t1(i4,a2)&
           -r1(i2,a3)*t2(i1,i3,a2,a4)*t1(i4,a1)&
           +r1(i2,a3)*t2(i1,i4,a1,a2)*t1(i3,a4)&
           -r1(i2,a3)*t2(i1,i4,a1,a4)*t1(i3,a2)&
           +r1(i2,a3)*t2(i1,i4,a2,a4)*t1(i3,a1)&
           -r1(i2,a3)*t2(i3,i4,a1,a2)*t1(i1,a4)&
           +r1(i2,a3)*t2(i3,i4,a1,a4)*t1(i1,a2)&
           -r1(i2,a3)*t2(i3,i4,a2,a4)*t1(i1,a1)&
           +r1(i2,a4)*t2(i1,i3,a1,a2)*t1(i4,a3)&
           -r1(i2,a4)*t2(i1,i3,a1,a3)*t1(i4,a2)&
           +r1(i2,a4)*t2(i1,i3,a2,a3)*t1(i4,a1)&
           -r1(i2,a4)*t2(i1,i4,a1,a2)*t1(i3,a3)&
           +r1(i2,a4)*t2(i1,i4,a1,a3)*t1(i3,a2)&
           -r1(i2,a4)*t2(i1,i4,a2,a3)*t1(i3,a1)&
           +r1(i2,a4)*t2(i3,i4,a1,a2)*t1(i1,a3)&
           -r1(i2,a4)*t2(i3,i4,a1,a3)*t1(i1,a2)&
           +r1(i2,a4)*t2(i3,i4,a2,a3)*t1(i1,a1)&
           +r1(i3,a1)*t2(i1,i2,a2,a3)*t1(i4,a4)&
           -r1(i3,a1)*t2(i1,i2,a2,a4)*t1(i4,a3)&
           +r1(i3,a1)*t2(i1,i2,a3,a4)*t1(i4,a2)&
           -r1(i3,a1)*t2(i1,i4,a2,a3)*t1(i2,a4)&
           +r1(i3,a1)*t2(i1,i4,a2,a4)*t1(i2,a3)&
           -r1(i3,a1)*t2(i1,i4,a3,a4)*t1(i2,a2)&
           +r1(i3,a1)*t2(i2,i4,a2,a3)*t1(i1,a4)&
           -r1(i3,a1)*t2(i2,i4,a2,a4)*t1(i1,a3)&
           +r1(i3,a1)*t2(i2,i4,a3,a4)*t1(i1,a2)&
           -r1(i3,a2)*t2(i1,i2,a1,a3)*t1(i4,a4)&
           +r1(i3,a2)*t2(i1,i2,a1,a4)*t1(i4,a3)&
           -r1(i3,a2)*t2(i1,i2,a3,a4)*t1(i4,a1)&
           +r1(i3,a2)*t2(i1,i4,a1,a3)*t1(i2,a4)&
           -r1(i3,a2)*t2(i1,i4,a1,a4)*t1(i2,a3)&
           +r1(i3,a2)*t2(i1,i4,a3,a4)*t1(i2,a1)&
           -r1(i3,a2)*t2(i2,i4,a1,a3)*t1(i1,a4)&
           +r1(i3,a2)*t2(i2,i4,a1,a4)*t1(i1,a3)&
           -r1(i3,a2)*t2(i2,i4,a3,a4)*t1(i1,a1)&
           +r1(i3,a3)*t2(i1,i2,a1,a2)*t1(i4,a4)&
           -r1(i3,a3)*t2(i1,i2,a1,a4)*t1(i4,a2)&
           +r1(i3,a3)*t2(i1,i2,a2,a4)*t1(i4,a1)&
           -r1(i3,a3)*t2(i1,i4,a1,a2)*t1(i2,a4)&
           +r1(i3,a3)*t2(i1,i4,a1,a4)*t1(i2,a2)&
           -r1(i3,a3)*t2(i1,i4,a2,a4)*t1(i2,a1)&
           +r1(i3,a3)*t2(i2,i4,a1,a2)*t1(i1,a4)&
           -r1(i3,a3)*t2(i2,i4,a1,a4)*t1(i1,a2)&
           +r1(i3,a3)*t2(i2,i4,a2,a4)*t1(i1,a1)&
           -r1(i3,a4)*t2(i1,i2,a1,a2)*t1(i4,a3)&
           +r1(i3,a4)*t2(i1,i2,a1,a3)*t1(i4,a2)&
           -r1(i3,a4)*t2(i1,i2,a2,a3)*t1(i4,a1)&
           +r1(i3,a4)*t2(i1,i4,a1,a2)*t1(i2,a3)&
           -r1(i3,a4)*t2(i1,i4,a1,a3)*t1(i2,a2)&
           +r1(i3,a4)*t2(i1,i4,a2,a3)*t1(i2,a1)&
           -r1(i3,a4)*t2(i2,i4,a1,a2)*t1(i1,a3)&
           +r1(i3,a4)*t2(i2,i4,a1,a3)*t1(i1,a2)&
           -r1(i3,a4)*t2(i2,i4,a2,a3)*t1(i1,a1)&
           -r1(i4,a1)*t2(i1,i2,a2,a3)*t1(i3,a4)&
           +r1(i4,a1)*t2(i1,i2,a2,a4)*t1(i3,a3)&
           -r1(i4,a1)*t2(i1,i2,a3,a4)*t1(i3,a2)&
           +r1(i4,a1)*t2(i1,i3,a2,a3)*t1(i2,a4)&
           -r1(i4,a1)*t2(i1,i3,a2,a4)*t1(i2,a3)&
           +r1(i4,a1)*t2(i1,i3,a3,a4)*t1(i2,a2)&
           -r1(i4,a1)*t2(i2,i3,a2,a3)*t1(i1,a4)&
           +r1(i4,a1)*t2(i2,i3,a2,a4)*t1(i1,a3)&
           -r1(i4,a1)*t2(i2,i3,a3,a4)*t1(i1,a2)&
           +r1(i4,a2)*t2(i1,i2,a1,a3)*t1(i3,a4)&
           -r1(i4,a2)*t2(i1,i2,a1,a4)*t1(i3,a3)&
           +r1(i4,a2)*t2(i1,i2,a3,a4)*t1(i3,a1)&
           -r1(i4,a2)*t2(i1,i3,a1,a3)*t1(i2,a4)&
           +r1(i4,a2)*t2(i1,i3,a1,a4)*t1(i2,a3)&
           -r1(i4,a2)*t2(i1,i3,a3,a4)*t1(i2,a1)&
           +r1(i4,a2)*t2(i2,i3,a1,a3)*t1(i1,a4)&
           -r1(i4,a2)*t2(i2,i3,a1,a4)*t1(i1,a3)&
           +r1(i4,a2)*t2(i2,i3,a3,a4)*t1(i1,a1)&
           -r1(i4,a3)*t2(i1,i2,a1,a2)*t1(i3,a4)&
           +r1(i4,a3)*t2(i1,i2,a1,a4)*t1(i3,a2)&
           -r1(i4,a3)*t2(i1,i2,a2,a4)*t1(i3,a1)&
           +r1(i4,a3)*t2(i1,i3,a1,a2)*t1(i2,a4)&
           -r1(i4,a3)*t2(i1,i3,a1,a4)*t1(i2,a2)&
           +r1(i4,a3)*t2(i1,i3,a2,a4)*t1(i2,a1)&
           -r1(i4,a3)*t2(i2,i3,a1,a2)*t1(i1,a4)&
           +r1(i4,a3)*t2(i2,i3,a1,a4)*t1(i1,a2)&
           -r1(i4,a3)*t2(i2,i3,a2,a4)*t1(i1,a1)&
           +r1(i4,a4)*t2(i1,i2,a1,a2)*t1(i3,a3)&
           -r1(i4,a4)*t2(i1,i2,a1,a3)*t1(i3,a2)&
           +r1(i4,a4)*t2(i1,i2,a2,a3)*t1(i3,a1)&
           -r1(i4,a4)*t2(i1,i3,a1,a2)*t1(i2,a3)&
           +r1(i4,a4)*t2(i1,i3,a1,a3)*t1(i2,a2)&
           -r1(i4,a4)*t2(i1,i3,a2,a3)*t1(i2,a1)&
           +r1(i4,a4)*t2(i2,i3,a1,a2)*t1(i1,a3)&
           -r1(i4,a4)*t2(i2,i3,a1,a3)*t1(i1,a2)&
           +r1(i4,a4)*t2(i2,i3,a2,a3)*t1(i1,a1)

    return
    end function lr_c4_r1t2t1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    real(kind=rglu) function lr_c4_r1t1t1t1(arr) result(ret)
    implicit none

    integer(kind=iglu), intent(in) :: arr(:)
    integer(kind=iglu)             :: i1,i2,i3,i4,a1,a2,a3,a4


    ret=0

    if (.NOT.pattern(1)) return

    i1=arr(1); i2=arr(2); i3=arr(3); i4=arr(4)
    a1=arr(5); a2=arr(6); a3=arr(7); a4=arr(8)

    ret=ret+r1(i1,a1)*t1(i2,a2)*t1(i3,a3)*t1(i4,a4)&
           -r1(i1,a1)*t1(i2,a2)*t1(i3,a4)*t1(i4,a3)&
           -r1(i1,a1)*t1(i2,a3)*t1(i3,a2)*t1(i4,a4)&
           +r1(i1,a1)*t1(i2,a3)*t1(i3,a4)*t1(i4,a2)&
           +r1(i1,a1)*t1(i2,a4)*t1(i3,a2)*t1(i4,a3)&
           -r1(i1,a1)*t1(i2,a4)*t1(i3,a3)*t1(i4,a2)&
           -r1(i1,a2)*t1(i2,a1)*t1(i3,a3)*t1(i4,a4)&
           +r1(i1,a2)*t1(i2,a1)*t1(i3,a4)*t1(i4,a3)&
           +r1(i1,a2)*t1(i2,a3)*t1(i3,a1)*t1(i4,a4)&
           -r1(i1,a2)*t1(i2,a3)*t1(i3,a4)*t1(i4,a1)&
           -r1(i1,a2)*t1(i2,a4)*t1(i3,a1)*t1(i4,a3)&
           +r1(i1,a2)*t1(i2,a4)*t1(i3,a3)*t1(i4,a1)&
           +r1(i1,a3)*t1(i2,a1)*t1(i3,a2)*t1(i4,a4)&
           -r1(i1,a3)*t1(i2,a1)*t1(i3,a4)*t1(i4,a2)&
           -r1(i1,a3)*t1(i2,a2)*t1(i3,a1)*t1(i4,a4)&
           +r1(i1,a3)*t1(i2,a2)*t1(i3,a4)*t1(i4,a1)&
           +r1(i1,a3)*t1(i2,a4)*t1(i3,a1)*t1(i4,a2)&
           -r1(i1,a3)*t1(i2,a4)*t1(i3,a2)*t1(i4,a1)&
           -r1(i1,a4)*t1(i2,a1)*t1(i3,a2)*t1(i4,a3)&
           +r1(i1,a4)*t1(i2,a1)*t1(i3,a3)*t1(i4,a2)&
           +r1(i1,a4)*t1(i2,a2)*t1(i3,a1)*t1(i4,a3)&
           -r1(i1,a4)*t1(i2,a2)*t1(i3,a3)*t1(i4,a1)&
           -r1(i1,a4)*t1(i2,a3)*t1(i3,a1)*t1(i4,a2)&
           +r1(i1,a4)*t1(i2,a3)*t1(i3,a2)*t1(i4,a1)&
           -r1(i2,a1)*t1(i1,a2)*t1(i3,a3)*t1(i4,a4)&
           +r1(i2,a1)*t1(i1,a2)*t1(i3,a4)*t1(i4,a3)&
           +r1(i2,a1)*t1(i1,a3)*t1(i3,a2)*t1(i4,a4)&
           -r1(i2,a1)*t1(i1,a3)*t1(i3,a4)*t1(i4,a2)&
           -r1(i2,a1)*t1(i1,a4)*t1(i3,a2)*t1(i4,a3)&
           +r1(i2,a1)*t1(i1,a4)*t1(i3,a3)*t1(i4,a2)&
           +r1(i2,a2)*t1(i1,a1)*t1(i3,a3)*t1(i4,a4)&
           -r1(i2,a2)*t1(i1,a1)*t1(i3,a4)*t1(i4,a3)&
           -r1(i2,a2)*t1(i1,a3)*t1(i3,a1)*t1(i4,a4)&
           +r1(i2,a2)*t1(i1,a3)*t1(i3,a4)*t1(i4,a1)&
           +r1(i2,a2)*t1(i1,a4)*t1(i3,a1)*t1(i4,a3)&
           -r1(i2,a2)*t1(i1,a4)*t1(i3,a3)*t1(i4,a1)&
           -r1(i2,a3)*t1(i1,a1)*t1(i3,a2)*t1(i4,a4)&
           +r1(i2,a3)*t1(i1,a1)*t1(i3,a4)*t1(i4,a2)&
           +r1(i2,a3)*t1(i1,a2)*t1(i3,a1)*t1(i4,a4)&
           -r1(i2,a3)*t1(i1,a2)*t1(i3,a4)*t1(i4,a1)&
           -r1(i2,a3)*t1(i1,a4)*t1(i3,a1)*t1(i4,a2)&
           +r1(i2,a3)*t1(i1,a4)*t1(i3,a2)*t1(i4,a1)&
           +r1(i2,a4)*t1(i1,a1)*t1(i3,a2)*t1(i4,a3)&
           -r1(i2,a4)*t1(i1,a1)*t1(i3,a3)*t1(i4,a2)&
           -r1(i2,a4)*t1(i1,a2)*t1(i3,a1)*t1(i4,a3)&
           +r1(i2,a4)*t1(i1,a2)*t1(i3,a3)*t1(i4,a1)&
           +r1(i2,a4)*t1(i1,a3)*t1(i3,a1)*t1(i4,a2)&
           -r1(i2,a4)*t1(i1,a3)*t1(i3,a2)*t1(i4,a1)&
           +r1(i3,a1)*t1(i1,a2)*t1(i2,a3)*t1(i4,a4)&
           -r1(i3,a1)*t1(i1,a2)*t1(i2,a4)*t1(i4,a3)&
           -r1(i3,a1)*t1(i1,a3)*t1(i2,a2)*t1(i4,a4)&
           +r1(i3,a1)*t1(i1,a3)*t1(i2,a4)*t1(i4,a2)&
           +r1(i3,a1)*t1(i1,a4)*t1(i2,a2)*t1(i4,a3)&
           -r1(i3,a1)*t1(i1,a4)*t1(i2,a3)*t1(i4,a2)&
           -r1(i3,a2)*t1(i1,a1)*t1(i2,a3)*t1(i4,a4)&
           +r1(i3,a2)*t1(i1,a1)*t1(i2,a4)*t1(i4,a3)&
           +r1(i3,a2)*t1(i1,a3)*t1(i2,a1)*t1(i4,a4)&
           -r1(i3,a2)*t1(i1,a3)*t1(i2,a4)*t1(i4,a1)&
           -r1(i3,a2)*t1(i1,a4)*t1(i2,a1)*t1(i4,a3)&
           +r1(i3,a2)*t1(i1,a4)*t1(i2,a3)*t1(i4,a1)&
           +r1(i3,a3)*t1(i1,a1)*t1(i2,a2)*t1(i4,a4)&
           -r1(i3,a3)*t1(i1,a1)*t1(i2,a4)*t1(i4,a2)&
           -r1(i3,a3)*t1(i1,a2)*t1(i2,a1)*t1(i4,a4)&
           +r1(i3,a3)*t1(i1,a2)*t1(i2,a4)*t1(i4,a1)&
           +r1(i3,a3)*t1(i1,a4)*t1(i2,a1)*t1(i4,a2)&
           -r1(i3,a3)*t1(i1,a4)*t1(i2,a2)*t1(i4,a1)&
           -r1(i3,a4)*t1(i1,a1)*t1(i2,a2)*t1(i4,a3)&
           +r1(i3,a4)*t1(i1,a1)*t1(i2,a3)*t1(i4,a2)&
           +r1(i3,a4)*t1(i1,a2)*t1(i2,a1)*t1(i4,a3)&
           -r1(i3,a4)*t1(i1,a2)*t1(i2,a3)*t1(i4,a1)&
           -r1(i3,a4)*t1(i1,a3)*t1(i2,a1)*t1(i4,a2)&
           +r1(i3,a4)*t1(i1,a3)*t1(i2,a2)*t1(i4,a1)&
           -r1(i4,a1)*t1(i1,a2)*t1(i2,a3)*t1(i3,a4)&
           +r1(i4,a1)*t1(i1,a2)*t1(i2,a4)*t1(i3,a3)&
           +r1(i4,a1)*t1(i1,a3)*t1(i2,a2)*t1(i3,a4)&
           -r1(i4,a1)*t1(i1,a3)*t1(i2,a4)*t1(i3,a2)&
           -r1(i4,a1)*t1(i1,a4)*t1(i2,a2)*t1(i3,a3)&
           +r1(i4,a1)*t1(i1,a4)*t1(i2,a3)*t1(i3,a2)&
           +r1(i4,a2)*t1(i1,a1)*t1(i2,a3)*t1(i3,a4)&
           -r1(i4,a2)*t1(i1,a1)*t1(i2,a4)*t1(i3,a3)&
           -r1(i4,a2)*t1(i1,a3)*t1(i2,a1)*t1(i3,a4)&
           +r1(i4,a2)*t1(i1,a3)*t1(i2,a4)*t1(i3,a1)&
           +r1(i4,a2)*t1(i1,a4)*t1(i2,a1)*t1(i3,a3)&
           -r1(i4,a2)*t1(i1,a4)*t1(i2,a3)*t1(i3,a1)&
           -r1(i4,a3)*t1(i1,a1)*t1(i2,a2)*t1(i3,a4)&
           +r1(i4,a3)*t1(i1,a1)*t1(i2,a4)*t1(i3,a2)&
           +r1(i4,a3)*t1(i1,a2)*t1(i2,a1)*t1(i3,a4)&
           -r1(i4,a3)*t1(i1,a2)*t1(i2,a4)*t1(i3,a1)&
           -r1(i4,a3)*t1(i1,a4)*t1(i2,a1)*t1(i3,a2)&
           +r1(i4,a3)*t1(i1,a4)*t1(i2,a2)*t1(i3,a1)&
           +r1(i4,a4)*t1(i1,a1)*t1(i2,a2)*t1(i3,a3)&
           -r1(i4,a4)*t1(i1,a1)*t1(i2,a3)*t1(i3,a2)&
           -r1(i4,a4)*t1(i1,a2)*t1(i2,a1)*t1(i3,a3)&
           +r1(i4,a4)*t1(i1,a2)*t1(i2,a3)*t1(i3,a1)&
           +r1(i4,a4)*t1(i1,a3)*t1(i2,a1)*t1(i3,a2)&
           -r1(i4,a4)*t1(i1,a3)*t1(i2,a2)*t1(i3,a1)

    return
    end function lr_c4_r1t1t1t1

