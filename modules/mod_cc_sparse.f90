    module coupledClusterSparse

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob, only: iglu,rglu,lglu,true,false

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    !   ~~~ CCSD ~~~   !
    real   (kind=rglu), dimension (:,:), allocatable :: pvd,t1,d1
    real   (kind=rglu), dimension   (:), allocatable :: vt,vd,vF

    integer(kind=iglu), dimension (:,:), allocatable :: t1mrEx,t1cIndex,t1Indexs
    integer(kind=iglu), dimension   (:), allocatable :: t1erow,t1numcol

    integer(kind=iglu), dimension (:,:), allocatable :: Indexs,Indexsbv,cIndex

    integer(kind=iglu), dimension   (:), allocatable :: erow,numcol,erowbv,numcolbv
    integer(kind=iglu), dimension (:,:), allocatable :: byExOrbbv,occEx,mrEx,byExOrb,ftfm

    integer(kind=iglu), dimension   (:), allocatable :: ferow,fnumcol,whOVf

    integer(kind=iglu), dimension (:,:), allocatable :: intersectOrbitals

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=iglu) :: Nue,t1Ne

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    public
    private :: iglu,rglu,lglu,true,false

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module coupledClusterSparse