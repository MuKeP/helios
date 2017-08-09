	module coupledClusterSparse

	use hdb, only: iglu,rglu,lglu,true,false

	!   ~~~ CCSD ~~~   !
	real   (kind=rglu), dimension (:,:), allocatable :: pvd,t1,d1
	real   (kind=rglu), dimension   (:), allocatable :: vt,vd,vF

	!   ~~~ Indexes ~~~   !
	integer(kind=iglu) :: Nue,t1Ne

	integer(kind=iglu), dimension (:,:), allocatable :: t1mrEx,t1cIndex,t1Indexs
	integer(kind=iglu), dimension   (:), allocatable :: t1erow,t1numcol

	integer(kind=iglu), dimension (:,:), allocatable :: Indexs,Indexsbv,iIndexs,cIndex

	integer(kind=iglu), dimension   (:), allocatable :: erow,numcol,erowbv,numcolbv
	integer(kind=iglu), dimension (:,:), allocatable :: byExOrbbv,occEx,mrEx,byExOrb,ftfm

	integer(kind=iglu), dimension   (:), allocatable :: ferow,fnumcol,whOVf

	integer(kind=iglu), dimension (:,:), allocatable :: intersectOrbitals

	end module coupledClusterSparse