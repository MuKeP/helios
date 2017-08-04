	subroutine iterator(dosub,maxiter,eps,oupattern)

	use glob

	implicit none

	interface
		subroutine dosub(iter,accuracy)
		use glob, only: iglu,rglu
		integer(kind=iglu), intent(in)  :: iter
		real   (kind=rglu), intent(out) :: accuracy(5)
		end subroutine dosub
	end interface

	integer(kind=iglu), intent(in) :: maxiter
	real   (kind=rglu), intent(in) :: eps
	character (len=*) , intent(in) :: oupattern









	return
	end subroutine iterator