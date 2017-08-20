	subroutine projection_lrccsd_spin_hf_intermediates_1

	use glob        , only: iglu,rglu
	use lrccsdModule, only: Nel,No,ou
	use lrccsdModule, only: R,F,t1,t2
	use lrccsdModule, only: ai1,ai2,ai3,ai4,ai5,ai6,ai7,Fab,Fij,Fia

	implicit none

	integer(kind=iglu) :: i,j,k,l,a,b,c,d
	real   (kind=rglu) :: sum,ss1,ss2,sum1,sum2,sum3,sum4,sum5,sum6,sum7


	!write (ou,'(/1X,A)') trim(ttime())//" Intermediates' calculation" 

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  1/10' 
	!$omp parallel default(shared) private(i,a,b,c,j,sum)
	!$omp do
	do i = 1,Nel
	do a = Nel+1,No
	do b = Nel+1,No
	do c = Nel+1,No
		sum=0
		do j = 1,Nel
			sum=sum+R(i,b,j,c)*t1(j,a)
		enddo
		ai7(i,a,b,c)=R(i,b,a,c)-sum
	enddo
	enddo
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  2/10' 
	!$omp parallel default(shared) private(i,j,k,a,c,sum)
	!$omp do
	do i = 1,Nel
	do j = 1,Nel
	do k = 1,Nel
	do a = Nel+1,No
		sum=0
		do c = Nel+1,No
			sum=sum+R(i,a,j,c)*t1(k,c)
		enddo
		ai6(i,j,k,a)=R(i,k,j,a)-sum
	enddo
	enddo
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  3/10'
	!$omp parallel default(shared) private(a,b,c,d,k,l,sum1,sum2,sum3)
	!$omp do
	do a = Nel+1,No
	do b = Nel+1,No
	do c = Nel+1,No
	do d = Nel+1,No
		sum1=0
		do k = 1,Nel
			sum1=sum1+t1(k,a)*R(k,c,b,d)-t1(k,b)*R(k,c,a,d)
		enddo

		sum2=0
		do k = 1,Nel
			do l = 1,Nel
				sum2=sum2+t1(k,a)*t1(l,b)*R(k,c,l,d)
			enddo
		enddo

		sum3=0
		do k = 1,Nel
			do l = 1,Nel
				sum3=sum3+t2(k,l,a,b)*R(k,c,l,d)
			enddo
		enddo
		 
		ai5(a,b,c,d)=(R(a,c,b,d)-sum1+sum2)/2+sum3/4
	enddo
	enddo
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  4/10' 
	!$omp parallel default(shared) private(i,j,k,l,c,d,sum1,sum2,sum3)
	!$omp do
	do i = 1,Nel
	do j = 1,Nel
	do k = 1,Nel
	do l = 1,Nel
		sum1=0
		do c = Nel+1,No
			sum1=sum1+t1(j,c)*R(k,i,l,c)-t1(i,c)*R(k,j,l,c)
		enddo

		sum2=0
		do c = Nel+1,No
			do d = Nel+1,No
				sum2=sum2+t1(i,c)*t1(j,d)*R(k,c,l,d)
			enddo
		enddo
		 
		sum3=0
		do c = Nel+1,No
			do d = Nel+1,No
				sum3=sum3+t2(i,j,c,d)*R(k,c,l,d)
			enddo
		enddo

		ai4(i,j,k,l)=(R(i,k,j,l)+sum1+sum2)/2+sum3/4
	enddo
	enddo
	enddo
	enddo
	!$omp end parallel
 
	!write (lrsu,'(1X,A)') trim(ttime())//' Part  5/10' 
	!$omp parallel default(shared) private(i,c,a,b,l,k,d,sum1,sum2,ss1,ss2,sum3,sum4,sum5,sum6,sum7)
	!$omp do
	do i = 1,Nel
	do c = Nel+1,No
	do a = Nel+1,No
	do b = Nel+1,No
		sum1=0
		do d = Nel+1,No
			sum1=sum1+t1(i,d)*ai5(a,b,c,d)
		enddo

		sum2=0
		do l = 1,Nel
			ss1=R(l,i,a,c)
			ss2=R(l,i,b,c)
			do k = 1,Nel
			do d = Nel+1,No
				ss1=ss1-t2(i,k,a,d)*R(l,c,k,d)
				ss2=ss2-t2(i,k,b,d)*R(l,c,k,d)
			enddo
			enddo
			sum2=sum2+t1(l,b)*ss1-t1(l,a)*ss2
		enddo

		sum3=0
		do l = 1,Nel
		do k = 1,Nel
			sum3=sum3+t1(l,a)*t1(k,b)*R(c,l,i,k)
		enddo
		enddo

		sum4=0
		do l = 1,Nel
		do k = 1,Nel
		do d = Nel+1,No
			sum4=sum4+t1(l,d)*t2(k,i,a,b)*R(c,k,d,l)
		enddo
		enddo
		enddo

		sum5=0
		do k = 1,Nel
			sum5=sum5+t2(i,k,a,b)*F(k,c)
		enddo

		sum6=0
		do k = 1,Nel
		do l = 1,Nel
			sum6=sum6+t2(k,l,a,b)*R(i,k,c,l)
		enddo
		enddo

		sum7=0
		do k = 1,Nel
		do d = Nel+1,No
			sum7=sum7+t2(i,k,a,d)*R(k,c,b,d)-t2(i,k,b,d)*R(k,c,a,d)	
		enddo
		enddo

		ai3(i,c,a,b)=-R(i,a,c,b)+2*sum1-sum2+sum3-sum4+sum5-sum6/2+sum7
	enddo
	enddo
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  6/10' 
	!$omp parallel default(shared) private(i,j,k,a,l,d,c,sum1,sum2,sum3,sum4,sum5,sum6,sum7,ss1,ss2)
	!$omp do		
	do i = 1,Nel
	do j = 1,Nel
	do k = 1,Nel
	do a = Nel+1,No

		sum1=0
		do l = 1,Nel
			sum1=sum1+t1(l,a)*ai4(i,j,k,l)
		enddo

		sum2=0
		do d = Nel+1,No
			ss1=R(j,k,d,a)
			ss2=R(i,k,d,a)
			do l = 1,Nel
			do c = Nel+1,No
				ss1=ss1-t2(j,l,a,c)*R(k,d,l,c)
				ss2=ss2-t2(i,l,a,c)*R(k,d,l,c)
			enddo
			enddo
			sum2=sum2+t1(i,d)*ss1-t1(j,d)*ss2
		enddo

		sum3=0
		do d = Nel+1,No
		do c = Nel+1,No
			sum3=sum3+t1(i,d)*t1(j,c)*R(k,d,a,c)
		enddo
		enddo

		sum4=0
		do l = 1,Nel
		do d = Nel+1,No
		do c = Nel+1,No
			sum4=sum4+t1(l,d)*t2(i,j,a,c)*R(k,d,l,c)
		enddo
		enddo
		enddo

		sum5=0
		do c = Nel+1,No
			sum5=sum5+t2(i,j,a,c)*F(k,c)
		enddo

		sum6=0
		do c = Nel+1,No
		do d = Nel+1,No
			sum6=sum6+t2(i,j,c,d)*R(k,c,a,d)
		enddo
		enddo

		sum7=0
		do l = 1,Nel
		do c = Nel+1,No
			sum7=sum7+t2(i,l,a,c)*R(j,k,c,l)
			sum7=sum7-t2(j,l,a,c)*R(i,k,c,l)
		enddo
		enddo

		ai2(i,j,k,a)=-R(i,k,j,a)+2*sum1+sum2-sum3-sum4+sum5-sum6/2+sum7
	enddo
	enddo
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  7/10'
	!$omp parallel default(shared) private(i,j,a,b,l,d,sum1,sum2,sum3)
	!$omp do
	do i = 1,Nel
	do j = 1,Nel
	do a = Nel+1,No
	do b = Nel+1,No

		sum1=0
		do l = 1,Nel
			sum1=sum1+t1(l,b)*R(j,i,l,a)
		enddo

		sum2=0
		do d = Nel+1,No
			sum2=sum2+t1(i,d)*R(j,a,b,d)
		enddo

		sum3=0
		do l = 1,Nel
		do d = Nel+1,No
			sum3=sum3+(t1(i,d)*t1(l,b)-t2(i,l,b,d))*R(j,a,l,d)
		enddo
		enddo

		ai1(i,a,j,b)=R(i,j,a,b)-sum1-sum2+sum3
	enddo
	enddo
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  8/10'
	!$omp parallel default(shared) private(a,b,c,l,k,sum1,sum2,sum3,sum4)
	!$omp do
	do a = Nel+1,No
	do b = Nel+1,No

		sum1=0
		do k = 1,Nel
			sum1=sum1+t1(k,a)*F(k,b)
		enddo
		
		sum2=0
		do k = 1,Nel
		do c = Nel+1,No
			sum2=sum2+t1(k,c)*R(k,b,a,c)
		enddo
		enddo
		
		sum3=0
		do k = 1,Nel
		do l = 1,Nel
		do c = Nel+1,No
			sum3=sum3+t1(k,c)*t1(l,a)*R(l,b,k,c)
		enddo
		enddo
		enddo
			
		sum4=0
		do k = 1,Nel
		do l = 1,Nel
		do c = Nel+1,No
			sum4=sum4+t2(l,k,a,c)*R(l,b,k,c)
		enddo
		enddo
		enddo

		Fab(a,b)=F(a,b)-sum1-sum2-sum3-sum4/2
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part  9/10'
	!$omp parallel default(shared) private(i,j,k,d,c,sum1,sum2,sum3,sum4)
	!$omp do
	do i = 1,Nel
	do j = 1,Nel

		sum1=0
		do c = Nel+1,No
			sum1=sum1+t1(i,c)*F(j,c)
		enddo
			
		sum2=0
		do k = 1,Nel
		do c = Nel+1,No
			sum2=sum2+t1(k,c)*R(j,i,k,c)
		enddo
		enddo
			
		sum3=0
		do k = 1,Nel
		do c = Nel+1,No
		do d = Nel+1,No
			sum3=sum3+t1(i,c)*t1(k,d)*R(j,c,k,d)
		enddo
		enddo
		enddo
			
		sum4=0
		do k = 1,Nel
		do d = Nel+1,No
		do c = Nel+1,No
			sum4=sum4+t2(i,k,d,c)*R(j,d,k,c)
		enddo
		enddo
		enddo

		Fij(i,j)=F(i,j)+sum1+sum2+sum3+sum4/2
	enddo
	enddo
	!$omp end parallel

	!write (lrsu,'(1X,A)') trim(ttime())//' Part 10/10'
	!$omp parallel default(shared) private(i,a,k,d,sum1)
	!$omp do
	do i = 1,Nel
	do a = Nel+1,No

		sum1=0
		do k = 1,Nel
		do d = Nel+1,No
			sum1=sum1+t1(k,d)*R(i,a,k,d)
		enddo
		enddo

		Fia(i,a)=F(i,a)+sum1
	enddo
	enddo
	!$omp end parallel		

	return
	end subroutine projection_lrccsd_spin_hf_intermediates_1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine projection_lrccsd_spin_hf_intermediates_2

	use glob        , only: iglu,rglu
	use lrccsdModule, only: Nel,No
	use lrccsdModule, only: R,F,t2,r1,r2
	use lrccsdModule, only: at1,at2,at3,at4,at5,at6,ai6,ai7

	implicit none

	integer(kind=iglu) :: i,j,k,l,a,b,c,d
	real   (kind=rglu) :: sum


	!$omp parallel default(shared) private(k,c,i,j,sum)
	!$omp do
	do i = 1,Nel
	do j = 1,Nel

		sum=0
		do k = 1,Nel
		do c = Nel+1,No
			sum=sum+r1(k,c)*ai6(j,k,i,c)
		enddo
		enddo

		at1(i,j)=sum
	enddo
	enddo
	!$omp end parallel
		
	!$omp parallel default(shared) private(k,c,a,b,sum)
	!$omp do
	do a = Nel+1,No
	do b = Nel+1,No

		sum=0
		do k = 1,Nel
		do c = Nel+1,No
			sum=sum+r1(k,c)*ai7(k,a,b,c)
		enddo
		enddo

		at2(a,b)=sum
	enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(k,c,i,j,d,sum)
	!$omp do
	do i = 1,Nel
	do j = 1,Nel

		sum=0
		do k = 1,Nel
		do c = Nel+1,No
		do d = Nel+1,No
			sum=sum+r2(i,k,c,d)*R(j,c,k,d)
		enddo
		enddo
		enddo

		at3(i,j)=sum/2
	enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(a,b,k,l,c,sum)
	!$omp do
	do a = Nel+1,No
	do b = Nel+1,No

		sum=0
		do k = 1,Nel
		do l = 1,Nel
		do c = Nel+1,No
			sum=sum+r2(k,l,a,c)*R(k,b,l,c)
		enddo
		enddo
		enddo

		at4(a,b)=sum/2
	enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(i,j,k,a,b,sum)
	!$omp do
	do i = 1,Nel
	do j = 1,Nel

		sum=0
		do k = 1,Nel
		do a = Nel+1,No
		do b = Nel+1,No
			sum=sum+r2(i,k,a,b)*t2(j,k,a,b)
		enddo
		enddo
		enddo

		at5(i,j)=sum/2
	enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(a,b,i,j,c,sum)
	!$omp do
	do a = Nel+1,No
	do b = Nel+1,No

		sum=0
		do i = 1,Nel
		do j = 1,Nel
		do c = Nel+1,No
			sum=sum+r2(i,j,a,c)*t2(i,j,b,c)
		enddo
		enddo
		enddo

		at6(a,b)=sum/2
	enddo
	enddo
	!$omp end parallel

	return
	end subroutine projection_lrccsd_spin_hf_intermediates_2