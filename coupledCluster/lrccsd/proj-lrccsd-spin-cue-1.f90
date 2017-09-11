	subroutine projection_lrccsd_singles_spin_cue

	use glob        , only: iglu,rglu
	use lrccsdModule, only: Fab,Fij,ai1,Fia,ai6,ai7
	use lrccsdModule, only: Nel,No,iapairs
	use lrccsdModule, only: r1,r2,d1
	use lrccsdModule, only: t1,t2,F,R

	implicit none

	integer(kind=iglu) :: i,j,k,l,a,b,c,d
	real   (kind=rglu) :: rez,sum,sum1,sum2,sum3,sum4


	d1=0
	!$omp parallel default(shared) private(i,a,j,b,k,c,rez,sum)
	!$omp do
	do i = 1,Nel
	do a = Nel+1,No
		if (btest(a,0).NE.btest(i,0)) cycle	

		rez=0

		sum=0
		do k = 1,Nel
			b=iapairs(k); c=a
			sum=sum+R(k,b,a,c)*(t1(k,c)*r1(i,b)-t1(k,b)*r1(i,c))
		enddo
		rez=rez-sum

		sum=0
		do b = Nel+1,No
			sum1=0
			do c = Nel+1,No
				l=iapairs(b); k=iapairs(c)
				sum1=sum1+R(b,l,c,k)*(t1(l,a)*t1(k,c)-t1(k,a)*t1(l,c)+t2(l,k,a,c))
			enddo

			sum2=0
			do k = 1,Nel
				sum2=sum2+t1(k,a)*F(k,b)
			enddo

			sum=sum+(-F(a,b)+sum1+sum2)*r1(i,b)
		enddo
		rez=rez-sum

		sum=0
		do b = Nel+1,No
			j=i; k=iapairs(b)
			sum=sum+R(i,j,b,k)*(t1(k,b)*r1(j,a)-t1(j,b)*r1(k,a))
		enddo
		rez=rez-sum

		sum=0
		do j = 1,Nel
			sum1=0
			do k = 1,Nel
				b=iapairs(j); c=iapairs(k)
				sum1=sum1+R(j,b,k,c)*(t1(i,b)*t1(k,c)-t1(i,c)*t1(k,b)-t2(i,k,c,b))
			enddo

			sum2=0
			do b = Nel+1,No
				sum2=sum2+t1(i,b)*F(j,b)
			enddo

			sum=sum+(F(i,j)+sum1+sum2)*r1(j,a)
		enddo
		rez=rez-sum

		! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		sum=0
		do j = 1,Nel
		do b = Nel+1,No
			sum=sum-ai1(i,b,j,a)*r1(j,b) 
		enddo
		enddo
		rez=rez+sum

		sum=0
		do j = 1,Nel
		do b = Nel+1,No
			sum=sum+Fia(j,b)*r2(i,j,a,b)
		enddo
		enddo
		rez=rez+sum

		sum=0
		do j = 1,Nel
		do k = 1,Nel
		do b = Nel+1,No
			sum=sum+ai6(j,k,i,b)*r2(j,k,a,b)
		enddo
		enddo
		enddo
		rez=rez-sum/2

		sum=0
		do j = 1,Nel
		do b = Nel+1,No
		do c = Nel+1,No
			sum=sum+ai7(j,a,b,c)*r2(i,j,b,c)
		enddo
		enddo
		enddo
		rez=rez-sum/2

		d1(i,a)=rez
	enddo
	enddo
	!$omp end parallel


	return
	end subroutine projection_lrccsd_singles_spin_cue