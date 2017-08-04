	subroutine projection_ccsd_singles_spin_cue_spare

	use coupledCluster      , only: Nel,No,Ne,R=>spin_cue_int,iapairs
	use coupledClusterSparse

	implicit none

	integer(kind=iglu) :: swit
	integer(kind=iglu) :: i,j,a,b,k,c,mm,nn,vv,sta1,sto1,sta2,sto2,vv2,vv3
	real   (kind=rglu) :: Ax,Bx,woint,wint,sss


	d1=0
	do mm = 1,Ne
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		b=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			c=Indexs(nn,2)

			Ax=vt(vv)/2

			d1(i,b)=d1(i,b)+R(b,b,j,c)*Ax
			d1(i,c)=d1(i,c)+R(c,b,j,c)*Ax
		enddo
	enddo

	do mm = 1,Ne
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		j=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			k=Indexs(nn,1)
			b=Indexs(nn,2)

			Ax=-vt(vv)/2

			d1(j,a)=d1(j,a)+R(j,j,k,b)*Ax
			d1(k,a)=d1(k,a)+R(j,k,k,b)*Ax
		enddo
	enddo

	do mm = 1,Ne
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			b=Indexs(nn,2)

			Ax=vt(vv)

			d1(i,a)=d1(i,a)+R(j,b,iapairs(b),iapairs(j))*Ax*t1(iapairs(b),iapairs(j))
			if (j.EQ.iapairs(b)) then
				do k = 1,Nel
					d1(i,a)=d1(i,a)+R(j,b,k,iapairs(k))*Ax*t1(k,iapairs(k))
				enddo
			endif
		enddo
	enddo

	do mm = 1,Ne
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		j=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			k=Indexs(nn,1)
			c=Indexs(nn,2)

			Ax=-vt(vv)/2

			b=iapairs(j)
			Bx=Ax*R(j,b,k,c)
			do i = 1,Nel
				d1(i,a)=d1(i,a)+Bx*t1(i,b)
			enddo
			b=iapairs(k)
			Bx=Ax*R(j,b,k,c)
			do i = 1,Nel
				d1(i,a)=d1(i,a)+Bx*t1(i,b)
			enddo
		enddo
	enddo

	do mm = 1,Ne
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		j=Indexs(mm,1)
		b=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			i=Indexs(nn,1)
			c=Indexs(nn,2)

			Ax=-vt(vv)/2

			k=iapairs(c)
			Bx=Ax*R(j,b,k,c)
			do a = Nel+1,No
				d1(i,a)=d1(i,a)+Bx*t1(k,a)
			enddo
			k=iapairs(b)
			Bx=Ax*R(j,b,k,c)
			do a = Nel+1,No
				d1(i,a)=d1(i,a)+Bx*t1(k,a)
			enddo
		enddo
	enddo

	do i = 1,Nel; do a = Nel+1,No
		d1(i,a)=d1(i,a)/4
	enddo; enddo

	do i = 1,Nel ! F(i,a)
		sta1=whOVf(i)
		sto1=ferow(i+1)-1
		do vv = sta1,sto1
			d1(i,fnumcol(vv))=d1(i,fnumcol(vv))+vF(vv)
		enddo
	enddo

	!$omp parallel default(shared) private(i,a,wint,woint,sss,sta1,sto1,vv,j,b,vv2,c,vv3,k,Ax)
	!$omp do
	do i = 1,Nel ! only t1
	do a = Nel+1,No

		wint=0; woint=0

		sss=0 ! F(a,b)*t1(i,b)
		sta1=whOVf(a)
		sto1=ferow(a+1)-1
		do vv = sta1,sto1
			sss=sss+vF(vv)*t1(i,fnumcol(vv))
		enddo
		woint=woint+sss

		sss=0 ! -F(i,j)*t1(j,a)
		sta1=ferow(i)
		sto1=whOVf(i)-1
		do vv = sta1,sto1
			sss=sss+vF(vv)*t1(fnumcol(vv),a)
		enddo
		woint=woint-sss

		sss=0
		sss=sss+R(i,i,a,a)*t1(i,a) ! -[ji||ab]*t1(j,b)
		if (iapairs(i).EQ.a) then
			do j = 1,i-1
				sss=sss+R(j,i,a,iapairs(j))*t1(j,iapairs(j))
			enddo
			do j = i+1,Nel
				sss=sss+R(j,i,a,iapairs(j))*t1(j,iapairs(j))
			enddo
		endif
		wint=wint-sss

		sss=0
		do j = 1,Nel ! -F(j,b)*t1(i,b)*t1(j,a)
			sta1=whOVf(j)
			sto1=ferow(j+1)-1
			do vv = sta1,sto1
				b=fnumcol(vv)
				sss=sss+vF(vv)*t1(i,b)*t1(j,a)
			enddo
		enddo
		woint=woint-sss

		sss=0
		do j = 1,Nel ! [ab||jc]*t1(i,b)*t1(j,c)
			sss=sss+R(a,a,j,iapairs(j))*t1(i,a)*t1(j,iapairs(j))
			sss=sss+R(a,iapairs(j),j,a)*t1(i,iapairs(j))*t1(j,a)
		enddo
		wint=wint+sss

		sss=0
		do j = 1,Nel ! -[ji||kb]*t1(j,a)*t1(k,b)
			sss=sss+R(j,i,i,iapairs(j))*t1(j,a)*t1(i,iapairs(j))
			sss=sss+R(i,i,j,iapairs(j))*t1(i,a)*t1(j,iapairs(j))
		enddo
		wint=wint-sss

		sss=0 ! -[jb||kc]*t1(j,b)*t1(i,c)*t1(k,a)
		do vv2 = 1,t1mrEx(i,0)
			c=t1mrEx(i,vv2)
			do vv3 = 1,t1mrEx(a,0)
				k=t1mrEx(a,vv3)

				Ax=t1(i,c)*t1(k,a)

				if (k.EQ.iapairs(c)) then
					do j = 1,Nel
						b=iapairs(j)
						sss=sss+R(j,b,k,c)*t1(j,b)*Ax
					enddo
				else
					j=iapairs(c)
					b=iapairs(k)
					sss=sss+R(j,b,k,c)*t1(j,b)*Ax
				endif
			enddo
		enddo
		wint=wint-sss

		d1(i,a)=d1(i,a)+woint+wint/4
	enddo
	enddo
	!$omp end parallel

	do j = 1,Nel
		sta1=whOVf(j)
		sto1=ferow(j+1)-1

		do vv = sta1,sto1
			b=fnumcol(vv)
			Ax=vF(vv)
			mm=cIndex(j,b)
			sta2=erow(mm)
			sto2=erow(mm+1)-1

			do vv2 = sta2,sto2
				nn=numcol(vv2)

				i=Indexs(nn,1)
				a=Indexs(nn,2)

				d1(i,a)=d1(i,a)+Ax*vt(vv2)
			enddo
		enddo
	enddo

	return
	end subroutine projection_ccsd_singles_spin_cue_spare