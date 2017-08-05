	subroutine projection_ccsd_doubles_spin_cue_spare

	use glob                , only: getThreadNumber
	use coupledCluster      , only: Nel,No,Ne,Nth,R=>spin_cue_int,iapairs
	use coupledClusterSparse

	implicit none

	integer(kind=iglu) :: ct,swit
	integer(kind=iglu) :: i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,sta2,sto2,sta3,sto3,mm2,nn2,vv2,vv3,vv4,ccc
	real   (kind=rglu) :: Ax,Bx,Cx,Dx,sum

	integer(kind=iglu) :: pi1,pi2

	integer(kind=iglu) :: intero,intere


	vd=0; pvd=0
	!$omp parallel default(shared) private(ct,mm,i,a,nn,j,b)
	!ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! <ij||ab> = [ia||jb] !done
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		j=iapairs(a); b=iapairs(i)

		nn=cIndex(j,b)
		vd(ftfm(mm,nn))=R(i,a,j,b)

		if (b.EQ.a) then
			do j = 1,Nel
				b=iapairs(j)
				nn=cIndex(j,b)
				vd(ftfm(mm,nn))=R(i,a,j,b)
			enddo
		endif
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,mm,sta1,sto1,i,a,vv,nn,j,b,vv2,k,c,Ax,vv3,xx1,xx2)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij|ab) [bj||kc]*t(i,k,a,c) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		a=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			b=Indexs(nn,2)

			if (j.EQ.iapairs(b)) then
				do vv2 = 1,occEx(mm,0)
					k=occEx(mm,vv2); c=iapairs(k)

					Ax=R(b,j,k,c)*vt( ftfm(cIndex(k,c),mm) )

					vv3=ftfm(mm,nn); pvd(vv3,ct)=pvd(vv3,ct)+Ax
					vv3=ftfm(nn,mm); pvd(vv3,ct)=pvd(vv3,ct)+Ax

					xx1=cIndex(i,b)
					xx2=cIndex(j,a)
					vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Ax
					vv3=ftfm(xx2,xx1); pvd(vv3,ct)=pvd(vv3,ct)-Ax
				enddo
			else
				k=j; c=b

				Ax=R(b,j,k,c)*vt( ftfm(cIndex(k,c),mm) )

				vv3=ftfm(mm,nn); pvd(vv3,ct)=pvd(vv3,ct)+Ax
				vv3=ftfm(nn,mm); pvd(vv3,ct)=pvd(vv3,ct)+Ax

				xx1=cIndex(i,b)
				xx2=cIndex(j,a)
				vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Ax
				vv3=ftfm(xx2,xx1); pvd(vv3,ct)=pvd(vv3,ct)-Ax
			endif
		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,mm,sta1,sto1,k,a,vv,nn,l,b)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! -0.5[ik||jl]*t(k,l,b,a) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1
		
		k=Indexs(mm,1)
		a=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)
			
			l=Indexs(nn,1)
			b=Indexs(nn,2)
			
			if (k.NE.l) pvd(vv,ct)=pvd(vv,ct)+vt(vv)*R(k,k,l,l) !k=i, l=j
			! simplification due to [kk||ll]-[lk||kl]=2[kk||ll]: _old 0.5d0*vt(vv)*(R(k,k,l,l)-R(l,k,k,l))
		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,mm,sta1,sto1,i,c,vv,nn,j,d)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! 0.5*[ac||bd]*t(i,j,c,d) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1
		
		i=Indexs(mm,1)
		c=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)
			
			j=Indexs(nn,1)
			d=Indexs(nn,2)

			if (c.NE.d) pvd(vv,ct)=pvd(vv,ct)+vt(vv)*R(c,c,d,d) ! _old 0.5d0*vt(vv)*(R(c,c,d,d)-R(d,c,c,d))
		enddo
	enddo				
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,k,l,a,b,c,d,mm,nn,sta1,sto1,sta2,sto2,sta3,sto3,vv,vv2,vv3,vv4,Ax,Bx,intero,pi1,pi2,xx1,xx2)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij) -0.5[kc||ld]*t(l,j,a,b)*t(k,i,c,d) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		l=Indexs(mm,1); d=iapairs(l)
		a=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			b=Indexs(nn,2)

			do vv2 = 1,mrEx(d,0)
				k=mrEx(d,vv2); c=iapairs(k)

				Ax=vt(vv)*R(k,c,l,d)

				sta2=byExOrb(nn,iapairs(a))         ; sto2=byExOrb(nn,iapairs(a)+1)-1
				sta3=byExOrb(cIndex(k,c),iapairs(d)); sto3=byExOrb(cIndex(k,c),iapairs(d)+1)-1

				intero=0
				do vv3 = sta2,sto2; pi1=Indexsbv(numcolbv(vv3),1)
				do vv4 = sta3,sto3; pi2=Indexsbv(numcolbv(vv4),1)
					if (pi1.EQ.pi2) then
						intero=intero+1; intersectOrbitals(intero,ct)=pi1
						exit
					endif
				enddo
				enddo

				do vv3 = 1,intero
					i=intersectOrbitals(vv3,ct)

					Bx=Ax*vt( ftfm(cIndex(k,c),cIndex(i,d)) )

					xx1=cIndex(i,a)
					xx2=cIndex(j,b)
					vv4=ftfm(xx1,xx2); pvd(vv4,ct)=pvd(vv4,ct)-Bx   !!!

					xx1=cIndex(j,a)
					xx2=cIndex(i,b)
					vv4=ftfm(xx1,xx2); pvd(vv4,ct)=pvd(vv4,ct)+Bx   !!!
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel
	
	!$omp parallel default(shared) private(ct,i,j,k,l,a,b,c,d,mm,nn,sta1,sto1,sta2,sto2,sta3,sto3,vv,vv2,vv3,vv4,Ax,Bx,intero,pi1,pi2,xx1,xx2)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ab) -0.5[kc||ld]*t(l,k,d,b)*t(j,i,c,a) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		a=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			c=Indexs(nn,2); k=iapairs(c)

			do vv2 = 1,mrEx(k,0)
				d=mrEx(k,vv2); l=iapairs(d)

				Ax=vt(vv)*R(k,c,l,d)

				sta2=byExOrb(mm,j)         ; sto2=byExOrb(mm,j+1)-1
				sta3=byExOrb(cIndex(l,d),k); sto3=byExOrb(cIndex(l,d),k+1)-1

				intero=0
				do vv3 = sta2,sto2; pi1=Indexs(numcol(vv3),2)
				do vv4 = sta3,sto3; pi2=Indexs(numcol(vv4),2)
					if (pi1.EQ.pi2) then
						intero=intero+1; intersectOrbitals(intero,ct)=pi1
						exit
					endif
				enddo
				enddo

				do vv3 = 1,intero
					b=intersectOrbitals(vv3,ct)

					Bx=Ax*vt( ftfm(cIndex(k,b),cIndex(l,d)) )

					xx1=cIndex(i,a)
					xx2=cIndex(j,b)
					vv4=ftfm(xx1,xx2); pvd(vv4,ct)=pvd(vv4,ct)-Bx

					xx1=cIndex(i,b)
					xx2=cIndex(j,a)
					vv4=ftfm(xx1,xx2); pvd(vv4,ct)=pvd(vv4,ct)+Bx
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,sta2,sto2,mm2,nn2,vv2,vv3,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij) -[kc||ld]*t(l,j,a,b)*t(k,c)*t(i,d)
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		j=Indexs(mm,1)
		b=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			l=Indexs(nn,1)
			a=Indexs(nn,2)

			Ax=vt(vv)

			do k = 1,Nel
				if (k.EQ.l) cycle

				d=iapairs(l); c=iapairs(k)
				Bx=Ax*R(k,c,l,d)

				sta2=byExOrbbv(mm,a)
				sto2=byExOrbbv(mm,a+1)-1
				do vv2 = sta2,sto2
					i=Indexsbv(numcolbv(vv2),1)

					Cx=Bx*(t1(i,d)*t1(k,c)-t1(i,c)*t1(k,d))

					xx1=cIndex(i,a)
					vv3=ftfm(xx1,mm); pvd(vv3,ct)=pvd(vv3,ct)-Cx

					xx1=cIndex(j,a)
					xx2=cIndex(i,b)
					vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Cx
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,sta1,sto1,sta2,sto2,mm2,nn2,vv2,vv3,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ab) -[kc||ld]*t(i,j,a,c)*t(l,d)*t(k,b) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			c=Indexs(nn,2)

			Ax=vt(vv)

			do d = Nel+1,No
				if (c.EQ.d) cycle

				k=iapairs(c); l=iapairs(d)
				Bx=Ax*R(k,c,l,d)

				sta2=byExOrb(mm,j)
				sto2=byExOrb(mm,j+1)-1
				do vv2 = sta2,sto2
					b=Indexs(numcol(vv2),2)

					Cx=Bx*(t1(k,b)*t1(l,d)-t1(l,b)*t1(k,d))

					vv3=ftfm(mm,cIndex(j,b)); pvd(vv3,ct)=pvd(vv3,ct)-Cx
					vv3=ftfm(cIndex(i,b),cIndex(j,a)); pvd(vv3,ct)=pvd(vv3,ct)+Cx
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,vv2,vv3,vv4,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij) -[kc||ld]*t(j,k,a,c)*t(l,i,d,b) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			b=Indexs(nn,2)

			do vv2 = 1,occEx(mm,0)
				k=occEx(mm,vv2); c=iapairs(k)

				do vv3 = 1,occEx(nn,0)
					l=occEx(nn,vv3); d=iapairs(l)

					if (k.EQ.l) cycle

					Ax= vt( ftfm( mm,cIndex(k,c) ) )*vt( ftfm( cIndex(l,d),nn ) )
					Bx=-vt( ftfm( mm,cIndex(k,d) ) )*vt( ftfm( cIndex(l,c),nn ) )

					Cx=R(k,c,l,d)*(Ax+Bx)

					pvd(vv,ct)=pvd(vv,ct)+Cx

					xx1=cIndex(j,a)
					xx2=cIndex(i,b)
					vv4=ftfm(xx1,xx2); pvd(vv4,ct)=pvd(vv4,ct)-Cx
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,sta1,sto1,sta2,sto2,vv2,vv3,vv4,Ax,intero,intere)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! 0.25[kc||ld]*t(i,j,c,d)*t(k,l,a,b) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		k=Indexs(mm,1); c=iapairs(k)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			l=Indexs(nn,1); d=iapairs(l)
			b=Indexs(nn,2)

			if (k.EQ.l) cycle

			Ax=R(k,c,l,d)*vt(vv)/2

			intero=0
			do vv2 = 1,occEx(cIndex(iapairs(c),d),0); pi1=occEx(cIndex(iapairs(c),d),vv2)
			do vv3 = 1,occEx(cIndex(iapairs(a),b),0); pi2=occEx(cIndex(iapairs(a),b),vv3)
				if (pi1.EQ.pi2) then
					intero=intero+1; intersectOrbitals(intero,ct)=pi1
					exit
				endif
			enddo
			enddo

			do vv2 = 1,intero
				i=intersectOrbitals(vv2,ct)

				do vv3 = 1,intero
					j=intersectOrbitals(vv3,ct)

					vv4=ftfm(cIndex(i,a),cIndex(j,b))

					Cx=vt( ftfm(cIndex(i,c),cIndex(j,d)) )
					pvd(vv4,ct)=pvd(vv4,ct)+Ax*Cx   !!!
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,vv2,vv3,vv4,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij|ab) [kc||bd]*t(k,i,c,a)*t(j,d) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		k=Indexs(mm,1)
		c=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			i=Indexs(nn,1)
			a=Indexs(nn,2)

			Ax=vt(vv)

			if (k.EQ.iapairs(c)) then
				do vv4 = 1,mrEx(i,0)
					b=mrEx(i,vv4); d=b
					Bx=Ax*R(k,c,b,d)

					do vv3 = byExOrbbv(nn,b),byExOrbbv(nn,b+1)-1
						j=Indexsbv(numcolbv(vv3),1)

						Cx=Bx*t1(j,d)

						xx2=cIndex(j,b)
						vv2=ftfm(nn,xx2); pvd(vv2,ct)=pvd(vv2,ct)+Cx
						vv2=ftfm(xx2,nn); pvd(vv2,ct)=pvd(vv2,ct)+Cx

						xx1=cIndex(j,a)
						xx2=cIndex(i,b)
						vv2=ftfm(xx1,xx2); pvd(vv2,ct)=pvd(vv2,ct)-Cx
						vv2=ftfm(xx2,xx1); pvd(vv2,ct)=pvd(vv2,ct)-Cx
					enddo
				enddo
			else
				b=c; d=iapairs(k)
				Bx=Ax*R(k,c,b,d)

				do vv3 = byExOrbbv(nn,b),byExOrbbv(nn,b+1)-1
					j=Indexsbv(numcolbv(vv3),1)

					Cx=Bx*t1(j,d)

					xx1=cIndex(i,a)
					xx2=cIndex(j,b)
					vv2=ftfm(nn,xx2); pvd(vv2,ct)=pvd(vv2,ct)+Cx
					vv2=ftfm(xx2,nn); pvd(vv2,ct)=pvd(vv2,ct)+Cx

					xx1=cIndex(j,a)
					xx2=cIndex(i,b)
					vv2=ftfm(xx1,xx2); pvd(vv2,ct)=pvd(vv2,ct)-Cx
					vv2=ftfm(xx2,xx1); pvd(vv2,ct)=pvd(vv2,ct)-Cx
				enddo
			endif
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,vv2,vv3,vv4,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij|ab) -[kc||lj]*t(k,i,c,a)*t(l,b) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		k=Indexs(mm,1)
		c=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			i=Indexs(nn,1)
			a=Indexs(nn,2)

			Ax=vt(vv)

			if (k.EQ.iapairs(c)) then
				do vv4 = 1,mrEx(a,0)
					j=mrEx(a,vv4); l=j
					Bx=Ax*R(k,c,l,j)

					do vv3 = byExOrb(nn,j),byExOrb(nn,j+1)-1
						b=Indexs(numcol(vv3),2)

						Cx=Bx*t1(l,b)

						xx2=cIndex(j,b)
						vv2=ftfm(nn,xx2); pvd(vv2,ct)=pvd(vv2,ct)-Cx
						vv2=ftfm(xx2,nn); pvd(vv2,ct)=pvd(vv2,ct)-Cx

						xx1=cIndex(j,a)
						xx2=cIndex(i,b)
						vv2=ftfm(xx1,xx2); pvd(vv2,ct)=pvd(vv2,ct)+Cx
						vv2=ftfm(xx2,xx1); pvd(vv2,ct)=pvd(vv2,ct)+Cx
					enddo
				enddo
			else
				l=iapairs(c); j=k
				Bx=Ax*R(k,c,l,j)

				do vv3 = byExOrb(nn,j),byExOrb(nn,j+1)-1
					b=Indexs(numcol(vv3),2)
					Cx=Bx*t1(l,b)

					xx2=cIndex(j,b)
					vv2=ftfm(nn,xx2); pvd(vv2,ct)=pvd(vv2,ct)-Cx
					vv2=ftfm(xx2,nn); pvd(vv2,ct)=pvd(vv2,ct)-Cx

					xx1=cIndex(j,a)
					xx2=cIndex(i,b)
					vv2=ftfm(xx1,xx2); pvd(vv2,ct)=pvd(vv2,ct)+Cx
					vv2=ftfm(xx2,xx1); pvd(vv2,ct)=pvd(vv2,ct)+Cx
				enddo
			endif
		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,vv2,vv3,Ax,Bx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij) -[kj||lc]*t(k,i,b,a)*t(l,c) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			k=Indexs(nn,1)
			b=Indexs(nn,2)

			Ax=vt(vv)
			do vv2 = byExOrbbv(mm,b),byExOrbbv(mm,b+1)-1
				j=Indexsbv(numcolbv(vv2),1)

				if (j.EQ.k) then
					do l = 1,Nel
						c=iapairs(l)

						Bx=Ax*R(k,j,l,c)*t1(l,c)

						xx2=cIndex(j,b)
						vv3=ftfm(mm,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Bx

						xx1=cIndex(j,a)
						xx2=cIndex(i,b)
						vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Bx
					enddo
				else
					l=j; c=iapairs(k)

					Bx=Ax*R(k,j,l,c)*t1(l,c)

					xx2=cIndex(j,b)
					vv3=ftfm(mm,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Bx

					xx1=cIndex(j,a)
					xx2=cIndex(i,b)
					vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Bx
				endif
			enddo
		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,vv2,vv3,Ax,Bx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ab) [cb||dk]*t(j,i,c,a)*t(k,d) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			c=Indexs(nn,2)

			Ax=vt(vv)
			do vv2 = byExOrb(mm,j),byExOrb(mm,j+1)-1
				b=Indexs(numcol(vv2),2)

				if (c.EQ.b) then
					do d = Nel+1,No
						k=iapairs(d)

						Bx=Ax*R(c,b,d,k)*t1(k,d)

						xx2=cIndex(j,b)
						vv3=ftfm(mm,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Bx

						xx1=cIndex(i,b)
						xx2=cIndex(j,a)
						vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Bx
					enddo
				else
					d=b; k=iapairs(c)

					Bx=Ax*R(c,b,d,k)*t1(k,d)

					xx2=cIndex(j,b)
					vv3=ftfm(mm,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Bx

					xx1=cIndex(i,b)
					xx2=cIndex(j,a)
					vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Bx
				endif
			enddo
		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,mm2,vv2,vv3,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij) 0.5[kc||lj]*t(k,l,a,b)*t(i,c) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		k=Indexs(mm,1)
		a=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			l=Indexs(nn,1)
			b=Indexs(nn,2)

			if (k.EQ.l) cycle

			Ax=vt(vv)/2

			c=iapairs(k); j=l
			Bx=R(k,c,l,j)*Ax
			mm2=cIndex(j,b)
			do vv2 = byExOrbbv(mm2,a),byExOrbbv(mm2,a+1)-1
				i=Indexsbv(numcolbv(vv2),1)
				Cx=Bx*t1(i,c)

				vv3=ftfm(cIndex(i,a),mm2); pvd(vv3,ct)=pvd(vv3,ct)+Cx

				xx1=cIndex(j,a)
				xx2=cIndex(i,b)
				vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Cx
			enddo

			c=iapairs(l); j=k
			Bx=R(k,c,l,j)*Ax
			mm2=cIndex(j,b)
			do vv2 = byExOrbbv(mm2,a),byExOrbbv(mm2,a+1)-1
				i=Indexsbv(numcolbv(vv2),1)
				Cx=Bx*t1(i,c)

				vv3=ftfm(cIndex(i,a),mm2); pvd(vv3,ct)=pvd(vv3,ct)+Cx

				xx1=cIndex(j,a)
				xx2=cIndex(i,b)
				vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Cx
			enddo

		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,mm2,vv2,vv3,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ab) -0.5[ca||dk]*t(i,j,c,d)*t(k,b) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		c=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			d=Indexs(nn,2)

			if (c.EQ.d) cycle

			Ax=vt(vv)/2

			a=c; k=iapairs(d)
			Bx=R(c,a,d,k)*Ax
			mm2=cIndex(i,a)
			do vv2 = byExOrb(mm2,j),byExOrb(mm2,j+1)-1
				b=Indexs(numcol(vv2),2)
				Cx=Bx*t1(k,b)

				vv3=ftfm(mm2,cIndex(j,b)); pvd(vv3,ct)=pvd(vv3,ct)-Cx

				xx1=cIndex(i,b)
				xx2=cIndex(j,a)
				vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Cx
			enddo

			a=d; k=iapairs(c)
			Bx=R(c,a,d,k)*Ax
			mm2=cIndex(i,a)
			do vv2 = byExOrb(mm2,j),byExOrb(mm2,j+1)-1
				b=Indexs(numcol(vv2),2)
				Cx=Bx*t1(k,b)

				vv3=ftfm(mm2,cIndex(j,b)); pvd(vv3,ct)=pvd(vv3,ct)-Cx

				xx1=cIndex(i,b)
				xx2=cIndex(j,a)
				vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Cx
			enddo

		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,xx1,xx2,sta1,sto1,nn2,vv2,vv3,Ax,Bx,Cx)
	ct=getThreadNumber()
	!$omp do
	do k = 1,Nel ! P(ij|ab) -[kc||ld]*t(l,j,d,b)*t(i,c)*t(k,a) !done
	do c = Nel+1,No

		if (k.EQ.iapairs(c)) then
			do l = 1,Nel
				d=iapairs(l)

				Ax=R(l,d,k,c)

				mm=cIndex(l,d)

				sta1=erow(mm)
				sto1=erow(mm+1)-1
				
				do vv = sta1,sto1
					nn=numcol(vv)

					j=Indexs(nn,1)
					b=Indexs(nn,2)

					Bx=Ax*vt(vv)

					do vv2 = erow(nn),erow(nn+1)-1
						nn2=numcol(vv2)

						i=Indexs(nn2,1)
						a=Indexs(nn2,2)

						Cx=Bx*t1(i,c)*t1(k,a)

						vv3=ftfm(nn2,nn); pvd(vv3,ct)=pvd(vv3,ct)-Cx
						vv3=ftfm(nn,nn2); pvd(vv3,ct)=pvd(vv3,ct)-Cx

						xx1=cIndex(j,a)
						xx2=cIndex(i,b)
						
						vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Cx
						vv3=ftfm(xx2,xx1); pvd(vv3,ct)=pvd(vv3,ct)+Cx
					enddo
				enddo
			enddo
		else
			l=iapairs(c); d=iapairs(k)

			Ax=R(l,d,k,c)

			mm=cIndex(l,d)

			sta1=erow(mm)
			sto1=erow(mm+1)-1

			do vv = sta1,sto1
				nn=numcol(vv)

				j=Indexs(nn,1)
				b=Indexs(nn,2)

				Bx=Ax*vt(vv)

				do vv2 = erow(nn),erow(nn+1)-1
					nn2=numcol(vv2)

					i=Indexs(nn2,1)
					a=Indexs(nn2,2)

					Cx=Bx*t1(i,c)*t1(k,a)

					vv3=ftfm(nn2,nn); pvd(vv3,ct)=pvd(vv3,ct)-Cx
					vv3=ftfm(nn,nn2); pvd(vv3,ct)=pvd(vv3,ct)-Cx

					xx1=cIndex(j,a)
					xx2=cIndex(i,b)
					
					vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Cx
					vv3=ftfm(xx2,xx1); pvd(vv3,ct)=pvd(vv3,ct)+Cx
				enddo
			enddo
		endif
	enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,sta1,sto1,mm2,vv2,vv3,vv4,Ax)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! 0.5[kc||ld]*t(k,l,a,b)*t(i,c)*t(j,b) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		k=Indexs(mm,1); c=iapairs(k)
		a=Indexs(mm,2)
		do vv = sta1,sto1
			nn=numcol(vv)

			l=Indexs(nn,1); d=iapairs(l)
			b=Indexs(nn,2)

			if (k.EQ.l) cycle

			Ax=R(k,c,l,d)*vt(vv)/2
			do vv2 = 1,mrEx(a,0)
				i=mrEx(a,vv2); mm2=cIndex(i,a)

				do vv3 = byExOrbbv(mm2,b),byExOrbbv(mm2,b+1)-1
					j=Indexsbv(numcolbv(vv3),1)

					vv4=ftfm(mm2,cIndex(j,b))
					pvd(vv4,ct)=pvd(vv4,ct)+Ax*(t1(i,c)*t1(j,d)-t1(i,d)*t1(j,c))
				enddo
			enddo

		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,i,j,a,b,k,l,c,d,mm,nn,vv,sta1,sto1,mm2,vv2,vv3,vv4,Ax)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! 0.5[kc||ld]*t(i,j,c,d)*t(k,a)*t(l,b) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1

		i=Indexs(mm,1)
		c=Indexs(mm,2); k=iapairs(c)
		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			d=Indexs(nn,2); l=iapairs(d)

			if (c.EQ.d) cycle

			Ax=R(k,c,l,d)*vt(vv)/2
			do vv2 = 1,mrEx(i,0)
				a=mrEx(i,vv2); mm2=cIndex(i,a)

				do vv3 = byExOrb(mm2,j),byExOrb(mm2,j+1)-1
					b=Indexs(numcol(vv3),2)

					vv4=ftfm(mm2,cIndex(j,b)); pvd(vv4,ct)=pvd(vv4,ct)+Ax*(t1(l,b)*t1(k,a)-t1(k,b)*t1(l,a))
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel	

	!$omp parallel default(shared) private(mm,sta1,sto1,i,a,vv,nn,j,b,sss,k,c,l,d,vv2,vv3)
	ct=getThreadNumber()
	!$omp do
	! Only t1
	do mm = 1,Ne
		sta1=erow(mm)
		sto1=erow(mm+1)-1
		
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			b=Indexs(nn,2)
			
			sum=0
			
			if (j.EQ.iapairs(b)) sum=sum+R(a,a,j,b)*t1(i,a)
			if (j.EQ.iapairs(a)) sum=sum+R(b,a,j,b)*t1(i,b)
			if (i.EQ.iapairs(b)) sum=sum-R(a,a,i,b)*t1(j,a)
			if (i.EQ.iapairs(a)) sum=sum-R(b,a,i,b)*t1(j,b)
			if (j.EQ.iapairs(b)) sum=sum-R(i,i,b,j)*t1(i,a)
			if (i.EQ.iapairs(b)) sum=sum-R(j,i,b,j)*t1(j,a)
			if (j.EQ.iapairs(a)) sum=sum+R(i,i,a,j)*t1(i,b)
			if (i.EQ.iapairs(a)) sum=sum+R(j,i,a,j)*t1(j,b)

			if (j.EQ.iapairs(a)) then
				do k = 1,Nel
					c=iapairs(k); sum=sum+R(a,j,k,c)*t1(i,c)*t1(k,b)
				enddo
			else
				k=j; c=a; sum=sum+R(a,j,k,c)*t1(i,c)*t1(k,b)
			endif

			if (j.EQ.iapairs(b)) then
				do k = 1,Nel
					c=iapairs(k); sum=sum-R(b,j,k,c)*t1(i,c)*t1(k,a)
				enddo
			else
				k=j; c=b; sum=sum-R(b,j,k,c)*t1(i,c)*t1(k,a)
			endif

			if (i.EQ.iapairs(a)) then
				do k = 1,Nel
					c=iapairs(k); sum=sum-R(a,i,k,c)*t1(j,c)*t1(k,b)
				enddo
			else
				k=i; c=a; sum=sum-R(a,i,k,c)*t1(j,c)*t1(k,b) 
			endif

			if (i.EQ.iapairs(b)) then
				do k = 1,Nel
					c=iapairs(k); sum=sum+R(b,i,k,c)*t1(j,c)*t1(k,a)
				enddo
			else
				k=i; c=b; sum=sum+R(b,i,k,c)*t1(j,c)*t1(k,a)
			endif

			k=i; l=j; sum=sum+R(i,k,j,l)*t1(k,a)*t1(l,b)
			k=j; l=i; sum=sum+R(i,k,j,l)*t1(k,a)*t1(l,b)
			c=a; d=b; sum=sum+R(a,c,b,d)*t1(i,c)*t1(j,d)
			c=b; d=a; sum=sum+R(a,c,b,d)*t1(i,c)*t1(j,d)

			do c = Nel+1,No
				k=iapairs(c); l=j; sum=sum+R(c,k,j,l)*t1(i,c)*t1(k,a)*t1(l,b)
				k=j; l=iapairs(c); sum=sum+R(c,k,j,l)*t1(i,c)*t1(k,a)*t1(l,b)
				k=iapairs(c); l=i; sum=sum-R(c,k,i,l)*t1(j,c)*t1(k,a)*t1(l,b)
				k=i; l=iapairs(c); sum=sum-R(c,k,i,l)*t1(j,c)*t1(k,a)*t1(l,b)
			enddo

			do k = 1,Nel
				c=a; d=iapairs(k); sum=sum-R(a,c,k,d)*t1(i,c)*t1(j,d)*t1(k,b)
				c=iapairs(k); d=a; sum=sum-R(a,c,k,d)*t1(i,c)*t1(j,d)*t1(k,b)
				c=b; d=iapairs(k); sum=sum+R(b,c,k,d)*t1(i,c)*t1(j,d)*t1(k,a)
				c=iapairs(k); d=b; sum=sum+R(b,c,k,d)*t1(i,c)*t1(j,d)*t1(k,a)
			enddo

			do vv2 = 1,t1mrEx(a,0) ! [kc||ld]*t(k,a)*t(l,b)*t(i,c)*t(j,d)
				k=t1mrEx(a,vv2)

				do vv3 = 1,t1mrEx(b,0)
					l=t1mrEx(b,vv3)

					if (k.NE.l) then
						c=iapairs(k); d=iapairs(l); sum=sum+R(k,c,l,d)*t1(k,a)*t1(l,b)*t1(i,c)*t1(j,d)
						c=iapairs(l); d=iapairs(k); sum=sum+R(k,c,l,d)*t1(k,a)*t1(l,b)*t1(i,c)*t1(j,d)
					endif
				enddo
			enddo
			
			vd(vv)=vd(vv)+sum
		enddo
	enddo
	!$omp end parallel

	! collect threads
	!$omp parallel default(shared) private(vv,vv2)
	!$omp do
	do vv = 1,Nue
		do vv2 = 1,Nth
			vd(vv)=vd(vv)+pvd(vv,vv2); pvd(vv,vv2)=0
		enddo
		vd(vv)=vd(vv)/4 ! 1/4 for integrals
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,mm,sta1,sto1,i,a,vv,nn,j,c,Ax,vv2,Bx,b,xx1,xx2,vv3)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ab) F(b,c)*t(i,j,a,c) !done
		sta1=erow(mm)
		sto1=erow(mm+1)-1
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			j=Indexs(nn,1)
			c=Indexs(nn,2)

			Ax=vt(vv)

			do vv2 = whOVf(c),ferow(c+1)-1
				b=fnumcol(vv2)

				Bx=vF(vv2)*Ax

				vv3=ftfm(mm,cIndex(j,b)); pvd(vv3,ct)=pvd(vv3,ct)+Bx

				xx1=cIndex(i,b)
				xx2=cIndex(j,a)
				vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)-Bx
			enddo
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,mm,sta1,sto1,i,a,vv,nn,k,b,Ax,vv2,Bx,j,xx1,xx2,vv3)
	ct=getThreadNumber()
	!$omp do
	do mm = 1,Ne ! P(ij) -F(k,j)*t(i,k,a,b)
		sta1=erow(mm)
		sto1=erow(mm+1)-1
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do vv = sta1,sto1
			nn=numcol(vv)

			k=Indexs(nn,1)
			b=Indexs(nn,2)

			Ax=vt(vv)

			do vv2 = ferow(k),whOVf(k)-1
				j=fnumcol(vv2)

				Bx=Ax*vF(vv2)

				vv3=ftfm(mm,cIndex(j,b)); pvd(vv3,ct)=pvd(vv3,ct)-Bx

				xx1=cIndex(j,a)
				xx2=cIndex(i,b)
				vv3=ftfm(xx1,xx2); pvd(vv3,ct)=pvd(vv3,ct)+Bx
			enddo
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,k,vv,c,Ax,a,mm2,vv2,nn2,j,b,Bx,vv3,i,Cx,xx1,xx2,vv4)
	ct=getThreadNumber()
	!$omp do
	do k = 1,Nel ! P(ij) -F(k,c)*t(k,j,a,b)*t(i,c)
		do vv = whOVf(k),ferow(k+1)-1
			c=fnumcol(vv)

			Ax=vF(vv)

			do a = Nel+1,No
				mm2=cIndex(k,a)

				do vv2 = erow(mm2),erow(mm2+1)-1
					nn2=numcol(vv2)

					j=Indexs(nn2,1)
					b=Indexs(nn2,2)

					Bx=Ax*vt(vv2)

					do vv3 = byExOrbbv(nn2,a),byExOrbbv(nn2,a+1)-1
						i=Indexsbv(numcolbv(vv3),1)

						Cx=Bx*t1(i,c)

						vv4=ftfm(cIndex(i,a),nn2); pvd(vv4,ct)=pvd(vv4,ct)-Cx

						xx1=cIndex(j,a)
						xx2=cIndex(i,b)
						vv4=ftfm(xx1,xx2); pvd(vv4,ct)=pvd(vv4,ct)+Cx
					enddo
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel

	!$omp parallel default(shared) private(ct,k,vv,c,Ax,j,mm2,vv2,nn2,i,a,Bx,vv3,b,Cx,xx1,xx2,vv4)
	ct=getThreadNumber()
	!$omp do
	do k = 1,Nel ! P(ab) -F(k,c)*t(j,i,c,a)*t(k,b)
		do vv = whOVf(k),ferow(k+1)-1
			c=fnumcol(vv)

			Ax=vF(vv)

			do j = 1,Nel
				mm2=cIndex(j,c)

				do vv2 = erow(mm2),erow(mm2+1)-1
					nn2=numcol(vv2)

					i=Indexs(nn2,1)
					a=Indexs(nn2,2)
					
					Bx=Ax*vt(vv2)

					do vv3 = byExOrb(nn2,j),byExOrb(nn2,j+1)-1
						b=Indexs(numcol(vv3),2)

						Cx=Bx*t1(k,b)

						vv4=ftfm(nn2,cIndex(j,b)); pvd(vv4,ct)=pvd(vv4,ct)-Cx

						xx1=cIndex(i,b)
						xx2=cIndex(j,a)
						vv4=ftfm(xx1,xx2); pvd(vv4,ct)=pvd(vv4,ct)+Cx
					enddo
				enddo
			enddo
		enddo
	enddo
	!$omp end parallel

	! collect threads
	!$omp parallel default(shared) private(vv,vv2)
	!$omp do
	do vv = 1,Nue
		do vv2=1,Nth
			vd(vv)=vd(vv)+pvd(vv,vv2); pvd(vv,vv2)=0
		enddo
	enddo
	!$omp end parallel

	return
	end subroutine projection_ccsd_doubles_spin_cue_spare