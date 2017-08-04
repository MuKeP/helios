	subroutine prepareSparseIndexInformation

	use coupledClusterSparse
	use coupledCluster      , only: Nel,No,Ne,Nth,cueDistance,mol,cuebd

	implicit none

	integer(kind=iglu) :: i,j,a,b,c,d,mm,nn,pp,vv,cnt,sta,sto
	real   (kind=rglu) :: thCentroid(4),rnt2,rpert2


!	mol%cueLevel(1)=150.d0
!	mol%cueLevel(2)=150.d0

!	write (*,*) cuebd%radius(1); stop

	Ne=0
	do i = 1,Nel
	do a = Nel+1,No
		if (cueDistance(i,a).GT.mol%cueLevel(2)) cycle !length of excitation is bigger than mol%cueLevel(2)
		Ne=Ne+1
	enddo
	enddo

	allocate ( Indexs(Ne,2), Indexsbv(Ne,2) )
	Indexs=0; Indexsbv=0

	allocate ( intersectOrbitals(Nel,Nth),intersectExcitations(Ne,Nth) )
	intersectOrbitals=0; intersectExcitations=0

!	Indexs contains orbital pairs by excitation.
!	sorted by occupied.
	mm=0
	do i = 1,Nel
	do a = Nel+1,No
		if (cueDistance(i,a).GT.mol%cueLevel(2)) cycle

		mm=mm+1
		Indexs(mm,1)=i
		Indexs(mm,2)=a
	enddo
	enddo

!	Indexsbv contains orbital pairs by excitation.
!	sorted by vacant. (name = Indexs by vacant)
	mm=0
	do a = Nel+1,No
	do i = 1,Nel
		if (cueDistance(i,a).GT.mol%cueLevel(2)) cycle

		mm=mm+1
		Indexsbv(mm,1)=i
		Indexsbv(mm,2)=a
	enddo
	enddo

!	cIndex contains number of excitations by orbitals.
!	unexistent excitations will get 0 number.
	allocate ( cIndex(No,No) )
	cIndex=0

	do mm = 1,Ne
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		cIndex(i,a)=mm
	enddo

!   ~~~~~ Prepare indexInfo for spare t1 ~~~~~   !
	t1Ne=0
	do i = 1,Nel
		do a = Nel+1,No
			if (cueDistance(i,a).GT.mol%cueLevel(1)) cycle
			t1Ne=t1Ne+1
		enddo
	enddo

	allocate ( t1erow(Nel+1),t1numcol(0:t1Ne),t1cIndex(1:Nel,Nel+1:No), t1Indexs(t1Ne,2) )
	t1erow=0; t1numcol=0; t1cIndex=0; t1Indexs=0

	mm=1
	do i = 1,Nel
		t1erow(i)=mm
		do a = Nel+1,No
			if (cueDistance(i,a).GT.mol%cueLevel(1)) cycle
			t1numcol(mm)=a
			t1cIndex(i,a)=mm

			t1Indexs(mm,1)=i
			t1Indexs(mm,2)=a

			mm=mm+1
		enddo
	enddo
	t1erow(Nel+1)=mm
!   ~~~~~~~~~ List of suitable excitations ~~~~~~~~~   !
	allocate ( t1mrEx(No,0:No) )
	! t1mrEx = in case l!=maxnei
	! t1mrEx(i,:) - list of excitations, that fits "i" by mol%cueLevel(1).
	t1mrEx=0
	do i = 1,Nel ! vacant by occupied.
		cnt=0
		do a = Nel+1,No
			if (cueDistance(i,a).GT.mol%cueLevel(1)) cycle
			cnt=cnt+1
			t1mrEx(i,cnt)=a
		enddo
		t1mrEx(i,0)=cnt
	enddo

	! t1mrEx(:,0) howmany spin-orbitals are in the list.

	do a = Nel+1,No ! occupied by vacant
		cnt=0
		do i = 1,Nel
			if (cueDistance(i,a).GT.mol%cueLevel(1)) cycle
			cnt=cnt+1
			t1mrEx(a,cnt)=i
		enddo
		t1mrEx(a,0)=cnt
	enddo
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!	byExOrb contains position of the first appearance for occupied MO
!	in current row (because all amplitudes are collected by rows and
!	"Indexs" array is sorted by occupied orbitals.
!	values will be defined later.

!	byExOrbbv is the same massive but for Indexsbv (when sorted by vacant orbitals)

	allocate ( byExOrb(Ne,Nel+1),byExOrbbv(Ne,Nel+1:No+1) )
	byExOrb=0; byExOrbbv=0

	allocate ( occEx(0:Ne,0:Nel) )
	occEx=0

!   ~~~~~~~~~ List of suitable excitations ~~~~~~~~~   !

!	occEx contains list of suitable occupied orbitals by excitation. (only for 1 diagram)
	do mm = 1,Ne
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		cnt=0
		do j = 1,Nel
			if (max(cueDistance(i,j),cueDistance(a,j)).GT.mol%cueLevel(2)) cycle
			cnt=cnt+1
			occEx(mm,cnt)=j
		enddo
		occEx(mm,0)=cnt
	enddo			

	allocate ( mrEx(No,0:No) )
	! mrEx = in case l!=maxnei
	! mrEx(i,:) - list of excitations, that fits "i" by cueDistance.
	mrEx=0
	do i = 1,Nel ! vacant by occupied.
		cnt=0
		do a = Nel+1,No
			if (cueDistance(i,a).GT.mol%cueLevel(2)) cycle
			cnt=cnt+1
			mrEx(i,cnt)=a
		enddo
		mrEx(i,0)=cnt
	enddo

	!mrEx(:,0) howmany spin-orbitals are in the list.
	do a = Nel+1,No ! occupied by vacant
		cnt=0
		do i = 1,Nel
			if (cueDistance(i,a).GT.mol%cueLevel(2)) cycle
			cnt=cnt+1
			mrEx(a,cnt)=i
		enddo
		mrEx(a,0)=cnt
	enddo
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	! size of the task.
	Nue=0
	do mm = 1,Ne
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		do nn = 1,Ne
			j=Indexs(nn,1)
			b=Indexs(nn,2)

			if (btest(i+j,0).NE.btest(a+b,0)) cycle

			thCentroid(1)=cueDistance(i,b)
			thCentroid(2)=cueDistance(j,a)
			thCentroid(3)=cueDistance(a,b)
			thCentroid(4)=cueDistance(i,j)
			if (maxval(thCentroid).GT.mol%cueLevel(2)) cycle

			Nue=Nue+1
		enddo
	enddo

	rnt2=float(Nel)*float(Nel)*float(Nel)*float(Nel); rpert2=100000.d0*float(Nue)/rnt2
	!write (*,*) rnt2,nue,rpert2

	!write (ou,300) Nel,Ne,Nue,rpert2
	!write (fu,300) Nel,Ne,Nue,rpert2
	!write (*,300) Nel,Ne,Nue,rpert2

	! amplitude arrays.
	allocate ( t1(1:Nel,Nel+1:No),d1(1:Nel,Nel+1:No) )
	allocate ( pvd(0:Nue,1:Nth) )
	allocate ( vt(0:Nue),vd(0:Nue) )

	! index information
	allocate ( numcol  (Nue),erow  (0:Ne+1) )
	allocate ( numcolbv(Nue),erowbv(0:Ne+1) )

	numcol=0; erow=0; numcolbv=0; erowbv=0

	! ftfm=find t for me
	! position in vvd&vvt vectors by two excitations.
	!! the most expensive by memory. !!

!	stop
	allocate ( ftfm(0:Ne,0:Ne) )

	! in this block we will define byExOrb
	! when sorted by occupied
	ftfm=0
	pp = 1
	erow(0)=pp
	do mm = 1,Ne
		i=Indexs(mm,1)
		a=Indexs(mm,2)

		erow(mm)=pp
		do nn = 1,Ne
			j=Indexs(nn,1)
			b=Indexs(nn,2)

			if (btest(i+j,0).NE.btest(a+b,0)) cycle

			thCentroid(1)=cueDistance(i,b)
			thCentroid(2)=cueDistance(j,a)
			thCentroid(3)=cueDistance(a,b)
			thCentroid(4)=cueDistance(i,j)
			if (maxval(thCentroid).GT.mol%cueLevel(2)) cycle
			numcol(pp)=nn

			ftfm(mm,nn)=pp
			pp=pp+1
		enddo
	enddo
	erow(Ne+1)=pp

	do mm = 1,Ne
		sta=erow(mm)
		sto=erow(mm+1)-1

		pp=Indexs(numcol(sta),1)
		byExOrb(mm,pp)=sta

		do vv = sta,sto
			j=Indexs(numcol(vv),1)

			if (j.NE.pp) then
				byExOrb(mm,j)=vv
				pp=j
			endif
		enddo
	enddo

	pp=Nue
	do mm = Ne,1,-1
		do j = Nel,1,-1
			if (byExOrb(mm,j).NE.0) then
				pp=byExOrb(mm,j)
			else
				byExOrb(mm,j)=pp
			endif
		enddo
	enddo

	do mm = 1,Ne-1
		byExOrb(mm,Nel+1)=byExOrb(mm+1,1)
	enddo
	byExOrb(Ne,Nel+1)=Nue

	!when sorted by vacant
	pp=1; erowbv(0)=pp
	do mm = 1,Ne
		i=Indexsbv(mm,1)
		a=Indexsbv(mm,2)

		erowbv(mm)=pp
		do nn = 1,Ne
			j=Indexsbv(nn,1)
			b=Indexsbv(nn,2)

			if (btest(i+j,0).NE.btest(a+b,0)) cycle

			thCentroid(1)=cueDistance(i,b)
			thCentroid(2)=cueDistance(j,a)
			thCentroid(3)=cueDistance(a,b)
			thCentroid(4)=cueDistance(i,j)
			if (maxval(thCentroid).GT.mol%cueLevel(2)) cycle
			numcolbv(pp)=nn

			pp=pp+1
		enddo
	enddo
	erowbv(Ne+1)=pp

	do mm = 1,Ne
		sta=erowbv(mm)
		sto=erowbv(mm+1)-1

		pp=Indexsbv(numcolbv(sta),2)
		byExOrbbv(mm,pp)=sta

		do vv = sta,sto
			b=Indexsbv(numcolbv(vv),2)

			if (b.NE.pp) then
				byExOrbbv(mm,b)=vv
				pp=b
			endif
		enddo
	enddo

	pp=Nue
	do mm = Ne,1,-1
		do b = No,Nel+1,-1
			if (byExOrbbv(mm,b).NE.0) then
				pp=byExOrbbv(mm,b)
			else
				byExOrbbv(mm,b)=pp
			endif
		enddo
	enddo

	do mm = 1,Ne-1
		byExOrbbv(mm,No+1)=byExOrbbv(mm+1,Nel+1)
	enddo
	byExOrbbv(Ne,No+1)=Nue

	! Fock spare information

	allocate (ferow(No+1),whOVf(No))

	!write (90,*) t1mrEx,t1cIndex,t1Indexs,t1erow,t1numcol
	!write (90,*) Indexs,Indexsbv ,cIndex ,erow,numcol,erowbv,numcolbv
	!write (90,*) byExOrbbv,occEx,mrEx,byExOrb,ftfm
	!stop

	!memy=memoryControl(30)+memoryControl(31)
	!write (su,200) memy; write (ou,200) memy; write (fu,200) memy; write (se,200) memy
200	format ( 4X,'All necessary memory:',2X,F12.4,' MB'/)
300 format (/4X,'Number of electrons:                ',2X,i10/&
&            4X,'Number of single excitations:       ',2X,i10/&
&            4X,'Number of accounted T2 amplitudes:  ',2X,i10/&
&            4X,'Promille of accounted T2 amplitudes:',2X,F10.5,'E-03%'/)

	return
	end subroutine prepareSparseIndexInformation

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	subroutine initSpareCC

	use coupledClusterSparse
	use coupledCluster      , only: Nel,No,gluCompare
	use coupledCluster      , only: F,NFnz

	implicit none

	integer(kind=iglu) :: i,j,k,l

	!no guess for new point.
	t1=0; vt=0
	d1=0; vd=0
	pvd=0

	!index information to use non-zero Fockian elements
	!"vector Fock"                vF     : is vectorized F with only non-zero elements
	!"Fock non-zero"              Fnz    : contains ia by position in vF vector
	!"Fock enter row"             ferow  : indicates "starts" of rows in vF
	!"Fock column number"         fnumcol: contains column number for every element in vF
	!"where is occ-vac separator" whOVf  : shows the start of occ block and the end of vac block for row

	! Usage

	! whOVf(c)   the start of vacant block in the "c" row of F.
	! whOVf(c)-1 the end of occupied block in the "c" row of F.

	! Occupied:
	!********* sta=ferow(c)
	!********* sto=whOVf(c)-1
	!********* do k = sta,sto
	!*********     Ax=vF(k); i=fnumcol(k)
	!********* enddo
	! is equal to
	!********* do i = 1,Nel
	!*********     Ax=F(c,i)
	!********* enddo
	!
	! Vacant:
	!********* sta=whOVf(c)
	!********* sto=ferow(c+1)-1
	!********* do k = sta,sto
	!*********     Ax=vF(k); b=fnumcol(k)
	!********* enddo
	! is equal to
	!********* do b = Nel+1,No
	!*********     Ax=F(c,b)
	!********* enddo

	if (allocated(vF)) deallocate (vF,fnumcol)

	allocate (vF(NFnz),fnumcol(NFnz))

	k=1
	do i = 1,No
		ferow(i)=k; l=0
		do j = 1,No
			if ((l.EQ.0).AND.(j.GT.Nel)) then
				whOVf(i)=k; l=1
			endif
			if (abs(F(i,j)).GT.gluCompare) then
				fnumcol(k)=j
				vF(k)=F(i,j)
				k=k+1
			endif
		enddo
		if (l.EQ.0) then
			whOVf(i)=k
		endif
	enddo
	ferow(No+1)=k

	!write (90,*) ferow,whOVf,fnumcol,vF
	!stop

	return
	end subroutine initSpareCC