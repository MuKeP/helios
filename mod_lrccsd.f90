	module lrccsdModule

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	use glob          , only: uch,uchGet,uchSet,rglu,iglu,lglu,true,false,rspu
	use coupledCluster, only: F,R,t1,t2,d1,d2
	use hdb           , only: mol,statesbd,lrbd,diisbd,ou

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	character (len=*), parameter :: lrVersion='3.000'       !5
	character (len=*), parameter :: lrDate   ='2017.08.13'  !10
	character (len=*), parameter :: lrAuthor ='Vladimir V. Ivanov' !18

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	real(kind=rglu), allocatable :: r1(:,:),r2(:,:,:,:)

	! intermediates
	real(kind=rglu), dimension (:,:,:,:), allocatable :: ai1,ai2,ai3,ai4,ai5,ai6,ai7
	real(kind=rglu), dimension (:,:)    , allocatable :: Fab,Fij,Fia,at1,at2,at3,at4,at5,at6

	real(kind=rglu)                                   :: accuracy(5)

	real   (kind=rglu), allocatable :: st1(:,:),st2(:,:),st3(:,:)
	real   (kind=rglu), allocatable :: sd1(:,:),sd2(:,:),sd3(:,:)

	real   (kind=rglu), allocatable :: diisVectors(:,:),diisValues(:)
	real   (kind=rglu), allocatable :: diisMatrix(:,:),diisCoefficients(:)

	! storage arrays
	real(kind=rglu),    allocatable :: lrHoldStateVectorU2(:,:,:,:,:)
	real(kind=rglu),    allocatable :: lrHoldStateVectorU1(:,:,:)
	real(kind=rglu),    allocatable :: lrHoldStateResults(:,:)
	real(kind=rglu),    allocatable :: lrHoldStateEnergy(:)
	real(kind=rglu),    allocatable :: lrGuessCofs(:,:),lrCueTransform(:,:)
	real(kind=rglu),    allocatable :: lrOrthogonality(:,:)	

	integer(kind=iglu), allocatable :: lrGuessConfNum(:),lrGuessConf(:,:,:)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	type(uch)          :: umethod
	integer(kind=iglu) :: maxiterLR,iount,spin,Nel,No,lrStateNum,lrsu,N,Nocc,currentState,Nd
	real   (kind=rglu) :: epsLR,lrOmega,currentEnergy
	logical(kind=lglu) :: lrdoOrthogonal

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	private

	public :: lrVersion,lrDate,lrAuthor
	public :: setLRParameters,initLR,energyLR,finalizeLR
	public :: lrGuessCofs,lrGuessConfNum,lrGuessConf,lrCueTransform

	public :: Nel,No,ou
	public :: ai1,ai2,ai3,ai4,ai5,ai6,ai7
	public :: Fab,Fij,Fia,at1,at2,at3,at4,at5,at6
	public :: r1,r2,d1,d2,t1,t2,R,F

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine setLRParameters(method)

		implicit none

		character (len=*)  :: method

		logical(kind=lglu) :: ulrdoOrthogonal



		umethod=uchSet(method)

		lrStateNum=statesbd%nStates
		spin=statesbd%spin
		lrdoOrthogonal=lrbd%orthogonalize
		N=mol%nAtoms
		Nel=mol%nEls
		No=2*N
		Nocc=Nel/2
		iount=ou
		Nd=diisbd%steps

		select case (uchGet(umethod))
			case ('1')
			case ('2')
			case ('3')
			case default

		end select

		call controlMemoryLR('general','allocate')

		return
		end subroutine setLRParameters

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine initLR

		implicit none

		integer(kind=iglu) :: k,i,j,a,b,vv


		!call controlMemoryLR('diis','allocate')
		call controlMemoryLR('intermediates','allocate')

		call projection_lrccsd_spin_hf_intermediates_1
		
		lrHoldStateVectorU1=0; lrHoldStateVectorU2=0
		if (lrStateNum.LE.1) lrdoOrthogonal=false

		do k = 1,lrStateNum
			currentState=k
			call guessState(k)

			call iterator(iterationLR,energylr,lrbd%maxiters,lrbd%accuracy*real(2**(k-1),rglu),false)

			lrHoldStateVectorU1(    :,:,k)=r1
			lrHoldStateVectorU2(:,:,:,:,k)=r2
			lrHoldStateEnergy  (        k)=lrOmega
		enddo

		call controlMemoryLR('intermediates','deallocate')

		call resortStates

		return
		end subroutine initLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine iterationLR(iteration,epsilon,saccuracy)

		implicit none

		integer(kind=iglu) :: iteration
		real   (kind=rglu) :: epsilon,saccuracy(5)


		select case (uchGet(umethod))
			case('spin-cue-ccsd')
			case('spin-u-ccsd','spin-r-ccsd')
		end select

		call projection_lrccsd_spin_hf_intermediates_2

		call projection_lrccsd_singles_spin_hf
		call projection_lrccsd_doubles_spin_hf

		call stepLR

		saccuracy=accuracy

		return

	!	if (maxval(accuracy).LT.epsilon) return

		if (diisbd%enabled) then
			call pushLRVectors

			if (iteration.GT.Nd) then
				call newLRVectors(Nd)
			else
				call newLRVectors(iteration)
			endif
		endif


		return
		end subroutine iterationLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine stepLR

		implicit none

		integer(kind=iglu) :: i,a,j,b
		real(kind=rglu)    :: Ax,sss,accur1,energy(5),rez,s1,s2


		s1=0
		s2=0
		!$omp parallel default(shared) private (i,a) reduction(+:s1,s2)
		!$omp do
		do i=1,Nel
			do a=Nel+1,No
				s1=s1+d1(i,a)*r1(i,a)
				s2=s2+r1(i,a)*r1(i,a)
			enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,j,a,b) reduction(+:s1,s2)
		!$omp do
		do i=1,Nel-1
			do j=i+1,Nel
				do a=Nel+1,No-1
					do b=a+1,No
						s1=s1+d2(i,j,a,b)*r2(i,j,a,b)
						s2=s2+r2(i,j,a,b)*r2(i,j,a,b)
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel

		lrOmega=s1/s2

		currentEnergy=lrOmega

		accuracy=-1
		!$omp parallel default(shared) private (i,a)
		!$omp do
		do i = 1,Nel
		do a = Nel+1,No
			d1(i,a)=d1(i,a)-lrOmega*r1(i,a)
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,a)
		!$omp do
		do i=1,Nel
		do a=Nel+1,No
			r1(i,a)=r1(i,a)-lrbd%iterStep*d1(i,a)
		enddo
		enddo
		!$omp end parallel
		
		accuracy(1)=maxval(abs(d1))

		!$omp parallel default(shared) private (i,j,a,b)
		!$omp do
		do i=1,Nel
		do a=Nel+1,No
		do j=1,Nel
		do b=Nel+1,No
			d2(i,j,a,b)=d2(i,j,a,b)-lrOmega*r2(i,j,a,b)
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,j,a,b)
		!$omp do
		do i=1,Nel
		do a=Nel+1,No
		do j=1,Nel
		do b=Nel+1,No
			r2(i,j,a,b)=r2(i,j,a,b)-lrbd%iterStep*d2(i,j,a,b)
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		accuracy(2)=maxval(abs(d2))

		if(spin.EQ.0) then
			!$omp parallel default(shared) private(i,a,rez)
			!$omp do		
			do i=1,Nel,2
				do a=Nel+1,No,2
					rez=( r1(i,a)+r1(i+1,a+1) )/2
					r1(i,a)    = rez
					r1(i+1,a+1)= rez
				enddo
			enddo
			!$omp end parallel
		else
			!$omp parallel default(shared) private(i,a,rez)
			!$omp do	
			do i=1,Nel,2
				do a=Nel+1,No,2
					rez=( r1(i,a)-r1(i+1,a+1) )/2
					r1(i,a)    = rez
					r1(i+1,a+1)=-rez
				enddo
			enddo
			!$omp end parallel
		endif

		call normu
		if (lrdoOrthogonal) call orthogonalization(currentState)

		return
		end subroutine stepLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine energyLR(energy)
		implicit none

		real(kind=rglu)    :: s1,s2,acm,alambda,c0,AEL,sss1,sss2,energy(5)
		integer(kind=iglu) :: i,j,a,b


		energy=0
		acm=real(8065.54d0,rglu)*lrOmega
		alambda=real(1.0d7,rglu)/acm

		s1=0
		do i=1,Nel     
		do a=Nel+1,No
			s1=s1+r1(i,a)**2
		enddo
		enddo

		s2=0
		do i=1, Nel-1
		do j=i+1, Nel
		do a=Nel+1,No-1
		do b=a+1, No
			s2=s2+r2(i,j,a,b)**2
		enddo
		enddo
		enddo
		enddo

		sss1=0
		sss2=0
		do  i=1,Nel
		do  a=Nel+1,No
			sss1=sss1+Fia(i,a)*r1(i,a)
			do  j=1,Nel
			do  b=Nel+1,No
				sss2=sss2+R(j,b,i,a)*r2(i,j,a,b)
			enddo
			enddo
		enddo
		enddo
      
		c0=(sss1+sss2/4)/lrOmega
		AEL=(s1+2*s2)/(s1+s2)


		energy(1)=currentEnergy
		energy(2)=c0
		energy(3)=AEL

		return
		end subroutine energyLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine finalizeLR

		implicit none


		call controlMemoryLR('general','deallocate')
		call controlMemoryLR('intermediates','deallocate')
		call controlMemoryLR('diis','deallocate')

		return
		end subroutine finalizeLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine pushLRVectors

		implicit none

		integer(kind=iglu) :: i,j,k,a,b,vv


		do k = Nd,2,-1
			do vv = 1,UBound(sd1,1)
				sd1(vv,k)=sd1(vv,k-1); st1(vv,k)=st1(vv,k-1)
			enddo

			do vv = 1,UBound(sd2,1)
				sd2(vv,k)=sd2(vv,k-1); st2(vv,k)=st2(vv,k-1)
			enddo
		enddo

		vv=0
		do i = 1,Nel
		do a = Nel+1,No
			if (btest(i,0).NE.btest(a,0)) cycle
			vv=vv+1; sd1(vv,1)=d1(i,a)-lrOmega*r1(i,a); st1(vv,1)=r1(i,a)
		enddo
		enddo

		vv=0
		do i = 1,Nel-1
		do j = i+1,Nel
		do a = Nel+1,No-1
		do b = a+1,No
			if (btest(i+j,0).NE.btest(a+b,0)) cycle
			vv=vv+1; sd2(vv,1)=d2(i,j,a,b)-lrOmega*r2(i,j,a,b); st2(vv,1)=r2(i,j,a,b)
		enddo
		enddo
		enddo
		enddo

		return
		end subroutine pushLRVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine newLRVectors(iteration)

		use math, only: gltred4

		implicit none

		integer(kind=iglu) :: iteration,k,l,vv,pp,i,j,a,b,c,mm,nn,ub
		real   (kind=rglu) :: sum


		do k = iteration,2,-1
		do l = iteration,2,-1
			diisMatrix(k,l)=diisMatrix(k-1,l-1)
		enddo
		enddo

		do k = 1,iteration
			call diisDotProduct(k,1)
		enddo

		if (iteration.EQ.1) return

		do k = 1,iteration+1
			diisMatrix(k,iteration+1)=1
			diisMatrix(iteration+1,k)=1
		enddo; diisMatrix(iteration+1,iteration+1)=0

		ub=iteration+1
		call gltred4(diisMatrix(1:ub,1:ub),diisVectors(1:ub,1:ub),diisValues(1:ub),ub,real(1.d-100,rglu),real(1.d-300,rglu))

		do k = 1,iteration
			diisCoefficients(k)=0
			do l = 1,iteration+1
				diisCoefficients(k)=diisCoefficients(k)+diisVectors(k,l)*diisVectors(iteration+1,l)/diisValues(l)
			enddo
		enddo

		vv=0
		do i = 1,Nel
		do a = Nel+1,No
			if (btest(i,0).NE.btest(a,0)) cycle
			vv=vv+1; sum=0
			do pp = 1,Nd
				sum=sum+diisCoefficients(pp)*st1(vv,pp)
			enddo
			r1(i,a)=sum
		enddo
		enddo

		vv=0
		do i = 1,Nel-1
		do j = i+1,Nel
		do a = Nel+1,No-1
		do b = a+1,No
			if (btest(i+j,0).NE.btest(a+b,0)) cycle
			vv=vv+1; sum=0
			do pp = 1,Nd
				sum=sum+diisCoefficients(pp)*st2(vv,pp)
			enddo
			r2(i,j,a,b)= sum
			r2(j,i,a,b)=-sum
			r2(i,j,b,a)=-sum
			r2(j,i,b,a)= sum
		enddo
		enddo
		enddo
		enddo

		return
		end subroutine newLRVectors

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine diisDotProduct(k1,k2)

		implicit none

		integer(kind=iglu) :: k1,k2,vv
		real   (kind=rglu) :: sum1,sum2,sum3


		sum1=0 ! singles
		!$omp parallel default(shared) private(vv) reduction(+:sum1)
		!$omp do
		do vv = 1,UBound(sd1,1)
			sum1=sum1+sd1(vv,k1)*sd1(vv,k2)
		enddo
		!$omp end parallel

		sum2=0 ! doubles
		!$omp parallel default(shared) private(vv) reduction(+:sum2)
		!$omp do
		do vv = 1,UBound(sd2,1)
			sum2=sum2+sd2(vv,k1)*sd2(vv,k2)
		enddo
		!$omp end parallel

		sum3=0 ! triples
		diisMatrix(k1,k2)=sum1+sum2+sum3; diisMatrix(k2,k1)=diisMatrix(k1,k2)

		return
		end subroutine diisDotProduct

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine resortStates

		implicit none

		integer(kind=iglu) :: glcount,swcount,i,j,a,b,k,l
		real(kind=rglu)    :: swEnergy,swVal,swRez(8)


		glcount=0
		do
			swcount=0
			do k = 1,lrStateNum-1
				l=k+1

				if (lrHoldStateEnergy(l).LT.lrHoldStateEnergy(k)) then
					swcount=swcount+1

					swRez=lrHoldStateResults(l,:)
					lrHoldStateResults(l,:)=lrHoldStateResults(k,:)
					lrHoldStateResults(k,:)=swRez

					swEnergy=lrHoldStateEnergy(l)
					lrHoldStateEnergy(l)=lrHoldStateEnergy(k)
					lrHoldStateEnergy(k)=swEnergy

					do i = 1,Nel
					do a = Nel+1,No
						swVal=lrHoldStateVectorU1(i,a,l)
						lrHoldStateVectorU1(i,a,l)=lrHoldStateVectorU1(i,a,k)
						lrHoldStateVectorU1(i,a,k)=swVal
						do j = 1,Nel
						do b = Nel+1,No
							swVal=lrHoldStateVectorU2(i,j,a,b,l)
							lrHoldStateVectorU2(i,j,a,b,l)=lrHoldStateVectorU2(i,j,a,b,k)
							lrHoldStateVectorU2(i,j,a,b,k)=swVal
						enddo
						enddo
					enddo
					enddo
				endif
			enddo
			glcount=glcount+swcount
			if (swcount.EQ.0) exit
		enddo

		if (glcount.NE.0) write (lrsu,100)

100		format (4X,'Incorrect order of excited states was fixed.')

		return
		end subroutine resortStates

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine guessState(istate)

		implicit none

		integer(kind=iglu) :: istate,k,i,j,a,b
		real(kind=rglu)    :: val,val2


		if (.NOT.allocated(lrguessCofs)) then
			do i = 1,Nocc
				r1(2*i-1,2*i-1+Nel)=real(0.1,rglu)
				if (spin.EQ.0) then
					r1(2*i  ,2*i+Nel)= real(0.1,rglu)
				else
					r1(2*i  ,2*i+Nel)=real(-0.1,rglu)
				endif
			enddo
			return
		endif

		!r1=0; r2=0
		!do k = 1,lrGuessConfNum(istate)
		!	i  =lrGuessConf(istate,k,1)
		!	a  =lrGuessConf(istate,k,2)
		!	val=lrGuessCofs(istate,k)
				
			!if (lrdoCUE) then
			!	do j = 1,Nocc
			!	do b = Nocc+1,N
			!		val2=r1(2*j-1,2*b-1-eno)
			!		val2=val2+val*lrCueTransform(j,i)*lrCueTransform(b,a)

			!		r1(2*j-1,2*b-1-eno)=val2
			!		if (spin.EQ.0) then
			!			r1(2*j  ,2*b-eno  )= val2
			!		else
			!			r1(2*j  ,2*b-eno  )=-val2
			!		endif
			!	enddo
			!	enddo
			!else
			!	if (spin.EQ.0) then
			!		r1(2*i-1,2*a-Nel-1)= val
			!	else
			!		r1(2*i  ,2*a-Nel  )=-val
			!	endif
			!endif
		!enddo

		return
		end subroutine guessState

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine orthogonalization(cstate)
		implicit none

		integer(kind=iglu) :: i,j,a,b,k,cstate,nstates
		real(kind=rglu)    :: sum1,sum2,normm


		nstates=cstate-1; if (nstates.LE.0) return

		lrOrthogonality=0
		do k = 1,nstates

			sum1=0; sum2=0
			do i = 1,Nel
			do a = Nel+1,No
				sum1=sum1+r1(i,a)*lrHoldStateVectorU1(i,a,k)
			enddo
			enddo

			do i = 1,Nel
			do j = 1,Nel
			do a = Nel+1,No
			do b = Nel+1,No
				sum2=sum2+r2(i,j,a,b)*lrHoldStateVectorU2(i,j,a,b,k)
			enddo
			enddo
			enddo
			enddo

			lrOrthogonality(1,k)=sum1+sum2
		enddo

		do i = 1,Nel
		do a = Nel+1,No
			sum1=0
			do k = 1,nstates
				sum1=sum1+lrHoldStateVectorU1(i,a,k)*lrOrthogonality(1,k)
			enddo
			r1(i,a)=r1(i,a)-sum1
		enddo
		enddo

		do i = 1,Nel
		do a = Nel+1,No
			do j = 1,Nel
			do b = Nel+1,No
				sum2=0
				do k = 1,nstates
					sum2=sum2+lrHoldStateVectorU2(i,j,a,b,k)*lrOrthogonality(1,k)
				enddo
				r2(i,j,a,b)=r2(i,j,a,b)-sum2
			enddo
			enddo
		enddo
		enddo

		return
		end subroutine orthogonalization

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !





!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !



!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine normu
		implicit none

		integer(kind=iglu) :: a,b,i,j
		real(kind=rglu)    :: s1,s2
      
      
		s1=0
		!$omp parallel default(shared) private (i,a) reduction(+:s1)
		!$omp do
		do i = 1,Nel
			do a = Nel+1,No
				s1=s1+r1(i,a)**2
			enddo
		enddo
		!$omp end parallel
            
		s2=0
		!$omp parallel default(shared) private (i,j,a,b) reduction(+:s2)
		!$omp do	
		do i = 1,Nel-1
			do j = i+1,Nel
				do a = Nel+1,No-1
					do b = a+1,No
						s2=s2+r2(i,j,a,b)**2
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel
      
		s1=real(1,rglu)/sqrt(s1+s2)

		!$omp parallel default(shared) private (i,a)
		!$omp do
		do i = 1,Nel
			do a = Nel+1,No
				r1(i,a)=r1(i,a)*s1
			enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,j,a,b,s2)
		!$omp do
		do i = 1,Nel-1
			do j = i+1,Nel
				do a = Nel+1,No-1
					do b = a+1,No
						s2=r2(i,j,a,b)*s1
						r2(i,j,a,b)= s2
						r2(j,i,a,b)=-s2
						r2(i,j,b,a)=-s2
						r2(j,i,b,a)= s2
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel

		return
		end subroutine normu

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		subroutine controlMemoryLR(section,action)
		implicit none

		character (len=*) :: section,action
		integer(kind=iglu) :: i,j,a,b,vv,err


		select case (section)

			case ('general')
				select case (action)
					case ('allocate')
						allocate (r1(Nel,Nel+1:No),r2(Nel,Nel,Nel+1:No,Nel+1:No))
						r1=0; r2=0

						allocate (lrHoldStateVectorU1(Nel,Nel,lrStateNum),&
						          lrHoldStateVectorU2(Nel,Nel,Nel,Nel,lrStateNum),&
						          lrHoldStateEnergy(lrStateNum),&
						          lrHoldStateResults(lrStateNum,8),&
						          lrOrthogonality(2,lrStateNum))

						lrHoldStateVectorU1=0; lrHoldStateVectorU2=0
						lrHoldStateEnergy=0; lrHoldStateResults=0
						lrOrthogonality=0

					case ('deallocate')
						deallocate (lrHoldStateEnergy,lrHoldStateVectorU1,lrHoldStateVectorU2,&
						            lrHoldStateResults,lrOrthogonality,stat=err)

				end select

			case ('intermediates')
				select case (action)
					case ('allocate')
						allocate (Fab(Nel+1:No,Nel+1:No),Fij(Nel,Nel),Fia(Nel,Nel+1:No))
						Fab=0; Fij=0; Fia=0

						allocate (ai1(Nel,Nel+1:No,Nel,Nel+1:No),ai2(Nel,Nel,Nel,Nel+1:No))
						allocate (ai3(Nel,Nel+1:No,Nel+1:No,Nel+1:No),ai4(Nel,Nel,Nel,Nel))
						allocate (ai5(Nel+1:No,Nel+1:No,Nel+1:No,Nel+1:No),ai6(Nel,Nel,Nel,Nel+1:No))
						allocate (ai7(Nel,Nel+1:No,Nel+1:No,Nel+1:No))
						ai1=0; ai2=0; ai3=0; ai4=0; ai5=0; ai6=0; ai7=0

						allocate (at1(Nel,Nel),at2(Nel+1:No,Nel+1:No),at3(Nel,Nel))
						allocate (at4(Nel+1:No,Nel+1:No),at5(Nel,Nel),at6(Nel+1:No,Nel+1:No))
						at1=0; at2=0; at3=0; at4=0; at5=0; at6=0

					case ('deallocate')
						deallocate (Fab,Fij,Fia,ai1,ai2,ai3,ai4,ai5,ai6,ai7,at1,at2,at3,at4,at5,at6, stat=err)

				end select

			case ('diis')
				select case (action)
					case ('allocate')
						vv=0
						do i = 1,Nel
						do a = Nel+1,No
							if (btest(a,0).NE.btest(i,0)) cycle
							vv=vv+1
						enddo
						enddo

						allocate (st1(vv,Nd),sd1(vv,Nd))

						vv=0
						do i = 1,Nel-1
						do a = Nel+1,No-1
						do j = i+1,Nel
						do b = a+1,No
							if (btest(a+b,0).NE.btest(i+j,0)) cycle

							vv=vv+1
						enddo
						enddo
						enddo
						enddo

						allocate (st2(vv,Nd),sd2(vv,Nd))

						allocate (diisVectors(Nd+1,Nd+1),diisValues(Nd+1))
						allocate (diisMatrix (Nd+1,Nd+1),diisCoefficients(Nd))

					case ('deallocate')
						deallocate (diisVectors,diisValues,diisMatrix,diisCoefficients, stat=err)
						deallocate (st1,sd1,st2,sd2, stat=err)

				end select

		end select

		return
		end subroutine controlMemoryLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

		function ttime()

		implicit none

		character (len=10) :: ttime
		integer(kind=iglu) :: now(3)


		call itime(now(1:3)); write (ttime,"(i2.2,':',i2.2,':',i2.2)") now

		return
		end function ttime

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end module lrccsdModule
