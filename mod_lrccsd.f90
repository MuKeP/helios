	!MS$if defined(__unix) ! 1=windows; 2=unix
		!MS$define OS=2
	!MS$else
		!MS$define OS=1
	!MS$endif
	
	!MS$if defined(_OPENMP) ! in case of openMP standard
		!MS$define opnmp=1
	!MS$else
		!MS$define opnmp=0
	!MS$endif

	module lrccsdModule

	use glob, only: rglu,iglu,lglu,true,false

	!MS$if(OS.EQ.2)
		!MS$if(opnmp.EQ.1)
			include "omp_lib.h"
		!MS$endif
	!MS$endif

	character (len=*), parameter :: lrVersion='2.203'       !5
	character (len=*), parameter :: lrDate   ='2015.08.07'  !10
	character (len=*), parameter :: lrAuthor ='Vladimir V. Ivanov' !18

	! glu = real for global use  (storage and non-accuracy-demanding procedures).
	! spu = real for special use (accuracy-demanding procedures).
	!integer*4, parameter :: rglu=r8kind, rspu=r8kind
	!integer*4, parameter :: iglu=i4kind, ispu=i8kind
	!integer*4, parameter :: lglu=l1kind

	!real(kind=rglu), parameter :: gluzero=real(0,rglu)
	!real(kind=rspu), parameter :: spuzero=real(0,rspu)

	! intermediates
	real(kind=rglu), dimension (:,:,:,:), allocatable :: ai1,ai2,ai3,ai4,ai5,ai6,ai7
	real(kind=rglu), dimension (:,:)    , allocatable :: Fab,Fij,Fia,at1,at2,at3,at4,at5,at6

	! work arrays
	real(kind=rglu), dimension (:,:)    , allocatable :: x1,lrt1,lru1,lrf
	real(kind=rglu), dimension (:,:,:,:), allocatable :: x2,lrt2,lru2,lrint

	! storage arrays
	real(kind=rglu),    allocatable :: lrHoldStateVectorU2(:,:,:,:,:)
	real(kind=rglu),    allocatable :: lrHoldStateVectorU1(:,:,:)
	real(kind=rglu),    allocatable :: lrHoldStateResults(:,:)
	real(kind=rglu),    allocatable :: lrHoldStateEnergy(:)
	real(kind=rglu),    allocatable :: lrGuessCofs(:,:),lrCueTransform(:,:)
	real(kind=rglu),    allocatable :: lrOrthogonality(:,:)	

	integer(kind=iglu), allocatable :: lrGuessConfNum(:),lrGuessConf(:,:,:)

	integer(kind=iglu) :: maxiterLR=10000,iount=10,spin=0
	real(kind=rglu)    :: epsLR=real(1d-10,rglu),stepLR=real(0.03d0,rglu)
	real(kind=rglu)    :: lrOmega

	integer(kind=iglu) :: eNo,eNv,lrStateNum,lrsu,N,Nocc
	
	logical(kind=lglu) :: lrdoCUE,lrdoOrthogonal

	private

	public :: lrVersion,lrDate,lrAuthor,lrt1,lrt2,lru1,lru2,lrf,lrint,eNo,eNv,lrsu
	public :: lrHoldStateEnergy,lrHoldStateResults,lrHoldStateVectorU1,lrHoldStateVectorU2
	public :: setLR_Params,doLR,lrFinalize
	public :: lrGuessCofs,lrGuessConfNum,lrGuessConf,lrCueTransform

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		function ttime()

		implicit none

		character (len=10) :: ttime
		integer(kind=iglu) :: now(3)


		call itime(now(1:3)); write (ttime,"(i2.2,':',i2.2,':',i2.2)") now

		return
		end function ttime

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer(kind=iglu) function setLR_Params(uN,accuracy,iterStep,sSpin,maxiter,unt,stateNumber,ulrdoOrthogonal)

		implicit none

		real(kind=rglu)    :: accuracy,iterStep
		integer(kind=iglu) :: maxiter,unt,stateNumber,sSpin,uN
		logical(kind=lglu) :: ulrdoOrthogonal


		N=uN
		eNo=N
		Nocc=N/2        !!!!!!!!!!!
		lrStateNum=stateNumber
		lrdoOrthogonal=ulrdoOrthogonal
		epsLR=accuracy; stepLR=iterStep; spin=sSpin; maxiterLR=maxiter; iount=unt

		if (allocated(lrHoldStateEnergy)) then
			deallocate (lrHoldStateEnergy,lrHoldStateVectorU1,lrHoldStateVectorU2,lrHoldStateResults,lrOrthogonality)		
		endif

		allocate (lrHoldStateVectorU1(eNo,eNo,lrStateNum),lrHoldStateVectorU2(eNo,eNo,eNo,eNo,lrStateNum),&
		&         lrHoldStateEnergy(lrStateNum),lrHoldStateResults(8,lrStateNum),lrOrthogonality(lrStateNum,2))

		setLR_Params=0; return
		end function setLR_Params

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine doLR(udoCUE)

		implicit none

		integer(kind=iglu) :: i,j,k,l,a,b,d,e,c,ii,iter,aa,iterProp,scalarM=10
		real(kind=rglu)    :: accur1,accur2,val,accur,alambda,s1
		real(kind=rglu)    :: ael,acm,s2,sss1,sss2,c0
		logical(kind=lglu) :: udoCUE


		lrdoCUE=udoCUE
		if (.NOT.allocated(lrHoldStateEnergy)) stop '  Necessary to call setLR_Params at first.'

		allocate (Fab(eNo,eNo),Fij(eNo,eNo),Fia(eNo,eNo),x1(eNo,eNo))
		Fab=0; Fij=0; Fia=0; x1=0

		allocate (ai1(eNo,eNo,eNo,eNo),ai2(eNo,eNo,eNo,eNo))
		allocate (ai3(eNo,eNo,eNo,eNo),ai4(eNo,eNo,eNo,eNo))
		allocate (ai5(eNo,eNo,eNo,eNo),ai6(eNo,eNo,eNo,eNo))
		allocate (ai7(eNo,eNo,eNo,eNo), x2(eNo,eNo,eNo,eNo))

		ai1=0; ai2=0; ai3=0; ai4=0
		ai5=0; ai6=0; ai7=0; x2 =0

		allocate (at1(eNo,eNo),at2(eNo,eNo),at3(eNo,eNo))
		allocate (at4(eNo,eNo),at5(eNo,eNo),at6(eNo,eNo))

		at1=0; at2=0; at3=0
		at4=0; at5=0; at6=0

		write (lrsu,100)

		call prepint1
		
		lrHoldStateVectorU1=0; lrHoldStateVectorU2=0
		if (lrStateNum.LE.1) lrdoOrthogonal=false
		do k = 1,lrStateNum

			call guessState(k)

			!call orthogonalization(k)

			accur1=0; accur2=0
			do iter=1,maxiterLR
				call prepint2
				call c1
				call c2
				call energylr
				call stepu1(accur1)
				call symm
				call stepu2(accur2)
				call normu

				if (lrdoOrthogonal) call orthogonalization(k)

				iterProp=iter-scalarM*(iter/scalarM)
				if (iterProp.EQ.1) write(lrsu, 101) trim(ttime()),iter,accur1,accur2,lrOmega
				if ((max(accur1,accur2).LE.epsLR*real(2**(k-1),rglu)).AND.(iter.GT.1)) exit

			enddo

			write(lrsu, 102) iter,accur1,accur2,lrOmega
			write(iount,102) iter,accur1,accur2,lrOmega

			acm=real(8065.54d0,rglu)*lrOmega
			alambda=real(1.0d7,rglu)/acm

			s1=0
			do ii=1,eNo     
			do aa=1,eNv
				s1=s1+lru1(ii,aa)**2
			enddo
			enddo

			s2=0
			do i=1, eNo-1
			do j=i+1, eNo
			do a=1, eNv-1
			do b=a+1, eNv
				s2=s2+lru2(i,j,a,b)**2
			enddo
			enddo
			enddo
			enddo

			sss1=0
			sss2=0
			do  i=1,eNo
			do  a=1,eNv
				sss1=sss1+Fia(i,a)*lru1(i,a)
				do  j=1,eNo
				do  b=1,eNv
					sss2=sss2+lrint(j,i,b+eNo,a+eNo)*lru2(i,j,a,b)
				enddo
				enddo
			enddo
			enddo
      
			c0=(sss1+sss2*real(0.25,rglu))/lrOmega
			AEL=(s1+real(2,rglu)*s2)/(s1+s2)

			lrHoldStateVectorU1(    :,:,k)=lru1
			lrHoldStateVectorU2(:,:,:,:,k)=lru2
			lrHoldStateEnergy  (        k)=lrOmega
			lrHoldStateResults (      1,k)=accur1
			lrHoldStateResults (      2,k)=accur2
			lrHoldStateResults (      3,k)=AEL
			lrHoldStateResults (      4,k)=c0
			lrHoldStateResults (      5,k)=dble(iter)

			write(lrsu,'(/3X,A,ES11.3)')  ' Reference State Weigth   ',c0
			write(lrsu,'( 3X,A,1X,F6.3)') ' Average Excitation Level ',AEL
		enddo

		call resortStates

		deallocate (Fab,Fij,Fia,x1,ai1,ai2,ai3,ai4,ai5,ai6,ai7,x2,at1,at2,at3,at4,at5,at6)

	100 format(//16X,'LR-CCSD by Vladimir V. Ivanov 29.06.2013'/)
	101 format(1X,A,' Iteration: ',i5,' Accuracy  ',ES10.3,1X,ES10.3,' Omega= ',F10.6)
	102 format(3X  ,' Completed. ',i5,' Accuracy  ',ES10.3,1X,ES10.3,' Omega= ',F10.6)

		return
		end subroutine doLR

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

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

					swRez=lrHoldStateResults(:,l)
					lrHoldStateResults(:,l)=lrHoldStateResults(:,k)
					lrHoldStateResults(:,k)=swRez

					swEnergy=lrHoldStateEnergy(l)
					lrHoldStateEnergy(l)=lrHoldStateEnergy(k)
					lrHoldStateEnergy(k)=swEnergy

					do i = 1,eNo
					do a = 1,eNv
						swVal=lrHoldStateVectorU1(i,a,l)
						lrHoldStateVectorU1(i,a,l)=lrHoldStateVectorU1(i,a,k)
						lrHoldStateVectorU1(i,a,k)=swVal
						do j = 1,eNo
						do b = 1,eNv
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

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine guessState(istate)

		implicit none

		integer(kind=iglu) :: istate,k,i,j,a,b
		real(kind=rglu)    :: val,val2


		if (.NOT.allocated(lrguessCofs)) then
			do i = 1,Nocc
				lru1(2*i-1,2*i-1)=real(0.1,rglu)
				if (spin.EQ.0) then
					lru1(2*i  ,2*i)= real(0.1,rglu)
				else
					lru1(2*i  ,2*i)=real(-0.1,rglu)
				endif
			enddo
			return
		endif

		lru1=0; lru2=0
		do k = 1,lrGuessConfNum(istate)
			i  =lrGuessConf(istate,k,1)
			a  =lrGuessConf(istate,k,2)
			val=lrGuessCofs(istate,k)
				
			if (lrdoCUE) then
				do j = 1,Nocc
				do b = Nocc+1,N
					val2=lru1(2*j-1,2*b-1-eno)
					val2=val2+val*lrCueTransform(j,i)*lrCueTransform(b,a)

					lru1(2*j-1,2*b-1-eno)=val2
					if (spin.EQ.0) then
						lru1(2*j  ,2*b-eno  )= val2
					else
						lru1(2*j  ,2*b-eno  )=-val2
					endif
				enddo
				enddo
			else
				if (spin.EQ.0) then
					lru1(2*i-1,2*a-eNo-1)= val
				else
					lru1(2*i  ,2*a-eNo  )=-val
				endif
			endif
		enddo

		return
		end subroutine guessState

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine lrFinalize

		implicit none


		if (allocated(lrHoldStateEnergy))   deallocate (lrHoldStateEnergy)
		if (allocated(lrHoldStateVectorU1)) deallocate (lrHoldStateVectorU1)
		if (allocated(lrHoldStateVectorU2)) deallocate (lrHoldStateVectorU2)
		if (allocated(lrHoldStateResults))  deallocate (lrHoldStateResults)
		if (allocated(lrOrthogonality))     deallocate (lrOrthogonality)

		return
		end subroutine lrFinalize

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine orthogonalization(cstate)

		implicit none

		integer(kind=iglu) :: i,j,a,b,k,cstate,nstates
		real(kind=rglu)    :: sss1,sss2,normm

		nstates=cstate-1; if (nstates.LE.0) return


!		sss1=0; sss2=0
!		do i = 1,eNo
!		do a = 1,eNo
!			sss1=sss1+lrHoldStateVectorU1(i,a,1)**2
!			do j = 1,eNo
!			do b = 1,eNo
!				sss2=sss2+lrHoldStateVectorU2(i,j,a,b,1)**2
!			enddo
!			enddo
!		enddo
!		enddo
		!write (*,*) 'prev vector norm',sss1+sss2
		!read (*,*)

!		sss1=0; sss2=0
!		do i = 1,eNo
!		do a = 1,eNo
!			sss1=sss1+lru1(i,a)**2
!			do j = 1,eNo
!			do b = 1,eNo
!				sss2=sss2+lru2(i,j,a,b)**2
!			enddo
!			enddo
!		enddo
!		enddo
		!write (*,*) 'cur vector norm',sss1+sss2
		!read (*,*)

		lrOrthogonality=0
		do k = 1,nstates

			sss1=0; sss2=0
			do i = 1,eNo
			do a = 1,eNo
				sss1=sss1+lru1(i,a)*lrHoldStateVectorU1(i,a,k)
			enddo
			enddo

			do i = 1,eNo
			do j = 1,eNo
			do a = 1,eNv
			do b = 1,eNv
				sss2=sss2+lru2(i,j,a,b)*lrHoldStateVectorU2(i,j,a,b,k)
			enddo
			enddo
			enddo
			enddo

			lrOrthogonality(k,1)=sss1 +sss2
			!lrOrthogonality(k,2)=sss2
		enddo

		!write (*,'(A,4(1X,ES10.3))') 'Orthogonality',lrOrthogonality(1:nstates,1)

		!if (maxval(abs(lrOrthogonality(:,1))).LT.epsLR) return

		do i = 1,eNo
		do a = 1,eNo
			sss1=0
			do k = 1,nstates
				sss1=sss1+lrHoldStateVectorU1(i,a,k)*lrOrthogonality(k,1)
			enddo
			lru1(i,a)=lru1(i,a)-sss1
		enddo
		enddo

	!	return

		do i = 1,eNo
		do a = 1,eNo
			do j = 1,eNo
			do b = 1,eNo
				sss2=0
				do k = 1,nstates
					sss2=sss2+lrHoldStateVectorU2(i,j,a,b,k)*lrOrthogonality(k,1)
				enddo
				lru2(i,j,a,b)=lru2(i,j,a,b)-sss2
			enddo
			enddo
		enddo
		enddo

!		sss1=0; sss2=0
!		do i = 1,eNo
!		do a = 1,eNo
!			sss1=sss1+lru1(i,a)*lrHoldStateVectorU1(i,a,1)
!			do j = 1,eNo
!			do b = 1,eNo
!				sss2=sss2+lru2(i,j,a,b)*lrHoldStateVectorU2(i,j,a,b,1)
!			enddo
!			enddo
!		enddo
!		enddo
		!write (*,*) 'orthogon',sss1+sss2
		!read (*,*)

		return
		end subroutine orthogonalization


!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine symm

		implicit none

		integer(kind=iglu) :: i,a
		real(kind=rglu)    :: sss


		if(spin.EQ.0) then
			!$omp parallel default(shared) private(i,a,sss)
			!$omp do		
			do i=1,eNo,2
				do a=1,eNv,2
					sss=real(0.5,rglu)*( lru1(i,a)+lru1(i+1,a+1) )
					lru1(i,a)    =sss
					lru1(i+1,a+1)=sss
				enddo
			enddo
			!$omp end parallel
		endif

		if(spin.EQ.1) then
			!$omp parallel default(shared) private(i,a,sss)
			!$omp do	
			do i=1,eNo,2
				do a=1,eNv,2
					sss=real(0.5,rglu)*( lru1(i,a)-lru1(i+1,a+1) )
					lru1(i,a)    = sss
					lru1(i+1,a+1)=-sss
				enddo
			enddo
			!$omp end parallel
		endif

		return
		end subroutine symm

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prepint1
		implicit none

		integer(kind=iglu) :: i,j,k,l,a,b,c,d
		real(kind=rglu)    :: sum,ss1,ss2,sum1,sum2,sum3,sum4,sum5,sum6,sum7


		write (lrsu,'(/1X,A)') trim(ttime())//" Intermediates' calculation" 

		write (lrsu,'(1X,A)') trim(ttime())//' Part  1/10' 
		!$omp parallel default(shared) private(i,a,b,c,j,sum)
		!$omp do
		do  i=1, eNo
		do  a=1, eNv
		do  b=1, eNv
		do  c=1, eNv
			sum=0
			do  j=1,eNo
				sum=sum+lrint(i,j,b+eNo,c+eNo)*lrt1(j,a)
			enddo
			ai7(i,a,b,c)=lrint(i,a+eNo,b+eNo,c+eNo)-sum
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part  2/10' 
		!$omp parallel default(shared) private(i,j,k,a,c,sum)
		!$omp do
		do  i=1, eNo
		do  j=1, eNo
		do  k=1, eNo
		do  a=1, eNv
			sum=0
			do  c=1,eNv
				sum=sum+lrint(i,j,a+eNo,c+eNo)*lrt1(k,c)
			enddo
			ai6(i,j,k,a)=lrint(i,j,k,a+eNo)-sum
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part  3/10'
		!$omp parallel default(shared) private(a,b,c,d,k,l,sum1,sum2,sum3)
		!$omp do
		do  a=1, eNv
		do  b=1, eNv
		do  c=1, eNv
		do  d=1, eNv
			sum1=0
			do  k=1, eNo
				sum1=sum1 +lrt1(k,a)*lrint(k,b+eNo,c+eNo,d+eNo)-lrt1(k,b)*lrint(k,a+eNo,c+eNo,d+eNo)
			enddo

			sum2=0
			do  k=1, eNo
				do  l=1, eNo
					sum2=sum2 +lrt1(k,a)*lrt1(l,b)*lrint(k,l,c+eNo,d+eNo)
				enddo
			enddo

			sum3=0
			do  k=1, eNo
				do  l=1, eNo
					sum3=sum3 +lrt2(k,l,a,b)*lrint(k,l,c+eNo,d+eNo)
				enddo
			enddo
		 
			ai5(a,b,c,d)=real(0.5,rglu)*( lrint(a+eNo,b+eNo,c+eNo,d+eNo)-sum1+sum2)+sum3*real(0.25,rglu)
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part  4/10' 
		!$omp parallel default(shared) private(i,j,k,l,c,d,sum1,sum2,sum3)
		!$omp do
		do  i=1, eNo
		do  j=1, eNo
		do  k=1, eNo
		do  l=1, eNo
			sum1=0
			do  c=1, eNv
				sum1=sum1 +lrt1(j,c)*lrint(k,l,i,c+eNo)-lrt1(i,c)*lrint(k,l,j,c+eNo)
			enddo

			sum2=0
			do  c=1, eNv
				do  d=1, eNv
					sum2=sum2 +lrt1(i,c)*lrt1(j,d)*lrint(k,l,c+eNo,d+eNo)
				enddo
			enddo
		 
			sum3=0
			do  c=1, eNv
				do  d=1, eNv
					sum3=sum3 +lrt2(i,j,c,d)*lrint(k,l,c+eNo,d+eNo)
				enddo
			enddo

			ai4(i,j,k,l)=real(0.5,rglu)*( lrint(i,j,k,l)+sum1+sum2)+sum3*real(0.25,rglu)
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel
 
		write (lrsu,'(1X,A)') trim(ttime())//' Part  5/10' 
		!$omp parallel default(shared) private(i,c,a,b,l,k,d,sum1,sum2,ss1,ss2,sum3,sum4,sum5,sum6,sum7)
		!$omp do
		do  i=1, eNo
		do  c=1, eNv
		do  a=1, eNv
		do  b=1, eNv
			sum1=0
			do  d=1, eNv
				sum1=sum1 +lrt1(i,d)*ai5(a,b,c,d)
			enddo

			sum2=00
			do  l=1, eNo
				ss1=lrint(l,a+eNo,i,c+eNo)
				ss2=lrint(l,b+eNo,i,c+eNo)
				do  k=1, eNo
					do  d=1, eNv
						ss1=ss1 -lrt2(i,k,a,d)*lrint(l,k,c+eNo,d+eNo)
						ss2=ss2 -lrt2(i,k,b,d)*lrint(l,k,c+eNo,d+eNo)
					enddo
				enddo
				sum2=sum2 +lrt1(l,b)*ss1-lrt1(l,a)*ss2
			enddo

			sum3=0
			do  l=1, eNo
				do  k=1, eNo
					sum3=sum3 +lrt1(l,a)*lrt1(k,b)*lrint(c+eNo,i,l,k)
				enddo
			enddo

			sum4=0
			do  l=1, eNo
			do  k=1, eNo
				do  d=1, eNv
					sum4=sum4 +lrt1(l,d)*lrt2(k,i,a,b)*lrint(c+eNo,d+eNo,k,l)
				enddo
			enddo
			enddo

			sum5=0
			do  k=1, eNo
				sum5=sum5 +lrt2(i,k,a,b)*lrf(k,c+eNo)
			enddo

			sum6=0
			do  k=1, eNo
				do  l=1, eNo
					sum6=sum6 +lrt2(k,l,a,b)*lrint(i,c+eNo,k,l)
				enddo
			enddo

			sum7=0
			do  k=1, eNo
				do  d=1, eNv
					sum7=sum7 +lrt2(i,k,a,d)*lrint(k,b+eNo,c+eNo,d+eNo)-lrt2(i,k,b,d)*lrint(k,a+eNo,c+eNo,d+eNo)	
				enddo
			enddo

			ai3(i,c,a,b)=-lrint(i,c+eNo,a+eNo,b+eNo)+sum1*real(2,rglu)-sum2+sum3-sum4+sum5-sum6*real(0.5,rglu)+sum7
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part  6/10' 
		!$omp parallel default(shared) private(i,j,k,a,l,d,c,sum1,sum2,sum3,sum4,sum5,sum6,sum7,ss1,ss2)
		!$omp do		
		do i=1, eNo
		do j=1, eNo
		do k=1, eNo
		do a=1, eNv

			sum1=0
			do  l=1, eNo
				sum1=sum1 +lrt1(l,a)*ai4(i,j,k,l)
			enddo

			sum2=0
			do  d=1, eNv
				ss1=lrint(j,d+eNo,k,a+eNo)
				ss2=lrint(i,d+eNo,k,a+eNo)
				do  l=1, eNo
					do  c=1, eNv
						ss1=ss1 -lrt2(j,l,a,c)*lrint(k,l,d+eNo,c+eNo)
						ss2=ss2 -lrt2(i,l,a,c)*lrint(k,l,d+eNo,c+eNo)
					enddo
				enddo
				sum2=sum2 +lrt1(i,d)*ss1-lrt1(j,d)*ss2
			enddo

			sum3=0
			do  d=1, eNv
				do  c=1, eNv
					sum3=sum3 +lrt1(i,d)*lrt1(j,c)*lrint(k,a+eNo,d+eNo,c+eNo)
				enddo
			enddo

			sum4=0
			do  l=1, eNo
			do  d=1, eNv
				do  c=1, eNv
					sum4=sum4 +lrt1(l,d)*lrt2(i,j,a,c)*lrint(k,l,d+eNo,c+eNo)
				enddo
			enddo
			enddo

			sum5=0
			do  c=1, eNv
				sum5=sum5 +lrt2(i,j,a,c)*lrf(k,c+eNo)
			enddo

			sum6=0
			do  c=1, eNv
				do  d=1, eNv
					sum6=sum6 +lrt2(i,j,c,d)*lrint(k,a+eNo,c+eNo,d+eNo)
				enddo
			enddo

			sum7=0
			do  l=1, eNo
				do  c=1, eNv
					sum7=sum7 +lrt2(i,l,a,c)*lrint(j,c+eNo,k,l)
					sum7=sum7 -lrt2(j,l,a,c)*lrint(i,c+eNo,k,l)
				enddo
			enddo

			ai2(i,j,k,a)= -lrint(i,j,k,a+eNo)+sum1*real(2,rglu)+sum2-sum3-sum4+sum5-sum6*real(0.5,rglu)+sum7
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part  7/10'
		!$omp parallel default(shared) private(i,j,a,b,l,d,sum1,sum2,sum3)
		!$omp do
		do  i=1, eNo
		do  a=1, eNv
		do  j=1, eNo
		do  b=1, eNv

			sum1=0
			do  l=1, eNo
				sum1=sum1 +lrt1(l,b)*lrint(j,l,i,a+eNo)
			enddo

			sum2=0
			do  d=1, eNv
				sum2=sum2 +lrt1(i,d)*lrint(j,b+eNo,a+eNo,d+eNo)
			enddo

			sum3=0
			do  l=1, eNo
			do  d=1, eNv
				sum3=sum3 +(lrt1(i,d)*lrt1(l,b)-lrt2(i,l,b,d))*lrint(j,l,a+eNo,d+eNo)
			enddo
			enddo

			ai1(i,a,j,b)=  lrint(i,a+eNo,j,b+eNo)-sum1-sum2+sum3
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part  8/10'
		!$omp parallel default(shared) private(a,b,c,l,k,sum1,sum2,sum3,sum4)
		!$omp do
		do a=1, eNv
		do b=1, eNv

			sum1=0
			do  k=1, eNo
				sum1=sum1 +lrt1(k,a)*lrf(k,b+eNo)
			enddo
		
			sum2=0
			do  k=1, eNo
			do  c=1, eNv
				sum2=sum2 +lrt1(k,c)*lrint(k,a+eNo,b+eNo,c+eNo)
			enddo
			enddo
		
			sum3=0
			do  k=1, eNo
			do  l=1, eNo
			do  c=1, eNv
				sum3=sum3 +lrt1(k,c)*lrt1(l,a)*lrint(l,k,b+eNo,c+eNo)
			enddo
			enddo
			enddo
			
			sum4=0
			do  k=1, eNo
			do  l=1, eNo
			do  c=1, eNv
				sum4=sum4 +lrt2(l,k,a,c)*lrint(l,k,b+eNo,c+eNo)
			enddo
			enddo
			enddo

			Fab(a,b)= lrf(a+eNo,b+eNo)-sum1-sum2-sum3-sum4*real(0.5,rglu)
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part  9/10'
		!$omp parallel default(shared) private(i,j,k,d,c,sum1,sum2,sum3,sum4)
		!$omp do
		do i=1, eNo
		do j=1, eNo

			sum1=0
			do  c=1, eNv
				sum1=sum1 +lrt1(i,c)*lrf(j,c+eNo)
			enddo
			
			sum2=0
			do  k=1, eNo
			do  c=1, eNv
				sum2=sum2 +lrt1(k,c)*lrint(j,k,i,c+eNo)
			enddo
			enddo
			
			sum3=0
			do  k=1, eNo
			do  c=1, eNv
			do  d=1, eNv
				sum3=sum3 +lrt1(i,c)*lrt1(k,d)*lrint(j,k,c+eNo,d+eNo)
			enddo
			enddo
			enddo
			
			sum4=0
			do  k=1, eNo
			do  d=1, eNv
			do  c=1, eNv
				sum4=sum4 +lrt2(i,k,d,c)*lrint(j,k,d+eNo,c+eNo)
			enddo
			enddo
			enddo

			Fij(i,j)=  lrf(i,j)+sum1+sum2+sum3+sum4*real(0.5,rglu)
		enddo
		enddo
		!$omp end parallel

		write (lrsu,'(1X,A)') trim(ttime())//' Part 10/10'
		!$omp parallel default(shared) private(i,a,k,d,sum1)
		!$omp do
		do  i=1, eNo
		do  a=1, eNv

			sum1=0
			do  k=1, eNo
			do  d=1, eNv
				sum1=sum1 +lrt1(k,d)*lrint(i,k,a+eNo,d+eNo)
			enddo
			enddo

			Fia(i,a)= lrf(i,a+eNo)+sum1
		enddo
		enddo
		!$omp end parallel		

		return
		end subroutine prepint1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine prepint2

		implicit none

		integer(kind=iglu) :: i,j,k,l,a,b,c,d
		real(kind=rglu)    :: sum


		!$omp parallel default(shared) private(k,c,i,j,sum)
		!$omp do
		do  i=1, eNo
		do  j=1, eNo

			sum=0
			do  k=1, eNo
			do  c=1, eNv
				sum=sum +lru1(k,c)*ai6(j,k,i,c)
			enddo
			enddo

			at1(i,j)=sum
		enddo
		enddo
		!$omp end parallel
		
		!$omp parallel default(shared) private(k,c,a,b,sum)
		!$omp do
		do  a=1, eNv
		do  b=1, eNv

			sum=0
			do  k=1, eNo
			do  c=1, eNv
				sum=sum +lru1(k,c)*ai7(k,a,b,c)
			enddo
			enddo

			at2(a,b)=sum
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private(k,c,i,j,d,sum)
		!$omp do
		do  i=1, eNo
		do  j=1, eNo

			sum=0
			do  k=1, eNo
			do  c=1, eNv
			do  d=1, eNv
				sum=sum +lru2(i,k,c,d)*lrint(j,k,c+eNo,d+eNo)
			enddo
			enddo
			enddo

			at3(i,j)= sum*real(0.5,rglu)
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private(a,b,k,l,c,sum)
		!$omp do
		do  a=1, eNv
		do  b=1, eNv

			sum=0
			do  k=1, eNo
			do  l=1, eNo
			do  c=1, eNv
				sum=sum +lru2(k,l,a,c)*lrint(k,l,b+eNo,c+eNo)
			enddo
			enddo
			enddo

			at4(a,b)= sum*real(0.5,rglu)
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private(i,j,k,a,b,sum)
		!$omp do
		do  i=1, eNo
		do  j=1, eNo

			sum=0
			do  k=1, eNo
			do  a=1, eNv
			do  b=1, eNv
				sum=sum +lru2(i,k,a,b)*lrt2(j,k,a,b)
			enddo
			enddo
			enddo

			at5(i,j)= sum*real(0.5,rglu)
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private(a,b,i,j,c,sum)
		!$omp do
		do  a=1, eNv
		do  b=1, eNv

			sum=0
			do  i=1, eNo
			do  j=1, eNo
			do  c=1, eNv
				sum=sum +lru2(i,j,a,c)*lrt2(i,j,b,c)
			enddo
			enddo
			enddo

			at6(a,b)= sum*real(0.5,rglu)
		enddo
		enddo
		!$omp end parallel

		return
		end subroutine prepint2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine c1
		implicit none

		integer(kind=iglu) :: i,j,k,a,b,c
		real(kind=rglu)    :: rez,sum1,sum2,sen1,sen2


		!$omp parallel default(shared) private(i,a,j,b,k,c,rez,sum1,sum2,sen1,sen2)
		!$omp do
		do i=1,eNo
		do a=1,eNv

			rez=0
      
			sum1=0; sum2=0
			sen1=0; sen2=0

			do  b=1, eNv
				sum1=sum1 +Fab(a,b)*lru1(i,b)
			enddo

			do  j=1, eNo
				sum1=sum1 -Fij(i,j)*lru1(j,a)
			enddo

			do  b=1, eNv
			do  j=1, eNo
				sum1=sum1 -ai1(i,b,j,a)*lru1(j,b) 
			enddo
			enddo

			do  b=1, eNv
			do  j=1, eNo
				sum2=sum2 +Fia(j,b)*lru2(i,j,a,b)
			enddo
			enddo

			do  b=1, eNv
			do  k=1, eNo
			do  j=1, eNo
				sen1=sen1 +ai6(j,k,i,b)*lru2(j,k,a,b)
			enddo
			enddo
			enddo

			do  b=1, eNv
			do  c=1, eNv
			do  j=1, eNo
				sen2=sen2 +ai7(j,a,b,c)*lru2(i,j,b,c)
			enddo
			enddo
			enddo
			sum2=sum2 -real(0.5,rglu)*(sen1+sen2)

			x1(i,a)= sum1 + sum2
		enddo
		enddo
		!$omp end parallel

		return
		end subroutine c1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine c2

		implicit none

		integer(kind=iglu) :: i,j,k,l,a,b,c,d,a1,b1,i1,j1
		real(kind=rglu)    :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,rez


		x2=0

		do i=1, eNo-1
		do j=i+1, eNo

		!$omp parallel default(shared) private (a,b,k,c,l,d,rez,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9)
		!$omp do
		do a=1, eNv-1
		do b=a+1, eNv
	
			rez=0
			if (btest(a+b,0).NE.btest(i+j,0)) cycle

			sum1=0; sum4=0
			do c=1, eNv
				sum1=sum1 +lru1(i,c)*ai3(j,c,a,b)&
		         &        -lru1(j,c)*ai3(i,c,a,b)
		         
				sum4=sum4 +fab(b,c)*lru2(i,j,a,c)&
		         &        -fab(a,c)*lru2(i,j,b,c)

			enddo

			do k=1, eNo
				sum1=sum1 +lru1(k,a)*ai2(i,j,k,b)&
		         &        -lru1(k,b)*ai2(i,j,k,a)

				sum4=sum4 +fij(i,k)*lru2(j,k,a,b)&
		         &        -fij(j,k)*lru2(i,k,a,b)
			enddo

			sum2=0; sum8=0
			do l=1, eNo
				sum2=sum2 +at1(i,l)*lrt2(j,l,a,b)&
		         &        -at1(j,l)*lrt2(i,l,a,b)

				sum8=sum8 +at3(i,l)*lrt2(j,l,a,b)-at3(j,l)*lrt2(i,l,a,b)
			enddo

			sum3=0; sum9=0
			do d=1, eNv
				sum3=sum3 +at2(a,d)*lrt2(i,j,b,d)&
		         &        -at2(b,d)*lrt2(i,j,a,d)
		         
				sum9=sum9 +at4(a,d)*lrt2(i,j,b,d)-at4(b,d)*lrt2(i,j,a,d)
			enddo

			sum5=0
			do k=1, eNo
			do c=1, eNv
				sum5=sum5 +lru2(j,k,a,c)*ai1(i,c,k,b)&
				&         -lru2(i,k,a,c)*ai1(j,c,k,b)&
		        &         +lru2(i,k,b,c)*ai1(j,c,k,a)&
		        &         -lru2(j,k,b,c)*ai1(i,c,k,a)
			enddo
			enddo

			sum6=0
			do k=1, eNo
			do l=1, eNo
				sum6=sum6 +lru2(k,l,a,b)*ai4(i,j,k,l)
			enddo
			enddo

			sum7=0
			do c=1, eNv
			do d=1, eNv
				sum7=sum7 +lru2(i,j,c,d)*ai5(a,b,c,d)
			enddo
			enddo

			rez=sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8+sum9

			x2(i,j,a,b)= rez
			x2(i,j,b,a)=-rez
			x2(j,i,a,b)=-rez
			x2(j,i,b,a)= rez
		enddo
		enddo

		!$omp end parallel		
		enddo
		enddo

		return
		end subroutine c2

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine normu
		implicit none

		integer(kind=iglu) :: a,b,i,j
		real(kind=rglu)    :: s1,s2
      
      
		s1=0
		!$omp parallel default(shared) private (i,a) reduction(+:s1)
		!$omp do
		do i = 1,eNo
			do a = 1,eNv
				s1=s1+lru1(i,a)**2
			enddo
		enddo
		!$omp end parallel
            
		s2=0
		!$omp parallel default(shared) private (i,j,a,b) reduction(+:s2)
		!$omp do	
		do i = 1,eNo-1
			do j = i+1,eNo
				do a = 1,eNv-1
					do b = a+1,eNv
						s2=s2+lru2(i,j,a,b)**2
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel
      
		s1=real(1,rglu)/sqrt(s1+s2)

		!$omp parallel default(shared) private (i,a)
		!$omp do
		do i = 1,eNo
			do a = 1,eNv
				lru1(i,a)=lru1(i,a)*s1
			enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,j,a,b,s2)
		!$omp do
		do i = 1,eNo-1
			do j = i+1,eNo
				do a = 1,eNv-1
					do b = a+1,eNv
						s2=lru2(i,j,a,b)*s1
						lru2(i,j,a,b)= s2
						lru2(j,i,a,b)=-s2
						lru2(i,j,b,a)=-s2
						lru2(j,i,b,a)= s2
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel

		return
		end subroutine normu

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine energylr
		implicit none

		integer(kind=iglu) :: i,j,a,b
		real(kind=rglu)    :: s1,s2


		s1=0
		s2=0
		!$omp parallel default(shared) private (i,a) reduction(+:s1,s2)
		!$omp do
		do i=1,eNo
			do a=1,eNv
				s1=s1+x1(i,a)*lru1(i,a)
				s2=s2+lru1(i,a)*lru1(i,a)
			enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,j,a,b) reduction(+:s1,s2)
		!$omp do
		do i=1,eNo-1
			do j=i+1,eNo
				do a=1,eNv-1
					do b=a+1,eNv
						s1=s1+x2(i,j,a,b)*lru2(i,j,a,b)
						s2=s2+lru2(i,j,a,b)*lru2(i,j,a,b)
					enddo
				enddo
			enddo
		enddo
		!$omp end parallel

		lrOmega=s1/s2
		return
		end subroutine energylr

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine stepu1(accur1)
		implicit none

		integer(kind=iglu) :: i,a
		real(kind=rglu)    :: Ax,sss,accur1


		!$omp parallel default(shared) private (i,a)
		!$omp do
		do i = 1,eNo
		do a = 1,eNo
			x1(i,a)=x1(i,a)-lrOmega*lru1(i,a)
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,a)
		!$omp do
		do i=1,eNo
		do a=1,eNv
			Ax=lru1(i,a)-stepLR*x1(i,a)
			!if (Ax*lru1(i,a).LT.0) then
			!	lru1(i,a)=0
			!else
				lru1(i,a)=Ax
			!endif
		enddo
		enddo
		!$omp end parallel
		
		accur1=maxval(abs(x1))

		return
		end subroutine stepu1

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine stepu2(accur2)

		implicit none

		integer(kind=iglu) :: i,j,a,b
		real(kind=rglu)    :: Ax,sss,accur2


		!$omp parallel default(shared) private (i,j,a,b)
		!$omp do
		do i=1,eNo
		do a=1,eNv
		do j=1,eNo
		do b=1,eNv
			x2(i,j,a,b)=x2(i,j,a,b)-lrOmega*lru2(i,j,a,b)
		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		!$omp parallel default(shared) private (i,j,a,b)
		!$omp do
		do i=1,eNo
		do a=1,eNv
		do j=1,eNo
		do b=1,eNv
			Ax=lru2(i,j,a,b)-stepLR*x2(i,j,a,b)

			!if (Ax*lru2(i,j,a,b).LT.0) then
			!	lru2(i,j,a,b)=0
			!else
				lru2(i,j,a,b)=Ax
			!endif

		enddo
		enddo
		enddo
		enddo
		!$omp end parallel

		accur2=maxval(abs(x2))

		return
		end subroutine stepu2

	end module lrccsdModule
