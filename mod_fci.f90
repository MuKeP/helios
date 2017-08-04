	module fciModule

	use glob    , only: true,false,rspu,rglu,iglu,lglu,void
	use math    , only: tred4
	use fcontrol

	private

	character (len=*), parameter :: ciVersion='4.400'       !5
	character (len=*), parameter :: ciDate   ='2017.07.07'  !10
	character (len=*), parameter :: ciAuthor ='Vladimir V. Ivanov'

	! glu = real for global use  (storage and non-accuracy-demanding procedures).
	! spu = real for special use (accuracy-demanding procedures).
	!integer*4, parameter :: rglu=r8kind, rspu=r8kind
	!integer*4, parameter :: iglu=i4kind, ispu=i8kind
	!integer*4, parameter :: lglu=l1kind

	!real(kind=rglu), parameter :: gluzero=real(0,rglu)
	!real(kind=rspu), parameter :: spuzero=real(0,rspu)

	! default parameter values
	integer(kind=iglu)     :: fciNSteps=7,fciMaxiter=1000,fciNStates=1,iount
	real(kind=rglu)        :: fciTol=real(1d-12,rglu),sparseThreshold=real(1d-12,rglu)

	real(kind=rglu)        :: fciEnergy
	integer(kind=iglu)     :: N,Nocca,Noccb,NEL,KF
	integer(kind=iglu)     :: NN,NN1

	real(kind=rglu), allocatable :: fciHoldStateEnergy(:)
	real(kind=rglu), allocatable :: fciHoldStateOrthog(:)
	real(kind=rglu), allocatable :: fciHoldStateVector(:,:)
	real(kind=rglu), allocatable :: fciHoldStateRDM   (:,:,:)

	real(kind=rglu), allocatable :: fH(:,:),fG(:,:) ! Hamiltonian, Gamma
	real(kind=rglu), allocatable :: X(:),Z(:)     ! D,GRD
	real(kind=rglu), allocatable :: GAM(:)        ! ~integrals
	
	! arrays for Knuth sparse algorythm
	real(kind=rglu)   , allocatable :: AN(:)
	integer(kind=iglu), allocatable :: JPA(:),NEXTR(:),JR(:),IA(:,:)
	real(kind=rglu)   , allocatable :: fR(:,:),TAU(:,:) ! RDM, supportive array

	! arrays for iteration procedure
	real(kind=rglu), allocatable :: eiMas(:,:),eiVec(:,:),eiVal(:)
	real(kind=rglu), allocatable :: OV(:),suppMass(:,:)
	logical(kind=lglu)           :: doRDM,doOffDiagonal

	public :: ciVersion,ciDate,ciAuthor,fH,fG,&
	&         setFCI_Params,ciFinalize,doFCI,fciHoldStateEnergy,fciHoldStateRDM

	contains

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer(kind=iglu) function setFCI_Params(uN,tol,iters,nsteps,nstates,udoRDM,udoOffDiagonal,uio)
		implicit none

		integer(kind=iglu)           :: iters,nsteps,nstates,uN
		real   (kind=rglu)           :: tol
		logical(kind=lglu), optional :: udoRDM,udoOffDiagonal
		integer(kind=iglu), optional :: uio


		doRDM=false;         if (present(udoRDM)) doRDM=udoRDM
		iount=6;             if (present(uio)) iount=uio
		doOffDiagonal=false; if (present(udoOffDiagonal)) doOffDiagonal=udoOffDiagonal
		N         =uN
		fciTol    =tol
		fciMaxiter=iters
		fciNSteps =nsteps
		fciNStates=nstates

		if (doOffDiagonal.AND.(.NOT.doRDM)) doOffDiagonal=false

		NOCCA=N/2; NOCCB=N/2   !!!!!!!!!!!!
		
		if (N.GE.16) then
			if (fciNSteps.GT.8) then
				write (iount,'(4X,A,i2,A)') 'Too much memory is necessary for',fciNSteps,' steps. It will be set 4.'
				fciNSteps=8
			endif
		endif

		if (allocated(fG)) deallocate (fG,fH)

		allocate (fG(N,N),fH(N,N))

		setFCI_Params=0; return
		end function setFCI_Params

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine doFCI

		implicit none

		integer(kind=iglu) :: I,J,IJ,II ! big vals
		integer(kind=iglu) :: iter,istate
		real   (kind=rglu) :: EEE,AN0,AMA,ecore,tbefore,tafter,totbefore,totafter


		if (N.GT.16) then; write (iount,'(A)') '   Too big size.'; stop; endif

		call MULTI

		IJ=0
		do I=1,NN
		do J=1,I
			IJ=IJ+1
		enddo
		enddo

		allocate (GAM(IJ),X(IJ),Z(IJ),TAU(N,N)); GAM=0; X=0; Z=0; TAU=0
		allocate (suppMass(NN1,fciNSteps),OV(NN1)); suppMass=0; OV=0
		allocate (eiMas(fciNSteps,fciNSteps),eiVec(fciNSteps,fciNSteps),eiVal(fciNSteps))

		eiMas=0; eiVec=0; eiVal=0

		if (.NOT.allocated(fciHoldStateEnergy)) then
			allocate(fciHoldStateEnergy(    fciNStates))
			allocate(fciHoldStateOrthog(    fciNStates))
			allocate(fciHoldStateVector(IJ ,fciNStates))
			allocate(fciHoldStateRDM   (N,N,fciNStates))
		endif

		fciHoldStateEnergy=0
		fciHoldStateOrthog=0
		fciHoldStateVector=0
		fciHoldStateRDM=   0

		ecore=0
		do i= 1,N-1
		do j= i+1,N
			ecore=ecore+fG(I,J)
		enddo
		enddo

		do istate = 1,fciNStates
			call SUPER
			call FSPARS

			IJ=0
			do I=1,NN
			do J=1,I
			   IJ=IJ+1
			   X(IJ)=0
			   if(I.EQ.J) X(IJ)=1
			enddo
			enddo

			if (istate.GE.3) then
				do I=1,NN
					II=index(i,i)
					X(II)=(-1.)**i
				enddo

			endif
			AN0=ANOR(X,NN,NN1)

			if (istate.EQ.1) then
				if (fciNStates.NE.1) write (iount,100)
			else
				write (iount,101) istate-1
			endif

			do iter = 1,fciMaxiter
				call operg(X,Z,EEE,AMA)
				if(AMA.LT.fciTol*real(2**(istate-1),rglu) ) exit
				call stepms(X,Z,EEE)
				call ortog(istate-1)

				write(iount,102) iter, AMA, EEE
			enddo

			write(iount,103) iter, AMA, EEE

			fciHoldStateEnergy(  istate)=EEE+ecore
			fciHoldStateVector(:,istate)=X

			if (doRDM) call density(istate)

			deallocate (AN,JPA,NEXTR,JR)
		enddo

		if (fciNStates.GT.2) call resortStates

100		format (/18X,'   ==========Ground state==========    '/)
101		format (/18X,'==========Excited state # ',i1,'=========='/)
102		format(4X,' Iteration ',I4,' Accuracy ',ES10.4,' Energy=',F20.14)
103		format(4X,' Completed ',I4,' Accuracy ',ES10.4,' Energy=',F20.14)

		deallocate (GAM,X,Z,IA,suppMass,eiMas,eiVec,OV,eiVal,TAU)

		return
		end subroutine doFCI

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine resortStates

		implicit none

		integer(kind=iglu) :: K
		integer(kind=iglu) :: i,j,l,d,swcount,glcount
		real   (kind=rglu) :: swEnergy,swRDM


		glcount=0
		do
			swcount=0
			do i = 2,fciNStates-1
				j=i+1
				if (fciHoldStateEnergy(j).LT.fciHoldStateEnergy(i)) then
					swcount=swcount+1

					swEnergy=fciHoldStateEnergy(j)
					fciHoldStateEnergy(j)=fciHoldStateEnergy(i)
					fciHoldStateEnergy(i)=swEnergy

					do d = 1,N
					do l = 1,N
						swRDM=fciHoldStateRDM(d,l,j)
						fciHoldStateRDM(d,l,j)=fciHoldStateRDM(d,l,i)
						fciHoldStateRDM(d,l,i)=swRDM
					enddo
					enddo

					do K = 1,UBound(X,1)
						X(K)=fciHoldStateVector(K,j)
						fciHoldStateVector(K,j)=fciHoldStateVector(K,i)
						fciHoldStateVector(K,i)=X(K)
					enddo
				endif
			enddo
			glcount=glcount+swcount
			if (swcount.EQ.0) exit
		enddo

		if (glcount.NE.0) write (iount,100)

100		format (4X,'Incorrect order of excited states was fixed.')

		return
		end subroutine resortStates

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine stepms(D,GRD,ALAM)

		implicit none

		integer(kind=iglu) :: ii,inom
		integer(kind=iglu) :: I,K
		real   (kind=rglu) :: D(:), GRD(:),ALAM, TT,AMA,AN0


		suppMass=0
		suppMass(:,1)=D

		eiMas=0
		eiMas(1,1)=ALAM
		eiMas(2,1)=ANOR( GRD, NN, NN1)
		eiMas(1,2)=eiMas(2,1)

		do K = 1,UBound(GRD,1)
			suppMass(K,2)=GRD(K)
		enddo
		do ii=2, fciNSteps
			inom=ii-1
			do K = 1,UBound(GRD,1)
				GRD(K)=suppMass(K,ii)
			enddo
			call OPERG(GRD,OV,TT,AMA)
			eiMas(ii,ii)=TT

			do K = 1,UBound(GRD,1)
				GRD(K)=suppMass(K,inom)
			enddo

			do I=1,NN1
				OV(I)=OV(I)-eiMas(inom,ii)*GRD(I)
			enddo

			if(ii.EQ.fciNSteps) exit

			eiMas(ii,ii+1)=ANOR(OV,NN,NN1)
			eiMas(ii+1,ii)=eiMas(ii,ii+1)
			do K = 1,UBound(OV,1)
				suppMass(K,ii+1)=OV(K)
			enddo
		enddo

		call tred4(eiMas,eiVec,eiVal,fciNSteps,real(1.d-100,rspu),real(1.d-300,rspu))

		D=0
		do ii = 1,fciNSteps
			OV=suppMass(:,ii)
			do I=1,NN1
				D(I)=D(I)+eiVec(ii,1)*OV(I)
			enddo
		enddo
		AN0=ANOR(D,NN,NN1)

		return
		end subroutine stepms

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine operg(D,GRD,EEE,AMA)

		implicit none

		real   (kind=rglu) :: D(:),GRD(:)
		real   (kind=rglu) :: EEE, AMA,SUM1,SUM2,SEN
		integer(kind=iglu) :: I,J,IJ,IE,NX,L,LJ,IL,KK


		IJ=0
		EEE=0
		do I=1,NN

			!$omp parallel default(shared) private(J,IJ,IE,SUM1,L,LJ,NX,SUM2,IL,SEN) reduction(+:EEE)
			!$omp do
			do J=1,I
				IJ=INDEX(I,J)
				
				IE=JR(I)
				SUM1=0
				do
					L=JPA(IE)
					LJ=INDEX(L,J)
					SUM1=SUM1+AN(IE)*D(LJ)
					NX=NEXTR(IE)
					if(NX.EQ.0) exit
					IE=NX
				enddo

				IE=JR(J)
				SUM2=0
				do
					L=JPA(IE)
					IL=INDEX(I,L)
					SUM2=SUM2+AN(IE)*D(IL)
					NX=NEXTR(IE)
					if(NX.EQ.0) exit
					IE=NX
				enddo

				GRD(IJ)=SUM1+SUM2+D(IJ)*GAM(IJ)

				SEN=GRD(IJ)*D(IJ)
				if(I.NE.J) SEN=SEN+SEN
				EEE=EEE+SEN
			enddo
			!$omp end parallel
		enddo

		AMA=0
		do i=1,NN1
			GRD(i)=GRD(i)-EEE*D(i)
			if (abs(GRD(i)).GT.AMA) AMA=abs(GRD(i))
		enddo

		return
		end subroutine operg

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine density(state)

		implicit none

		integer(kind=iglu) :: I,J
		integer(kind=iglu) :: state,kk,ll
		real   (kind=rglu) :: T

		fciHoldStateRDM(:,:,state)=0
		write (iount,'(/4X,A,L/)') 'RDM preparation. Do off-diagonal elemnts:',doOffDiagonal
		do I=1,NN
			if (mod(I,200).EQ.0) write(iount,'(1X,F5.1\)') 100.*I/NN
			do J=I,NN
				T=VSPS(I,J)
				call TAY(I,J)
				if (doOffDiagonal) then
					do kk=1,N
					do ll = kk,N
						fciHoldStateRDM(kk,ll,state)=fciHoldStateRDM(kk,ll,state)+real(2,rglu)*T*TAU(kk,ll)
					enddo
					enddo
					do kk = 1,N
					do ll = kk,N
						fciHoldStateRDM(ll,kk,state)=fciHoldStateRDM(kk,ll,state)
					enddo
					enddo
				else
					do kk=1,N
						fciHoldStateRDM(kk,kk,state)=fciHoldStateRDM(kk,kk,state)+real(2,rglu)*T*TAU(kk,kk)
					enddo
				endif
			enddo
		enddo
		write (iount,*)

		return 
		end subroutine density

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine TAY(I,J)

		implicit none

		integer(kind=iglu) :: K,L,I,J,ISI,ISJ,IP,JP,IAA,IBB,IIA,IIB,KOL,M


		TAU=0; KOL=0
		DO K=1,KF
		DO L=1,KF
			IF(IA(I,K).EQ.IA(J,L)) KOL=KOL+1
		enddo
		enddo
 
		IF((KF-KOL).EQ.1) then
			ISI=0; ISJ=0
			IP=1;  JP=1
			DO K=1,KF
				ISI=ISI+IA(I,K)
				ISJ=ISJ+IA(J,K)
				IP=IP*IA(I,K)
				JP=JP*IA(J,K)
			enddo
			IAA=int(IP,kind=8)*(ISI-ISJ)/(IP-JP)
			IBB=int(JP,kind=8)*(ISI-ISJ)/(IP-JP)
			DO K=1,KF
				IF(IAA.EQ.IA(I,K)) IIA=K
				IF(IBB.EQ.IA(J,K)) IIB=K
			enddo
			TAU(IAA,IBB)=-1.**(IIA+IIB)
		endif
		if(KOL.EQ.KF) then
 			do M=1,KF
				K=IA(I,M)
				TAU(K,K)=1
			enddo
		endif

		return
		end subroutine TAY

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		real*8 function ANOR(D,NN,NN1)

		implicit none

		real   (kind=rglu) :: D(:)
		real   (kind=rglu) :: SS,SEN
		integer(kind=iglu) :: IJ,I,J,NN,NN1

		SS=0
		IJ=0
		do I=1,NN
			do J=1,I
				IJ=IJ+1
				SEN=D(IJ)*D(IJ)
				if(I.NE.J) SEN=SEN+SEN
				SS=SS+SEN
			enddo
		enddo

		ANOR=SQRT(SS)
		SS=1./ANOR

		do I=1,NN1
			D(I)=SS*D(I)
		enddo

		return
		end function ANOR

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine fspars

		implicit none

		real   (kind=rglu) :: PRSNT
		integer(kind=iglu) :: I,J,IJ,NJ


		do I=1,NN
			NJ=0
			do J=1,NN
				IJ=INDEX(I,J)
				if( abs(Z(IJ)).LT.sparseThreshold) cycle
				NEL=NEL+1
			enddo
		enddo

		allocate(AN(NEL),JPA(NEL),NEXTR(NEL),JR(NN))


		AN=0; JPA=0; NEXTR=0; JR=0

      
		NEL=0
		do I=1,NN
			NJ=0
			do J=1,NN
				IJ=INDEX(I,J)
				if( abs(Z(IJ)).LT.sparseThreshold) cycle
				NEL=NEL+1
				NJ=NJ+1
				AN(NEL)=Z(IJ)

				JPA(NEL)=J
				if(NJ.NE.1) then
					NEXTR(NEL-1)=NEL
					cycle
				endif
				JR(I)=NEL
			enddo
		enddo
		PRSNT=real(100,rglu)-real(100,rglu)*NEL/NN*NN

		return
		end subroutine fspars

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine super

		implicit none

		real   (kind=rglu) :: GAM1,GAM2
		integer(kind=iglu) :: I,J,K,L,M,IJ,ISS,JSS,LAA,LBB,K1,KOL,ILA,ILB,LL,KK,M1,L1
		integer(kind=iglu) :: IP,JP


		GAM=0; X=0; Z=0

		K1=KF-1
		IJ=0
		do I=1,NN
		do J=1,I
			IJ=IJ+1
			Z(IJ)=0
			KOL=0
			do M=1,KF
			do L=1,KF
				if (IA(I,M).EQ.IA(J,L)) KOL=KOL+1
			enddo
			enddo
 
			if (KOL+1.EQ.KF) then

				ISS=0; JSS=0; IP=1; JP=1

				do L=1,KF
					ISS=ISS+IA(I,L)
					JSS=JSS+IA(J,L)
					IP=IP*IA(I,L)
					JP=JP*IA(J,L)
				enddo

				LAA=int(IP,kind=8)*(ISS-JSS)/(IP-JP)
				LBB=int(JP,kind=8)*(ISS-JSS)/(IP-JP)

				if (dabs(fH(LAA,LBB)).GE.sparseThreshold) then
					do L=1,KF
						if(LAA.eq.IA(I,L)) ILA=L
						if(LBB.eq.IA(J,L)) ILB=L
					enddo
					LL=ILA+ILB
					Z(IJ)=-1.**LL*fH(LAA,LBB)
				endif
			endif
      
			if(KOL.EQ.KF) then
				do K=1,KF
					KK=IA(I,K)
					Z(IJ)=Z(IJ)+fH(KK,KK)
				enddo
			endif
 
			GAM1=0
			if(K1.NE.0) then
				do M=1,K1
				do L=M+1,KF
					GAM1=GAM1+fG(IA(I,M),IA(I,L))+fG(IA(J,M),IA(J,L))
				enddo
				enddo
			endif

			GAM2=0
			do M=1,KF
			do L=1,KF
			   M1=IA(I,M)
			   L1=IA(J,L)
			   GAM2=GAM2+fG(M1,L1)
			enddo
			enddo

			GAM(IJ)= GAM1+GAM2
		enddo
		enddo

		return
		end subroutine super

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine multi

		implicit none

		integer(kind=iglu)              :: I,K1,L,P
		integer(kind=iglu), allocatable :: IND(:)
		integer(kind=iglu)              :: fileid


		KF=(NOCCA+NOCCB)/2

		allocate (IND(KF))

		fileid=fcNewID()

		open (fileid,status='scratch',form='binary')
		K1=KF+1; NN=1
		do I=1,KF
			IND(I)=I
			write (fileid) int(1,kind=iglu),I,I !IA(1,I)=I
		enddo

		P=KF
		do
			if(P.LT.1) exit
			if(IND(KF).EQ.N) then
				P=P-1
			else
				P=KF
			endif

			if(P.LT.1) exit

			I=K1
			do
				I=I-1
				if(I.LT.P) exit
				IND(I)=IND(P)+I-P+1
			enddo
 			NN=NN+1
			do L=1,KF
				write (fileid) NN,L,IND(L) !IA(NN,L)=IND(L)
 			enddo
		enddo

		if (allocated(IA)) deallocate (IA)
		allocate (IA(NN,KF))
		rewind(fileid)

		do
			if (eof(fileid)) exit
			read (fileid) NN,L,K1
			IA(NN,L)=K1
		enddo
		close (fileid)
		void=fcNullID(fileid)

		NN1=NN*(NN+1)/2

		deallocate (IND)
      
		return
		end subroutine multi

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		integer(kind=iglu) function index(I,J)
		implicit none
		integer(kind=iglu) :: I,J


		if(I.GT.J) then
			index=I*(I-1)/2+J
		else
			index=J*(J-1)/2+I
		endif

		return
		end function index

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		real(kind=rglu) function VSPS(I,J)
		implicit none
		real(kind=rglu)    :: T
		integer(kind=iglu) :: I,J,K,IK,JK


		T=0
		do K=1,NN
			IK=INDEX(I,K)
			JK=INDEX(J,K)
			T=T+X(IK)*X(JK)
		enddo
		VSPS=T; return
		end function VSPS

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine ortog(state)

		implicit none

		real   (kind=rglu) :: spr,sen
		integer(kind=iglu) :: state,k
		integer(kind=iglu) :: I,J,IJ


		if (state.EQ.0) return

		do k = 1,state
			IJ=0; spr=0
			do I=1,NN
			do J=1,I
				IJ=IJ+1
				sen=X(IJ)*fciHoldStateVector(IJ,k)
				if (I.NE.J) sen=sen+sen
				spr=spr+sen
			enddo
			enddo
			fciHoldStateOrthog(k)=spr
		enddo

!		spr=maxval(abs(fciHoldStateOrthog(1:state))); if (spr.LT.fciTol) return
!		write (iount,"(5X,'Ortogonality=',1X,ES12.4)") spr

		do k = 1,state
			do I=1,NN1
			   X(I)=X(I)-fciHoldStateOrthog(k)*fciHoldStateVector(I,k)
			enddo
		enddo

        return
        end subroutine ortog

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

		subroutine ciFinalize

		implicit none


		if (allocated(fciHoldStateEnergy)) deallocate (fciHoldStateEnergy)
		if (allocated(fciHoldStateRDM))    deallocate (fciHoldStateRDM)
		if (allocated(fciHoldStateVector)) deallocate (fciHoldStateVector)
		if (allocated(fciHoldStateOrthog)) deallocate (fciHoldStateOrthog)

		return
		end subroutine ciFinalize

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

	end module fciModule
