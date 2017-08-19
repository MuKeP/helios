	subroutine showPolarizability(method,ou,width)

	use glob     , only: rspu,r16kind,iglu,rglu,lglu,void,true,false,uch,uchGet
	use txtParser, only: tpLowerCase,tpUpperCase,tpFill,tpAdjustc
	use printmod , only: prStrByVal
	use hdb      , only: HartreeEnergy,BohrRadius,dipoleToDeby,gammaToesu
	use hdb      , only: Et,MEt,MMEt,mol,polarizbd
	use derivat  , only: deShareParams,deLagDeriv,deLSMDeriv,deFinalize

	implicit none

	character (len=*) , intent(in) :: method
	integer(kind=iglu), intent(in) :: ou,width

	real(kind=rglu), parameter :: sigmaBondAlpha=3.48d0,CCsigmaBondGamma=220,CHsigmaBondGamma=240

	character (len=1)  ::    sX,   sY
	character (len=2)  ::   sXX,  sYY
	character (len=3)  ::  sXXX, sYYY, sXYY, sYXX
	character (len=4)  :: sXXXX,sYYYY,sXXYY,sYYXX
	character          :: sq*2,scale*3

	integer(kind=iglu) :: i,j,k,Nel,Nat,pp,ppp,sta,sto

	real   (kind=rglu) :: ProEt(-(polarizbd%nPoints-1)/2:(polarizbd%nPoints-1)/2)

	real   (kind=rglu) :: ddx(4),aax(4),bbx(3),ggx(3),bstDeriv(4,10)
	real   (kind=rglu) :: dx,dy,dz,dmod
	real   (kind=rglu) :: axx,ayy,azz,avA
	real   (kind=rglu) :: bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod
	real   (kind=rglu) :: gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averG

	real   (kind=rglu) :: sigmaContributionAlpha,sigmaContributionGamma

	real   (kind=rglu) :: esu,PDip,PPol,PHpl1,PHpl2
	real   (kind=rglu) :: transition(0:8)


	Nel=mol%nEls; Nat=mol%uniqueAtoms; scale=uchGet(polarizbd%scales)
	select case ( len_trim(  uchGet(polarizbd%scales)  ) )
		case (1); sq='ss'
		case (2); sq='ds'
		case (3); sq='ts'
	end select

	sta=-(polarizbd%nPoints-1)/2
	sto= (polarizbd%nPoints-1)/2
	sigmaContributionAlpha=  sigmaBondAlpha*(mol%uniqueBonds+mol%chBonds)
	sigmaContributionGamma=CCsigmaBondGamma*mol%uniqueBonds+CHsigmaBondGamma*mol%chBonds

	ppp=5; pp=len(method)+2+ppp*2

	transition(0)=1; transition(1)=-dipoleToDeby
	do i = 2,UBound(transition,1)
		transition(i)=-HartreeEnergy**(i-1)/BohrRadius**(i)
	enddo
	PDip=-dipoleToDeby; PPol=transition(2); PHpl1=transition(3); PHpl2=transition(4)
	esu=gammaToesu

	void=deShareParams(polarizbd%nPoints,polarizbd%derivStep,MMEt)

	write (ou,111) tpFill(pp,'~'),tpFill(ppp,'~'),method,tpFill(ppp,'~'),tpFill(pp,'~')

	select case (sq)

		case ('ss')
			sX=scale(1:1); sXX=sX//sX; sXXX=sX//sXX; sXXXX=sXX//sXX
			ProEt=Et
			do i = sta,sto
				write (ou,101) i,Et(i)
			enddo

			select case (polarizbd%nPoints)

				case (9)
					ddx(1)=deLagDeriv(3,tpLowerCase(sX))*PDip
					ddx(2)=deLagDeriv(5,tpLowerCase(sX))*PDip
					ddx(3)=deLagDeriv(7,tpLowerCase(sX))*PDip
					ddx(4)=deLagDeriv(9,tpLowerCase(sX))*PDip

					aax(1)=deLagDeriv(3,tpLowerCase(sXX))*PPol
					aax(2)=deLagDeriv(5,tpLowerCase(sXX))*PPol
					aax(3)=deLagDeriv(7,tpLowerCase(sXX))*PPol
					aax(4)=deLagDeriv(9,tpLowerCase(sXX))*PPol

					bbx(1)=deLagDeriv(5,tpLowerCase(sXXX))*PHpl1
					bbx(2)=deLagDeriv(7,tpLowerCase(sXXX))*PHpl1
					bbx(3)=deLagDeriv(9,tpLowerCase(sXXX))*PHpl1

					ggx(1)=deLagDeriv(5,tpLowerCase(sXXXX))*PHpl2
					ggx(2)=deLagDeriv(7,tpLowerCase(sXXXX))*PHpl2
					ggx(3)=deLagDeriv(9,tpLowerCase(sXXXX))*PHpl2

					write (ou,200) trim(method),sX
					write (ou,201) 
					write (ou,211) ddx(1),ddx(2),ddx(3),ddx(4)
					write (ou,212) aax(1),aax(2),aax(3),aax(4)
					write (ou,281) aax(4)/3,aax(4)/(3*Nel)
					if (Nat.NE.Nel) write (ou,282) aax(4)/(3*Nat)
					write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
					write (ou,213) bbx(1),bbx(2),bbx(3)
					write (ou,214) ggx(1),ggx(2),ggx(3),ggx(1)*esu,ggx(2)*esu,ggx(3)*esu
					write (ou,283) ggx(3)/5,ggx(3)/(5*Nel)
					if (Nat.NE.Nel) write (ou,284) ggx(3)/(5*Nat)
					write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
					void=lsmInterpolation(1)
				
				case (7)
					ddx(1)=deLagDeriv(3,tpLowerCase(sX))*PDip
					ddx(2)=deLagDeriv(5,tpLowerCase(sX))*PDip
					ddx(3)=deLagDeriv(7,tpLowerCase(sX))*PDip

					aax(1)=deLagDeriv(3,tpLowerCase(sXX))*PPol
					aax(2)=deLagDeriv(5,tpLowerCase(sXX))*PPol
					aax(3)=deLagDeriv(7,tpLowerCase(sXX))*PPol

					bbx(1)=deLagDeriv(5,tpLowerCase(sXXX))*PHpl1
					bbx(2)=deLagDeriv(7,tpLowerCase(sXXX))*PHpl1

					ggx(1)=deLagDeriv(5,tpLowerCase(sXXXX))*PHpl2
					ggx(2)=deLagDeriv(7,tpLowerCase(sXXXX))*PHpl2

					write (ou,200) trim(method),sX
					write (ou,221) 
					write (ou,231) ddx(1),ddx(2),ddx(3)
					write (ou,232) aax(1),aax(2),aax(3)
					write (ou,281) aax(3)/3,aax(3)/(3*Nel)
					if (Nat.NE.Nel) write (ou,282) aax(3)/(3*Nat)
					write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
					write (ou,233) bbx(1),bbx(2)
					write (ou,234) ggx(1),ggx(2),ggx(1)*esu,ggx(2)*esu
					write (ou,283) ggx(2)/5,ggx(2)/(5*Nel)
					if (Nat.NE.Nel) write (ou,284) ggx(2)/(5*Nat)
					write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
					void=lsmInterpolation(1)

				case (5)
					ddx(1)=deLagDeriv(3,tpLowerCase(sX))*PDip
					ddx(2)=deLagDeriv(5,tpLowerCase(sX))*PDip

					aax(1)=deLagDeriv(3,tpLowerCase(sXX))*PPol
					aax(2)=deLagDeriv(5,tpLowerCase(sXX))*PPol

					bbx(1)=deLagDeriv(5,tpLowerCase(sXXX))*PHpl1

					ggx(1)=deLagDeriv(5,tpLowerCase(sXXXX))*PHpl2

					write (ou,200) trim(method),sX
					write (ou,241) 
					write (ou,251) ddx(1),ddx(2)
					write (ou,252) aax(1),aax(2)
					write (ou,281) aax(2)/3,aax(2)/(3*Nel)
					if (Nat.NE.Nel) write (ou,282) aax(2)/(3*Nat)
					write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
					write (ou,253) bbx(1)
					write (ou,254) ggx(1),ggx(1)*esu
					write (ou,283) ggx(1)/5,ggx(1)/(5*Nel)
					if (Nat.NE.Nel) write (ou,284) ggx(1)/(5*Nat)
					write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
					void=lsmInterpolation(1)

				case (3)
					ddx(1)=deLagDeriv(3,tpLowerCase(sX))*PDip
					aax(1)=deLagDeriv(3,tpLowerCase(sXX))*PPol

					write (ou,200) trim(method),sX
					write (ou,261) 
					write (ou,271) ddx(1)
					write (ou,272) aax(1)
					write (ou,281) aax(1)/3,aax(1)/(3*Nel)
					if (Nat.NE.Nel) write (ou,282) aax(1)/(3*Nat)
					write (ou,500) sigmaContributionAlpha,sigmaBondAlpha

			end select


		case ('ds')
			scale=tpUpperCase(scale)
			sX=Scale(1:1)
			sY=Scale(2:2)
			sXX=sX//sX; sYY=sY//sY
			sXXX=sX//sXX; sYYY=sY//sYY
			sXYY=sX//sYY; sYXX=sY//sXX
			sXXXX=sXX//sXX; sYYYY=sYY//sYY
			sXXYY=sXX//sYY; sYYXX=sYY//sXX
			scale=tpLowerCase(scale)
		
			do i = sta,sto
				do j = sta,sto
					write (ou,102) i,j,MEt(i,j)
				enddo
				write (ou,*)
			enddo

			do ! ~~~~~~~~~ Dipole moment ~~~~~~~~~ !
				if (polarizbd%nPoints.LT.3) exit
				write (ou,301)

				dx=deLagDeriv(3,tpLowerCase(sX))*PDip
				dy=deLagDeriv(3,tpLowerCase(sY))*PDip
				dmod=sqrt(dx**2+dy**2)

				write (ou,311) sX,dx,sY,dy,dmod

				bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,4)=dmod

				if (polarizbd%nPoints.EQ.3) exit

				dx=deLagDeriv(5,tpLowerCase(sX))*PDip
				dy=deLagDeriv(5,tpLowerCase(sY))*PDip
				dmod=sqrt(dx**2+dy**2)

				write (ou,312) sX,dx,sY,dy,dmod

				bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,4)=dmod

				if (polarizbd%nPoints.EQ.5) exit

				dx=deLagDeriv(7,tpLowerCase(sX))*PDip
				dy=deLagDeriv(7,tpLowerCase(sY))*PDip
				dmod=sqrt(dx**2+dy**2)

				write (ou,313) sX,dx,sY,dy,dmod

				bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,4)=dmod

				if (polarizbd%nPoints.EQ.7) exit

				dx=deLagDeriv(9,tpLowerCase(sX))*PDip
				dy=deLagDeriv(9,tpLowerCase(sY))*PDip
				dmod=sqrt(dx**2+dy**2)

				write (ou,314) sX,dx,sY,dy,dmod

				bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,4)=dmod
				exit
			enddo

			do ! ~~~~~~~~~ Polarizability ~~~~~~~~~ !
				if (polarizbd%nPoints.LT.3) exit
				write (ou,321)

				axx = deLagDeriv(3,tpLowerCase(sXX))*PPol
				ayy = deLagDeriv(3,tpLowerCase(sYY))*PPol
				avA = (axx+ayy)/3

				write (ou,331) sXX,axx,sYY,ayy,avA

				bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,4)=avA

				if (polarizbd%nPoints.EQ.3) then
					if (Nel.NE.Nat) then
						write (ou,391) avA/Nat
					endif

					write (ou,392) avA/Nel,myCont(avA,axx,ayy,3)
					write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
					exit
				endif

				axx = deLagDeriv(5,tpLowerCase(sXX))*PPol
				ayy = deLagDeriv(5,tpLowerCase(sYY))*PPol
				avA = (axx+ayy)/3

				write (ou,332) sXX,axx,sYY,ayy,avA

				bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,4)=avA

				if (polarizbd%nPoints.EQ.5) then
					if (Nel.NE.Nat) then
						write (ou,391) avA/Nat
					endif
					write (ou,392) avA/Nel,myCont(avA,axx,ayy,3)
					write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
					exit
				endif

				axx = deLagDeriv(7,tpLowerCase(sXX))*PPol
				ayy = deLagDeriv(7,tpLowerCase(sYY))*PPol
				avA = (axx+ayy)/3

				write (ou,333) sXX,axx,sYY,ayy,avA

				bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,4)=avA

				if (polarizbd%nPoints.EQ.7) then
					if (Nel.NE.Nat) then
						write (ou,391) avA/Nat
					endif

					write (ou,392) avA/Nel,myCont(avA,axx,ayy,3)
					write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
					exit
				endif

				axx = deLagDeriv(9,tpLowerCase(sXX))*PPol
				ayy = deLagDeriv(9,tpLowerCase(sYY))*PPol
				avA = (axx+ayy)/3

				write (ou,334) sXX,axx,sYY,ayy,avA

				bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,4)=avA

				if (Nel.NE.Nat) then
					write (ou,391) avA/Nat
				endif
				write (ou,392) avA/Nel,myCont(avA,axx,ayy,3)
				write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
				exit
			enddo

			do ! ~~~~~~~~~ 1st Hyperpolarizability ~~~~~~~~~ !
				if (polarizbd%nPoints.LT.5) exit
				write (ou,341)

				bxxx = deLagDeriv(5,tpLowerCase(sXXX))*PHpl1
				byyy = deLagDeriv(5,tpLowerCase(sYYY))*PHpl1
				bxyy = deLagDeriv(5,tpLowerCase(sXYY))*PHpl1
				byxx = deLagDeriv(5,tpLowerCase(sYXX))*PHpl1
				bmod=sqrt((bxxx+bxyy)**2+(byyy+byxx)**2)

				bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,4)=bxyy
				bstDeriv(3,5)=byxx; bstDeriv(3,10)=bmod

				write (ou,351) sXXX,bxxx,sYYY,byyy,sXYY,bxyy,sYXX,byxx,bmod

				if (polarizbd%nPoints.EQ.5) exit

				bxxx = deLagDeriv(7,tpLowerCase(sXXX))*PHpl1
				byyy = deLagDeriv(7,tpLowerCase(sYYY))*PHpl1
				bxyy = deLagDeriv(7,tpLowerCase(sXYY))*PHpl1
				byxx = deLagDeriv(7,tpLowerCase(sYXX))*PHpl1
				bmod=sqrt((bxxx+bxyy)**2+(byyy+byxx)**2)

				bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,4)=bxyy
				bstDeriv(3,5)=byxx; bstDeriv(3,10)=bmod

				write (ou,352) sXXX,bxxx,sYYY,byyy,sXYY,bxyy,sYXX,byxx,bmod

				if (polarizbd%nPoints.EQ.7) exit

				bxxx = deLagDeriv(9,tpLowerCase(sXXX))*PHpl1
				byyy = deLagDeriv(9,tpLowerCase(sYYY))*PHpl1
				bxyy = deLagDeriv(9,tpLowerCase(sXYY))*PHpl1
				byxx = deLagDeriv(9,tpLowerCase(sYXX))*PHpl1
				bmod=sqrt((bxxx+bxyy)**2+(byyy+byxx)**2)

				bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,4)=bxyy
				bstDeriv(3,5)=byxx; bstDeriv(3,10)=bmod

				write (ou,353) sXXX,bxxx,sYYY,byyy,sXYY,bxyy,sYXX,byxx,bmod
				exit
			enddo

			do ! ~~~~~~~~~ 2nd Hyperpolarizability ~~~~~~~~~ !
				if (polarizbd%nPoints.LT.5) exit
				write (ou,361)

				gxxxx = deLagDeriv(5,tpLowerCase(sXXXX))*PHpl2
				gyyyy = deLagDeriv(5,tpLowerCase(sYYYY))*PHpl2
				gxxyy = deLagDeriv(5,tpLowerCase(sXXYY))*PHpl2
				gyyxx = deLagDeriv(5,tpLowerCase(sYYXX))*PHpl2
				averg = (gxxxx+gyyyy+gxxyy+gyyxx)/5

				bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,4)=gxxyy
				bstDeriv(4,5)=gyyxx; bstDeriv(4,10)=averg
					
				write (ou,371) sXXXX,gxxxx,sYYYY,gyyyy,sXXYY,gxxyy,sYYXX,gyyxx,averg,averg*esu

				if (polarizbd%nPoints.EQ.5) then
					if (Nel.NE.Nat) then
						write (ou,393) averg/Nat
					endif
					write (ou,394) averg/Nel,myCont(averg,gxxxx,gyyyy,5)
					write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
					ProEt = MEt(sta:sto,0); void=lsmInterpolation(2)
					exit
				endif

				gxxxx = deLagDeriv(7,tpLowerCase(sXXXX))*PHpl2
				gyyyy = deLagDeriv(7,tpLowerCase(sYYYY))*PHpl2
				gxxyy = deLagDeriv(7,tpLowerCase(sXXYY))*PHpl2
				gyyxx = deLagDeriv(7,tpLowerCase(sYYXX))*PHpl2
				averg = (gxxxx+gyyyy+gxxyy+gyyxx)/5

				bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,4)=gxxyy
				bstDeriv(4,5)=gyyxx; bstDeriv(4,10)=averg

				write (ou,372) sXXXX,gxxxx,sYYYY,gyyyy,sXXYY,gxxyy,sYYXX,gyyxx,averg,averg*esu

				if (polarizbd%nPoints.EQ.7) then
					if (Nel.NE.Nat) then
						write (ou,393) averg/Nat
					endif
					write (ou,394) averg/Nel,myCont(averg,gxxxx,gyyyy,5)

					write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
					ProEt = MEt(sta:sto,0); void=lsmInterpolation(2)
					exit
				endif

				gxxxx = deLagDeriv(9,tpLowerCase(sXXXX))*PHpl2
				gyyyy = deLagDeriv(9,tpLowerCase(sYYYY))*PHpl2
				gxxyy = deLagDeriv(9,tpLowerCase(sXXYY))*PHpl2
				gyyxx = deLagDeriv(9,tpLowerCase(sYYXX))*PHpl2
				averg = (gxxxx+gyyyy+gxxyy+gyyxx)/5

				bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,4)=gxxyy
				bstDeriv(4,5)=gyyxx; bstDeriv(4,10)=averg
					
				write (ou,373) sXXXX,gxxxx,sYYYY,gyyyy,sXXYY,gxxyy,sYYXX,gyyxx,averg,averg*esu

				if (Nel.NE.Nat) then
					write (ou,393) averg/Nat
				endif

				write (ou,394) averg/Nel,myCont(averg,gxxxx,gyyyy,5)
				write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
				ProEt = MEt(sta:sto,0); void=lsmInterpolation(2)
				exit
			enddo


		case ('ts')

		do k = sta,sto
			do i = sta,sto
				do j = sta,sto
					if ((i.NE.0).AND.(j.NE.0).AND.(k.NE.0)) cycle
					write (ou,104) i,j,k,MMEt(i,j,k)
				enddo
			enddo
		enddo
		write (ou,105)

		do ! ~~~~~~~~~ Dipole moment ~~~~~~~~~ !
			if (polarizbd%nPoints.LT.3) exit
			write (ou,401)

			dx=deLagDeriv(3,'x')*PDip
			dy=deLagDeriv(3,'y')*PDip
			dz=deLagDeriv(3,'z')*PDip
			dmod=sqrt(dx**2+dy**2+dz**2)

			write (ou,411) dx,dy,dz,dmod

			bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,3)=dz; bstDeriv(1,4)=dmod

			if (polarizbd%nPoints.EQ.3) exit

			dx=deLagDeriv(5,'x')*PDip
			dy=deLagDeriv(5,'y')*PDip
			dz=deLagDeriv(5,'z')*PDip
			dmod=sqrt(dx**2+dy**2+dz**2)

			write (ou,412) dx,dy,dz,dmod

			bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,3)=dz; bstDeriv(1,4)=dmod

			if (polarizbd%nPoints.EQ.5) exit

			dx=deLagDeriv(7,'x')*PDip
			dy=deLagDeriv(7,'y')*PDip
			dz=deLagDeriv(7,'z')*PDip
			dmod=sqrt(dx**2+dy**2+dz**2)

			write (ou,413) dx,dy,dz,dmod

			bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,3)=dz; bstDeriv(1,4)=dmod

			if (polarizbd%nPoints.EQ.7) exit

			dx=deLagDeriv(9,'x')*PDip
			dy=deLagDeriv(9,'y')*PDip
			dz=deLagDeriv(9,'z')*PDip
			dmod=sqrt(dx**2+dy**2+dz**2)

			write (ou,414) dx,dy,dz,dmod

			bstDeriv(1,1)=dx; bstDeriv(1,2)=dy; bstDeriv(1,3)=dz; bstDeriv(1,4)=dmod
			exit
		enddo

		do ! ~~~~~~~~~ Polarizability ~~~~~~~~~ !
			if (polarizbd%nPoints.LT.3) exit
			write (ou,421)

			axx = deLagDeriv(3,'xx')*PPol
			ayy = deLagDeriv(3,'yy')*PPol
			azz = deLagDeriv(3,'zz')*PPol

			avA = (axx+ayy+azz)/3

			write (ou,431) axx,ayy,azz,avA

			bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,3)=azz; bstDeriv(2,4)=avA

			if (polarizbd%nPoints.EQ.3) then
				if (Nel.NE.Nat) then
					write (ou,491) avA/Nat
				endif

				write (ou,492) avA/Nel,myCont(avA,axx,ayy,3)
				write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
				exit
			endif

			axx = deLagDeriv(5,'xx')*PPol
			ayy = deLagDeriv(5,'yy')*PPol
			azz = deLagDeriv(5,'zz')*PPol

			avA = (axx+ayy+azz)/3

			write (ou,432) axx,ayy,azz,avA

			bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,3)=azz; bstDeriv(2,4)=avA

			if (polarizbd%nPoints.EQ.5) then
				if (Nel.NE.Nat) then
					write (ou,491) avA/Nat
				endif
				write (ou,492) avA/Nel,myCont(avA,axx,ayy,3)
				write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
				exit
			endif

			axx = deLagDeriv(7,'xx')*PPol
			ayy = deLagDeriv(7,'yy')*PPol
			azz = deLagDeriv(7,'zz')*PPol

			avA = (axx+ayy+azz)/3

			write (ou,433) axx,ayy,azz,avA

			bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,3)=azz; bstDeriv(2,4)=avA

			if (polarizbd%nPoints.EQ.7) then
				if (Nel.NE.Nat) then
					write (ou,491) avA/Nat
				endif

				write (ou,492) avA/Nel,myCont(avA,axx,ayy,3)
				write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
				exit
			endif

			axx = deLagDeriv(9,'xx')*PPol
			ayy = deLagDeriv(9,'yy')*PPol
			azz = deLagDeriv(9,'zz')*PPol

			avA = (axx+ayy+azz)/3

			write (ou,434) axx,ayy,azz,avA

			bstDeriv(2,1)=axx; bstDeriv(2,2)=ayy; bstDeriv(2,3)=azz; bstDeriv(2,4)=avA

			if (Nel.NE.Nat) then
				write (ou,491) avA/Nat
			endif
			write (ou,492) avA/Nel,myCont(avA,axx,ayy,3)
			write (ou,500) sigmaContributionAlpha,sigmaBondAlpha
			exit
		enddo

		do ! ~~~~~~~~~ 1st Hyperpolarizability ~~~~~~~~~ !
			if (polarizbd%nPoints.LT.5) exit
			write (ou,441)
				
			bxxx = deLagDeriv(5,'xxx')*PHpl1
			byyy = deLagDeriv(5,'yyy')*PHpl1
			bzzz = deLagDeriv(5,'zzz')*PHpl1
			bxyy = deLagDeriv(5,'xyy')*PHpl1
			byxx = deLagDeriv(5,'yxx')*PHpl1
			bxzz = deLagDeriv(5,'xzz')*PHpl1
			bzxx = deLagDeriv(5,'zxx')*PHpl1
			byzz = deLagDeriv(5,'yzz')*PHpl1
			bzyy = deLagDeriv(5,'zyy')*PHpl1
			bmod=sqrt((bxxx+bxyy+bxzz)**2+(byyy+byxx+byzz)**2+(bzzz+bzxx+bzyy)**2)

			write (ou,451) bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod

			bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,3)=bzzz; bstDeriv(3,4)=bxyy
			bstDeriv(3,5)=byxx; bstDeriv(3,6)=bxzz; bstDeriv(3,7)=bzxx; bstDeriv(3,8)=byzz
			bstDeriv(3,9)=bzyy; bstDeriv(3,10)=bmod

			if (polarizbd%nPoints.EQ.5) exit

			bxxx = deLagDeriv(7,'xxx')*PHpl1
			byyy = deLagDeriv(7,'yyy')*PHpl1
			bzzz = deLagDeriv(7,'zzz')*PHpl1
			bxyy = deLagDeriv(7,'xyy')*PHpl1
			byxx = deLagDeriv(7,'yxx')*PHpl1
			bxzz = deLagDeriv(7,'xzz')*PHpl1
			bzxx = deLagDeriv(7,'zxx')*PHpl1
			byzz = deLagDeriv(7,'yzz')*PHpl1
			bzyy = deLagDeriv(7,'zyy')*PHpl1
			bmod=sqrt((bxxx+bxyy+bxzz)**2+(byyy+byxx+byzz)**2+(bzzz+bzxx+bzyy)**2)

			bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,3)=bzzz; bstDeriv(3,4)=bxyy
			bstDeriv(3,5)=byxx; bstDeriv(3,6)=bxzz; bstDeriv(3,7)=bzxx; bstDeriv(3,8)=byzz
			bstDeriv(3,9)=bzyy; bstDeriv(3,10)=bmod

			write (ou,452) bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod

			if (polarizbd%nPoints.EQ.7) exit

			bxxx = deLagDeriv(9,'xxx')*PHpl1
			byyy = deLagDeriv(9,'yyy')*PHpl1
			bzzz = deLagDeriv(9,'zzz')*PHpl1
			bxyy = deLagDeriv(9,'xyy')*PHpl1
			byxx = deLagDeriv(9,'yxx')*PHpl1
			bxzz = deLagDeriv(9,'xzz')*PHpl1
			bzxx = deLagDeriv(9,'zxx')*PHpl1
			byzz = deLagDeriv(9,'yzz')*PHpl1
			bzyy = deLagDeriv(9,'zyy')*PHpl1
			bmod=sqrt((bxxx+bxyy+bxzz)**2+(byyy+byxx+byzz)**2+(bzzz+bzxx+bzyy)**2)

			bstDeriv(3,1)=bxxx; bstDeriv(3,2)=byyy; bstDeriv(3,3)=bzzz; bstDeriv(3,4)=bxyy
			bstDeriv(3,5)=byxx; bstDeriv(3,6)=bxzz; bstDeriv(3,7)=bzxx; bstDeriv(3,8)=byzz
			bstDeriv(3,9)=bzyy; bstDeriv(3,10)=bmod

			write (ou,453) bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod
			exit
		enddo

		do ! ~~~~~~~~~ 2nd Hyperpolarizability ~~~~~~~~~ !
			if (polarizbd%nPoints.LT.5) exit
			write (ou,461)

			gxxxx = deLagDeriv(5,'xxxx')*PHpl2
			gyyyy = deLagDeriv(5,'yyyy')*PHpl2
			gzzzz = deLagDeriv(5,'zzzz')*PHpl2
			gxxyy = deLagDeriv(5,'xxyy')*PHpl2
			gyyxx = deLagDeriv(5,'yyxx')*PHpl2
			gxxzz = deLagDeriv(5,'xxzz')*PHpl2
			gzzxx = deLagDeriv(5,'zzxx')*PHpl2
			gyyzz = deLagDeriv(5,'yyzz')*PHpl2
			gzzyy = deLagDeriv(5,'zzyy')*PHpl2
			averg = (gxxxx+gyyyy+gzzzz+gxxyy+gyyxx+gxxzz+gzzxx+gyyzz+gzzyy)/5

			bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,3)=gzzzz; bstDeriv(4,4)=gxxyy
			bstDeriv(4,5)=gyyxx; bstDeriv(4,6)=gxxzz; bstDeriv(4,7)=gzzxx; bstDeriv(4,8)=gyyzz
			bstDeriv(4,9)=gzzyy; bstDeriv(4,10)=averg
				
			write (ou,471) gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averg,averg*esu

			if (polarizbd%nPoints.EQ.5) then
				if (Nel.NE.Nat) then
					write (ou,493) averg/Nat
				endif
				write (ou,494) averg/Nel,myCont(averg,gxxxx,gyyyy,5)

				write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
				ProEt = MMEt(sta:sto,0,0); void=lsmInterpolation(3)
				exit
			endif

			gxxxx = deLagDeriv(7,'xxxx')*PHpl2
			gyyyy = deLagDeriv(7,'yyyy')*PHpl2
			gzzzz = deLagDeriv(7,'zzzz')*PHpl2
			gxxyy = deLagDeriv(7,'xxyy')*PHpl2
			gyyxx = deLagDeriv(7,'yyxx')*PHpl2
			gxxzz = deLagDeriv(7,'xxzz')*PHpl2
			gzzxx = deLagDeriv(7,'zzxx')*PHpl2
			gyyzz = deLagDeriv(7,'yyzz')*PHpl2
			gzzyy = deLagDeriv(7,'zzyy')*PHpl2
			averg = (gxxxx+gyyyy+gzzzz+gxxyy+gyyxx+gxxzz+gzzxx+gyyzz+gzzyy)/5

			bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,3)=gzzzz; bstDeriv(4,4)=gxxyy
			bstDeriv(4,5)=gyyxx; bstDeriv(4,6)=gxxzz; bstDeriv(4,7)=gzzxx; bstDeriv(4,8)=gyyzz
			bstDeriv(4,9)=gzzyy; bstDeriv(4,10)=averg

			write (ou,472) gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averg,averg*esu

			if (polarizbd%nPoints.EQ.7) then
				if (Nel.NE.Nat) then
					write (ou,493) averg/Nat
				endif
				write (ou,494) averg/Nel,myCont(averg,gxxxx,gyyyy,5)

				write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
				ProEt = MMEt(sta:sto,0,0); void=lsmInterpolation(3)
				exit
			endif

			gxxxx = deLagDeriv(9,'xxxx')*PHpl2
			gyyyy = deLagDeriv(9,'yyyy')*PHpl2
			gzzzz = deLagDeriv(9,'zzzz')*PHpl2
			gxxyy = deLagDeriv(9,'xxyy')*PHpl2
			gyyxx = deLagDeriv(9,'yyxx')*PHpl2
			gxxzz = deLagDeriv(9,'xxzz')*PHpl2
			gzzxx = deLagDeriv(9,'zzxx')*PHpl2
			gyyzz = deLagDeriv(9,'yyzz')*PHpl2
			gzzyy = deLagDeriv(9,'zzyy')*PHpl2
			averg = (gxxxx+gyyyy+gzzzz+gxxyy+gyyxx+gxxzz+gzzxx+gyyzz+gzzyy)/5

			bstDeriv(4,1)=gxxxx; bstDeriv(4,2)=gyyyy; bstDeriv(4,3)=gzzzz; bstDeriv(4,4)=gxxyy
			bstDeriv(4,5)=gyyxx; bstDeriv(4,6)=gxxzz; bstDeriv(4,7)=gzzxx; bstDeriv(4,8)=gyyzz
			bstDeriv(4,9)=gzzyy; bstDeriv(4,10)=averg
				
			write (ou,473) gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,averg,averg*esu
			if (Nel.NE.Nat) then
				write (ou,493) averg/Nat
			endif

			write (ou,494) averg/Nel,myCont(averg,gxxxx,gyyyy,5)

			write (ou,501) sigmaContributionGamma,CCsigmaBondGamma,CHsigmaBondGamma
			ProEt = MMEt(sta:sto,0,0); void=lsmInterpolation(3)
			exit
		enddo

	end select

	if (polarizbd%nPoints.GE.3) write (ou,502)
	if (polarizbd%nPoints.GE.5) write (ou,503)

200 format (///4X,'~~~~~~~~~~',1X,A,1X,4X,'Scale = 'A1,1X,'~~~~~~~~~~'/)

201 format (17X,'3 points',15X,'5 points',19X,'7 points',19X,'9 points'/)
211 format ('Dip Mom', 5X,F11.4, 1X,'Deby',7X,F11.4, 1X,'Deby',11X,F11.4, 1X,'Deby',11X,F11.4, 1X,'Deby'/)
212 format ('Polariz',    F16.4, 1X,'a.u.',2X,F16.4, 1X,'a.u.', 6X,F16.4, 1X,'a.u.', 6X,F16.4, 1X,'a.u.'/)
213 format ('Hyper-1',24X,                    F15.7, 1X,'a.u.', 7X,F15.7, 1X,'a.u.', 7X,F15.7, 1X,'a.u.'/)
214 format ('Hyper-2',26X,                    ES13.6,1X,'a.u.', 9X,ES13.6,1X,'a.u.', 9X,ES13.6,1X,'a.u.'/&
&				      24X,                    F15.1,'x10^-36 esu',1X,F15.1,'x10^-36 esu',1X,F15.1,'x10^-36 esu')

221 format (17X,'3 points',15X,'5 points',19X,'7 points'/)
231 format ('Dip Mom', 5X,F11.4, 1X,'Deby',7X,F11.4, 1X,'Deby',11X,F11.4, 1X,'Deby'/)
232 format ('Polariz',    F16.4, 1X,'a.u.',2X,F16.4, 1X,'a.u.', 6X,F16.4, 1X,'a.u.'/)
233 format ('Hyper-1',24X,                    F15.7, 1X,'a.u.', 7X,F15.7, 1X,'a.u.'/)
234 format ('Hyper-2',26X,                    ES13.6,1X,'a.u.', 9X,ES13.6,1X,'a.u.'/&
&				      24X,                    F15.1,'x10^-36 esu',1X,F15.1,'x10^-36 esu')

241 format (17X,'3 points',15X,'5 points'/)
251 format ('Dip Mom', 5X,F11.4, 1X,'Deby',7X,F11.4, 1X,'Deby'/)
252 format ('Polariz',    F16.4, 1X,'a.u.',2X,F16.4, 1X,'a.u.'/)
253 format ('Hyper-1',24X,                    F15.7, 1X,'a.u.'/)
254 format ('Hyper-2',26X,                    ES13.6,1X,'a.u.'/&
&				      24X,                    F15.1,'x10^-36 esu')

261 format (17X,'3 points'/)
271 format ('Dip Mom', 5X,F11.4, 1X,'Deby'/)
272 format ('Polariz',    F16.4, 1X,'a.u.')

281 format (/'         Average Polarizability                =',1X,F13.4,1X,'a.u.'/&
&			 'Specific Average Polarizability (per electron) =',1X,F13.4,1X,'a.u.' )
282 format ( 'Specific Average Polarizability (per atom)     =',1X,F13.4,1X,'a.u.'/)

283 format (/'         Average HyperPol - 2                =',1X,ES14.7,1X,'a.u.'/&
&			 'Specific Average HyperPol - 2 (per electron) =',1X,ES14.7,1X,'a.u.')
284 format ( 'Specific Average HyperPol - 2 (per atom)     =',1X,ES14.7,1X,'a.u.')

101 format ( 4X,'F',1X,i2,1X,F25.16)
102 format ( 4X,'F',1X,i2,1X,i2,1X,F25.16)

103 format (/4X,A,1X,'energies'/)
104 format ( 4X,'F',1X,i2,1X,i2,1X,i2,F25.16)
105 format (/4X,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'/)

111	format (//2X,A/&
&			  2X,A,1X,A,1X,A/&
&			  2X,A//)
112 format ( /2X,A/)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

301 format (/ 7X,'Dipole moment (Deby)'/)
311 format (  2X,'3 points',3X,A1,' =',1X,F7.3,4X,A1,' =',1X,F7.3,4X,'|D| =',1X,F7.3)
312 format (  2X,'5 points',3X,A1,' =',1X,F7.3,4X,A1,' =',1X,F7.3,4X,'|D| =',1X,F7.3)
313 format (  2X,'7 points',3X,A1,' =',1X,F7.3,4X,A1,' =',1X,F7.3,4X,'|D| =',1X,F7.3)
314 format (  2X,'9 points',3X,A1,' =',1X,F7.3,4X,A1,' =',1X,F7.3,4X,'|D| =',1X,F7.3)

321 format (/ 7X,'Polarizability (a.u.)'/)
331 format (  2X,'3 points',2X,A2,' =',1X,F13.4,2X,A2,' =',1X,F13.4,2X,'<alpha> =',1X,F13.4)
332 format (  2X,'5 points',2X,A2,' =',1X,F13.4,2X,A2,' =',1X,F13.4,2X,'<alpha> =',1X,F13.4)
333 format (  2X,'7 points',2X,A2,' =',1X,F13.4,2X,A2,' =',1X,F13.4,2X,'<alpha> =',1X,F13.4)
334 format (  2X,'9 points',2X,A2,' =',1X,F13.4,2X,A2,' =',1X,F13.4,2X,'<alpha> =',1X,F13.4)

391 format (/2X,'Specific Average Polarizability (per atom)     =',1X,F13.4,1X,'a.u.')
392 format (/2X,'Specific Average Polarizability (per electron) =',1X,F13.4,1X,'a.u.'/&
&			 2X,'Contribution of General Component =',3X,F11.4,'%')

393 format (/2X,'Specific Average HyperPol - 2 (per atom)     =',9X,ES14.7,1X,'a.u.')
394 format (/2X,'Specific Average HyperPol - 2 (per electron) =',9X,ES14.7,1X,'a.u.'/&
&			 2X,'Contribution of General Component =',3X,F11.4,'%')

341 format (//7X,'HyperPol - 1 (a.u.)'/)
351 format (  2X,'5 points',3X,A3,' =',1X,F14.6,3X,A3,' =',1X,F14.6/&
&			  2X,          11X,A3,' =',1X,F14.6,3X,A3,' =',1X,F14.6,3X,'|B| =',1X,F14.6)
352 format (  2X,'7 points',3X,A3,' =',1X,F14.6,3X,A3,' =',1X,F14.6/&
&			  2X,          11X,A3,' =',1X,F14.6,3X,A3,' =',1X,F14.6,3X,'|B| =',1X,F14.6)
353 format (  2X,'9 points',3X,A3,' =',1X,F14.6,3X,A3,' =',1X,F14.6/&
&			  2X,          11X,A3,' =',1X,F14.6,3X,A3,' =',1X,F14.6,3X,'|B| =',1X,F14.6)

361 format (//7X,'HyperPol - 2 (a.u.)'/)
371 format (  2X,'5 points',2X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
&                          12X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
&             2X,'<HyperPol - 2> =',   2X,ES14.7,3X,'=',2X,F17.7,'E-36 esu')
372 format (  2X,'7 points',2X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
&                          12X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
&             2X,'<HyperPol - 2> =',   2X,ES14.7,3X,'=',2X,F17.7,'E-36 esu')
373 format (  2X,'9 points',2X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
&                          12X,A4,' =',2X,ES14.7,5X,A4,' =',2X,ES14.7/&
&             2X,'<HyperPol - 2> =',   2X,ES14.7,3X,'=',2X,F17.7,'E-36 esu')

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

401 format (/ 7X,'Dipole moment (Deby)'/)
411 format (  2X,'3 points',3X,'X =',1X,F7.3,4X,'Y =',1X,F7.3,4X,'Z =',1X,F7.3,4X,'|D| =',1X,F7.3)
412 format (  2X,'5 points',3X,'X =',1X,F7.3,4X,'Y =',1X,F7.3,4X,'Z =',1X,F7.3,4X,'|D| =',1X,F7.3)
413 format (  2X,'7 points',3X,'X =',1X,F7.3,4X,'Y =',1X,F7.3,4X,'Z =',1X,F7.3,4X,'|D| =',1X,F7.3)
414 format (  2X,'9 points',3X,'X =',1X,F7.3,4X,'Y =',1X,F7.3,4X,'Z =',1X,F7.3,4X,'|D| =',1X,F7.3)

421 format (/ 7X,'Polarizability (a.u.)'/)
431 format (  2X,'3 points',2X,'XX =',1X,F13.4,2X,'YY =',1X,F13.4,2X,'ZZ =',1X,F13.4,2X,'<alpha> =',1X,F13.4)
432 format (  2X,'5 points',2X,'XX =',1X,F13.4,2X,'YY =',1X,F13.4,2X,'ZZ =',1X,F13.4,2X,'<alpha> =',1X,F13.4)
433 format (  2X,'7 points',2X,'XX =',1X,F13.4,2X,'YY =',1X,F13.4,2X,'ZZ =',1X,F13.4,2X,'<alpha> =',1X,F13.4)
434 format (  2X,'9 points',2X,'XX =',1X,F13.4,2X,'YY =',1X,F13.4,2X,'ZZ =',1X,F13.4,2X,'<alpha> =',1X,F13.4)

491 format (/2X,'Specific Average Polarizability (per atom)     =',1X,F13.4,1X,'a.u.')
492 format (/2X,'Specific Average Polarizability (per electron) =',1X,F13.4,1X,'a.u.'/&
&			 2X,'Contribution of General Component =',3X,F11.4,'%')

441 format (//7X,'HyperPol - 1 (a.u.)'/)
451 format (  2X,'5 points',3X,'XXX =',1X,F14.6,3X,'YYY =',1X,F14.6,3X,'ZZZ =',1X,F14.6/&
&			  2X,          11X,'XYY =',1X,F14.6,3X,'YXX =',1X,F14.6,/&
&			  2X,          11X,'XZZ =',1X,F14.6,3X,'ZXX =',1X,F14.6,/&
&			  2X,          11X,'YZZ =',1X,F14.6,3X,'ZYY =',1X,F14.6,3X,'|B| =',1X,F14.6)

452 format (  2X,'7 points',3X,'XXX =',1X,F14.6,3X,'YYY =',1X,F14.6,3X,'ZZZ =',1X,F14.6/&
&			  2X,          11X,'XYY =',1X,F14.6,3X,'YXX =',1X,F14.6,/&
&			  2X,          11X,'XZZ =',1X,F14.6,3X,'ZXX =',1X,F14.6,/&
&			  2X,          11X,'YZZ =',1X,F14.6,3X,'ZYY =',1X,F14.6,3X,'|B| =',1X,F14.6)

453 format (  2X,'9 points',3X,'XXX =',1X,F14.6,3X,'YYY =',1X,F14.6,3X,'ZZZ =',1X,F14.6/&
&			  2X,          11X,'XYY =',1X,F14.6,3X,'YXX =',1X,F14.6,/&
&			  2X,          11X,'XZZ =',1X,F14.6,3X,'ZXX =',1X,F14.6,/&
&			  2X,          11X,'YZZ =',1X,F14.6,3X,'ZYY =',1X,F14.6,3X,'|B| =',1X,F14.6)

461 format (//7X,'HyperPol - 2 (a.u.)'/)
471 format (  2X,'5 points',2X,'XXXX =',2X,ES14.7,5X,'YYYY =',2X,ES14.7,5X,'ZZZZ =',2X,ES14.7/&
&                          12X,'XXYY =',2X,ES14.7,5X,'YYXX =',2X,ES14.7/&
&                          12X,'XXZZ =',2X,ES14.7,5X,'ZZXX =',2X,ES14.7/&
&                          12X,'YYZZ =',2X,ES14.7,5X,'ZZYY =',2X,ES14.7/&
&             2X,'<HyperPol - 2> =',    2X,ES14.7,3X,'=',2X,F17.7,'E-36 esu')

472 format (  2X,'7 points',2X,'XXXX =',2X,ES14.7,5X,'YYYY =',2X,ES14.7,5X,'ZZZZ =',2X,ES14.7/&
&                          12X,'XXYY =',2X,ES14.7,5X,'YYXX =',2X,ES14.7/&
&                          12X,'XXZZ =',2X,ES14.7,5X,'ZZXX =',2X,ES14.7/&
&                          12X,'YYZZ =',2X,ES14.7,5X,'ZZYY =',2X,ES14.7/&
&             2X,'<HyperPol - 2> =',    2X,ES14.7,3X,'=',2X,F17.7,'E-36 esu')

473 format (  2X,'9 points',2X,'XXXX =',2X,ES14.7,5X,'YYYY =',2X,ES14.7,5X,'ZZZZ =',2X,ES14.7/&
&                          12X,'XXYY =',2X,ES14.7,5X,'YYXX =',2X,ES14.7/&
&                          12X,'XXZZ =',2X,ES14.7,5X,'ZZXX =',2X,ES14.7/&
&                          12X,'YYZZ =',2X,ES14.7,5X,'ZZYY =',2X,ES14.7/&
&             2X,'<HyperPol - 2> =',    2X,ES14.7,3X,'=',2X,F17.7,'E-36 esu')

493 format (/2X,'Specific Average HyperPol - 2 (per atom)     =',9X,ES14.7,1X,'a.u.')
494 format (/2X,'Specific Average HyperPol - 2 (per electron) =',9X,ES14.7,1X,'a.u.'/&
&			 2X,'Contribution of General Component =',3X,F11.4,'%')

500 format (/2X,'Sigma contribution to average polarizability =',1X,F11.2/&
&            2X,F5.3,' a.u. per C-C and C-H bonds [1]')
501 format (/2X,'Sigma contribution to average 2nd hyperpolarizability =',1X,ES10.3/&
&            2X,F5.0,' a.u. per C-C bond [2]'/&
&            2X,F5.0,' a.u. per C-H bond [2]')
502 format (/2X,'________________________________________________________'/&
&			 2X,'[1] Yu.B. Visotski and B.C. Briantsev, Quantum Chemistry'/&
&			 2X,'    of Radicals and Ion-Radicals with Conjugated Bonds'/&
&			 2X,'    (Donetsk, DonGuet, 2004), p. 207, in Russian.')
503 format ( 2X,'[2] B.M. Pierce, J. Chem. Phys. 91, 791 (1989).')

	call deFinalize; return

	contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	integer(kind=iglu) function lsmInterpolation(scl) result(ret)

	implicit none

	integer(kind=iglu), intent(in) :: scl
	character                      :: line*510,value*100,pFormat*9,tmplt*4,crt*3

	integer(kind=iglu)             :: k,l,m,sc1,sc2,sc3,pAccuracy,overAllLen,cLen
	real   (kind=rglu)             :: cofs(0:4),avar	


	ret=0; if (rspu.NE.r16kind) return


	write (ou,'(/A)') tpAdjustc('LSM INSET TO CONTROL NUMERICAL STEADINESS',width,'=')

	cofs(0)=1; cofs(1)=PDip; cofs(2)=PPol; cofs(3)=PHpl1; cofs(4)=PHpl2
		
	pAccuracy=5; overAllLen=4+4+pAccuracy
	pFormat='(ES'//prStrByVal(overAllLen)//'.'//prStrByVal(pAccuracy)//')'
	overAllLen=overAllLen+1
	crt='xyz'
	select case (scl)

		case(1)
			sc1=index( crt,scale(1:1) )
			line=tpFill(line)
			do k = 0,4
				line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
			enddo; cLen=5+5*overAllLen
			write (ou,'(A)') line(1:cLen)

			do k = 0,4
				avar=transition(k)*deLSMDeriv(repeat(crt(sc1:sc1),k))
				if (abs(avar).LT.1D-5) avar=0
				value=tpFill(value); write (value,fmt=pFormat) avar
				line(5+1+k*overAllLen:5+(k+1)*overAllLen)=value(1:overAllLen)
			enddo
			write (ou,'(A)') line(1:cLen)

		case(2)
			sc1=index( crt,scale(1:1) )
			sc2=index( crt,scale(2:2) )

			select case (scale(1:2))
				case ('xy'); sc3=3
				case ('xz'); sc3=2
				case ('yz'); sc3=1
			end select

			line=tpFill(line)

			do k = 0,4
				line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
			enddo; cLen=5+5*overAllLen

			write (ou,'(A)') line(1:cLen)

			do l = 0,4
				line=tpFill(line); tmplt=tpFill(tmplt)
				tmplt=tpFill(l,crt(sc2:sc2))

				line(2:5)=adjustr( tmplt )
				do k = 0,4
					if (k+l.LE.4) then
						avar=transition(k+l)*deLSMDeriv(repeat(scale(1:1),k)//repeat(scale(2:2),l))
						if (abs(avar).LT.1D-5) avar=0.d0
						value=tpFill(value); write (value,fmt=pFormat) avar
						line(5+1+k*overAllLen:5+(k+1)*overAllLen)=value(1:overAllLen)
					else
						line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpFill(overAllLen)
					endif
				enddo

				write (ou,'(A)') line(1:cLen)
			enddo

		case(3)
			!xy
			sc1=1; sc2=2; sc3=3; line=tpFill(line)

			do k = 0,4
				line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
			enddo; cLen=5+5*overAllLen

			write (ou,'(A)') line(1:cLen)

			m=0
			do l = 0,4
				line=tpFill(line); tmplt=tpFill(tmplt)
				tmplt=tpFill(l,crt(sc2:sc2))

				line(2:5)=adjustr( tmplt )
				do k = 0,4
					if (k+l.LE.4) then
						avar=transition(k+l)*deLSMDeriv(repeat(crt(1:1),k)//repeat(crt(2:2),l)//repeat(crt(3:3),m))
						if (abs(avar).LT.1D-5) avar=0.d0
						value=tpFill(value); write (value,fmt=pFormat) avar
						line(5+1+k*overAllLen:5+(k+1)*overAllLen)=value(1:overAllLen)
					else
						line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpFill(overAllLen)
					endif
				enddo

				write (ou,'(A)') line(1:cLen)
			enddo

			!xz
			sc1=1; sc2=3; sc3=2; line=tpFill(line)

			do k = 0,4
				line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
			enddo; cLen=5+5*overAllLen

			write (ou,'(/A)') line(1:cLen)

			l=0
			do m = 0,4
				line=tpFill(line); tmplt=tpFill(tmplt)
				tmplt=tpFill(m,crt(sc2:sc2))

				line(2:5)=adjustr( tmplt )
				do k = 0,4
					if (k+m.LE.4) then
						avar=transition(k+m)*deLSMDeriv(repeat(crt(1:1),k)//repeat(crt(2:2),l)//repeat(crt(3:3),m))
						if (abs(avar).LT.1D-5) avar=0.d0
						value=tpFill(value); write (value,fmt=pFormat) avar
						line(5+1+k*overAllLen:5+(k+1)*overAllLen)=value(1:overAllLen)
					else
						line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpFill(overAllLen)
					endif
				enddo

				write (ou,'(A)') line(1:cLen)
			enddo

			!yz
			sc1=2; sc2=3; sc3=1; line=tpFill(line)

			do k = 0,4
				line(5+1+k*overAllLen:5+(k+1)*overAllLen)=tpAdjustc( tpFill(k,crt(sc1:sc1)),overAllLen )
			enddo; cLen=5+5*overAllLen

			write (ou,'(/A)') line(1:cLen)

			k=0
			do m = 0,4
				line=tpFill(line); tmplt=tpFill(tmplt)
				tmplt=tpFill(m,crt(sc2:sc2))

				line(2:5)=adjustr( tmplt )
				do l = 0,4
					if (l+m.LE.4) then
						avar=transition(l+m)*deLSMDeriv(repeat(crt(1:1),k)//repeat(crt(2:2),l)//repeat(crt(3:3),m))
						if (abs(avar).LT.1D-5) avar=0.d0
						value=tpFill(value); write (value,fmt=pFormat) avar
						line(5+1+l*overAllLen:5+(l+1)*overAllLen)=value(1:overAllLen)
					else
						line(5+1+l*overAllLen:5+(l+1)*overAllLen)=tpFill(overAllLen)
					endif
				enddo

				write (ou,'(A)') line(1:cLen)
			enddo
	end select

	write (ou,'(A)') tpFill(width,'=')

	return
	end function lsmInterpolation

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	real(kind=rglu) function myCont(aver,fcomp,scomp,div)
	implicit none

	real   (kind=rglu) :: fcomp,scomp,maxgg,aver
	integer(kind=iglu) :: div


	if (abs(fcomp).GT.abs(scomp)) then
		maxgg=abs(fcomp)/div
	else
		maxgg=abs(scomp)/div
	endif
	myCont=100.*maxgg/aver; return
	end function myCont

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

	end subroutine showPolarizability