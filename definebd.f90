	subroutine definebd

	use hdb
	use datablock, only: bdShareVariable,bdCollect

	implicit none

	character (len=*), parameter :: bdstart='%name{',bdstop='}',bdcomment='#',&
	                                bdaccord='=',bdseparator=','

	integer(kind=iglu) :: err
	integer(kind=8)    :: addr(10)


	err=fcBanID(in,true); if (err.LE.0) then; write (eu,*) 'Procedure definebd. Internal i/o error.'; stop; endif

	addr(1)=bdShareVariable(generalbd%methods        ,'method'          ,opt=false,several=true,&
	expect='list(huckel+hf+cue-ccs+r-ccd+u-ccd+mp2+mp3+cue-ccsd+r-ccsd+u-ccsd+r-ccsd(t)+cue-ccsdt+u-ccsdt+r-ccsdt+fci)')

	addr(2)=bdShareVariable(generalbd%task           ,'jobtype'         ,opt=false,several=true,&
	expect='list(polarizability+density+coulson+hypercharges+energy+wf-analize)')

	addr(3)=bdShareVariable(generalbd%coulombType,'coulomb-integrals'   ,opt=true ,several=true,def='ohno-klopman',&
	expect='list(ohno-klopman+mataga-nishimoto+hubbard)')

	addr(4)=bdShareVariable(generalbd%alternation,'bond-alternation'    ,opt=true ,def=gluZero,expect='range(0,1)')

	void=bdCollect('general',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(systembd%memory,         'memory'          ,opt=true,def=real(1024,kind=rglu))
	addr(2)=bdShareVariable(systembd%nNodes,         'nproc'           ,opt=true,def=1,expect='range(1,128)')
	addr(3)=bdShareVariable(systembd%allowRestart,   'restart'         ,opt=true,def=false)
	addr(4)=bdShareVariable(systembd%verboselvl,     'verbose'         ,opt=true,def=2,expect='range(1,6)')
	addr(5)=bdShareVariable(systembd%allowMutt,      'email-send'      ,opt=true,def=false)
	addr(6)=bdShareVariable(systembd%muttSendTared,  'email-send-tared',opt=true,def=true)
	addr(7)=bdShareVariable(systembd%muttDestination,'email-reciver'   ,opt=true,def='quant@mail.com')
	addr(8)=bdShareVariable(systembd%ignoreSIGHUP,   'ignore-sighup'   ,opt=true,def=true)
	void=bdCollect('system',addr(1:8),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(geometrybd%symmetryAccount      ,'symmetry-account'   ,opt=true,def=true)
	addr(2)=bdShareVariable(geometrybd%searchPlanar(1)      ,'search-planarity'   ,opt=true,def=true)
	addr(3)=bdShareVariable(geometrybd%searchLinear(1)      ,'search-linearity'   ,opt=true,def=true)
	addr(4)=bdShareVariable(geometrybd%symmetryTolerance    ,'symmetry-threshold' ,opt=true,def=real(12,kind=rglu),expect='range(1,15)')
	addr(5)=bdShareVariable(geometrybd%searchTolerance      ,'search-threshold'   ,opt=true,def=real(5,kind=rglu) ,expect='range(1,15)')
	addr(6)=bdShareVariable(geometrybd%randomDisplacement(1),'random-displacement',opt=true,def=gluZero ,expect='range(0,5)')
	void=bdCollect('geometry',addr(1:6),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(statesbd%nStates,'excited-states-number',opt=true,def=0,expect='range(0,5)')
	addr(2)=bdShareVariable(statesbd%spin   ,'spin'                 ,opt=true,def=0,expect='range(0,10)')
	void=bdCollect('states',addr(1:2),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(polarizbd%scales,'components' ,opt=true,def='xyz',expect='list(x,y,z,xy,yx,xz,zx,yz,zy,xyz)')
	addr(2)=bdShareVariable(polarizbd%nPoints,'points'    ,opt=true,def=5,expect='list(1,3,5,7,9)')
	addr(3)=bdShareVariable(polarizbd%maxPower,'max-power',opt=true,def=4,expect='range(1,4)')
	addr(4)=bdShareVariable(polarizbd%derivStep,'step'    ,opt=true,def=real(0.025d0,kind=rglu),expect='range(0,1)')
	void=bdCollect('polariz',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(densitybd%dtype,'elements'             ,opt=true,def='all',expect='list(none,charges,orders,all)')
	addr(2)=bdShareVariable(densitybd%nPoints,'derive-points'      ,opt=true,def=3,expect='list(1,3,5,7,9)')
	addr(3)=bdShareVariable(densitybd%prntAccuracy,'print-accuracy',opt=true,def=5)
	addr(4)=bdShareVariable(densitybd%derivStep,'derive-step'      ,opt=true,def=real(0.025d0,kind=rglu),expect='range(0,1)')
	addr(5)=bdShareVariable(densitybd%NOAnalize,'natural-orbitals' ,opt=true,def=false)
	void=bdCollect('rdm',addr(1:5),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(coulsonbd%ctype,'elements'      ,opt=true,def='atom-atom',expect='list(atom-atom,atom-bond,bond-bond)')
	addr(2)=bdShareVariable(coulsonbd%nPoints,'points'      ,opt=true,def=3,expect='list(1,3,5,7,9)')
	addr(3)=bdShareVariable(coulsonbd%derivPower,'max-power',opt=true,def=2,expect='range(1,4)')
	addr(4)=bdShareVariable(coulsonbd%derivStep,'step'      ,opt=true,def=real(0.02d0,kind=rglu),expect='range(0,1)')
	void=bdCollect('coulson',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(hyperchargesbd%scales,'components'           ,opt=true,def='xyz',expect='list(x,y,z,xy,yx,xz,zx,yz,zy,xyz)')
	addr(2)=bdShareVariable(hyperchargesbd%naPoints,'derive-alpha-points',opt=true,def=3,expect='list(1,3,5,7,9)')
	addr(3)=bdShareVariable(hyperchargesbd%nfPoints,'derive-field-points',opt=true,def=5,expect='list(1,3,5,7,9)')
	addr(4)=bdShareVariable(hyperchargesbd%derivPower,'max-power'        ,opt=true,def=4,expect='range(1,4)')
	addr(5)=bdShareVariable(hyperchargesbd%derivaStep,'derive-alpha-step',opt=true,def=real(0.025d0,kind=rglu),expect='range(0,1)')
	addr(6)=bdShareVariable(hyperchargesbd%derivfStep,'derive-field-step',opt=true,def=real(0.025d0,kind=rglu),expect='range(0,1)')
	void=bdCollect('hypercharges',addr(1:6),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(cuebd%radius(1),'t1-limit'    ,opt=true,def=-1)
	addr(2)=bdShareVariable(cuebd%radius(2),'t2-limit'    ,opt=true,def=-2)
	addr(3)=bdShareVariable(cuebd%radius(3),'t3-limit'    ,opt=true,def=-3)
	addr(4)=bdShareVariable(cuebd%sparse   ,'spare-algorithm' ,opt=true,def=false)
	addr(5)=bdShareVariable(cuebd%showBasis,'export-basis-xyz',opt=true,def=false)
	void=bdCollect('cue',addr(1:5),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(fcibd%nSteps,'steps'                ,opt=true,def=5,expect='range(1,10)')
	addr(2)=bdShareVariable(fcibd%maxiters,'max-iterations'     ,opt=true,def=200)
	addr(3)=bdShareVariable(fcibd%accuracy,'accuracy'           ,opt=true,def=real(12,kind=rglu),expect='range(1,15)')
	addr(4)=bdShareVariable(fcibd%zeroThreshold,'zero-threshold',opt=true,def=real(14,kind=rglu),expect='range(1,15)')
	void=bdCollect('fci',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(scfbd%maxiters,'max-iterations',opt=true,def=400)
	addr(2)=bdShareVariable(scfbd%accuracy,'accuracy'      ,opt=true,def=real(12,kind=rglu),expect='range(1,15)')
	addr(3)=bdShareVariable(scfbd%iterStep,'iter-step'     ,opt=true,def=real(0.05d0,kind=rglu),expect='range(0,1)')
	addr(4)=bdShareVariable(scfbd%guess,   'density-guess' ,opt=true,def='huckel',expect='list(huckel,unitmatrix)')
	void=bdCollect('scf',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(pipekbd%maxiters,'max-iterations',opt=true,def=800)
	addr(2)=bdShareVariable(pipekbd%accuracy,'accuracy'      ,opt=true,def=real(10,kind=rglu),expect='range(1,15)')
	addr(3)=bdShareVariable(pipekbd%enabled ,'enabled'       ,opt=true,def=false)
	void=bdCollect('local',addr(1:3),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(ccbd%projType,'projection-type',opt=true,def='2-1',expect='list(1,2-1)')
	addr(2)=bdShareVariable(ccbd%maxiters,'max-iterations' ,opt=true,def=2000)
	addr(3)=bdShareVariable(ccbd%accuracy,'accuracy'       ,opt=true,def=real(12,kind=rglu),expect='range(1,15)')
	addr(4)=bdShareVariable(ccbd%iterStep(1),'t1-step'     ,opt=true,def=gluUnity,expect='range(0,1)')
	addr(5)=bdShareVariable(ccbd%iterStep(2),'t2-step'     ,opt=true,def=gluUnity,expect='range(0,1)')
	addr(6)=bdShareVariable(ccbd%iterStep(3),'t3-step'     ,opt=true,def=gluUnity,expect='range(0,1)')
	addr(7)=bdShareVariable(ccbd%forceSpin  ,'force-spin'  ,opt=true,def=false)
	void=bdCollect('coupled-cluster',addr(1:7),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(lrbd%guess   ,'guess'         ,opt=true,def='cis',expect='list(cis,rpa)')
	addr(2)=bdShareVariable(lrbd%maxiters,'max-iterations',opt=true,def=10000)
	addr(3)=bdShareVariable(lrbd%accuracy,'accuracy'      ,opt=true,def=real(6,kind=rglu),expect='range(1,15)')
	addr(4)=bdShareVariable(lrbd%iterStep,'iter-step'     ,opt=true,def=real(0.03d0,kind=rglu),expect='range(0,1)')
	void=bdCollect('linear-responce',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	addr(1)=bdShareVariable(diisbd%storage,'storage',opt=true,def='ram',expect='list(ram,hdd)')
	addr(2)=bdShareVariable(diisbd%steps  ,'steps'    ,opt=true,def=20,expect='range(2,50)')
	addr(3)=bdShareVariable(diisbd%enabled,'enabled',opt=true,def=false)
	void=bdCollect('diis',addr(1:3),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

	void=bdCollect('molecule',addr(1:1),bdstart,bdstop,bdcomment,bdaccord,bdseparator,true,in)

	return
	end subroutine definebd