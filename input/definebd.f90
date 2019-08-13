    subroutine definebd

    use glob,      only: assignment(=)
    use glob,      only: rglu,iglu,true,void,false
    use hdb,       only: in,eu,lrg,scfg,ierror
    use hdb,       only: generalbd,systembd,geometrybd,statesbd,polarizbd,iterationbd
    use hdb,       only: densitybd,coulsonbd,hyperchargesbd,cuebd,fcibd,scfbd,lrbd,ccbd
    use hdb,       only: localbd,fieldbd
    use fcontrol,  only: fcNewID,fcBanID
    use datablock, only: bdShareVariable,bdCollect

    implicit none

    character (len=*), parameter :: bdstart='%name{',bdstop='}',bdcomment='#',&
                                    bdaccord='=',bdseparator=','

    integer(kind=iglu) :: err
    integer(kind=8)    :: addr(15)


    err=fcBanID(in,true)
    if (err.LE.0) then
        ierror='Procedure definebd. Internal i/o error. position #1. contact developer.'
        call primaryInformation('error')
    endif

    addr(1)=bdShareVariable(generalbd%methods        ,'method'          ,opt=false,several=true,&
    expect='list(all,huckel,rhf,cue-ccs,r-ccd,u-ccd,mp2,mp3,cue-ccsd,r-ccsd,u-ccsd,r-ccsd(t),cue-ccsdt,u-ccsdt,r-ccsdt,fci)')
    addr(2)=bdShareVariable(generalbd%task           ,'jobtype'         ,opt=false,several=true,&
    expect='list(polarizability,density,coulson,hypercharges,energy,wf-analysis)')
    addr(3)=bdShareVariable(generalbd%coulombType,'coulomb-integrals'   ,opt=true ,several=false,def='ohno-klopman',&
    expect='list(ohno-klopman,mataga-nishimoto,hubbard)')
    addr(4)=bdShareVariable(generalbd%alternation,'bond-alternation'    ,opt=true ,def=0._rglu,expect='range(-1:1)')
    void=bdCollect('general',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr( 1)=bdShareVariable(systembd%memory,          'memory'                ,opt=true,def=1024._rglu)
    addr( 2)=bdShareVariable(systembd%nNodes,          'nproc'                 ,opt=true,def=1,expect='range(1:128)')
    addr( 3)=bdShareVariable(systembd%allowRestart,    'restart'               ,opt=true,def=false)                        !TODO
    addr( 4)=bdShareVariable(systembd%verboselvl,      'verbose'               ,opt=true,def=2,expect='range(1:6)')        !TODO
    addr( 5)=bdShareVariable(systembd%ignoreSIGHUP,    'ignore-sighup'         ,opt=true,def=true)
    addr( 6)=bdShareVariable(systembd%harvest,         'harvest'               ,opt=true,def=false)
    addr( 7)=bdShareVariable(systembd%memoryReport,    'memory-report'         ,opt=true,def=false)
    addr( 8)=bdShareVariable(systembd%memoryUnits,     'memory-units'          ,opt=true,def='Mb')                         !CHCK
    addr( 9)=bdShareVariable(systembd%memoryThreshold, 'memory-print-threshold',opt=true,def=8._rglu,expect='range(0:1024)')
    addr(10)=bdShareVariable(systembd%throughHeader,   'through-header'        ,opt=true,def='== Header of through file ==',register=true)
    addr(11)=bdShareVariable(systembd%throughFile,     'through-file'          ,opt=true,def='',register=true)
    addr(12)=bdShareVariable(systembd%throughEnable(1),'through-enable'        ,opt=true,def=false)
    addr(13)=bdShareVariable(systembd%throughPrefix,   'through-line-prefix'   ,opt=true,def='',register=true)
    !addr(10)=bdShareVariable(systembd%allowMutt,      'email-send'      ,opt=true,def=false)
    !addr(11)=bdShareVariable(systembd%muttSendTared,  'email-send-tared',opt=true,def=true)
    !addr(12)=bdShareVariable(systembd%muttDestination,'email-reciver'   ,opt=true,def='quant@mail.com')
    void=bdCollect('system',addr(1:13),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(iterationbd%chkStagnation      ,'enable-stagnation-check',opt=true,def=true)                 !CHCK
    addr(2)=bdShareVariable(iterationbd%chkDivergence      ,'enable-divergence-check',opt=true,def=true)
    addr(3)=bdShareVariable(iterationbd%chkStopIteration(1),'enable-stop-iteration'  ,opt=true,def=true)
    addr(4)=bdShareVariable(iterationbd%feelDivergence     ,'divergence-threshold'   ,opt=true,def=3._rglu  ,expect='range(2:8)',potentiate='+')
    addr(5)=bdShareVariable(iterationbd%feelStagnation     ,'stagnation-threshold'   ,opt=true,def=60._rglu ,expect='range(20:70)')
    addr(6)=bdShareVariable(iterationbd%thresholdStagnation,'stagnation-check-after' ,opt=true,def=200._rglu,expect='range(50:)')
    addr(7)=bdShareVariable(iterationbd%printFrequency     ,'print-frequency'        ,opt=true,def=0.5_rglu ,expect='range(0:)')
    addr(8)=bdShareVariable(iterationbd%printNotRearly     ,'print-upper-threshold'  ,opt=true,def=20._rglu ,expect='range(1:)')
    void=bdCollect('iteration',addr(1:8),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(geometrybd%symmetryAccount      ,'symmetry-account'   ,opt=true,def=true)
    addr(2)=bdShareVariable(geometrybd%searchPlanar(1)      ,'search-planarity'   ,opt=true,def=true)
    addr(3)=bdShareVariable(geometrybd%searchLinear(1)      ,'search-linearity'   ,opt=true,def=true)
    addr(4)=bdShareVariable(geometrybd%symmetryTolerance    ,'symmetry-threshold' ,opt=true,def=10._rglu,expect='range(1:15)',potentiate='-')
    addr(5)=bdShareVariable(geometrybd%searchTolerance      ,'search-threshold'   ,opt=true,def=5._rglu ,expect='range(1:15)',potentiate='-')
    addr(6)=bdShareVariable(geometrybd%randomDisplacement(1),'random-displacement',opt=true,def=0._rglu ,expect='range(0:5)')
    void=bdCollect('geometry',addr(1:6),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(statesbd%nStates,'excited-states-number',opt=true,def=0,expect='range(0:3)')
    addr(2)=bdShareVariable(statesbd%spin   ,'spin'                 ,opt=true,def=0,expect='range(0:1)')
    void=bdCollect('states',addr(1:2),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(polarizbd%scales   ,'components' ,opt=true,def='x',expect='list(x,y,z,xy,yx,xz,zx,yz,zy,xyz)')
    addr(2)=bdShareVariable(polarizbd%nPoints  ,'points'     ,opt=true,def=5,expect='range(1:9|odd)')
    addr(3)=bdShareVariable(polarizbd%maxPower ,'max-power'  ,opt=true,def=4,expect='range(1:4)')
    addr(4)=bdShareVariable(polarizbd%derivStep,'derive-step',opt=true,def=0.020_rglu,expect='range(0:1)')
    void=bdCollect('polariz',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(densitybd%dtype         ,'elements'        ,opt=true,def='none',expect='list(none,charges,orders,scharges,sorders,all)')
    addr(2)=bdShareVariable(densitybd%nPoints       ,'derive-points'   ,opt=true,def=3,expect='range(1:9|odd)')
    addr(3)=bdShareVariable(densitybd%prntAccuracy  ,'print-accuracy'  ,opt=true,def=5)
    addr(4)=bdShareVariable(densitybd%derivStep     ,'derive-step'     ,opt=true,def=0.025_rglu,expect='range(0:1)')
    addr(5)=bdShareVariable(densitybd%NOAnalize     ,'natural-orbitals',opt=true,def=false)
    addr(6)=bdShareVariable(densitybd%scharges      ,'selected-charges',opt=true,def='none',expect='any')
    addr(7)=bdShareVariable(densitybd%sorders       ,'selected-orders' ,opt=true,def='none',expect='any')
    addr(8)=bdShareVariable(densitybd%gcharges      ,'collect-charges' ,opt=true,def='none',expect='any')
    addr(9)=bdShareVariable(densitybd%forceNumerical,'force-numerical' ,opt=true,def=false)
    void=bdCollect('rdm',addr(1:9),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(coulsonbd%ctype         ,'elements'             ,opt=true,def='atom-atom',expect='list(none,atom-atom,atom-bond,bond-bond)')
    addr(2)=bdShareVariable(coulsonbd%selected      ,'selected-elements'    ,opt=true,def='none',expect='list(none,diagonals,offdiagonals)')
    addr(3)=bdShareVariable(coulsonbd%nPoints       ,'points'               ,opt=true,def=3,expect='range(1:9|odd)')
    addr(4)=bdShareVariable(coulsonbd%derivPower    ,'max-power'            ,opt=true,def=2,expect='range(1:4|even)')
    addr(5)=bdShareVariable(coulsonbd%derivStep     ,'derive-step'          ,opt=true,def=0.02_rglu,expect='range(0:1)')
    addr(6)=bdShareVariable(coulsonbd%prntAccuracy  ,'print-accuracy'       ,opt=true,def=5)
    addr(7)=bdShareVariable(coulsonbd%sdiagonals    ,'selected-diagonals'   ,opt=true,def='none',expect='any')
    addr(8)=bdShareVariable(coulsonbd%soffdiagonals ,'selected-offdiagonals',opt=true,def='none',expect='any')
    void=bdCollect('coulson',addr(1:8),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(hyperchargesbd%scales    ,'components'         ,opt=true,def='xyz',expect='list(x,y,z,xy,yx,xz,zx,yz,zy,xyz)')
    addr(2)=bdShareVariable(hyperchargesbd%naPoints  ,'derive-alpha-points',opt=true,def=5,expect='range(1:9|odd)')
    addr(3)=bdShareVariable(hyperchargesbd%nfPoints  ,'derive-field-points',opt=true,def=5,expect='range(1:9|odd)')
    addr(4)=bdShareVariable(hyperchargesbd%derivPower,'max-power'          ,opt=true,def=4,expect='range(1:4)')
    addr(5)=bdShareVariable(hyperchargesbd%derivaStep,'derive-alpha-step'  ,opt=true,def=0.025_rglu,expect='range(0:1)')
    addr(6)=bdShareVariable(hyperchargesbd%derivfStep,'derive-field-step'  ,opt=true,def=0.025_rglu,expect='range(0:1)')
    addr(7)=bdShareVariable(hyperchargesbd%scharges  ,'selected-charges'   ,opt=true,def='none',expect='any')
    addr(8)=bdShareVariable(hyperchargesbd%gcharges  ,'collect-charges'    ,opt=true,def='none',expect='any')
    addr(9)=bdShareVariable(hyperchargesbd%dtype     ,'elements'           ,opt=true,def='all',expect='list(none,all,scharges)')
    void=bdCollect('hypercharges',addr(1:9),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(cuebd%radius(0),'t1-limit-only-t1',opt=true,def=3)
    addr(2)=bdShareVariable(cuebd%radius(1),'t1-limit'        ,opt=true,def=-1)
    addr(3)=bdShareVariable(cuebd%radius(2),'t2-limit'        ,opt=true,def=-2)
    addr(4)=bdShareVariable(cuebd%radius(3),'t3-limit'        ,opt=true,def=-3)
    addr(5)=bdShareVariable(cuebd%sparse   ,'spare-algorithm' ,opt=true,def=false)
    addr(6)=bdShareVariable(cuebd%showBasis,'export-basis-xyz',opt=true,def=false)                                       !TODO
    void=bdCollect('cue',addr(1:6),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(fcibd%nSteps       ,'steps'         ,opt=true,def=5,expect='range(1:10)')
    addr(2)=bdShareVariable(fcibd%maxiters     ,'max-iterations',opt=true,def=200)
    addr(3)=bdShareVariable(fcibd%accuracy     ,'accuracy'      ,opt=true,def=12._rglu,expect='range(1:15)',potentiate='-')
    addr(4)=bdShareVariable(fcibd%zeroThreshold,'zero-threshold',opt=true,def=14._rglu,expect='range(1:15)',potentiate='-')
    void=bdCollect('fci',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(scfbd%maxiters         ,'max-iterations'       ,opt=true,def=400)
    addr(2)=bdShareVariable(scfbd%accuracy         ,'accuracy'             ,opt=true,def=12._rglu,expect='range(1:15)',potentiate='-')
    addr(3)=bdShareVariable(scfbd%iterStep         ,'iter-step'            ,opt=true,def=0.05_rglu,expect='range(0:1)')
    addr(4)=bdShareVariable(scfbd%guess            ,'density-guess'        ,opt=true,def='huckel',expect='list(huckel,unt,manual)')
    addr(5)=bdShareVariable(scfbd%keep             ,'keep-density'         ,opt=true,def=true)
    addr(6)=bdShareVariable(scfbd%exctype          ,'excited-states-method',opt=true,def='cis',expect='list(cis,rpa)')
    addr(7)=bdShareVariable(scfbd%achieveSolution  ,'force-locate-solution',opt=true,def=false)
    addr(8)=bdShareVariable(scfbd%iterStepVariation,'iter-step-variation'  ,opt=true,def=0.10_rglu,expect='range(0:1)')
    addr(9)=bdShareVariable(scfbd%iterStepChange   ,'variation-rate'       ,opt=true,def=0.005_rglu,expect='range(0:1)')
    void=bdCollect('scf',addr(1:9),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(localbd%maxiters,'max-iterations',opt=true,def=400)
    addr(2)=bdShareVariable(localbd%accuracy,'accuracy'      ,opt=true,def=8._rglu,expect='range(1:15)',potentiate='-')
    addr(3)=bdShareVariable(localbd%enabled ,'enabled'       ,opt=true,def=false)
    addr(4)=bdShareVariable(localbd%method  ,'procedure'     ,opt=true,def='pipek-mezey',expect='any')
    void=bdCollect('local',addr(1:4),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr( 1)=bdShareVariable(ccbd%projType      ,'projection-type'           ,opt=true,def='2-1',expect='list(1,2-1)')
    addr( 2)=bdShareVariable(ccbd%maxiters      ,'max-iterations'            ,opt=true,def=2000)
    addr( 3)=bdShareVariable(ccbd%accuracy      ,'accuracy'                  ,opt=true,def=12._rglu,expect='range(1:15)',potentiate='-')
    addr( 4)=bdShareVariable(ccbd%iterStep(1)   ,'t1-step'                   ,opt=true,def=1._rglu,expect='range(0:1)')
    addr( 5)=bdShareVariable(ccbd%iterStep(2)   ,'t2-step'                   ,opt=true,def=1._rglu,expect='range(0:1)')
    addr( 6)=bdShareVariable(ccbd%iterStep(3)   ,'t3-step'                   ,opt=true,def=1._rglu,expect='range(0:1)')
    addr( 7)=bdShareVariable(ccbd%forceSpin     ,'force-spin'                ,opt=true,def=false)
    addr( 8)=bdShareVariable(ccbd%diisEnabled   ,'diis-enabled'              ,opt=true,def=true)
    addr( 9)=bdShareVariable(ccbd%diisSteps     ,'diis-steps'                ,opt=true,def=20,expect='range(2:50)')
    addr(10)=bdShareVariable(ccbd%diisStorage   ,'diis-storage'              ,opt=true,def='ram',expect='list(ram,hdd)') !TODO
    addr(11)=bdShareVariable(ccbd%storeIntegrals,'store-integrals'           ,opt=true,def=true)                         !X3
    addr(12)=bdShareVariable(ccbd%wfSwitches    ,'wf-switches'               ,opt=true,def=int(2#111110, kind=iglu),expect='range(1:255)')
    addr(13)=bdShareVariable(ccbd%printThreshold,'amplitudes-print-threshold',opt=true,def=2._rglu,expect='range(1:15)',potentiate='-')
    void=bdCollect('coupled-cluster',addr(1:13),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr( 1)=bdShareVariable(lrbd%guess                 ,'guess'                   ,opt=true,def='cis',expect='list(unt,cis,rpa,manual)')
    addr( 2)=bdShareVariable(lrbd%guessThreshold        ,'guess-threshold'         ,opt=true,def=0.1_rglu,expect='range(0:1)')
    addr( 3)=bdShareVariable(lrbd%maxiters              ,'max-iterations'          ,opt=true,def=10000)
    addr( 4)=bdShareVariable(lrbd%accuracy              ,'accuracy'                ,opt=true,def=6._rglu,expect='range(1:15)',potentiate='-')
    addr( 5)=bdShareVariable(lrbd%iterStep(1)           ,'r1-step'                 ,opt=true,def=0.03_rglu,expect='range(0:1)')
    addr( 6)=bdShareVariable(lrbd%iterStep(2)           ,'r2-step'                 ,opt=true,def=0.03_rglu,expect='range(0:1)')
    addr( 7)=bdShareVariable(lrbd%orthogonalize         ,'orthogonalization'       ,opt=true,def=false)
    addr( 8)=bdShareVariable(lrbd%diisEnabled           ,'diis-enabled'            ,opt=true,def=true)
    addr( 9)=bdShareVariable(lrbd%diisSteps             ,'diis-steps'              ,opt=true,def=20,expect='range(2:50)')
    addr(10)=bdShareVariable(lrbd%diisStorage           ,'diis-storage'            ,opt=true,def='ram',expect='list(ram,hdd)')          !TODO
    addr(11)=bdShareVariable(lrbd%storeSolution(1)      ,'store-solution'          ,opt=true,def=true)
    addr(12)=bdShareVariable(lrbd%storeSolutionThreshold,'store-solution-threshold',opt=true,def=2._rglu,expect='range(1:4)',potentiate='-')
    addr(13)=bdShareVariable(lrbd%storeSolutionMode     ,'store-solution-mode'     ,opt=true,def='r1r2',expect='list(r1,r2,r1r2)')
    void=bdCollect('linear-response',addr(1:13),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    addr(1)=bdShareVariable(fieldbd%strength(1),'along-x',opt=true,def=0._rglu,expect='any')
    addr(2)=bdShareVariable(fieldbd%strength(2),'along-y',opt=true,def=0._rglu,expect='any')
    addr(3)=bdShareVariable(fieldbd%strength(3),'along-z',opt=true,def=0._rglu,expect='any')
    void=bdCollect('field',addr(1:3),bdstart,bdstop,bdcomment,bdaccord,bdseparator,false,0)

    void=bdCollect('molecule',addr(1:1),bdstart,bdstop,bdcomment,bdaccord,bdseparator,true,in)
    void=bdCollect('lrguess' ,addr(1:1),bdstart,bdstop,bdcomment,bdaccord,bdseparator,true,lrg)
    void=bdCollect('scfguess',addr(1:1),bdstart,bdstop,bdcomment,bdaccord,bdseparator,true,scfg)

    call definebdDescription

    return
    end subroutine definebd
