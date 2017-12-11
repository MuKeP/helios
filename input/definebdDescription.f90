    subroutine definebdDescription

    use glob,      only: iglu,void
    use hdb,       only: generalbd,systembd,iterationbd,geometrybd,statesbd,polarizbd,diisbd
    use hdb,       only: densitybd,coulsonbd,hyperchargesbd,cuebd,fcibd,scfbd,lrbd,ccbd,pipekbd
    use datablock, only: bdAddDescription,bdVariableAddDescription
    use txtParser, only: tpTranslateEscapes

    implicit none


    void=bdAddDescription('general','contains general semi-empirical parameters and calculation settings.')
    void=bdVariableAddDescription(loc(generalbd%methods),&
                                  tpTranslateEscapes('specifies method(s) to be used for calculation:\n'//&
                                                     'huckel           - Molecular orbitals Huckel method\n'//&
                                                     'rhf              - Restricted Hartree-Fock self-consistent field\n'//&
                                                     'fci              - Full configuration interaction method\n'//&
                                                     'mp2              - Moller-Plesset 2nd order perturbation theory\n'//&
                                                     'mp3              - Moller-Plesset 3rd order perturbation theory\n'//&
                                                     'u-ccd            - unrelaxed CC with doubles\n'//&
                                                     'r-ccd            - relaxed CC with doubles\n'//&
                                                     'u-ccsd           - unrelaxed CC with singles and doubles\n'//&
                                                     'r-ccsd           - relaxed CC with singles and doubles\n'//&
                                                     'r-ccsd(t)        - relaxed CCSD with non-iterative triples correction\n'//&
                                                     'u-ccsdt          - unrelaxed CC with singles, doubles and triples\n'//&
                                                     'r-ccsdt          - relaxed CC with singles, doubles and triples\n'//&
                                                     'cue-ccs          - cue CC with singles (relaxation of cue reference)\n'//&
                                                     'cue-ccsd         - cue CC with singles and doubles\n'//&
                                                     'cue-ccsdt        - cue CC with singles, doubles and triples') )

    void=bdVariableAddDescription(loc(generalbd%task),&
                                  tpTranslateEscapes('specifies property to be calculated:\n'//&
                                                     'energy           - energy of the system\n'//&
                                                     'wf-analysis      - perform wave-function analysis (in case of CC methods)\n'//&
                                                     'density          - RDM1 elements calculation\n'//&
                                                     'polarizability   - polarizabilities and hyperpolarizabilities\n'//&
                                                     'hypercharges     - hypercharges and hypermoments (according to Bredas\n'//&
                                                     '                   real-space finite field approach)\n'//&
                                                     'coulson          - Coulson atom-atom, atom-bond, bond-bond polarizabilities') )

    void=bdVariableAddDescription(loc(generalbd%coulombType),&
                                  tpTranslateEscapes('specifies empirical formula to describe two-center Coulomb integrals:\n'//&
                                                     'ohno-klopman     - Ohno-Klopman formula\n'//&
                                                     'mataga-nishimoto - Mataga-Nishimoto formula\n'//&
                                                     'hubbard          - use Hubbard approximation (only diagonal elements)') )

    void=bdVariableAddDescription(loc(generalbd%alternation),'specifies the values of single/double bond resonant integrals alternation')

    void=bdAddDescription('system','contains global system settings for calculation.')
    void=bdVariableAddDescription(loc(systembd%memory)      ,'specifies maximum amount of shared RAM (in MB, i.e. x1024*1024 bytes)')
    void=bdVariableAddDescription(loc(systembd%nNodes)      ,'specifies maximum number of cores to use for parallel execution')
    void=bdVariableAddDescription(loc(systembd%allowRestart),'store intermediate data to restart in case of unexpected termination')
    void=bdVariableAddDescription(loc(systembd%verboselvl)  ,'specifies admissible verbose level for output')
    void=bdVariableAddDescription(loc(systembd%ignoreSIGHUP),'specifies action on sighup recieved (only POSIX systems)')

    void=bdAddDescription('iteration','contains general settings for iteration procedures performed during calculation')
    void=bdVariableAddDescription(loc(iterationbd%chkStagnation)      ,'enables iteration procedure stagnation determination')
    void=bdVariableAddDescription(loc(iterationbd%chkDivergence)      ,'enables iteration procedure divergence determination')
    void=bdVariableAddDescription(loc(iterationbd%chkStopIteration(1)),'enables manual termination of iterations with signal (only POSIX systems)')
    void=bdVariableAddDescription(loc(iterationbd%feelDivergence)     ,'specifies accuracy threshold to determine divergence')
    void=bdVariableAddDescription(loc(iterationbd%feelStagnation)     ,'specifies iteration part of successful iterations to determine stagnation')
    void=bdVariableAddDescription(loc(iterationbd%thresholdStagnation),'specifies iteration number to start stagnation determination after')
    void=bdVariableAddDescription(loc(iterationbd%printFrequency)     ,'specifies print frequency during iterations (in seconds)')
    void=bdVariableAddDescription(loc(iterationbd%printNotRearly)     ,'specifies iteration frequency to be print not rearly than (in iterations)')

    void=bdAddDescription('geometry','contains settings for geometry manipulations')
    void=bdVariableAddDescription(loc(geometrybd%symmetryAccount)      ,'enables account of symmetry')
    void=bdVariableAddDescription(loc(geometrybd%searchPlanar(1))      ,'enables planarity determination')
    void=bdVariableAddDescription(loc(geometrybd%searchLinear(1))      ,'enables linearity determination')
    void=bdVariableAddDescription(loc(geometrybd%symmetryTolerance)    ,'specifies symmetry threshold (for comparison of eigenvectors)')
    void=bdVariableAddDescription(loc(geometrybd%searchTolerance)      ,'specifies symmerty search threshold (for planarity/linearity determination)')
    void=bdVariableAddDescription(loc(geometrybd%randomDisplacement(1)),'specifies the standard deviation of displacement made for coordinates')

    void=bdAddDescription('states','contains settings for excited states computation')
    void=bdVariableAddDescription(loc(statesbd%nStates),'specifies the number of excited states to calculate')
    void=bdVariableAddDescription(loc(statesbd%spin)   ,'specifies total spin of the state')

    void=bdAddDescription('polariz','contains Finite Field settings for polarizability calculation')
    void=bdVariableAddDescription(loc(polarizbd%scales)   ,'specifies components of (hyper)polarizabilities to be calculated')
    void=bdVariableAddDescription(loc(polarizbd%nPoints)  ,'specifies number of points to be used for derivation')
    void=bdVariableAddDescription(loc(polarizbd%maxPower) ,'specifies maximum power of derivative to be obtained')
    void=bdVariableAddDescription(loc(polarizbd%derivStep),'specifies derivation step for Finite Field method')

    void=bdAddDescription('rdm','contains settings for RDM1 preparation')
    void=bdVariableAddDescription(loc(densitybd%dtype)       ,'specifies RDM1 elements to be calculated')
    void=bdVariableAddDescription(loc(densitybd%nPoints)     ,'specifies number of points to be used for derivation')
    void=bdVariableAddDescription(loc(densitybd%prntAccuracy),'specifies print accuracy for RDM1')
    void=bdVariableAddDescription(loc(densitybd%derivStep)   ,'specifies derivation step')
    void=bdVariableAddDescription(loc(densitybd%NOAnalize)   ,'enables Natural Orbitals production (requires elements=all)')

    void=bdAddDescription('coulson','contains settings for Coulson polarizabilities calculation')
    void=bdVariableAddDescription(loc(coulsonbd%ctype)     ,'specifies type of polarizability to be calculated')
    void=bdVariableAddDescription(loc(coulsonbd%nPoints)   ,'specifies number of points to be used for derivation')
    void=bdVariableAddDescription(loc(coulsonbd%derivPower),'specifies maximum power of derivative to be obtained')
    void=bdVariableAddDescription(loc(coulsonbd%derivStep) ,'specifies derivation step for Finite Field method')

    void=bdAddDescription('hypercharges','contains settings for Bredas real-space Finite Field approach')
    void=bdVariableAddDescription(loc(hyperchargesbd%scales)    ,'specifies components of hypercharges to be calculated')
    void=bdVariableAddDescription(loc(hyperchargesbd%naPoints)  ,'specifies number of points for derivative by Coulomb integral')
    void=bdVariableAddDescription(loc(hyperchargesbd%nfPoints)  ,'specifies number of points for derivative by applied field')
    void=bdVariableAddDescription(loc(hyperchargesbd%derivPower),'specifies maximum power of derivative to be obtained')
    void=bdVariableAddDescription(loc(hyperchargesbd%derivaStep),'specifies derivation step for derivative by Coulomb integral')
    void=bdVariableAddDescription(loc(hyperchargesbd%derivfStep),'specifies derivation step for derivative by applied field')

    void=bdAddDescription('cue','contains settings of the cue basis')
    void=bdVariableAddDescription(loc(cuebd%radius(0)),'specifies limitation for t1^3 and t1^4 diagrams')
    void=bdVariableAddDescription(loc(cuebd%radius(1)),'specifies limitation for all t1 amplitudes')
    void=bdVariableAddDescription(loc(cuebd%radius(2)),'specifies limitation for all t2 amplitudes')
    void=bdVariableAddDescription(loc(cuebd%radius(3)),'specifies limitation for all t3 amplitudes')
    void=bdVariableAddDescription(loc(cuebd%sparse)   ,'enables spare algorithm (for cue-ccsd only)')
    void=bdVariableAddDescription(loc(cuebd%showBasis),'enables generation of xyz file with cue MO centroids')

    void=bdAddDescription('fci','contains settings for FCI calculation')
    void=bdVariableAddDescription(loc(fcibd%nSteps)       ,'specifies number of points to construct Krylov space')
    void=bdVariableAddDescription(loc(fcibd%maxiters)     ,'specifies maximum number of iterations to be performed')
    void=bdVariableAddDescription(loc(fcibd%accuracy)     ,'specifies accuracy to be reached')
    void=bdVariableAddDescription(loc(fcibd%zeroThreshold),'specifies threshold for zero-elements detection')

    void=bdAddDescription('scf','contains settings for SCF procedure')
    void=bdVariableAddDescription(loc(scfbd%maxiters),'specifies maximum number of iterations to be performed')
    void=bdVariableAddDescription(loc(scfbd%accuracy),'specifies accuracy to be reached')
    void=bdVariableAddDescription(loc(scfbd%iterStep),'specifies iteration step')
    void=bdVariableAddDescription(loc(scfbd%guess)   ,tpTranslateEscapes('specifies guess type for SCF:\n'//&
                                                                         'huckel - start with huckel density\n'//&
                                                                         'unt    - start with unitary matrix\n'//&
                                                                         'manual - manually defined guess (see scfguess block)') )
    void=bdVariableAddDescription(loc(scfbd%exctype) ,tpTranslateEscapes('specifies excited states calculation method:\n'//&
                                                                         'cis    - configuration interaction singles\n'//&
                                                                         'rpa    - random phase appriximation') )

    void=bdAddDescription('local','contains settings for localization procedure')
    void=bdVariableAddDescription(loc(pipekbd%enabled) ,'enables localization procedure')
    void=bdVariableAddDescription(loc(pipekbd%maxiters),'specifies maximum number of iterations to be performed')
    void=bdVariableAddDescription(loc(pipekbd%accuracy),'specifies accuracy to be reached')

    void=bdAddDescription('coupled-cluster','contains settings for coupled-cluster methods')
    void=bdVariableAddDescription(loc(ccbd%projType),tpTranslateEscapes('specifies projection type onto doubles (for spatial orbitals):\n'//&
                                                                        '1   - standard approach  [delta(ijab)]\n'//&
                                                                        '2-1 - biorthogonal basis [2*delta(ijab)+delta(ijba)]') )
    void=bdVariableAddDescription(loc(ccbd%maxiters)      ,'specifies maximum number of iterations to be performed')
    void=bdVariableAddDescription(loc(ccbd%accuracy)      ,'specifies accuracy to be reached')
    void=bdVariableAddDescription(loc(ccbd%iterStep(1))   ,'specifies iteration step for t1 amplitudes (for non-diis calculations)')
    void=bdVariableAddDescription(loc(ccbd%iterStep(2))   ,'specifies iteration step for t2 amplitudes (for non-diis calculations)')
    void=bdVariableAddDescription(loc(ccbd%iterStep(3))   ,'specifies iteration step for t3 amplitudes (for non-diis calculations)')
    void=bdVariableAddDescription(loc(ccbd%forceSpin)     ,'compels to use spin-orbitals even when it is not necessary')
    void=bdVariableAddDescription(loc(ccbd%storeIntegrals),'specifies whether to prepare integrals once or calculate it on demand (non-spare)')

    void=bdAddDescription('linear-response','contains settings for linear-response coupled-cluster methods')
    void=bdVariableAddDescription(loc(lrbd%guess),tpTranslateEscapes('specifies guess type for linear response procedure:\n'//&
                                                                     'unt    - t1(HOMO,LUMO) amplitude set to be unity\n'//&
                                                                     'cis    - use CIS solution as guess\n'//&
                                                                     'rpa    - use RPA solution as guess\n'//&
                                                                     'manual - manually defined guess (see lrguess block)') )
    void=bdVariableAddDescription(loc(lrbd%guessThreshold),'specifies guess threshold for selection (guess={cis,rpa})')
    void=bdVariableAddDescription(loc(lrbd%maxiters)      ,'specifies maximum number of iterations to be performed')
    void=bdVariableAddDescription(loc(lrbd%accuracy)      ,'specifies accuracy to be reached')
    void=bdVariableAddDescription(loc(lrbd%iterStep(1))   ,'specifies iteration step for r1 amplitudes (for non-diis calculations)')
    void=bdVariableAddDescription(loc(lrbd%iterStep(2))   ,'specifies iteration step for r2 amplitudes (for non-diis calculations)')
    void=bdVariableAddDescription(loc(lrbd%orthogonalize) ,'enables orthogonalization of states')

    void=bdAddDescription('diis','contains settings for DIIS procedure')
    void=bdVariableAddDescription(loc(diisbd%storage),tpTranslateEscapes('specifies where to store amplitudes and corrections vectors:\n'//&
                                                                         'ram - RAM\n'//&
                                                                         'hdd - hard disk (slow)') )
    void=bdVariableAddDescription(loc(diisbd%steps)  ,'specifies steps number for interpolation')
    void=bdVariableAddDescription(loc(diisbd%enabled),'enables DIIS procedure')

    void=bdAddDescription('molecule',tpTranslateEscapes('contains molecular information in the following format:\n\n'//&
                                                        'MOLECULE NAME\n'//&
                                                        'NUMBER OF ATOMS\n'//&
                                                        'NUMBER OF BONDS\n'//&
                                                        'NUMBER OF ELECTRONS IN ALPHA SHELL\n'//&
                                                        'NUMBER OF ELECTRONS IN BETA SHELL\n\n'//&
                                                        'AT1  X1  Y1  Z1  IP1  CI1  EL1\n'//&
                                                        'AT2  X2  Y2  Z2  IP2  CI2  EL2\n'//&
                                                        '...\n\n'//&
                                                        'AT = atom type (C,N,O,S,...)\n'//&
                                                        'X,Y,Z = Cartesian coordinates (angstroms)\n'//&
                                                        'IP = ionization potential (eV)\n'//&
                                                        'CI = one-center Coulomb integral (eV)\n'//&
                                                        'EL = number of electrons provided to the pi-shell by atom\n'//&
                                                        '     i.e. one from C, one or electron pair from N, etc.\n\n'//&
                                                        'FA1  SA1  RI1  K1  D1\n'//&
                                                        'FA2  SA2  RI2  K2  D2\n'//&
                                                        '...\n\n'//&
                                                        'FA = first  atom of the bond\n'//&
                                                        'SA = second atom of the bond\n'//&
                                                        'RI = resonant integral (eV)\n'//&
                                                        'K = kind of the bond (single/double)\n'//&
                                                        'D = bond length (angstroms)\n\n'//&
                                                        '* any blank lines/spaces/tabs are allowed\n'//&
                                                       '* every atom and every bond must be written on the new line') )

    void=bdAddDescription('lrguess',tpTranslateEscapes('contains guess information for LR methods in the following format:\n\n'//&
                                                       'R1 amplitudes        <= key word\n'//&
                                                       'i1  a1  r1(i1,a1)\n'//&
                                                       'i2  a2  r1(i2,a2)\n'//&
                                                       '...\n\n'//&
                                                       'i = occupied orbital\n'//&
                                                       'a = vacant orbital\n'//&
                                                       'r1(i,a) = value of corresponding amplitude\n\n'//&
                                                       'R2 section        <= key word\n'//&
                                                       'i1  j1  a1  b1  r2(i1,j1,a1,b1)\n'//&
                                                       'i2  j2  a2  b2  r2(i2,j2,a2,b2)\n'//&
                                                       '...\n\n'//&
                                                       'i, j = occupied orbitals\n'//&
                                                       'a, b = vacant orbitals\n'//&
                                                       'r2(i,j,a,b) = value of corresponding amplitude\n\n'//&
                                                       '* key words must start every section\n'//&
                                                       '* violation of ranges will lead to error during read (N - number of atoms)\n'//&
                                                       '  - occupied indices must be in the range [1,N]\n'//&
                                                       '  - vacant indices must be in the range [N+1,2N]\n'//&
                                                       '* any sections could be omitted') )

    void=bdAddDescription('scfguess',tpTranslateEscapes('contains guess information for SCF procedure in the following format:\n\n'//&
                                                        'mu1  nu1  P(mu1,nu1)\n'//&
                                                        'mu2  nu2  P(mu2,nu2)\n'//&
                                                        '...\n\n'//&
                                                        'mu, nu = indices of RDM1 element to be set\n'//&
                                                        'P(mu,nu) = RDM1 element value to be set\n'//&
                                                        '* violation of ranges will lead to error during read (N - number of atoms)\n'//&
                                                        '  - indices must be in the range [1..N]') )

    return
    end subroutine definebdDescription