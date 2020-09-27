# 27-Sep-2020 Changed directory structure

0) Moved to src directory
1) Removed several useless files

# 27-Sep-2020 Intermediate state for further work. Project requires cleaning.

0) MINOR: Fine print
   * cue/cueStructureAnalysis
   * input/parseInput
   * modules/mod_huckel - printHuckelSolution
   * modules/mod_localize
   * modules/mod_scf - printSCFSolution
1) MINOR: plugs
   * modules/mod_cc - printCCSolution
   * modules/mod_fci - printFCISolution
2) MINOR: Remove restriction on bond kind - 'nonset' available
   * general/symmetrySettings
   * input/readMoleculeInformation
3) MINOR: Prepare for global structure rework
4) MINOR: Wave-function analysis
5) REFACTOR: Input - through, parametrization
   * input/checkInput
   * input/definebd
   * input/definebdDescription
   * input/parseInput
   * input/readMoleculeInformation
   * modules/mod_helios
   * modules/mod_property
6) REFACTOR: configurations moved to include files
   * coupledCluster/includes/cc_configurations
   * coupledCluster/includes/lr_configurations
7) FIX: Polarizability output
8) TYPOS
   * calculation/gIterator
   * modules/mod_global
   * modules/mod_helios
9) TODO: Tests
