# generated automatically with command line:
# fmakefile.py --appname=helios --encoding=cp1251 --ignore-path=_tmp 
# paltform: Windows

NAME=helios.exe
COM=ifort
PFLAGS=/O1 /traceback /fpp /Qdiag-disable:7000,7954,8290,8291 /nologo
SFLAGS=/Qopenmp

OBJS= \
.\modules\mod_filecontrol.obj .\modules\mod_global.obj \
.\modules\mod_cc_sparse.obj .\modules\mod_math.obj .\modules\mod_orient.obj \
.\modules\mod_painter.obj .\modules\mod_sort.obj .\modules\mod_txtParser.obj \
.\modules\mod_approximate.obj .\modules\mod_derivative.obj \
.\modules\mod_parse.obj .\modules\mod_print.obj .\modules\mod_blockdata.obj \
.\modules\mod_helios.obj .\calculation\gIterator.obj \
.\calculation\putPoint.obj .\calculation\showPolarizability.obj \
.\cue\cueOrbitals.obj .\cue\cueStructureAnalysis.obj \
.\general\primaryInformation.obj .\general\symmertySettings.obj \
.\input\defineArguments.obj .\input\definebd.obj \
.\input\definebdDescription.obj .\input\inputHelp.obj .\input\parseInput.obj \
.\input\planarityCheck.obj .\input\readMoleculeInformation.obj \
.\modules\mod_excited.obj .\modules\mod_fci.obj .\modules\mod_huckel.obj \
.\modules\mod_scf.obj .\modules\mod_cc.obj \
.\coupledCluster\sparse-cue-ccsd.obj \
.\coupledCluster\ccd\proj-ccd-spat-hf-2.obj \
.\coupledCluster\ccd\proj-ccd-spin-hf-2.obj \
.\coupledCluster\ccs\proj-ccs-spat-cue-1.obj \
.\coupledCluster\ccs\proj-ccs-spin-cue-1.obj \
.\coupledCluster\ccsd\proj-ccsd-spat-cue-1.obj \
.\coupledCluster\ccsd\proj-ccsd-spat-cue-2.obj \
.\coupledCluster\ccsd\proj-ccsd-spat-hf-1.obj \
.\coupledCluster\ccsd\proj-ccsd-spat-hf-2.obj \
.\coupledCluster\ccsd\proj-ccsd-spin-cue-1-spare.obj \
.\coupledCluster\ccsd\proj-ccsd-spin-cue-1.obj \
.\coupledCluster\ccsd\proj-ccsd-spin-cue-2-spare.obj \
.\coupledCluster\ccsd\proj-ccsd-spin-cue-2.obj \
.\coupledCluster\ccsd\proj-ccsd-spin-hf-1.obj \
.\coupledCluster\ccsd\proj-ccsd-spin-hf-2.obj \
.\coupledCluster\ccsdt\proj-ccsdt-spin-cue-1.obj \
.\coupledCluster\ccsdt\proj-ccsdt-spin-cue-2.obj \
.\coupledCluster\ccsdt\proj-ccsdt-spin-cue-3.obj \
.\coupledCluster\ccsdt\proj-ccsdt-spin-hf-1.obj \
.\coupledCluster\ccsdt\proj-ccsdt-spin-hf-2.obj \
.\coupledCluster\ccsdt\proj-ccsdt-spin-hf-3.obj .\modules\mod_lrccsd.obj \
.\coupledCluster\lrccsd\proj-lrccsd-spin-hf-1.obj \
.\coupledCluster\lrccsd\proj-lrccsd-spin-hf-2.obj \
.\coupledCluster\lrccsd\proj-lrccsd-spin-hf-intermediates.obj \
.\general\versionControl.obj .\modules\mod_mbpt.obj .\modules\mod_property.obj \
.\main.obj


MODS= \
fcontrol.mod glob.mod coupledclustersparse.mod math.mod orientation.mod \
painter.mod sorts.mod txtparser.mod approx.mod derivat.mod argsparser.mod \
printmod.mod datablock.mod hdb.mod excitedstates.mod fci.mod huckel.mod \
scf.mod coupledcluster.mod lrccsdmodule.mod mbpt.mod property.mod


$(NAME): $(OBJS)
	$(COM) $(OBJS) $(SFLAGS) -o $(NAME)

.\modules\mod_filecontrol.obj: .\modules\mod_filecontrol.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_filecontrol.f90 -o .\modules\mod_filecontrol.obj
.\modules\mod_global.obj: .\modules\mod_global.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_global.f90 -o .\modules\mod_global.obj
.\modules\mod_cc_sparse.obj: .\modules\mod_global.obj .\modules\mod_cc_sparse.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_cc_sparse.f90 -o .\modules\mod_cc_sparse.obj
.\modules\mod_math.obj: .\modules\mod_global.obj .\modules\mod_math.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_math.f90 -o .\modules\mod_math.obj
.\modules\mod_orient.obj: .\modules\mod_global.obj .\modules\mod_orient.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_orient.f90 -o .\modules\mod_orient.obj
.\modules\mod_painter.obj: .\modules\mod_global.obj .\modules\mod_painter.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_painter.f90 -o .\modules\mod_painter.obj
.\modules\mod_sort.obj: .\modules\mod_global.obj .\modules\mod_sort.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_sort.f90 -o .\modules\mod_sort.obj
.\modules\mod_txtParser.obj: .\modules\mod_global.obj .\modules\mod_filecontrol.obj .\modules\mod_txtParser.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_txtParser.f90 -o .\modules\mod_txtParser.obj
.\modules\mod_approximate.obj: .\modules\mod_math.obj .\modules\mod_global.obj .\modules\mod_txtParser.obj .\modules\mod_approximate.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_approximate.f90 -o .\modules\mod_approximate.obj
.\modules\mod_derivative.obj: .\modules\mod_math.obj .\modules\mod_global.obj .\modules\mod_txtParser.obj .\modules\mod_derivative.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_derivative.f90 -o .\modules\mod_derivative.obj
.\modules\mod_parse.obj: .\modules\mod_global.obj .\modules\mod_txtParser.obj .\modules\mod_parse.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_parse.f90 -o .\modules\mod_parse.obj
.\modules\mod_print.obj: .\modules\mod_txtParser.obj .\modules\mod_sort.obj .\modules\mod_global.obj .\modules\mod_filecontrol.obj .\modules\mod_print.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_print.f90 -o .\modules\mod_print.obj
.\modules\mod_blockdata.obj: .\modules\mod_print.obj .\modules\mod_txtParser.obj .\modules\mod_global.obj .\modules\mod_filecontrol.obj .\modules\mod_blockdata.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_blockdata.f90 -o .\modules\mod_blockdata.obj
.\modules\mod_helios.obj: .\modules\mod_blockdata.obj .\modules\mod_txtParser.obj .\modules\mod_filecontrol.obj .\modules\mod_print.obj .\modules\mod_parse.obj .\modules\mod_global.obj .\modules\mod_helios.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_helios.f90 -o .\modules\mod_helios.obj
.\calculation\gIterator.obj: .\modules\mod_print.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\modules\mod_txtParser.obj .\calculation\gIterator.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\calculation\gIterator.f90 -o .\calculation\gIterator.obj
.\calculation\putPoint.obj: .\modules\mod_global.obj .\modules\mod_helios.obj .\calculation\putPoint.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\calculation\putPoint.f90 -o .\calculation\putPoint.obj
.\calculation\showPolarizability.obj: .\modules\mod_txtParser.obj .\modules\mod_print.obj .\modules\mod_derivative.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\calculation\showPolarizability.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\calculation\showPolarizability.f90 -o .\calculation\showPolarizability.obj
.\cue\cueOrbitals.obj: .\modules\mod_print.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\cue\cueOrbitals.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\cue\cueOrbitals.f90 -o .\cue\cueOrbitals.obj
.\cue\cueStructureAnalysis.obj: .\modules\mod_print.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\modules\mod_txtParser.obj .\cue\cueStructureAnalysis.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\cue\cueStructureAnalysis.f90 -o .\cue\cueStructureAnalysis.obj
.\general\primaryInformation.obj: .\modules\mod_txtParser.obj .\modules\mod_print.obj .\modules\mod_parse.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\general\primaryInformation.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\general\primaryInformation.f90 -o .\general\primaryInformation.obj
.\general\symmertySettings.obj: .\modules\mod_math.obj .\modules\mod_txtParser.obj .\modules\mod_print.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\general\symmertySettings.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\general\symmertySettings.f90 -o .\general\symmertySettings.obj
.\input\defineArguments.obj: .\modules\mod_parse.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\input\defineArguments.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\input\defineArguments.f90 -o .\input\defineArguments.obj
.\input\definebd.obj: .\modules\mod_global.obj .\modules\mod_blockdata.obj .\modules\mod_helios.obj .\modules\mod_filecontrol.obj .\input\definebd.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\input\definebd.f90 -o .\input\definebd.obj
.\input\definebdDescription.obj: .\modules\mod_global.obj .\modules\mod_blockdata.obj .\modules\mod_helios.obj .\modules\mod_txtParser.obj .\input\definebdDescription.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\input\definebdDescription.f90 -o .\input\definebdDescription.obj
.\input\inputHelp.obj: .\modules\mod_global.obj .\modules\mod_blockdata.obj .\modules\mod_helios.obj .\modules\mod_txtParser.obj .\input\inputHelp.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\input\inputHelp.f90 -o .\input\inputHelp.obj
.\input\parseInput.obj: .\modules\mod_blockdata.obj .\modules\mod_txtParser.obj .\modules\mod_print.obj .\modules\mod_parse.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\input\parseInput.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\input\parseInput.f90 -o .\input\parseInput.obj
.\input\planarityCheck.obj: .\modules\mod_orient.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\input\planarityCheck.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\input\planarityCheck.f90 -o .\input\planarityCheck.obj
.\input\readMoleculeInformation.obj: .\modules\mod_print.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\modules\mod_txtParser.obj .\input\readMoleculeInformation.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\input\readMoleculeInformation.f90 -o .\input\readMoleculeInformation.obj
.\modules\mod_excited.obj: .\modules\mod_math.obj .\modules\mod_global.obj .\modules\mod_print.obj .\modules\mod_helios.obj .\modules\mod_excited.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_excited.f90 -o .\modules\mod_excited.obj
.\modules\mod_fci.obj: .\modules\mod_math.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\modules\mod_filecontrol.obj .\modules\mod_fci.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_fci.f90 -o .\modules\mod_fci.obj
.\modules\mod_huckel.obj: .\modules\mod_math.obj .\modules\mod_global.obj .\modules\mod_print.obj .\modules\mod_helios.obj .\modules\mod_huckel.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_huckel.f90 -o .\modules\mod_huckel.obj
.\modules\mod_scf.obj: .\modules\mod_math.obj .\modules\mod_txtParser.obj .\modules\mod_print.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\modules\mod_scf.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_scf.f90 -o .\modules\mod_scf.obj
.\modules\mod_cc.obj: .\modules\mod_math.obj .\modules\mod_cc_sparse.obj .\modules\mod_scf.obj .\modules\mod_txtParser.obj .\modules\mod_print.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\modules\mod_cc.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_cc.f90 -o .\modules\mod_cc.obj
.\coupledCluster\sparse-cue-ccsd.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\modules\mod_cc_sparse.obj .\modules\mod_helios.obj .\coupledCluster\sparse-cue-ccsd.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\sparse-cue-ccsd.f90 -o .\coupledCluster\sparse-cue-ccsd.obj
.\coupledCluster\ccd\proj-ccd-spat-hf-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccd\proj-ccd-spat-hf-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccd\proj-ccd-spat-hf-2.f90 -o .\coupledCluster\ccd\proj-ccd-spat-hf-2.obj
.\coupledCluster\ccd\proj-ccd-spin-hf-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccd\proj-ccd-spin-hf-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccd\proj-ccd-spin-hf-2.f90 -o .\coupledCluster\ccd\proj-ccd-spin-hf-2.obj
.\coupledCluster\ccs\proj-ccs-spat-cue-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccs\proj-ccs-spat-cue-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccs\proj-ccs-spat-cue-1.f90 -o .\coupledCluster\ccs\proj-ccs-spat-cue-1.obj
.\coupledCluster\ccs\proj-ccs-spin-cue-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccs\proj-ccs-spin-cue-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccs\proj-ccs-spin-cue-1.f90 -o .\coupledCluster\ccs\proj-ccs-spin-cue-1.obj
.\coupledCluster\ccsd\proj-ccsd-spat-cue-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spat-cue-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spat-cue-1.f90 -o .\coupledCluster\ccsd\proj-ccsd-spat-cue-1.obj
.\coupledCluster\ccsd\proj-ccsd-spat-cue-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spat-cue-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spat-cue-2.f90 -o .\coupledCluster\ccsd\proj-ccsd-spat-cue-2.obj
.\coupledCluster\ccsd\proj-ccsd-spat-hf-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spat-hf-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spat-hf-1.f90 -o .\coupledCluster\ccsd\proj-ccsd-spat-hf-1.obj
.\coupledCluster\ccsd\proj-ccsd-spat-hf-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spat-hf-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spat-hf-2.f90 -o .\coupledCluster\ccsd\proj-ccsd-spat-hf-2.obj
.\coupledCluster\ccsd\proj-ccsd-spin-cue-1-spare.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\modules\mod_cc_sparse.obj .\coupledCluster\ccsd\proj-ccsd-spin-cue-1-spare.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spin-cue-1-spare.f90 -o .\coupledCluster\ccsd\proj-ccsd-spin-cue-1-spare.obj
.\coupledCluster\ccsd\proj-ccsd-spin-cue-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spin-cue-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spin-cue-1.f90 -o .\coupledCluster\ccsd\proj-ccsd-spin-cue-1.obj
.\coupledCluster\ccsd\proj-ccsd-spin-cue-2-spare.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\modules\mod_cc_sparse.obj .\coupledCluster\ccsd\proj-ccsd-spin-cue-2-spare.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spin-cue-2-spare.f90 -o .\coupledCluster\ccsd\proj-ccsd-spin-cue-2-spare.obj
.\coupledCluster\ccsd\proj-ccsd-spin-cue-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spin-cue-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spin-cue-2.f90 -o .\coupledCluster\ccsd\proj-ccsd-spin-cue-2.obj
.\coupledCluster\ccsd\proj-ccsd-spin-hf-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spin-hf-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spin-hf-1.f90 -o .\coupledCluster\ccsd\proj-ccsd-spin-hf-1.obj
.\coupledCluster\ccsd\proj-ccsd-spin-hf-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsd\proj-ccsd-spin-hf-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsd\proj-ccsd-spin-hf-2.f90 -o .\coupledCluster\ccsd\proj-ccsd-spin-hf-2.obj
.\coupledCluster\ccsdt\proj-ccsdt-spin-cue-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-1.f90 -o .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-1.obj
.\coupledCluster\ccsdt\proj-ccsdt-spin-cue-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-2.f90 -o .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-2.obj
.\coupledCluster\ccsdt\proj-ccsdt-spin-cue-3.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-3.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-3.f90 -o .\coupledCluster\ccsdt\proj-ccsdt-spin-cue-3.obj
.\coupledCluster\ccsdt\proj-ccsdt-spin-hf-1.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-1.f90 -o .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-1.obj
.\coupledCluster\ccsdt\proj-ccsdt-spin-hf-2.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-2.f90 -o .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-2.obj
.\coupledCluster\ccsdt\proj-ccsdt-spin-hf-3.obj: .\modules\mod_cc.obj .\modules\mod_global.obj .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-3.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-3.f90 -o .\coupledCluster\ccsdt\proj-ccsdt-spin-hf-3.obj
.\modules\mod_lrccsd.obj: .\modules\mod_math.obj .\modules\mod_cc.obj .\modules\mod_txtParser.obj .\modules\mod_scf.obj .\modules\mod_print.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\modules\mod_lrccsd.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_lrccsd.f90 -o .\modules\mod_lrccsd.obj
.\coupledCluster\lrccsd\proj-lrccsd-spin-hf-1.obj: .\modules\mod_lrccsd.obj .\modules\mod_global.obj .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-1.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-1.f90 -o .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-1.obj
.\coupledCluster\lrccsd\proj-lrccsd-spin-hf-2.obj: .\modules\mod_lrccsd.obj .\modules\mod_global.obj .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-2.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-2.f90 -o .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-2.obj
.\coupledCluster\lrccsd\proj-lrccsd-spin-hf-intermediates.obj: .\modules\mod_lrccsd.obj .\modules\mod_global.obj .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-intermediates.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-intermediates.f90 -o .\coupledCluster\lrccsd\proj-lrccsd-spin-hf-intermediates.obj
.\general\versionControl.obj: .\modules\mod_math.obj .\modules\mod_lrccsd.obj .\modules\mod_blockdata.obj .\modules\mod_filecontrol.obj .\modules\mod_txtParser.obj .\modules\mod_orient.obj .\modules\mod_derivative.obj .\modules\mod_painter.obj .\modules\mod_fci.obj .\modules\mod_parse.obj .\modules\mod_helios.obj .\modules\mod_print.obj .\modules\mod_approximate.obj .\modules\mod_global.obj .\general\versionControl.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\general\versionControl.f90 -o .\general\versionControl.obj
.\modules\mod_mbpt.obj: .\modules\mod_print.obj .\modules\mod_global.obj .\modules\mod_helios.obj .\modules\mod_scf.obj .\modules\mod_mbpt.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_mbpt.f90 -o .\modules\mod_mbpt.obj
.\modules\mod_property.obj: .\modules\mod_cc.obj .\modules\mod_huckel.obj .\modules\mod_lrccsd.obj .\modules\mod_txtParser.obj .\modules\mod_scf.obj .\modules\mod_print.obj .\modules\mod_mbpt.obj .\modules\mod_fci.obj .\modules\mod_helios.obj .\modules\mod_excited.obj .\modules\mod_global.obj .\modules\mod_property.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\modules\mod_property.f90 -o .\modules\mod_property.obj
.\main.obj: .\modules\mod_cc.obj .\modules\mod_lrccsd.obj .\modules\mod_blockdata.obj .\modules\mod_txtParser.obj .\modules\mod_print.obj .\modules\mod_sort.obj .\modules\mod_property.obj .\modules\mod_helios.obj .\modules\mod_global.obj .\main.f90
	$(COM) -c $(PFLAGS) $(SFLAGS) .\main.f90 -o .\main.obj

.PHONY: rm_objs rm_mods rm_app clean cleanall remake build set_env

rm_objs:
	rm -f $(OBJS)

rm_mods:
	rm -f $(MODS)

rm_app:
	rm -f $(NAME)

clean:
	$(MAKE) rm_objs
	$(MAKE) rm_mods

cleanall:
	$(MAKE) clean
	$(MAKE) rm_app

remake:
	$(MAKE) cleanall
	$(MAKE)

build:
	$(MAKE) cleanall
	$(MAKE)
	$(MAKE) clean

