# Microsoft Developer Studio Project File - Name="helios" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=helios - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "helios.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "helios.mak" CFG="helios - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "helios - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "helios - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "helios - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# SUBTRACT F90 /warn:unused
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "helios - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# SUBTRACT F90 /warn:unused
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "helios - Win32 Release"
# Name "helios - Win32 Debug"
# Begin Group "modules"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\mod_approximate.f90
DEP_F90_MOD_A=\
	".\Release\glob.mod"\
	".\Release\math.mod"\
	".\Release\txtParser.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_blockdata.f90
DEP_F90_MOD_B=\
	".\Release\fcontrol.mod"\
	".\Release\glob.mod"\
	".\Release\printmod.mod"\
	".\Release\txtParser.mod"\
	

!IF  "$(CFG)" == "helios - Win32 Release"

!ELSEIF  "$(CFG)" == "helios - Win32 Debug"

# SUBTRACT F90 /check:bounds

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mod_cc.f90
DEP_F90_MOD_C=\
	".\Release\coupledClusterSparse.mod"\
	".\Release\hdb.mod"\
	".\Release\math.mod"\
	".\Release\scf.mod"\
	

!IF  "$(CFG)" == "helios - Win32 Release"

# ADD F90 /warn:unused

!ELSEIF  "$(CFG)" == "helios - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mod_cc_sparse.f90
DEP_F90_MOD_CC=\
	".\Release\hdb.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_derivative.f90
DEP_F90_MOD_D=\
	".\Release\glob.mod"\
	".\Release\math.mod"\
	".\Release\txtParser.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_fci.f90
DEP_F90_MOD_F=\
	".\Release\fcontrol.mod"\
	".\Release\glob.mod"\
	".\Release\math.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_filecontrol.f90
# End Source File
# Begin Source File

SOURCE=.\mod_global.f90
NODEP_F90_MOD_G=\
	".\Release\ifport.mod"\
	".\Release\omp_lib.h"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_helios.f90
DEP_F90_MOD_H=\
	".\Release\argsParser.mod"\
	".\Release\datablock.mod"\
	".\Release\fcontrol.mod"\
	".\Release\glob.mod"\
	".\Release\printmod.mod"\
	".\Release\txtParser.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_lrccsd.f90
DEP_F90_MOD_L=\
	".\Release\glob.mod"\
	
NODEP_F90_MOD_L=\
	".\Release\omp_lib.h"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_math.f90
DEP_F90_MOD_M=\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_orient.f90
DEP_F90_MOD_O=\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_painter.f90
DEP_F90_MOD_P=\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_parse.f90
DEP_F90_MOD_PA=\
	".\Release\glob.mod"\
	".\Release\txtParser.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_print.f90
DEP_F90_MOD_PR=\
	".\Release\fcontrol.mod"\
	".\Release\glob.mod"\
	".\Release\txtParser.mod"\
	

!IF  "$(CFG)" == "helios - Win32 Release"

!ELSEIF  "$(CFG)" == "helios - Win32 Debug"

# ADD F90 /warn:unused

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\mod_scf.f90
DEP_F90_MOD_S=\
	".\Release\hdb.mod"\
	".\Release\math.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\mod_txtParser.f90
DEP_F90_MOD_T=\
	".\Release\fcontrol.mod"\
	".\Release\glob.mod"\
	
# End Source File
# End Group
# Begin Group "general"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\main.f90
DEP_F90_MAIN_=\
	".\Release\coupledCluster.mod"\
	".\Release\hdb.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\primaryInformation.f90
DEP_F90_PRIMA=\
	".\Release\hdb.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\symmertySettings.f90
DEP_F90_SYMME=\
	".\Release\hdb.mod"\
	".\Release\math.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\versionControl.f90
DEP_F90_VERSI=\
	".\Release\approx.mod"\
	".\Release\argsParser.mod"\
	".\Release\datablock.mod"\
	".\Release\derivat.mod"\
	".\Release\fciModule.mod"\
	".\Release\fcontrol.mod"\
	".\Release\glob.mod"\
	".\Release\hdb.mod"\
	".\Release\lrccsdModule.mod"\
	".\Release\math.mod"\
	".\Release\orientation.mod"\
	".\Release\painter.mod"\
	".\Release\printmod.mod"\
	".\Release\txtParser.mod"\
	
# End Source File
# End Group
# Begin Group "input"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\defineArguments.f90
DEP_F90_DEFIN=\
	".\Release\argsParser.mod"\
	".\Release\hdb.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\definebd.f90
DEP_F90_DEFINE=\
	".\Release\datablock.mod"\
	".\Release\hdb.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\parseInput.f90
DEP_F90_PARSE=\
	".\Release\argsParser.mod"\
	".\Release\datablock.mod"\
	".\Release\hdb.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\planarityCheck.f90
DEP_F90_PLANA=\
	".\Release\hdb.mod"\
	".\Release\orientation.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\readMoleculeInformation.f90
DEP_F90_READM=\
	".\Release\hdb.mod"\
	
# End Source File
# End Group
# Begin Group "calculation"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\gIterator.f90
DEP_F90_GITER=\
	".\Release\glob.mod"\
	".\Release\hdb.mod"\
	
# End Source File
# End Group
# Begin Group "cue"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\cueOrbitals.f90
DEP_F90_CUEOR=\
	".\Release\hdb.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\cueStructureAnalysis.f90
DEP_F90_CUEST=\
	".\Release\hdb.mod"\
	
# End Source File
# End Group
# Begin Group "coupledCluster"

# PROP Default_Filter ""
# Begin Group "ccs"

# PROP Default_Filter ""
# Begin Source File

SOURCE=".\proj-ccs-spat-cue-1.f90"
DEP_F90_PROJ_=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccs-spin-cue-1.f90"
DEP_F90_PROJ_C=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# End Group
# Begin Group "ccd"

# PROP Default_Filter ""
# Begin Source File

SOURCE=".\proj-ccd-spat-hf-2.f90"
DEP_F90_PROJ_CC=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccd-spin-hf-2.f90"
DEP_F90_PROJ_CCD=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# End Group
# Begin Group "ccsd"

# PROP Default_Filter ""
# Begin Source File

SOURCE=".\proj-ccsd-spat-cue-1.f90"
DEP_F90_PROJ_CCS=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spat-cue-2.f90"
DEP_F90_PROJ_CCSD=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spat-hf-1.f90"
DEP_F90_PROJ_CCSD_=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spat-hf-2.f90"
DEP_F90_PROJ_CCSD_S=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spin-cue-1-spare.f90"
DEP_F90_PROJ_CCSD_SP=\
	".\Release\coupledCluster.mod"\
	".\Release\coupledClusterSparse.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spin-cue-1.f90"
DEP_F90_PROJ_CCSD_SPI=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spin-cue-2-spare.f90"
DEP_F90_PROJ_CCSD_SPIN=\
	".\Release\coupledCluster.mod"\
	".\Release\coupledClusterSparse.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spin-cue-2.f90"
DEP_F90_PROJ_CCSD_SPIN_=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spin-hf-1.f90"
DEP_F90_PROJ_CCSD_SPIN_H=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsd-spin-hf-2.f90"
DEP_F90_PROJ_CCSD_SPIN_HF=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# End Group
# Begin Group "ccsdt"

# PROP Default_Filter ""
# Begin Source File

SOURCE=".\proj-ccsdt-spin-cue-1.f90"
DEP_F90_PROJ_CCSDT=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsdt-spin-cue-2.f90"
DEP_F90_PROJ_CCSDT_=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsdt-spin-cue-3.f90"
DEP_F90_PROJ_CCSDT_S=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsdt-spin-hf-1.f90"
DEP_F90_PROJ_CCSDT_SP=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsdt-spin-hf-2.f90"
DEP_F90_PROJ_CCSDT_SPI=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\proj-ccsdt-spin-hf-3.f90"
DEP_F90_PROJ_CCSDT_SPIN=\
	".\Release\coupledCluster.mod"\
	".\Release\glob.mod"\
	
# End Source File
# End Group
# Begin Source File

SOURCE=".\prepare-ccsd-sparse.f90"
DEP_F90_PREPA=\
	".\Release\coupledCluster.mod"\
	".\Release\coupledClusterSparse.mod"\
	
# End Source File
# End Group
# End Target
# End Project
