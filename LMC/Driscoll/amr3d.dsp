# Microsoft Developer Studio Project File - Name="amr3d" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=amr3d - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "amr3d.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "amr3d.mak" CFG="amr3d - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "amr3d - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "amr3d - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /assume:noaccuracy_sensitive /compile_only /debug:none /iface:nomixed_str_len_arg /iface:cref /include:"Release/" /math_library:fast /nologo /threads /tune:k7 /warn:nofileopt /unroll:4
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GR /GX /O2 /I "." /I ".." /I "..\..\iamrlib\InflowForce" /I "..\..\LMC\\" /I "..\..\iamrlib\\" /I "..\..\tensorMG\\" /I "..\..\hgproj\\" /I "..\..\mglib\\" /I "..\..\amrlib\\" /I "..\..\bndrylib\\" /I "..\..\BoxLib\\" /I "C:\WMPI\include" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D BL_PRVERSION=5 /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_USE_CHEM" /D "BL_DO_FLCT" /D for="if(0);else for" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib dformt.lib /nologo /subsystem:console /machine:I386 /include:"__matherr"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /dbglibs /debug:full /iface:nomixed_str_len_arg /iface:cref /include:"Debug/" /nologo /threads /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MTd /W3 /Gm /GR /GX /ZI /Od /I "." /I ".." /I "..\..\iamrlib\InflowForce" /I "..\..\LMC\\" /I "..\..\iamrlib\\" /I "..\..\tensorMG\\" /I "..\..\hgproj\\" /I "..\..\mglib\\" /I "..\..\amrlib\\" /I "..\..\bndrylib\\" /I "..\..\BoxLib\\" /I "C:\WMPI\include" /D "_CONSOLE" /D "_MBCS" /D "_DEBUG" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D BL_PRVERSION=5 /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_USE_CHEM" /D "BL_DO_FLCT" /D for="if(0);else for" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib dformt.lib /nologo /subsystem:console /debug /machine:I386 /include:"__matherr" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "amr3d - Win32 Release"
# Name "amr3d - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\..\mglib\ABecLaplacian.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Amr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\amr_multi.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\AmrLevel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Arena.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BArena.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BaseFab.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\BCRec.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLProfiler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLThread.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLWorkQueue.cpp
# End Source File
# Begin Source File

SOURCE=..\..\mglib\BndryData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\BndryRegister.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\boundary.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Box.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxArray.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxDomain.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxLib.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxList.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\CArena.cpp
# End Source File
# Begin Source File

SOURCE=..\..\mglib\CGSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\LMC\ChemDriver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Cluster.cpp
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\CoordSys.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Derive.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\Diffusion.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\DistributionMapping.cpp
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\DivVis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\LMC\drm19.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\ErrorList.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\FabArray.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\FabConv.cpp
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\FabSet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\FArrayBox.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\fill_patch.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\FluxRegister.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\FPC.cpp
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\Geometry.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\Godunov.cpp
# End Source File
# Begin Source File

SOURCE=..\..\LMC\HeatTransfer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_multi1.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_multi2.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_multi3.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_projector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hgparallel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\LMC\HT_setup.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\IndexType.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\interface.cpp
# End Source File
# Begin Source File

SOURCE=..\..\mglib\InterpBndryData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Interpolater.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\interpolator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\IntVect.cpp
# End Source File
# Begin Source File

SOURCE=..\..\mglib\Laplacian.cpp
# End Source File
# Begin Source File

SOURCE=..\..\mglib\LinOp.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacBndry.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacOperator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacOutFlowBC.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacProj.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\main.cpp
# End Source File
# Begin Source File

SOURCE=..\..\mglib\Mask.cpp
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCCGSolver.cpp
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCInterpBndryData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCLinOp.cpp
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCMultiGrid.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\MultiFab.cpp
# End Source File
# Begin Source File

SOURCE=..\..\mglib\MultiGrid.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\NavierStokes.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Orientation.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\OutFlowBC.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\ParallelDescriptor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\ParmParse.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\Projection.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\ProjOutFlowBC.cpp
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\RealBox.cpp
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\restrictor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\SlabStat.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\StateData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\StateDescriptor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\StationData.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\SyncRegister.cpp
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\TagBox.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\UseCount.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Utility.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\ViscBndry.cpp
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\ViscBndryTensor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\VisMF.cpp
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "H"
# Begin Source File

SOURCE=..\..\mglib\ABec_F.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\ABecLaplacian.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Amr.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\amr_defs.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\amr_multi.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\AmrLevel.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Arena.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\ArrayLim.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\ArrayViewEXT.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BArena.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\BC_TYPES.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\BCRec.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLassert.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLFort.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLVERSION.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\BndryData.H
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\BndryRegister.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\boundary.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\BoundCond.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Box.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxArray.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxDomain.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxLib.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BoxList.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\CArena.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\cdwrk.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\CG_F.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\CGSolver.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\ChemDriver.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\ChemDriver_F.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Cluster.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\CONSTANTS.H
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\CoordSys.H
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\COORDSYS_F.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Derive.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\DERIVE_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\Diffusion.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\DIFFUSION_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\DistributionMapping.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\DivVis.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\DivVis_F.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\ErrorList.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\FabConv.H
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\FabSet.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\FArrayBox.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\fill_patch.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\FLUXREG_F.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\FluxRegister.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\FPC.H
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\Geometry.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\GODCOMM_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\Godunov.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\GODUNOV_F.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\HeatTransfer.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\HEATTRANSFER_F.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_multi.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_projector.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hgparallel.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\htdata.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\IndexType.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\interface.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\INTERP_F.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\InterpBndryData.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\INTERPBNDRYDATA_F.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\Interpolater.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\interpolator.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\IntVect.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\Laplacian.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\LevelBld.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\LinOp.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\LO_BCTYPES.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\LO_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Looping.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\LP_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacBndry.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacOperator.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MACOPERATOR_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacOpMacDrivers.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacOpProjDrivers.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacOutFlowBC.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MACOUTFLOWBC_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MacProj.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MACPROJ_F.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\MAKESLICE_F.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\Mask.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCCGSolver.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCInterpBndryData.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCINTERPBNDRYDATA_F.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCLinOp.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCLO_F.H
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCMultiGrid.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\MG_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\MultiFab.H
# End Source File
# Begin Source File

SOURCE=..\..\mglib\MultiGrid.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\NavierStokes.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\NAVIERSTOKES_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Orientation.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\OutFlowBC.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\ParallelDescriptor.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\ParmParse.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\PROB_AMR_F.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\PROB_F.H
# End Source File
# Begin Source File

SOURCE=probdata.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Profiler.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\Projection.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\PROJECTION_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\ProjOutFlowBC.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\PROJOUTFLOWBC_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\REAL.H
# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\RealBox.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\RegType.H
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\restrictor.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\SlabStat.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\SLABSTAT_F.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\SLABSTAT_HT_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\SLABSTAT_NS_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\SPACE.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\SPACE_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\SPECIALIZE_F.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\StateData.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\StateDescriptor.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\StationData.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\SYNCREG_F.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\SyncRegister.H
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\TagBox.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Thread.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\UseCount.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\Utility.H
# End Source File
# Begin Source File

SOURCE=..\..\LMC\visc.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\ViscBndry.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\ViscBndryTensor.H
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\VISCOPERATOR_F.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\VisMF.H
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\WorkQueue.H
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter ""
# Begin Group "TEMPS-Debug"

# PROP Default_Filter "FOR"
# Begin Source File

SOURCE=Debug\Fort\ABec_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\ABec_UTIL.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\ARRAYLIM_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\CG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\ChemDriver_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\ChemDriver_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\COORDSYS_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\DERIVE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\DERIVE_HT_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\DIFFUSION_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\DV_3D1.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\DV_3D2.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\DV_3D3.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\DV_3D4.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\FILCC_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\FLUCTFILE.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\FLUXREG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\GODUNOV_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\GODUNOV_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\HEATTRANSFER_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\HEATTRANSFER_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\INFL_FORCE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\INFL_FORCE_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\INTERP_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\INTERPBNDRYDATA_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\LO_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\LO_UTIL.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\LP_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\MACOPERATOR_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\MACOUTFLOWBC_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\MACPROJ_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\MAKESLICE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\MCINTERPBNDRYDATA_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\MCLO_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\MG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\NAVIERSTOKES_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\PROB_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\PROJECTION_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\PROJOUTFLOWBC_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\SLABSTAT_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\SLABSTAT_HT_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\SLABSTAT_NS_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\SPECIALIZE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\SYNCREG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Debug\Fort\VISCOPERATOR_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

!ENDIF 

# End Source File
# End Group
# Begin Group "TEMPS-Release"

# PROP Default_Filter "FOR"
# Begin Source File

SOURCE=Release\Fort\ABec_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\ABec_UTIL.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\ARRAYLIM_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\CG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\ChemDriver_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\ChemDriver_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\COORDSYS_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\DERIVE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\DERIVE_HT_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\DIFFUSION_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\DV_3D1.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\DV_3D2.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\DV_3D3.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\DV_3D4.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\FILCC_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\FLUCTFILE.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\FLUXREG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\GODUNOV_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\GODUNOV_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\HEATTRANSFER_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\HEATTRANSFER_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\INFL_FORCE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\INFL_FORCE_F.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\INTERP_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\INTERPBNDRYDATA_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\LO_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\LO_UTIL.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\LP_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\MACOPERATOR_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\MACOUTFLOWBC_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\MACPROJ_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\MAKESLICE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\MCINTERPBNDRYDATA_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\MCLO_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\MG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\NAVIERSTOKES_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\PROB_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\PROJECTION_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\PROJOUTFLOWBC_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\SLABSTAT_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\SLABSTAT_HT_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\SLABSTAT_NS_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\SPECIALIZE_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\SYNCREG_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=Release\Fort\VISCOPERATOR_3D.FOR

!IF  "$(CFG)" == "amr3d - Win32 Release"

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# End Group
# Begin Source File

SOURCE=..\..\mglib\ABec_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\ABec_3D.F
InputName=ABec_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\ABec_3D.F
InputName=ABec_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\mglib\ABec_UTIL.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\ABec_UTIL.F
InputName=ABec_UTIL

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\ABec_UTIL.F
InputName=ABec_UTIL

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\hgproj\amr_real3d.f
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\ARRAYLIM_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\ARRAYLIM_3D.F
InputName=ARRAYLIM_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\ARRAYLIM_3D.F
InputName=ARRAYLIM_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLBoxLib_F.f
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLParmParse_F.f
# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\BLutil_F.f
# End Source File
# Begin Source File

SOURCE=..\..\mglib\CG_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\CG_3D.F
InputName=CG_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\CG_3D.F
InputName=CG_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\ChemDriver_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\ChemDriver_3D.F
InputName=ChemDriver_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\ChemDriver_3D.F
InputName=ChemDriver_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\ChemDriver_F.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\ChemDriver_F.F
InputName=ChemDriver_F

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\ChemDriver_F.F
InputName=ChemDriver_F

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\bndrylib\COORDSYS_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\bndrylib\COORDSYS_3D.F
InputName=COORDSYS_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\bndrylib\COORDSYS_3D.F
InputName=COORDSYS_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\DERIVE_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\DERIVE_3D.F
InputName=DERIVE_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=.\DERIVE_3D.F
InputName=DERIVE_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\DERIVE_HT_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\DERIVE_HT_3D.F
InputName=DERIVE_HT_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\DERIVE_HT_3D.F
InputName=DERIVE_HT_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\DIFFUSION_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\DIFFUSION_3D.F
InputName=DIFFUSION_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\DIFFUSION_3D.F
InputName=DIFFUSION_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\drm19_100.f
# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\DV_3D1.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D1.F
InputName=DV_3D1

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D1.F
InputName=DV_3D1

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\DV_3D2.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D2.F
InputName=DV_3D2

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D2.F
InputName=DV_3D2

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\DV_3D3.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D3.F
InputName=DV_3D3

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D3.F
InputName=DV_3D3

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\DV_3D4.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D4.F
InputName=DV_3D4

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\DV_3D4.F
InputName=DV_3D4

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\EGini.f
DEP_F90_EGINI=\
	"..\eg.cmn"\
	
# End Source File
# Begin Source File

SOURCE=..\..\LMC\EGSlib.f
DEP_F90_EGSLI=\
	"..\eg.cmn"\
	
# End Source File
# Begin Source File

SOURCE=..\..\amrlib\FILCC_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\FILCC_3D.F
InputName=FILCC_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\FILCC_3D.F
InputName=FILCC_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\InflowForce\FLUCTFILE.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\InflowForce\FLUCTFILE.F
InputName=FLUCTFILE

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\InflowForce\FLUCTFILE.F
InputName=FLUCTFILE

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\amrlib\FLUXREG_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\FLUXREG_3D.F
InputName=FLUXREG_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\FLUXREG_3D.F
InputName=FLUXREG_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=GODUNOV_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=GODUNOV_3D.F
InputName=GODUNOV_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=GODUNOV_3D.F
InputName=GODUNOV_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\GODUNOV_F.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\GODUNOV_F.F
InputName=GODUNOV_F

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\GODUNOV_F.F
InputName=GODUNOV_F

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\HEATTRANSFER_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\HEATTRANSFER_3D.F
InputName=HEATTRANSFER_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\HEATTRANSFER_3D.F
InputName=HEATTRANSFER_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\HEATTRANSFER_F.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\HEATTRANSFER_F.F
InputName=HEATTRANSFER_F

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\HEATTRANSFER_F.F
InputName=HEATTRANSFER_F

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_avg3d.f
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_multi3d.f
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_multi3d_full.f
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_multi3d_terrain.f
# End Source File
# Begin Source File

SOURCE=..\..\hgproj\hg_proj3d.f
# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\InflowForce\INFL_FORCE_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\InflowForce\INFL_FORCE_3D.F
InputName=INFL_FORCE_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\InflowForce\INFL_FORCE_3D.F
InputName=INFL_FORCE_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\InflowForce\INFL_FORCE_F.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\InflowForce\INFL_FORCE_F.F
InputName=INFL_FORCE_F

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\InflowForce\INFL_FORCE_F.F
InputName=INFL_FORCE_F

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\amrlib\INTERP_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\INTERP_3D.F
InputName=INTERP_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\INTERP_3D.F
InputName=INTERP_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\mglib\INTERPBNDRYDATA_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\INTERPBNDRYDATA_3D.F
InputName=INTERPBNDRYDATA_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\INTERPBNDRYDATA_3D.F
InputName=INTERPBNDRYDATA_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\mglib\LO_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\LO_3D.F
InputName=LO_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\LO_3D.F
InputName=LO_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\mglib\LO_UTIL.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\LO_UTIL.F
InputName=LO_UTIL

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\LO_UTIL.F
InputName=LO_UTIL

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\mglib\LP_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\LP_3D.F
InputName=LP_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\LP_3D.F
InputName=LP_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MACOPERATOR_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\MACOPERATOR_3D.F
InputName=MACOPERATOR_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\MACOPERATOR_3D.F
InputName=MACOPERATOR_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MACOUTFLOWBC_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\MACOUTFLOWBC_3D.F
InputName=MACOUTFLOWBC_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\MACOUTFLOWBC_3D.F
InputName=MACOUTFLOWBC_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\MACPROJ_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\MACPROJ_3D.F
InputName=MACPROJ_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\MACPROJ_3D.F
InputName=MACPROJ_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\amrlib\MAKESLICE_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\MAKESLICE_3D.F
InputName=MAKESLICE_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\MAKESLICE_3D.F
InputName=MAKESLICE_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCINTERPBNDRYDATA_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\MCINTERPBNDRYDATA_3D.F
InputName=MCINTERPBNDRYDATA_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\MCINTERPBNDRYDATA_3D.F
InputName=MCINTERPBNDRYDATA_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\tensorMG\MCLO_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\MCLO_3D.F
InputName=MCLO_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\tensorMG\MCLO_3D.F
InputName=MCLO_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\mglib\MG_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\MG_3D.F
InputName=MG_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\mglib\MG_3D.F
InputName=MG_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\NAVIERSTOKES_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\NAVIERSTOKES_3D.F
InputName=NAVIERSTOKES_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\NAVIERSTOKES_3D.F
InputName=NAVIERSTOKES_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=PROB_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=PROB_3D.F
InputName=PROB_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=PROB_3D.F
InputName=PROB_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\PROJECTION_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\PROJECTION_3D.F
InputName=PROJECTION_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\PROJECTION_3D.F
InputName=PROJECTION_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\PROJOUTFLOWBC_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\PROJOUTFLOWBC_3D.F
InputName=PROJOUTFLOWBC_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\PROJOUTFLOWBC_3D.F
InputName=PROJOUTFLOWBC_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\amrlib\SLABSTAT_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\SLABSTAT_3D.F
InputName=SLABSTAT_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\amrlib\SLABSTAT_3D.F
InputName=SLABSTAT_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\SLABSTAT_HT_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\SLABSTAT_HT_3D.F
InputName=SLABSTAT_HT_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\LMC\SLABSTAT_HT_3D.F
InputName=SLABSTAT_HT_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\SLABSTAT_NS_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\SLABSTAT_NS_3D.F
InputName=SLABSTAT_NS_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\SLABSTAT_NS_3D.F
InputName=SLABSTAT_NS_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\BoxLib\SPECIALIZE_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\BoxLib\SPECIALIZE_3D.F
InputName=SPECIALIZE_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\BoxLib\SPECIALIZE_3D.F
InputName=SPECIALIZE_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\SYNCREG_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\SYNCREG_3D.F
InputName=SYNCREG_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\SYNCREG_3D.F
InputName=SYNCREG_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\iamrlib\VISCOPERATOR_3D.F

!IF  "$(CFG)" == "amr3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\VISCOPERATOR_3D.F
InputName=VISCOPERATOR_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "amr3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\..\iamrlib\VISCOPERATOR_3D.F
InputName=VISCOPERATOR_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S.. /S..\..\iamrlib\InflowForce /S..\..\LMC\ /S..\..\iamrlib\ /S..\..\tensorMG\ /S..\..\hgproj\ /S..\..\mglib\ /S..\..\amrlib\ /S..\..\bndrylib\ /S..\..\BoxLib\ /DBL_LANG_FORT  /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM /DBL_DO_FLCT  /DBL_PRVERSION=5 $(InputPath) | perl ..\..\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\LMC\vode.f
DEP_F90_VODE_=\
	"..\conp.H"\
	
# End Source File
# End Group
# End Target
# End Project
