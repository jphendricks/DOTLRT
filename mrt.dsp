# Microsoft Developer Studio Project File - Name="mrt" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=mrt - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "mrt.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "mrt.mak" CFG="mrt - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "mrt - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "mrt - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "mrt - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /stack:0x5f5e100,0x5f5e100 /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "mrt - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /stack:0x5f5e100,0x5f5e100 /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "mrt - Win32 Release"
# Name "mrt - Win32 Debug"
# Begin Source File

SOURCE=.\absh2o.f90
# End Source File
# Begin Source File

SOURCE=.\absn2.f90
# End Source File
# Begin Source File

SOURCE=.\calc_fbw_temp_weight_scat.f90
# End Source File
# Begin Source File

SOURCE=.\calc_mon_temp_weight_scat.f90
# End Source File
# Begin Source File

SOURCE=.\calc_passband_freq.f90
# End Source File
# Begin Source File

SOURCE=.\calc_tot_ext.f90
# End Source File
# Begin Source File

SOURCE=.\calcprofile.f90
# End Source File
# Begin Source File

SOURCE=.\configure.f90
# End Source File
# Begin Source File

SOURCE=.\core95.f90
# End Source File
# Begin Source File

SOURCE=.\d3lec.f90
# End Source File
# Begin Source File

SOURCE=.\do_tb_gvh94.f90
# End Source File
# Begin Source File

SOURCE=.\Dtb94.f90
# End Source File
# Begin Source File

SOURCE=.\fresnel_refl.f90
# End Source File
# Begin Source File

SOURCE=.\gaussj.f90
# End Source File
# Begin Source File

SOURCE=.\gaussq.f90
# End Source File
# Begin Source File

SOURCE=.\GeoJacobian.f90
# End Source File
# Begin Source File

SOURCE=.\get_instr_spec.f90
# End Source File
# Begin Source File

SOURCE=.\GetProfile.f90
# End Source File
# Begin Source File

SOURCE=.\Hg_phmat.f90
# End Source File
# Begin Source File

SOURCE=.\hydro_master_derivatives.f90
# End Source File
# Begin Source File

SOURCE=.\InterpSR.f90
# End Source File
# Begin Source File

SOURCE=.\jacobi.f90
# End Source File
# Begin Source File

SOURCE=.\Jacobian.f90
# End Source File
# Begin Source File

SOURCE=.\kirchoff_ocean_refl.f90
# End Source File
# Begin Source File

SOURCE=.\mex_matlab_mrt.f90
# End Source File
# Begin Source File

SOURCE=.\mrt.f90
# End Source File
# Begin Source File

SOURCE=.\o2abs.f90
# End Source File
# Begin Source File

SOURCE=.\ocean_refl.f90
# End Source File
# Begin Source File

SOURCE=.\rd.f90
# End Source File
# Begin Source File

SOURCE=.\rf.f90
# End Source File
# Begin Source File

SOURCE=.\RT_Jacobian.f90
# End Source File
# Begin Source File

SOURCE=.\Tb94.f90
# End Source File
# Begin Source File

SOURCE=.\variables.f90
# End Source File
# End Target
# End Project
