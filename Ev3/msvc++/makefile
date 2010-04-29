# Microsoft Developer Studio Generated NMAKE File, Based on Ev3.dsp
!IF "$(CFG)" == ""
CFG=Ev3 - Win32 Release
!MESSAGE No configuration specified. Defaulting to Ev3 - Win32 Release.
!ENDIF 

!IF "$(CFG)" != "Ev3 - Win32 Release" && "$(CFG)" != "Ev3 - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Ev3.mak" CFG="Ev3 - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Ev3 - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Ev3 - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Ev3 - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\Ev3.exe"


CLEAN :
	-@erase "$(INTDIR)\examplemain.obj"
	-@erase "$(INTDIR)\expression.obj"
	-@erase "$(INTDIR)\parser.obj"
	-@erase "$(OUTDIR)\Ev3.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MD /W3 /GX /O2 /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "SUNWIN" /D "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)\Ev3.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\Ev3.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=/nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\Ev3.pdb" /machine:I386 /out:"$(OUTDIR)\Ev3.exe" 
LINK32_OBJS= \
	"$(INTDIR)\expression.obj" \
	"$(INTDIR)\examplemain.obj" \
	"$(INTDIR)\parser.obj"

"$(OUTDIR)\Ev3.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "Ev3 - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\Ev3.exe"


CLEAN :
	-@erase "$(INTDIR)\examplemain.obj"
	-@erase "$(INTDIR)\expression.obj"
	-@erase "$(INTDIR)\parser.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\Ev3.exe"
	-@erase "$(OUTDIR)\Ev3.ilk"
	-@erase "$(OUTDIR)\Ev3.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MDd /W3 /Gm /GX /ZI /Od /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "SUNWIN" /D "_AFXDLL" /D "_MBCS" /Fp"$(INTDIR)\Ev3.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\Ev3.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=/nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\Ev3.pdb" /debug /machine:I386 /out:"$(OUTDIR)\Ev3.exe" 
LINK32_OBJS= \
	"$(INTDIR)\expression.obj" \
	"$(INTDIR)\examplemain.obj" \
	"$(INTDIR)\parser.obj"

"$(OUTDIR)\Ev3.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("Ev3.dep")
!INCLUDE "Ev3.dep"
!ELSE 
!MESSAGE Warning: cannot find "Ev3.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "Ev3 - Win32 Release" || "$(CFG)" == "Ev3 - Win32 Debug"
SOURCE="\liberti\phd\expression\v3\examplemain.cxx"

"$(INTDIR)\examplemain.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE="\liberti\phd\expression\v3\expression.cxx"

"$(INTDIR)\expression.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE="\liberti\phd\expression\v3\parser.cxx"

"$(INTDIR)\parser.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)

!ENDIF 

