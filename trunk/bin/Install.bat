@ECHO OFF

REM
REM	This script is used to copy the .dll in this folder to the system folder.
REM

REM	ADD-BY-LEETEN 11/09/2011-BEGIN
REM	Change the path to the folder of this script.
CD %~dp0	
ECHO Current Directory is %CD%

REM	Copy the .dlls to the system folder.
REM	ADD-BY-LEETEN 11/09/2011-END

ECHO Copy the *.dlls to %windir%\system.
COPY *.dll %windir%\system\

PAUSE

REM
REM	$Log: not supported by cvs2svn $
REM	Revision 1.2  2009-01-27 05:08:38  leeten
REM
REM	[2009/01/26]
REM	1. [CHANGE] Pause  the script at the end to verify the result.
REM
REM	Revision 1.1  2008/08/22 13:48:04  leeten
REM	
REM	[2008/08/22]
REM	1. [FIRST TIME CHECKIN] This script is used to copy the dll in this folder to the system folder.
REM	
REM	