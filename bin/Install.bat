@ECHO OFF

REM
REM	This script is used to copy the .dll in this folder to the system folder.
REM

ECHO Copy the *.dlls to %windir%\system.

COPY *.dll %windir%\system\

PAUSE

REM
REM	$Log: not supported by cvs2svn $
REM	Revision 1.1  2008/08/22 13:48:04  leeten
REM	
REM	[2008/08/22]
REM	1. [FIRST TIME CHECKIN] This script is used to copy the dll in this folder to the system folder.
REM	
REM	