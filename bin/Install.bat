REM
REM	This script is used to copy the .dll in this folder to the system folder.
REM

@ECHO OFF

ECHO Copy the *.dlls to %windir%\system.

COPY *.dll %windir%\system\

REM
REM	$Log: not supported by cvs2svn $
REM	