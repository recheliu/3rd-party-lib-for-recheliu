/*

This files define static variables or other global variables 
declated in the headers.
*/

#include "../include/arerror.h"
#include "../include/arpackf.h"

extern "C" 
{
	CDebug F77NAME(debug);
}

ArpackError::ErrorCode ArpackError::code;