cmake_minimum_required (VERSION 2.8)
project (lib)

set(IS_LIB_STATIC ON CACHE BOOL "Flag to decide whether static libraries are built.")
set(WITH_NRRD OFF CACHE BOOL "Flag to decide whether NrrdIO is built. NrrdIO so far can only work for 32-bit architecture.")
set(WITH_OPENGL OFF CACHE BOOL "Flag to decide whether OpenGL-related librarties are built. ") # ADD-BY-LEETEN 04/21/2012
if( IS_LIB_STATIC )
	set(IS_STATIC			ON	CACHE	BOOL	"Build static FREEGLUT." 	FORCE)	# freeglut
	set(BUILD_SHARED_LIBS	OFF	CACHE	BOOL	"Build static zlib." 		FORCE)	# zlib
	set(GLUIDLL				OFF	CACHE	BOOL	"Build static glui." 		FORCE)	# glui
	set(IS_LIB3DS_STATIC	ON	CACHE	BOOL	"Build static lib3ds." 		FORCE)	# lib3ds
endif()

add_subdirectory(zlib-1.2.6)			# zlib	
if( WITH_OPENGL )	# ADD-BY-LEETEN 04/21/2012-BEGIN
	add_subdirectory(freeglut-2.8.0) 		# freeglut
	add_subdirectory(glui-2.36/src)			# glui
# ADD-BY-LEETEN 04/21/2012-BEGIN
	add_dependencies(glui freeglut)
endif()	
# ADD-BY-LEETEN 04/21/2012-END

# add_subdirectory(NrrdIO)				# ITKNrrdIO
# ADD-BY-LEETEN 04/20/2012-BEGIN
if( WITH_NRRD )
	add_subdirectory(NrrdIO)			
endif()
# ADD-BY-LEETEN 04/20/2012-END

add_subdirectory(lib3ds-1.2.0)			# lib3ds
# ADD-BY-LEETEN 04/20/2012-BEGIN
# MOD-BY-LEETEN 04/21/2012-FROM:	if( DEFINED(WIN32) )
# DEL-BY-LEETEN 07/09/2012-BEGIN
# if( DEFINED WIN32 )
# # MOD-BY-LEETEN 04/21/2012-END
	# add_subdirectory(pthreads-w32-2-8-0-release)	
# endif()
# DEL-BY-LEETEN 07/09/2012-END
# ADD-BY-LEETEN 04/20/2012-END

# ADD-BY-LEETEN 04/21/2012-BEGIN
add_subdirectory(isotable)
# ADD-BY-LEETEN 04/21/2012-END

# DE:-BY-LEETEN 04/21/2012:	add_dependencies(glui freeglut)
# add_dependencies(ITKNrrdIO zlib)	

# GLEW

# LIB3DS

#
# $Log$
#
