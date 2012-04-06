cmake_minimum_required (VERSION 2.8)
project (lib)

set(IS_LIB_STATIC ON CACHE BOOL "Flag to decide whether static libraries are built.")
if( IS_LIB_STATIC )
	set(IS_STATIC			ON	CACHE	BOOL	"Build static FREEGLUT." 	FORCE)	# freeglut
	set(BUILD_SHARED_LIBS	OFF	CACHE	BOOL	"Build static zlib." 		FORCE)	# zlib
	set(GLUIDLL				OFF	CACHE	BOOL	"Build static glui." 		FORCE)	# glui
	set(IS_LIB3DS_STATIC	ON	CACHE	BOOL	"Build static lib3ds." 		FORCE)	# lib3ds
endif()

add_subdirectory(zlib-1.2.6)			# zlib	
add_subdirectory(freeglut-2.8.0) 		# freeglut
add_subdirectory(glui-2.36/src)			# glui
# add_subdirectory(NrrdIO)				# ITKNrrdIO
add_subdirectory(lib3ds-1.2.0)			# lib3ds

add_dependencies(glui freeglut)
# add_dependencies(ITKNrrdIO zlib)	

# GLEW

# LIB3DS

#
# $Log$
#
