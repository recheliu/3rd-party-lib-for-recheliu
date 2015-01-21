cmake_minimum_required (VERSION 2.8)
project(ThirdPartyLib)

set(IS_LIB_STATIC ON CACHE BOOL "Flag to decide whether static libraries are built.")
set(WITH_NRRD ON CACHE BOOL "Flag to decide whether NrrdIO is built. NrrdIO so far can only work for 32-bit architecture." FORCE)
set(WITH_OPENGL OFF CACHE BOOL "Flag to decide whether OpenGL-related librarties are built. ") # ADD-BY-LEETEN 04/21/2012
if( IS_LIB_STATIC )
	set(IS_STATIC			ON	CACHE	BOOL	"Build static FREEGLUT." 	FORCE)	# freeglut
	set(BUILD_SHARED_LIBS	OFF	CACHE	BOOL	"Build static zlib." 		FORCE)	# zlib
	set(GLUIDLL				OFF	CACHE	BOOL	"Build static glui." 		FORCE)	# glui
	set(IS_LIB3DS_STATIC	ON	CACHE	BOOL	"Build static lib3ds." 		FORCE)	# lib3ds
endif()

add_subdirectory(lodepng)				
# ADD-BY-LEETEN 2015/01/20-BEGIN
add_subdirectory(arpack)				
# ADD-BY-LEETEN 2015/01/20-END

add_subdirectory(zlib-1.2.6)			# zlib	
include_directories(zlib-1.2.6)
include_directories(${CMAKE_BINARY_DIR}/zlib-1.2.6)
if( WITH_OPENGL )	
	add_subdirectory(freeglut-2.8.0) 		# freeglut
	add_subdirectory(glui-2.36/src)			# glui
	add_dependencies(glui freeglut)
	set(GLUT_DIR ${PROJECT_SOURCE_DIR}/freeglut-2.8.0 CACHE PATH "Path of GLUT." FORCE)
endif()	

# add_subdirectory(NrrdIO)				# ITKNrrdIO
if( WITH_NRRD )
	add_subdirectory(NrrdIO)			
endif()

add_subdirectory(lib3ds-1.2.0)			# lib3ds
if( DEFINED WIN32 ) 
	include_directories(lib3ds-1.2.0) 
endif() # ADD-BY-LEETEN 10/29/2012

add_subdirectory(isotable)

set(IJK_DIR					${PROJECT_BINARY_DIR}	CACHE	PATH	"" 	FORCE)	
set(LIBRARY_OUTPUT_PATH		${PROJECT_BINARY_DIR}/lib	CACHE	PATH	"" 	FORCE)	

# GLM
# The CMake file in GLM is desigend to build the test programs.
# Since only the headers are needed, the install rule is manually added here 
install(
	DIRECTORY	 glm-0.9.4.1/glm/
	DESTINATION "include/glm")

# add_dependencies(ITKNrrdIO zlib)	

# GLEW

# LIB3DS

configure_file (
  "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  @ONLY
  )

install(
	FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
	DESTINATION ".")

#
# $Log$
#
