#	Required cache entries:
# 		FREEGLUT
#		THIRDPARTY_DIR
#		OPENCV_DIR
#		OPENCV_VER

cmake_minimum_required (VERSION 2.8)
project (glui)

# The required packages
find_package(OpenGL)
find_package(GLUT)

# The required variables
set(FREEGLUT OFF CACHE BOOL "Flag to decide whether FREEGLUT is used.")
set(GLUT_DIR "~/lib" CACHE PATH "Path of GLUT.")
set(GLUIDLL OFF CACHE BOOL "Flag to decide whether the libraries are compiled to .dll.")

set(libfilename		 			glui32)
set(libfilename_debug		 	glui32d)

if( GLUIDLL )
	set(libfilename_postfix		 	${libfilename_postfix}dll)
	add_definitions(
		-DGLUIDLL
	)
endif()

if( FREEGLUT )
	set(libfilename_postfix		 	${libfilename_postfix}_freeglut)
	add_definitions(
		-DGLUI_FREEGLUT
	)
	if( NOT GLUIDLL )
		add_definitions(
			-DFREEGLUT_STATIC
		)
	endif()
endif()

add_definitions(
	-DGLUI_BUILDING_LIB
)

include_directories(
	${GLUT_DIR}/include
	"include"
)

link_directories(
	${GLUT_DIR}/lib
)

file(GLOB HEADER_FILES "*.h" )
file(GLOB GL_HEADER_FILES "include/GL/*.h" )
file(GLOB SOURCE_FILES "*.cpp" )

if( GLUIDLL )
	add_library ( ${PROJECT_NAME} SHARED
		${HEADER_FILES}
		${GL_HEADER_FILES}
		${SOURCE_FILES}
		)
else()
	add_library ( ${PROJECT_NAME} STATIC
		${HEADER_FILES}
		${GL_HEADER_FILES}
		${SOURCE_FILES}
		)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME "${libfilename_debug}${libfilename_postfix}"
	RELEASE_OUTPUT_NAME "${libfilename}${libfilename_postfix}"
)

install(
	TARGETS ${PROJECT_NAME}
	RUNTIME	DESTINATION bin
	ARCHIVE	DESTINATION	lib	)	

install(
	FILES ${GL_HEADER_FILES}
	DESTINATION "include/GL/")
	
#
# $Log$
#
