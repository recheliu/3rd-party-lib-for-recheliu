#	Required cache entries:
# 		FREEGLUT
#		THIRDPARTY_DIR
#		OPENCV_DIR
#		OPENCV_VER

cmake_minimum_required (VERSION 2.8)
project (freeglut)

# The required packages
find_package(OpenGL)

# The required variables
set(IS_STATIC OFF CACHE BOOL "Flag to decide whether the libraries are compiled to .dll.")

set(libfilename		 			freeglut)

if( IS_STATIC )
	set(libfilename		 	${libfilename}_static)

	add_definitions(
		-DFREEGLUT_STATIC
	)
endif()

include_directories(
	"include"
)

file(GLOB HEADER_FILES "src/*.h" )
file(GLOB GL_HEADER_FILES "include/GL/*.h" )
file(GLOB SOURCE_FILES "src/*.c" )

if( IS_STATIC )
	add_library ( ${PROJECT_NAME} STATIC
		${HEADER_FILES}
		${GL_HEADER_FILES}
		${SOURCE_FILES}
		)
else()
	add_library ( ${PROJECT_NAME} SHARED
		${HEADER_FILES}
		${GL_HEADER_FILES}
		${SOURCE_FILES}
		)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES 
	OUTPUT_NAME "${libfilename}"
)

install(
	TARGETS ${PROJECT_NAME}
	ARCHIVE	DESTINATION	lib	)	

install(
	FILES ${GL_HEADER_FILES}
	DESTINATION "include/GL/")
	
#
# $Log$
#
