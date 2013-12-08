cmake_minimum_required (VERSION 2.8)
project (lib3ds)

# The required packages
find_package(OpenGL)
find_package(GLUT)

# The required variables
set(IS_LIB3DS_STATIC OFF CACHE BOOL "Flag to decide whether LIB3DS is static.")

set(libfilename		 			lib3ds-120)

if( IS_LIB3DS_STATIC )
	set(libfilename		 	${libfilename}s)
else()
	add_definitions(
		-D_LIB3DS_DLL
	)
endif()

file(GLOB HEADER_FILES "lib3ds/*.h" )
file(GLOB SOURCE_FILES "lib3ds/*.c" )

include_directories(".")

if( IS_LIB3DS_STATIC )
	add_library ( ${PROJECT_NAME} STATIC
		${HEADER_FILES}
		${SOURCE_FILES}
		)
else()
	add_library ( ${PROJECT_NAME} SHARED
		${HEADER_FILES}
		${SOURCE_FILES}
		)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES 
  OUTPUT_NAME "${libfilename}" 
	DEBUG_OUTPUT_NAME 	"${libfilename}d"
	RELEASE_OUTPUT_NAME "${libfilename}"
)

install(
	TARGETS ${PROJECT_NAME}
	RUNTIME	DESTINATION bin
	ARCHIVE	DESTINATION	lib	)	

install(
	FILES ${HEADER_FILES}
	DESTINATION "include/lib3ds/")
	
#
# $Log$
#
