cmake_minimum_required (VERSION 2.8)
project (pthread)

set(libfilename	${PROJECT_NAME})

file(GLOB HEADER_FILES "*.h" )
file(GLOB SOURCE_FILES "*.c" )

if( IS_LIB_STATIC )
	add_definitions( -DPTW32_STATIC_LIB )
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

include_directories(
	"."
)

set_target_properties(${PROJECT_NAME} PROPERTIES 
	# MOD-BY-LEETEN 04/21/2012-FROM:	OUTPUT_NAME "${libfilename}"
	DEBUG_OUTPUT_NAME 	"${libfilename}d"
	RELEASE_OUTPUT_NAME "${libfilename}"
	# MOD-BY-LEETEN 04/21/2012-END
)

install(
	TARGETS ${PROJECT_NAME}
	RUNTIME	DESTINATION bin
	ARCHIVE	DESTINATION	lib
	)	

install(
	FILES "pthread.h"
	DESTINATION "include"
	)
	
#
# $Log$
#
