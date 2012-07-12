PROJECT(player)

find_package(OpenGL)

add_definitions(
	-DFREEGLUT_STATIC
	)

ADD_EXECUTABLE(${PROJECT_NAME}
	player.c
	)

if( DEFINED WIN32 )
	target_link_libraries(${PROJECT_NAME}
		debug	lib3ds-120sd	optimized	lib3ds-120s
		freeglut_static.lib
		)
else()
	target_link_libraries(${PROJECT_NAME}
		lib3ds-120s
		)
endif()

install(
	TARGETS ${PROJECT_NAME}
	RUNTIME	DESTINATION "bin"
	ARCHIVE	DESTINATION	"lib"	)	

