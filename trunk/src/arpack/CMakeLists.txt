CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(arpack)

file(GLOB HEADER_FILES 	"include/*.h" )
file(GLOB LIBRARY_FILES "lib/*" )
file(GLOB BINARY_FILES 	"bin/*" )
	
add_custom_target (${PROJECT_NAME}
	SOURCES	${HEADER_FILES} 
)

install(
	FILES ${HEADER_FILES}
	DESTINATION "include/arpack")

install(
	FILES ${LIBRARY_FILES}
	DESTINATION "lib")

install(
	FILES ${BINARY_FILES}
	DESTINATION "bin")

