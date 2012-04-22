PROJECT(IJKMCUBE)

#---------------------------------------------------------

IF (NOT DEFINED ${IJK_DIR})
  SET(IJK_DIR "../.." CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMCUBE_DIR ijkmcube)
SET(NRRD_LIBDIR "${IJK_DIR}/lib")

#---------------------------------------------------------


INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
# TEST-DEL:		LINK_LIBRARIES(ITKNrrdIO z)

# TEST-MOD-FROM:	ADD_EXECUTABLE(ijkmcube ijkmcube.cxx isosurfaceD.cxx ijktable.cxx)
LIST(APPEND HEADER_FILES
	isosurfaceD.h
	ijktable.h
	)

ADD_LIBRARY(ijktable 
	isosurfaceD.cxx 
	ijktable.cxx
	${HEADER_FILES}
	)
set_target_properties(ijktable PROPERTIES 
	DEBUG_OUTPUT_NAME 	"ijktabled"
	RELEASE_OUTPUT_NAME "ijktable"
)

ADD_EXECUTABLE(ijkmcube ijkmcube.cxx)
if( DEFINED WIN32 )
	target_link_libraries(ijkmcube
		debug	ITKNrrdIOd	optimized	ITKNrrdIO
		debug	zlibd		optimized	zlib
		)
else()
	target_link_libraries(ijkmcube
		ITKNrrdIO 
		z
		)
endif()
target_link_libraries(ijkmcube
	ijktable
	)
add_dependencies(ijkmcube ijktable)
# TEST-MOD-END

SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
# TEST-MOD-FROM:
	# INSTALL(TARGETS ijkmcube DESTINATION "bin")
# TEST-MOD-TO:
install(
	TARGETS ijktable
	RUNTIME	DESTINATION "bin"
	ARCHIVE	DESTINATION	"lib"	)	

install(
	TARGETS ijkmcube
	RUNTIME	DESTINATION "bin"
	ARCHIVE	DESTINATION	"lib"	)	

install(
	FILES ${HEADER_FILES}
	DESTINATION "include")

# TEST-MOD-END


# TEST-DEL:	ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY .. COMMAND tar cvf ${IJKMCUBE_DIR}/ijkmcube.tar ${IJKMCUBE_DIR}/README ${IJKMCUBE_DIR}/*.cxx ${IJKMCUBE_DIR}/*.h ${IJKMCUBE_DIR}/*.txx ${IJKMCUBE_DIR}/CMakeLists.txt )

