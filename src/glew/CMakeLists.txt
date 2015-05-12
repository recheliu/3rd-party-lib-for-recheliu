cmake_minimum_required (VERSION 2.8)
project(GLEW)

set(GLEW_ARCHIVE_PREFIX     "glew-1.12.0")
set(GLEW_DIR                "${CMAKE_CURRENT_SOURCE_DIR}/${GLEW_ARCHIVE_PREFIX}")
set(PYTHON_EXEC "/usr/bin/python" CACHE FILEPATH "Executable to the python executable.")

if(NOT EXISTS ${PYTHON_EXEC})
    message(FATAL_ERROR "Cannot locate the python executable in ${}.")
endif()

if( DEFINED WIN32 ) 
    set(GLEW_WIN32_ARCHIVE_FILENAME   "${GLEW_ARCHIVE_PREFIX}-win32.zip")
    
    # TODO: Download the pre-compiled .zip?
    # Since the link can be not-trivial,
    # a copy is mirrored to my repo.
    
    # Use the python script to extract the downloaded zip.
    execute_process(
        COMMAND ${PYTHON_EXEC} ${CMAKE_SOURCE_DIR}/python/unzip.py ${GLEW_WIN32_ARCHIVE_FILENAME}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    install(
        DIRECTORY "${GLEW_DIR}/include/GL"
        DESTINATION "include")

    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(PLATFORM_DIRNAME    "x64")
    else()
        set(PLATFORM_DIRNAME    "Win32")
    endif()	

	file(GLOB RUNTIME_FILES "${GLEW_DIR}/bin/Release/${PLATFORM_DIRNAME}/*" )
    install(
        FILES ${RUNTIME_FILES}
        DESTINATION "bin")
        
	file(GLOB ARCHIVE_FILES "${GLEW_DIR}/lib/Release/${PLATFORM_DIRNAME}/*" )
    install(
        FILES ${ARCHIVE_FILES} 
        DESTINATION "lib")
    # Setup the installation rule.
else()
    # TODO: Donwload the archive for source codes.
    message(WARNING "Currently this repo cannot build/install glew on non-Windows platform yet.")
endif() 

