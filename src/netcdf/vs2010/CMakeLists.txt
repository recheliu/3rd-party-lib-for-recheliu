cmake_minimum_required (VERSION 2.8)
project(NETCDF4)

set(PYTHON_EXEC "/usr/bin/python" CACHE FILEPATH "Executable to the python executable.")

if(NOT EXISTS ${PYTHON_EXEC})
    message(FATAL_ERROR "Cannot locate the python executable in ${PYTHON_EXEC}.")
endif()

FUNCTION(INSTALL_ARCHIVE ARCHIVE_FILENAME ARCHIVE_PREFIX)
    execute_process(
        COMMAND ${PYTHON_EXEC} ${CMAKE_SOURCE_DIR}/python/unzip.py ${ARCHIVE_FILENAME}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set(ARCHIVE_DIRPATH "${CMAKE_CURRENT_SOURCE_DIR}/${ARCHIVE_PREFIX}")

    foreach(SUBDIR_PREFIX "include" "bin" "lib")
        SET(SUBDIR_DIRPATH "${ARCHIVE_DIRPATH}/${SUBDIR_PREFIX}")
        if( EXISTS ${SUBDIR_DIRPATH} ) 
            install(
            DIRECTORY "${SUBDIR_DIRPATH}"
            DESTINATION ".")
        endif()
    endforeach()
ENDFUNCTION()

if( DEFINED WIN32 AND CMAKE_SIZEOF_VOID_P EQUAL 8 ) 
    INSTALL_ARCHIVE("hdf5.vs2010_x64.zip"   "hdf5")
    INSTALL_ARCHIVE("netcdf-4.2.vs2010_x64.zip"   "netcdf-4.2")
else()
    # TODO: Donwload the archive for source codes.
    message(WARNING "This script can only install NetCDF 4.2/HDF5 for Win64 platform.")
endif() 

