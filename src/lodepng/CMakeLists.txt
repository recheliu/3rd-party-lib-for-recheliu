# The files here are copied from http://lodev.org/lodepng/. The author of this library is Lode Vandevenne. I just add a CMakefile for my own usage.

cmake_minimum_required(VERSION 2.8)
PROJECT(lodepng)

INCLUDE_REGULAR_EXPRESSION("^.*.h$")

add_library(${PROJECT_NAME}
	lodepng.h
	lodepng.cpp
)

set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME "${PROJECT_NAME}d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}"
)

INSTALL(TARGETS ${PROJECT_NAME} DESTINATION lib)
INSTALL(FILES "lodepng.h" 		DESTINATION include)

