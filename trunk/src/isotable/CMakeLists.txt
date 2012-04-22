cmake_minimum_required (VERSION 2.8)
project (isotable)

add_subdirectory(ijkmcube.v0.1/ijkmcube)	
install(
	DIRECTORY "isosurface_table.v3"
	DESTINATION ".")

#
# $Log$
#
