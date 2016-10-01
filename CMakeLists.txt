project(MachUp Fortran)
cmake_minimum_required(VERSION 3.5)

# command-line options
option(LTO "Turn on for Link-Time Optimization (default:OFF)")
set(ndv 0 CACHE STRING "Number of DNAD Derivatives")


# setup build type information
if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE "Release" CACHE STRING
		"No build type set! Valid types: Debug Release RelWithDebInfo MinSizeRel" 
	FORCE )
endif()

# input lists
include(common/CMakeLists.txt)
include(src/CMakeLists.txt)

# set the executable name
# There is a separate executable for non-Release versions
set( EXEC "Machup" )
if( NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
	set( EXEC "${EXEC}_${CMAKE_BUILD_TYPE}" )
endif()

# executable information
set( EXEC "${EXEC}.out")
set( EXEC_DIR "bin/")

add_executable(${EXEC} ${CommonFiles} ${SrcFiles})
set_target_properties(${EXEC} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${EXEC_DIR}")

# using the gfortran compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.9)
		message(FATAL_ERROR "gfortran version must be 4.9 or higher")
	endif()

	# if link-time optimization has been requested
	if( LTO )
		set( LTO_FLAGS "-flto" )
	endif()

	# necessary compilation flags with gfortran
	set( GNU_FLAGS "-cpp -ffree-line-length-none -fdefault-real-8 ${LTO_FLAGS}" )
	set_target_properties(${EXEC} PROPERTIES COMPILE_FLAGS "${GNU_FLAGS}")
	set_target_properties(${EXEC} PROPERTIES LINK_FLAGS "${LTO_FLAGS}")

	# release time optimization flags
	set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -march=native" )

	# debug flags
	set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -fcheck=all")
	set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -Wall -fcheck=all")

# else, error, because I don't know the right flags for this compiler!
else()
	message(FATAL_ERROR "Build script currently only works with gfortran!")
endif()

# DNAD definitions
if( ndv )
	add_definitions(-Ddnad -Dndv=${ndv})
endif()


