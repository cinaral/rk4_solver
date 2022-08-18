#****************************************#
#* Cmake, language and compiler options *#
#****************************************#
cmake_minimum_required(VERSION 3.13)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(
	-fdiagnostics-color=always #* colored debug output
	-fmessage-length=0 #* disable line wrapping (default=72 char) 
	-Wall #* all warnings
	-Wextra #* more warnings
	-Wpedantic #* more warnings
	-fno-math-errno #* disable errno
	-fno-rtti #* disable RTTI 
	-fno-exceptions #* disable exceptions 
	-fno-non-call-exceptions 
	-fno-use-cxa-atexit #* disable destructors
	-fno-common #* declare globals once
	-ffreestanding #* standard library may be unavailable
	-D__USE_SINGLE_PRECISION__
	#-fmax-errors=5 #* max number of errors displayed
	#-Wno-unused-variable
	#-Wno-unused-parameter
	#-Wno-unused-but-set-parameter
	#-fstack-protector-strong #* stack guard, costs one write and read per function, needs canary
)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
            "Default build type: Debug" FORCE)
endif()

#* where to look for the project header and source files
set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

#* set up directories for binary, library and data files 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/dat)

#* set the number of CPUs 
include(ProcessorCount)
ProcessorCount(N)

if(NOT N EQUAL 0)
	set(CTEST_BUILD_FLAGS -j${N})
	set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${N})
endif()
message("> Number of cores set to ${N}")

#* disable annoying deprecation warnings
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

#* enable testing
enable_testing()

#***********#
#* Project *#
#***********#

project(rk4_solver VERSION 1.0.0 LANGUAGES CXX)

add_executable(rk4_solver-benchmark ${SRC_DIR}/rk4_solver-benchmark.cpp)

target_include_directories(rk4_solver-benchmark
	PUBLIC 
		${INCLUDE_DIR}	
)

target_compile_options(rk4_solver-benchmark
	PUBLIC
		-O3
		-m64
		-mavx2
		#-ffast-math #* feeling brave?
)

 
#***********# 
#* Testing *#
#***********#

set(TEST_NAMES
	test-rk4_solver-sine
	test-rk4_solver-motor
	#test-rk4_solver-bouncing_ball
)

#* loops through all tests
foreach(ELEMENT ${TEST_NAMES})
	add_executable(${ELEMENT} test/${ELEMENT}.cpp)

	target_include_directories( ${ELEMENT}
		PUBLIC 
			${INCLUDE_DIR} )

	add_test(NAME              ${ELEMENT} 
		     COMMAND           ${ELEMENT}
		     WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endforeach(ELEMENT ${TEST_NAMES})

#************#
#* Examples *#
#************#

set(EXAMPLE_NAMES
	example_step
	example_loop
	#example_event
)

foreach(ELEMENT ${EXAMPLE_NAMES})
	add_executable(${ELEMENT} examples/${ELEMENT}.cpp)

	target_include_directories(${ELEMENT}
		PUBLIC 
			${INCLUDE_DIR}	
	)
endforeach(ELEMENT ${EXAMPLE_NAMES})
