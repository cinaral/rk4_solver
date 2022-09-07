#***********#
#* Project *#
#***********#
cmake_minimum_required(VERSION 3.13)

set(rk4_solver_VERSION 1.0.2)
set(matrix_VERSION 1.0.0)
set(matrix_rw_VERSION 1.0.1)

project(rk4_solver VERSION ${rk4_solver_VERSION} LANGUAGES CXX)

include(FetchContent)

#* matrix, used for basic matrix operations
FetchContent_Declare(matrix URL https://github.com/cinaral/matrix/releases/download/v${matrix_VERSION}/matrix-src.zip)
FetchContent_MakeAvailable(matrix)
set(matrix_INCLUDE_DIR ${matrix_SOURCE_DIR}/include)

#* only if we're in the root in order to not populate FetchContent()
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
	option(BUILD_TESTS "Build tests" ON)
	option(BUILD_EXAMPLES "Build examples" ON)
	option(BUILD_BENCHMARK "Build tests" ON)

	#****************************************#
	#* CMake, language and compiler options *#
	#****************************************#
	##* [Optional]:
	#set(CMAKE_CXX_STANDARD 17)
	#set(CMAKE_CXX_STANDARD_REQUIRED ON)
	#set(CMAKE_CXX_EXTENSIONS OFF)

	#* project compilation flags
	add_compile_options(
		-fdiagnostics-color=always #* colored debug output
		-fmessage-length=0 #* disable line wrapping (default=72 char) 
		-Wall #* all warnings
		-Wextra #* more warnings
		-Wpedantic #* more warnings
		#-D__USE_SINGLE_PRECISION__ #* use single precision floats over double precision
		#-fno-math-errno #* disable errno 
		##* [Optional]: for embedded systems the following flags may be enabled
		#-fno-rtti #* disable RTTI 
		#-fno-exceptions #* disable exceptions  
		#-fno-non-call-exceptions  
		#-fno-use-cxa-atexit #* disable destructors 
		#-fno-common #* declare globals once 
		#-ffreestanding #* standard library may be unavailable 
	)

	#* where to look for the project header and source files
	set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
	set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

	#* set up directories for binary, library and data files 
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
	file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/dat)
	
	#***********# 
	#* Testing *#
	#***********#
	if(BUILD_TESTS)
		#* matrix_rw, used for testing
		FetchContent_Declare(matrix_rw URL https://github.com/cinaral/matrix_rw/releases/download/v${matrix_rw_VERSION}/matrix_rw-src.zip)
		FetchContent_MakeAvailable(matrix_rw)
		set(matrix_rw_INCLUDE_DIR ${matrix_rw_SOURCE_DIR}/include)

		#* enable testing
		enable_testing()

		#* list of the tests
		set(TEST_NAMES
			test-rk4_solver-sine
			test-rk4_solver-motor
			test-rk4_solver-ball
		)

		#* loops through all tests
		foreach(ELEMENT ${TEST_NAMES})
			add_executable(${ELEMENT} test/${ELEMENT}.cpp)

			target_include_directories( ${ELEMENT}
				PUBLIC 
					${INCLUDE_DIR}
					${matrix_INCLUDE_DIR}
					${matrix_rw_INCLUDE_DIR} )

			add_test(NAME              ${ELEMENT} 
					COMMAND           ${ELEMENT}
					WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
		endforeach(ELEMENT ${TEST_NAMES})
	endif()

	#************#
	#* Examples *#
	#************#
	if(BUILD_EXAMPLES)
		#* list of the examples
		set(EXAMPLE_NAMES
			example_step
			example_loop
			example_event
		)

		#* loops through all examples
		foreach(ELEMENT ${EXAMPLE_NAMES})
			add_executable(${ELEMENT} examples/${ELEMENT}.cpp)

			target_include_directories(${ELEMENT}
				PUBLIC 
					${INCLUDE_DIR}	
					${matrix_INCLUDE_DIR}
			)
		endforeach(ELEMENT ${EXAMPLE_NAMES})
	endif()

	#*************#
	#* Benchmark *#
	#*************#
	if(BUILD_BENCHMARK)
		add_executable(rk4_solver-benchmark ${SRC_DIR}/rk4_solver-benchmark.cpp)

		target_include_directories(rk4_solver-benchmark
			PUBLIC 
				${INCLUDE_DIR}	
				${matrix_INCLUDE_DIR}
		)

		target_compile_options(rk4_solver-benchmark
			PUBLIC
				-O3
				-m64
				-mavx2
				#-ffast-math #* feeling brave?
		)
	endif()

	#************************#
	#* Packaged source code *#
	#************************#
	set(PACKAGE_FILES include/matrix.hpp include/rk4_solver.hpp CMakeLists.txt LICENSE)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip
		COMMAND ${CMAKE_COMMAND} -E tar c ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip --format=zip -- ${PACKAGE_FILES}
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		DEPENDS ${PACKAGE_FILES})
	add_custom_target(${PROJECT_NAME}_package DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip)

endif()