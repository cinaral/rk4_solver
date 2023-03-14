cmake_minimum_required(VERSION 3.13)

project(
	rk4_solver
	VERSION 1.0.11
	LANGUAGES CXX
)

#* dependencies
include(FetchContent)

set(matrix_op_VERSION 1.0.7) #* matrix_op, used for basic matrix operations
set(matrix_rw_VERSION 1.0.5) #* matrix_rw, used for testing

FetchContent_Declare(
	matrix_op 
	URL https://github.com/cinaral/matrix_op/releases/download/v${matrix_op_VERSION}/src.zip
)
FetchContent_MakeAvailable(matrix_op)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
	FetchContent_Declare(
		matrix_rw 
		URL https://github.com/cinaral/matrix_rw/releases/download/v${matrix_rw_VERSION}/src.zip
	)
	FetchContent_MakeAvailable(matrix_rw)
endif()



if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
	option(DEBUG_BUILD "Build everything in debug mode" ON)
	option(BUILD_TESTS "Build tests" ON)
	option(BUILD_EXAMPLES "Build examples" ON)
	option(BUILD_BENCHMARKS "Build benchmarks" ON)
	option(SINGLE_PRECISION "Single precision" OFF)
	option(HEAP_USAGE "Heap usage" OFF)

	#* where to look for the project header and source files
	set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
	set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
	set(TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/test)
	set(EXAMPLE_DIR ${CMAKE_CURRENT_LIST_DIR}/examples)
	set(BENCHMARK_DIR ${CMAKE_CURRENT_LIST_DIR}/benchmarks)

	#* tests, examples, benchmarks to compile
	set(TEST_NAMES
		sine-test
		motor-test
		ball-test
		first_order-test
	)
	set(EXAMPLE_NAMES
		step-example
		loop-example
		event-example
		first_order-example
	)
	set(BENCHMARK_NAMES
		benchmark
	)

	#* files to package
	set(PACKAGE_FILES 
		include/${PROJECT_NAME}/
		include/${PROJECT_NAME}.hpp 
		CMakeLists.txt 
		LICENSE
	)

	#* set up output directories
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
	file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/dat)

	#* project compile options
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	set(CMAKE_CXX_EXTENSIONS OFF)

	add_compile_options(
		-fdiagnostics-color=always #* colored debug output
		-fmessage-length=0 #* disable line wrapping (default=72 char) 
		-Wall #* all warnings
		-Wextra #* more warnings
		-Wpedantic #* more warnings
		-Wno-array-bounds #* disable false array bounds warning
	)
	
	if(DEBUG_BUILD)
		add_compile_options(
			-g #* debug symbols
			-O0 #* no optimization, -Og may work better with some debuggers
		)
	else()
		add_compile_options(
			-O3 #* release optimization
		)
	endif()

	if(SINGLE_PRECISION)
		add_compile_options(-DUSE_SINGLE_PRECISION) #* use single precision floats
	endif()

	if(HEAP_USAGE)
		add_compile_options(-DDO_USE_HEAP)
	endif()
	
	#***********#
	#* Testing *#
	#***********#
	if(BUILD_TESTS)
		enable_testing()

		foreach(ELEMENT ${TEST_NAMES})
			add_executable(
				${ELEMENT}
				${TEST_DIR}/${ELEMENT}.cpp						
			)
			target_include_directories(
				${ELEMENT} PRIVATE 
				${INCLUDE_DIR}
				${matrix_op_SOURCE_DIR}/include
				${matrix_rw_SOURCE_DIR}/include
			)
			target_compile_options(
				${ELEMENT} PRIVATE 
			)
			add_test(
				NAME ${ELEMENT} 
				COMMAND ${ELEMENT}
				WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
			)
		endforeach(ELEMENT ${TEST_NAMES})
	endif()

	#************#
	#* Examples *#
	#************#
	if(BUILD_EXAMPLES)
		foreach(ELEMENT ${EXAMPLE_NAMES})
			add_executable(
				${ELEMENT} 
				${EXAMPLE_DIR}/${ELEMENT}.cpp
			)
			target_include_directories(
				${ELEMENT} PRIVATE 
				${INCLUDE_DIR}	
				${matrix_op_SOURCE_DIR}/include
			)
		endforeach(ELEMENT ${EXAMPLE_NAMES})
	endif()

	#*************#
	#* Benchmark *#
	#*************#
	if(BUILD_BENCHMARKS)
		foreach(ELEMENT ${BENCHMARK_NAMES})
			add_executable(
				${ELEMENT} 
				${BENCHMARK_DIR}/${ELEMENT}.cpp
			)
			target_include_directories(
				${ELEMENT} PRIVATE 
				${INCLUDE_DIR}	
				${matrix_op_SOURCE_DIR}/include
			)
			target_compile_options(
				${ELEMENT} PRIVATE
				-O3 #* optimization level
				-m64 #* x64
				-mavx2 #* enable avx2
				-fno-exceptions #* disable exceptions  
				#-fno-math-errno #* disable errno 
				#-ffast-math #* feeling brave? (may not improve performance)
			)
		endforeach(ELEMENT ${BENCHMARK_NAMES})
	endif()

	#************************#
	#* Packaged source code *#
	#************************#
	add_custom_command(
		OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/src.zip
		COMMAND ${CMAKE_COMMAND} -E tar c ${CMAKE_CURRENT_BINARY_DIR}/src.zip --format=zip -- ${PACKAGE_FILES}
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		DEPENDS ${PACKAGE_FILES}
	)
	add_custom_target(
		src_package 
		DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/src.zip
	)
endif()