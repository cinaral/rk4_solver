#***********#
#* Project *#
#***********#
cmake_minimum_required(VERSION 3.13)

project(
	rk4_solver
	VERSION 1.0.3 
	LANGUAGES CXX
)

include(FetchContent)

#* matrix_op, used for basic matrix operations
set(matrix_op_VERSION 1.0.1)

FetchContent_Declare(matrix_op URL https://github.com/cinaral/matrix_op/releases/download/v${matrix_op_VERSION}/matrix_op-src.zip)
FetchContent_MakeAvailable(matrix_op)
set(matrix_op_INCLUDE_DIR ${matrix_op_SOURCE_DIR}/include)

#********#
#* Self *#
#********#
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
	option(SINGLE_PRECISION "Single precision" OFF)
	option(HEAP_USAGE "Single precision" OFF)
	option(BUILD_TESTS "Build tests" ON)
	option(BUILD_EXAMPLES "Build examples" ON)
	option(BUILD_BENCHMARKS "Build bencmarks" ON)

	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	set(CMAKE_CXX_EXTENSIONS OFF)

	add_compile_options(
		-fdiagnostics-color=always #* colored debug output
		-fmessage-length=0 #* disable line wrapping (default=72 char) 
		-Wall #* all warnings
		-Wextra #* more warnings
		-Wpedantic #* more warnings
	)

	if(SINGLE_PRECISION)
		add_compile_options(-D__USE_SINGLE_PRECISION__) #* use single precision floats
	endif()

	if(HEAP_USAGE)
		add_compile_options(-D__DO_USE_HEAP__) #* use single precision floats
	endif()

	#* where to look for the project header and source files
	set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
	set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

	#* set up directories for binary, library, data and parameter files 
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
	file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/dat)
	
	#***********#
	#* Testing *#
	#***********#
	if(BUILD_TESTS)
		set(TEST_NAMES
			test-rk4_solver-sine
			test-rk4_solver-motor
			test-rk4_solver-ball
		)
		enable_testing()

		#* matrix_rw, used for testing
		set(matrix_rw_VERSION 1.0.2)
		
		FetchContent_Declare(
			matrix_rw 
			URL https://github.com/cinaral/matrix_rw/releases/download/v${matrix_rw_VERSION}/matrix_rw-src.zip
		)
		FetchContent_MakeAvailable(matrix_rw)
		set(matrix_rw_INCLUDE_DIR ${matrix_rw_SOURCE_DIR}/include)

		foreach(ELEMENT ${TEST_NAMES})
			add_executable(
				${ELEMENT} 
				test/${ELEMENT}.cpp						
			)

			target_include_directories(
				${ELEMENT} PRIVATE 
				${INCLUDE_DIR}
				${matrix_op_INCLUDE_DIR}
				${matrix_rw_INCLUDE_DIR}
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
		set(EXAMPLE_NAMES
			example_step
			example_loop
			example_event
		)

		foreach(ELEMENT ${EXAMPLE_NAMES})
			add_executable(
				${ELEMENT} 
				examples/${ELEMENT}.cpp)

			target_include_directories(
				${ELEMENT} PRIVATE 
				${INCLUDE_DIR}	
				${matrix_op_INCLUDE_DIR}
			)
		endforeach(ELEMENT ${EXAMPLE_NAMES})
	endif()

	#*************#
	#* Benchmark *#
	#*************#
	if(BUILD_BENCHMARKS)
		set(BENCHMARK_NAMES
			rk4_solver-benchmark
		)

		foreach(ELEMENT ${BENCHMARK_NAMES})
			add_executable(
				${ELEMENT} 
				src/${ELEMENT}.cpp
			)

			target_include_directories(
				${ELEMENT} PRIVATE 
				${INCLUDE_DIR}	
				${matrix_op_INCLUDE_DIR}
			)

			target_compile_options(
				${ELEMENT} PUBLIC
				-O3 #* optimization level
				-m64 #* x64
				-mavx2 #* enable avx2
				#-fno-math-errno #* disable errno 
				#-ffast-math #* feeling brave?
			)
		endforeach(ELEMENT ${BENCHMARK_NAMES})
	endif()

	#************************#
	#* Packaged source code *#
	#************************#
	set(PACKAGE_FILES 
		include/${PROJECT_NAME}/
		include/${PROJECT_NAME}.hpp 
		CMakeLists.txt 
		LICENSE
	)

	add_custom_command(
		OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip
		COMMAND ${CMAKE_COMMAND} -E tar c ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip --format=zip -- ${PACKAGE_FILES}
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		DEPENDS ${PACKAGE_FILES}
	)
	add_custom_target(${PROJECT_NAME}_package DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip)

endif()