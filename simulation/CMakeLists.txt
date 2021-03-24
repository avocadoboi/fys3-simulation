cmake_minimum_required(VERSION 3.18.0)
project(Simulation VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin/)

#-----------------------------
# Find/add dependencies

add_subdirectory(external/matplotplusplus)
find_package(mp-units CONFIG REQUIRED)
find_package(fmt REQUIRED)

#-----------------------------

add_executable(pingisboll source/pingisboll.cpp)
add_executable(oscillation source/svängning.cpp)
add_executable(golfboll source/golfboll.cpp)

set(TARGETS pingisboll oscillation golfboll)

foreach (TARGET ${TARGETS})
	target_include_directories(${TARGET} PRIVATE include)
	target_precompile_headers(${TARGET} PRIVATE include/common.hpp)

	target_compile_features(${TARGET} PUBLIC cxx_std_20)
	set_target_properties(${TARGET} PROPERTIES CXX_EXTENSIONS off)

	#-----------------------------
	# Warnings for all targets.

	if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
		target_compile_options(${TARGET} PRIVATE
			/experimental:external
			/external:anglebrackets
			/external:W0 
			/WX /Wall
		)
	else ()
		target_compile_options(${TARGET} PRIVATE
			-Werror -Wall -Wpedantic -Wextra 
			-Wimplicit-fallthrough=5 
			-Wduplicated-branches -Wduplicated-cond 
			-Wcast-qual -Wconversion
			
			-Wno-parentheses
			-Wno-missing-field-initializers
		)
	endif ()

	#-----------------------------
	# Link any dependencies...

	target_link_libraries(${TARGET} PRIVATE matplot)
	
	target_link_libraries(${TARGET} PRIVATE mp-units::mp-units)

	target_link_libraries(${TARGET} PRIVATE fmt::fmt)
endforeach(TARGET)