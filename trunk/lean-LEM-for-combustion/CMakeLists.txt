cmake_minimum_required (VERSION 3.5)
project (LEMlibrary)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(Threads REQUIRED)

include_directories($ENV{CANTERAINCL})
link_directories($ENV{CANTERALIB})

#add_executable(demo demo.cpp)

set(source_dir "${PROJECT_SOURCE_DIR}/")
file(GLOB source_files "${source_dir}/*.cpp")
list(REMOVE_ITEM source_files "lemPremixedFlameTest.cpp")

add_library(LEMlib STATIC ${source_files})

target_link_libraries(LEMlib cantera_shared fmt Threads::Threads)
#add_executable(test1  lemPremixedFlameTest.cpp)

#target_link_libraries(test1 LEMlib cantera_shared fmt Threads::Threads)

install(TARGETS LEMlib DESTINATION ${source_dir}/lib)
