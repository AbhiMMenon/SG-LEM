cmake_minimum_required (VERSION 3.5)
project (splicTesteer)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(Threads REQUIRED)
find_package(yaml-cpp REQUIRED)

include_directories($ENV{CANTERAINCL} $ENV{LEMINCL})
link_directories($ENV{CANTERALIB} $ENV{LEMLIB})


set(source_dir "${PROJECT_SOURCE_DIR}/")
file(GLOB source_files "${source_dir}/*.cpp")

add_executable(spliceFlame  ${source_files})
target_link_libraries(spliceFlame PRIVATE LEM cantera_shared fmt Threads::Threads )

install(TARGETS spliceFlame DESTINATION ${source_dir})
