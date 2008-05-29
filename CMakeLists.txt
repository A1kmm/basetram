SET(CMAKE_CXX_FLAGS_PROFILE "-O3 -DNDEBUG -pg -ggdb" CACHE STRING "Profiling compiler settings")
SET(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg" CACHE STRING "Profiling linker settings")

FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ../parsegenbank)

ADD_EXECUTABLE(basetram BaSeTraM.cpp)
TARGET_LINK_LIBRARIES(basetram boost_serialization boost_program_options boost_iostreams boost_regex boost_filesystem boost_mpi boost_system GenBankParser bz2 z mpi)
