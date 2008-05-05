FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ../parsegenbank)

ADD_EXECUTABLE(basetram BaSeTraM.cpp)
TARGET_LINK_LIBRARIES(basetram boost_program_options boost_iostreams boost_regex boost_filesystem GenBankParser)
