project("sptial ancestry analysis")

get_filename_component(ProjectName ${CMAKE_CURRENT_LIST_DIR} NAME)

FILE(GLOB SourceFiles *.h *.cpp)
list(REMOVE_ITEM SourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

find_package(Boost 1.60.0 REQUIRED COMPONENTS thread program_options)

add_library(${ProjectName}_Lib STATIC ${SourceFiles})

include_directories(${Boost_INCLUDE_DIRS})

add_executable(${ProjectName} main.cpp)

target_link_libraries(${ProjectName} ${ProjectName}_Lib)
target_link_libraries(${ProjectName} ${Boost_LIBRARIES})

install(TARGETS ${ProjectName} DESTINATION bin)
