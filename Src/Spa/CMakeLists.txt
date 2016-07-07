project("sptial ancestry analysis")

get_filename_component(ProjectName ${CMAKE_CURRENT_LIST_DIR} NAME)

FILE(GLOB SourceFiles *.h *.cpp)
FILE(GLOB TestSourceFiles Tests/*.h Tests/*.cpp)
list(REMOVE_ITEM SourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

find_package(Boost 1.60.0 REQUIRED COMPONENTS thread program_options filesystem unit_test_framework)

add_library(${ProjectName}Lib STATIC ${SourceFiles})

include_directories(${Boost_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_LIST_DIR})

add_executable(${ProjectName} main.cpp)
add_executable(${ProjectName}Tests ${TestSourceFiles})

target_link_libraries(${ProjectName} ${ProjectName}Lib)
target_link_libraries(${ProjectName} ${Boost_LIBRARIES})

target_link_libraries(${ProjectName}Tests ${ProjectName}Lib)
target_link_libraries(${ProjectName}Tests ${Boost_LIBRARIES})

add_test(NAME ${ProjectName}Tests COMMAND ${ProjectName}Tests --home ${CMAKE_CURRENT_LIST_DIR})

install(TARGETS ${ProjectName} DESTINATION bin)
