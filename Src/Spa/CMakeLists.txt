project("sptial ancestry analysis")

get_filename_component(ProjectName ${CMAKE_CURRENT_LIST_DIR} NAME)

FILE(GLOB SourceFiles *.h *.cpp)
FILE(GLOB TestSourceFiles Tests/*.h Tests/*.cpp)
list(REMOVE_ITEM SourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

find_package(Boost 1.60.0 REQUIRED COMPONENTS thread program_options unit_test_framework)

add_library(${ProjectName}_Lib STATIC ${SourceFiles})

include_directories(${Boost_INCLUDE_DIRS})

add_executable(${ProjectName} main.cpp)
add_executable(${ProjectName}_Tests ${TestSourceFiles})

target_link_libraries(${ProjectName} ${ProjectName}_Lib)
target_link_libraries(${ProjectName} ${Boost_LIBRARIES})

target_link_libraries(${ProjectName}_Tests ${ProjectName}_Lib)
target_link_libraries(${ProjectName}_Tests ${Boost_LIBRARIES})

add_test(NAME ${ProjectName}_Tests COMMAND ${ProjectName}_Tests --home ${CMAKE_CURRENT_LIST_DIR})

install(TARGETS ${ProjectName} DESTINATION bin)
