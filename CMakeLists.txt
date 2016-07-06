cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")

if (APPLE)
	# On OSX, brew now installs Boost 1.60 by default, which is
	# incompatible with our existing test code. Instead, we
	# install 1.55 and expand the include and link paths.
	set(ENV{CMAKE_PREFIX_PATH} "$ENV{CMAKE_PREFIX_PATH}:/usr/local/opt/boost")
endif()

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

enable_testing()

add_subdirectory(src)
