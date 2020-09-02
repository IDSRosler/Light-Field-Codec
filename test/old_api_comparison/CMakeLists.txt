cmake_minimum_required(VERSION 3.15)

file(GLOB SOURCES *.cpp)
get_filename_component(CURRENT_TEST ${CMAKE_CURRENT_SOURCE_DIR} NAME) 
add_executable(${CURRENT_TEST} ${SOURCES})
target_include_directories(${CURRENT_TEST} PUBLIC ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(${CURRENT_TEST} PUBLIC ${libs_to_link})

add_test(
    NAME ${CURRENT_TEST}
    COMMAND ${CURRENT_TEST}
)