cmake_minimum_required(VERSION 3.10)
project(LF_Codec)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g3;-ggdb>"
  "$<$<CONFIG:RELEASE>:-O3>"
)

if(NOT "$ENV{TRANSF_QUANT}" STREQUAL "")
    message(STATUS "Using TRANSF_QUANT=$ENV{TRANSF_QUANT}")
    add_compile_definitions(TRANSF_QUANT=$ENV{TRANSF_QUANT})
endif()
if(NOT "$ENV{QUANTIZATION}" STREQUAL "")
    message(STATUS "Using QUANTIZATION=$ENV{QUANTIZATION}")
    add_compile_definitions(QUANTIZATION=$ENV{QUANTIZATION})
endif()
if(NOT "$ENV{DPCM_DC}" STREQUAL "")
    message(STATUS "Using DPCM_DC=$ENV{DPCM_DC}")
    add_compile_definitions(DPCM_DC=$ENV{DPCM_DC})
endif()
if(NOT "$ENV{STATISTICS_LOCAL}" STREQUAL "")
    message(STATUS "Using STATISTICS_LOCAL=$ENV{STATISTICS_LOCAL}")
    add_compile_definitions(STATISTICS_LOCAL=$ENV{STATISTICS_LOCAL})
endif()
if(NOT "$ENV{STATISTICS_GLOBAL}" STREQUAL "")
    message(STATUS "Using STATISTICS_GLOBAL=$ENV{STATISTICS_GLOBAL}")
    add_compile_definitions(STATISTICS_GLOBAL=$ENV{STATISTICS_GLOBAL})
endif()
if(NOT "$ENV{STATISTICS_TIME}" STREQUAL "")
    message(STATUS "Using STATISTICS_TIME=$ENV{STATISTICS_TIME}")
    add_compile_definitions(STATISTICS_TIME=$ENV{STATISTICS_TIME})
endif()
if(NOT "$ENV{TRACE_TRANSF}" STREQUAL "")
    message(STATUS "Using TRACE_TRANSF=$ENV{TRACE_TRANSF}")
    add_compile_definitions(TRACE_TRANSF=$ENV{TRACE_TRANSF})
endif()
if(NOT "$ENV{TRACE_QUANT}" STREQUAL "")
    message(STATUS "Using TRACE_QUANT=$ENV{TRACE_QUANT}")
    add_compile_definitions(TRACE_QUANT=$ENV{TRACE_QUANT})
endif()
if(NOT "$ENV{TRACE_LRE}" STREQUAL "")
    message(STATUS "Using TRACE_LRE=$ENV{TRACE_LRE}")
    add_compile_definitions(TRACE_LRE=$ENV{TRACE_LRE})
endif()
if(NOT "$ENV{LFCODEC_USE_PREDICTION}" STREQUAL "")
    message(STATUS "Using LFCODEC_USE_PREDICTION=$ENV{LFCODEC_USE_PREDICTION}")
    add_compile_definitions(LFCODEC_USE_PREDICTION=$ENV{LFCODEC_USE_PREDICTION})
endif()
if(NOT "$ENV{LFCODEC_FORCE_DCT_NON_LUMA}" STREQUAL "")
    message(STATUS "Using LFCODEC_FORCE_DCT_NON_LUMA=$ENV{LFCODEC_FORCE_DCT_NON_LUMA}")
    add_compile_definitions(LFCODEC_FORCE_DCT_NON_LUMA=$ENV{LFCODEC_FORCE_DCT_NON_LUMA})
endif()
if(NOT "$ENV{LFCODEC_TRANSFORM_HISTOGRAM}" STREQUAL "")
    message(STATUS "Using LFCODEC_TRANSFORM_HISTOGRAM=$ENV{LFCODEC_TRANSFORM_HISTOGRAM}")
    add_compile_definitions(LFCODEC_TRANSFORM_HISTOGRAM=$ENV{LFCODEC_TRANSFORM_HISTOGRAM})
endif()
if(NOT "$ENV{LFCODEC_QUANTIZATION_EXTRA_VOLUMES}" STREQUAL "")
    message(STATUS "Using LFCODEC_QUANTIZATION_EXTRA_VOLUMES=$ENV{LFCODEC_QUANTIZATION_EXTRA_VOLUMES}")
    add_compile_definitions(LFCODEC_QUANTIZATION_EXTRA_VOLUMES=$ENV{LFCODEC_QUANTIZATION_EXTRA_VOLUMES})
endif()



find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()


enable_testing()
add_subdirectory(extern/googletest)
add_subdirectory(extern/spdlog)
add_subdirectory(test)

set(SOURCES
    src/EncoderParameters.cpp
    src/Point4D.cpp
    src/Transform.cpp
    src/Quantization.cpp
    src/LightField.cpp
    src/Statistics.cpp
    src/LRE.cpp
    src/DpcmDC.cpp
    src/Timer.cpp
    src/CAVLC.cpp
    src/EncBitstreamWriter.cpp
    src/EncBitstreamVLC.cpp
    src/EncBitstreamBuffer.cpp
    src/Prediction.cpp
    src/deprecated/OldTransform.cpp
    src/utils.cpp
    src/ScanOrder.cpp
    src/EntropyEncoder.cpp
    src/EncSymbol.cpp
    src/ArithmeticEncoder.cpp
    src/Tree.cpp
        src/SubpartitionModel.cpp src/SubpartitionModel.h src/ArithmeticEntropyEncoder.cpp src/ArithmeticEntropyEncoder.h)


add_executable(LF_Codec src/main.cpp)

add_library(lfcodec ${SOURCES})
target_include_directories(lfcodec PRIVATE src)
target_link_libraries(LF_Codec PUBLIC
    lfcodec
    spdlog::spdlog
)

find_package(OpenCV QUIET)
add_compile_definitions(LFCODEC_USE_OPENCV=${OPENCV_FOUND})

if(OPENCV_FOUND)
    message(STATUS "OpenCV Found")
    target_link_libraries(lfcodec PUBLIC ${OpenCV_LIBS})
endif()
