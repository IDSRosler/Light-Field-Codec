cmake_minimum_required(VERSION 3.0)
project(Light_Field_Codec)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(Light_Field_Codec
        src/CAVLC.cpp
        src/CAVLC.h
        src/clip.h
        src/DpcmDC.cpp
        src/DpcmDC.h
        src/EncBitstreamBuffer.cpp
        src/EncBitstreamBuffer.h
        src/EncBitstreamVLC.cpp
        src/EncBitstreamVLC.h
        src/EncBitstreamWriter.cpp
        src/EncBitstreamWriter.h
        src/EncoderParameters.cpp
        src/EncoderParameters.h
        src/LightField.cpp
        src/LightField.h
        src/LRE.cpp
        src/LRE.h
        src/main.cpp
        src/Point4D.cpp
        src/Point4D.h
        src/Prediction.cpp
        src/Prediction.h
        src/Quantization.cpp
        src/Quantization.h
        src/Statistics.cpp
        src/Statistics.h
        src/Time.cpp
        src/Time.h
        src/Transform.cpp
        src/Transform.h
        src/Typedef.h)
