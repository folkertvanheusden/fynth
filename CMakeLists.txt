project(fynth)
cmake_minimum_required(VERSION 3.2)
add_definitions("-std=gnu++17 -ggdb3")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fstack-protector-strong -fno-inline")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer-fno-inline")

add_executable(
  fynth
  error.cpp
  fft.cpp
  filter.cpp
  main.cpp
  sample.cpp
  sf2.cpp
  terminal.cpp
  utils.cpp
)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads)
target_link_libraries(fynth Threads::Threads)

include(FindPkgConfig)

pkg_check_modules(ALSA REQUIRED alsa)
target_link_libraries(fynth ${ALSA_LIBRARIES})
target_include_directories(fynth PUBLIC ${ALSA_INCLUDE_DIRS})
target_compile_options(fynth PUBLIC ${ALSA_CFLAGS_OTHER})

pkg_check_modules(SNDFILE REQUIRED sndfile)
target_link_libraries(fynth ${SNDFILE_LIBRARIES})
target_include_directories(fynth PUBLIC ${SNDFILE_INCLUDE_DIRS})
target_compile_options(fynth PUBLIC ${SNDFILE_CFLAGS_OTHER})

find_package(Curses REQUIRED)
pkg_check_modules(NCURSES libpipewire-0.3)
target_link_libraries(fynth ${CURSES_LIBRARIES})
target_include_directories(fynth PUBLIC ${CURSES_INCLUDE_DIRS})
target_compile_options(fynth PUBLIC ${CURSES_CFLAGS_OTHER})

find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)
target_link_libraries(fynth fftw3)

# pipewire
find_package(PipeWire)
pkg_check_modules(PIPEWIRE libpipewire-0.3)
target_link_libraries(fynth ${PIPEWIRE_LIBRARIES})
target_include_directories(fynth PUBLIC ${PIPEWIRE_INCLUDE_DIRS})
target_compile_options(fynth PUBLIC ${PIPEWIRE_CFLAGS_OTHER})

set_target_properties(fynth PROPERTIES OUTPUT_NAME fynth)

install(TARGETS fynth
    COMPONENT fynth
    RUNTIME DESTINATION "/usr/bin/"
)

#install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE" "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
#        DESTINATION "/usr/share/fynth"
#        COMPONENT "fynth")

set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Folkert van Heusden <mail@vanheusden.com>")
set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_DEBIAN_COMPRESSION_TYPE "lzma")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_STRIP_FILES OFF)
set(CPACK_DEB_AUTO_GENERATED_FILE_NAME ON)
set(CPACK_DEBIAN_DEBUGINFO_PACKAGE ON)
set(CPACK_DEBIAN_PACKAGE_DEBUG ON)
set(CPACK_DEBIAN_RUNTIME_DEBUGINFO_PACKAGE ON)
set(CPACK_DEBIAN_fynth_DEBUGINFO_PACKAGE ON)

include(CPack)
