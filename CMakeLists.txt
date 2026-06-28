cmake_minimum_required(VERSION 3.10)

if(NOT DEFINED CMAKE_C_COMPILER)
  set(CMAKE_C_COMPILER gcc)
endif()

project(MEDS9923 C)

set(LIBS -lssl -lcrypto)

if(DEBUG)
  set(COMPILE_OPTIONS -g -Wall)
  set(COMPILE_DEFINITIONS DEBUG)
  set(OBJDIR debug)
else()
  set(COMPILE_OPTIONS -O3 -Wall)
  set(OBJDIR build)
endif()

if(LOG)
  list(APPEND COMPILE_DEFINITIONS DEBUG)
endif()

set(OBJDIR_PATH "${CMAKE_SOURCE_DIR}/${OBJDIR}")

set(COMMON_SOURCES
  meds.c
  util.c
  osfreq.c
  fips202.c
  field.c
  matrixmod.c
  matrixelim.c
  triform.c
  canonical.c
  bitstream.c
  randombytes.c
)

set(EXES bench PQCgenKAT_sign)

if(EXISTS "${CMAKE_SOURCE_DIR}/test_meds2end_seedless.c")
  list(APPEND EXES test_meds2end_seedless)
endif()

file(GLOB KAT_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/PQCsignKAT_*")

foreach(EXE IN LISTS EXES)
  add_executable(${EXE} ${EXE}.c ${COMMON_SOURCES})
  target_compile_options(${EXE} PRIVATE ${COMPILE_OPTIONS})
  target_compile_definitions(${EXE} PRIVATE ${COMPILE_DEFINITIONS})
  target_include_directories(${EXE} PRIVATE
    "${OBJDIR_PATH}"
    "${CMAKE_SOURCE_DIR}/NIST"
  )
  target_link_libraries(${EXE} PRIVATE ${LIBS})
  set_target_properties(${EXE} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${OBJDIR_PATH}"
  )
endforeach()

add_custom_target(default ALL DEPENDS bench PQCgenKAT_sign)

add_executable(test_speed
  test_speed/test_speed.c
  test_speed/cpucycles.c
  test_speed/speed_print.c
  ${COMMON_SOURCES}
)
target_compile_options(test_speed PRIVATE ${COMPILE_OPTIONS})
target_compile_definitions(test_speed PRIVATE ${COMPILE_DEFINITIONS})
target_include_directories(test_speed PRIVATE
  "${OBJDIR_PATH}"
  "${CMAKE_SOURCE_DIR}/NIST"
  "${CMAKE_SOURCE_DIR}/test_speed"
)
target_link_libraries(test_speed PRIVATE ${LIBS})
set_target_properties(test_speed PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${OBJDIR_PATH}"
)

add_custom_target(run_kat
  COMMAND "$<TARGET_FILE:PQCgenKAT_sign>"
  DEPENDS PQCgenKAT_sign
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)

add_custom_target(run_bench
  COMMAND "$<TARGET_FILE:bench>"
  DEPENDS bench
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)

add_custom_target(run_speed
  COMMAND "$<TARGET_FILE:test_speed>"
  DEPENDS test_speed
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES
  "${CMAKE_SOURCE_DIR}/build"
  "${CMAKE_SOURCE_DIR}/debug"
  ${KAT_FILES}
)
