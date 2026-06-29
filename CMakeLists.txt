cmake_minimum_required(VERSION 3.10)

if(NOT DEFINED CMAKE_C_COMPILER)
  set(CMAKE_C_COMPILER gcc)
endif()

project(MEDS2endGen C)

find_package(OpenSSL REQUIRED)

if(DEBUG)
  set(COMMON_COMPILE_OPTIONS -g -Wall)
  set(COMMON_COMPILE_DEFINITIONS DEBUG)
else()
  set(COMMON_COMPILE_OPTIONS -O3 -Wall)
  set(COMMON_COMPILE_DEFINITIONS)
endif()

if(LOG)
  list(APPEND COMMON_COMPILE_DEFINITIONS DEBUG)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

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
  corank1.c
  trine_expand.c
  trine_codec.c
  bitstream.c
  randombytes.c
)

set(TRINE_LIB_TARGETS)
set(TRINE_KAT_TARGETS)
set(TRINE_SPEED_TARGETS)

function(add_trine_parameter_set suffix param_id)
  set(lib_target "meds2endgen_${suffix}")
  set(kat_target "kat_${suffix}")
  set(speed_target "speed_${suffix}")

  add_library(${lib_target} STATIC ${COMMON_SOURCES})
  target_compile_options(${lib_target} PRIVATE ${COMMON_COMPILE_OPTIONS})
  target_compile_definitions(${lib_target}
    PUBLIC PARAMS=${param_id}
    PRIVATE ${COMMON_COMPILE_DEFINITIONS}
  )
  target_include_directories(${lib_target} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/NIST"
  )
  target_link_libraries(${lib_target} PUBLIC OpenSSL::SSL OpenSSL::Crypto)

  add_executable(${kat_target} PQCgenKAT_sign.c)
  target_compile_options(${kat_target} PRIVATE ${COMMON_COMPILE_OPTIONS})
  target_compile_definitions(${kat_target} PRIVATE ${COMMON_COMPILE_DEFINITIONS})
  target_link_libraries(${kat_target} PRIVATE ${lib_target})

  add_executable(${speed_target}
    test_speed/test_speed.c
    test_speed/cpucycles.c
    test_speed/speed_print.c
  )
  target_compile_options(${speed_target} PRIVATE ${COMMON_COMPILE_OPTIONS})
  target_compile_definitions(${speed_target} PRIVATE ${COMMON_COMPILE_DEFINITIONS})
  target_include_directories(${speed_target} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/test_speed"
  )
  target_link_libraries(${speed_target} PRIVATE ${lib_target})

  add_custom_target("run_kat_${suffix}"
    COMMAND "$<TARGET_FILE:${kat_target}>"
    DEPENDS ${kat_target}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  )

  add_custom_target("run_speed_${suffix}"
    COMMAND "$<TARGET_FILE:${speed_target}>" 2 2
    DEPENDS ${speed_target}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  )

  set(TRINE_LIB_TARGETS ${TRINE_LIB_TARGETS} ${lib_target} PARENT_SCOPE)
  set(TRINE_KAT_TARGETS ${TRINE_KAT_TARGETS} ${kat_target} PARENT_SCOPE)
  set(TRINE_SPEED_TARGETS ${TRINE_SPEED_TARGETS} ${speed_target} PARENT_SCOPE)
endfunction()

add_trine_parameter_set(balanced_i 1)
add_trine_parameter_set(balanced_iii 2)
add_trine_parameter_set(balanced_v 3)
add_trine_parameter_set(shortsig_i 4)
add_trine_parameter_set(shortsig_iii 5)
add_trine_parameter_set(shortsig_v 6)

add_custom_target(all_parameter_sets
  DEPENDS
    ${TRINE_LIB_TARGETS}
    ${TRINE_KAT_TARGETS}
    ${TRINE_SPEED_TARGETS}
)

add_custom_target(all_kat_targets DEPENDS ${TRINE_KAT_TARGETS})
add_custom_target(all_speed_targets DEPENDS ${TRINE_SPEED_TARGETS})

add_custom_target(default ALL
  DEPENDS kat_balanced_i speed_balanced_i
)

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES
  "${CMAKE_BINARY_DIR}/PQCsignKAT_*"
)
