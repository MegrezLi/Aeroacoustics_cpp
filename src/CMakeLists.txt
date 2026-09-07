# Shared source list keeps static and shared acoustic libraries in sync.
set(acoustic_sources
  acoustics/empirical_models.cpp
  acoustics/tno.cpp
  acoustics/spectrum.cpp
  acoustics/geometry.cpp
  acoustics/driver.cpp
  numerics/quadrature.cpp
  numerics/backend.cpp
  interfaces/c_api.cpp
)
add_library(aeroacoustics STATIC ${acoustic_sources})
add_library(aeroacoustics_shared SHARED ${acoustic_sources})
foreach(library aeroacoustics aeroacoustics_shared)
  target_include_directories(${library} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/aeroacoustics>)
  if(AEROACOUSTICS_USE_MKL)
    target_compile_definitions(${library} PUBLIC AERO_USE_MKL)
    target_link_libraries(${library} PUBLIC MKL::MKL)
  endif()
endforeach()
set_target_properties(aeroacoustics PROPERTIES POSITION_INDEPENDENT_CODE ON)
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${acoustic_sources})

if(AEROACOUSTICS_BUILD_TURBINE)
  add_subdirectory(turbine)
endif()
add_subdirectory(apps)
