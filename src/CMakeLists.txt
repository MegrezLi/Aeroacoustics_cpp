# Shared source list keeps static and shared acoustic libraries in sync.
set(acoustic_sources
  acoustics/boundary_layer.cpp
  acoustics/directivity.cpp
  acoustics/bpm_spectral_shapes.cpp
  acoustics/bpm_trailing_edge.cpp
  acoustics/bpm_other_sources.cpp
  acoustics/inflow_noise.cpp
  acoustics/kernel_compat.cpp
  acoustics/tno.cpp
  acoustics/howe_chase.cpp
  acoustics/trailing_edge_options.cpp
  acoustics/quantities.cpp
  propagation/outdoor.cpp
  acoustics/spectrum.cpp
  acoustics/workspace.cpp
  acoustics/geometry.cpp
  acoustics/driver.cpp
  numerics/quadrature.cpp
  numerics/backend.cpp
  interfaces/c_api.cpp
)
add_library(acoustic_objects OBJECT ${acoustic_sources})
set_target_properties(acoustic_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(acoustic_objects PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/aeroacoustics>)
if(AEROACOUSTICS_USE_MKL)
  target_compile_definitions(acoustic_objects PUBLIC AERO_USE_MKL)
  target_link_libraries(acoustic_objects PUBLIC MKL::MKL)
endif()
add_library(aeroacoustics STATIC)
add_library(aeroacoustics_shared SHARED)
foreach(library aeroacoustics aeroacoustics_shared)
  target_link_libraries(${library} PUBLIC acoustic_objects)
endforeach()
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${acoustic_sources})

if(AEROACOUSTICS_BUILD_TURBINE)
  add_subdirectory(turbine)
endif()
add_subdirectory(apps)
