set(turbine_sources
  io/case_input.cpp
  io/engineering.cpp
  io/surface.cpp
  io/metrics_input.cpp
  analysis/statistics.cpp
  analysis/metrics.cpp
  inflow/grid_wind.cpp
  control/controller.cpp
  aerodynamics/bem.cpp
  aerodynamics/unsteady.cpp
  structure/blade_dynamics.cpp
  coupling/mesh.cpp
  coupling/rotor.cpp
  coupling/solver.cpp
  coupling/integrator.cpp
  coupling/modules.cpp
  coupling/structural_adapter.cpp
  simulation/acoustic_adapter.cpp
  simulation/results.cpp
  simulation/file_output.cpp
  simulation/simulation.cpp
  simulation/batch.cpp
)
add_library(turbine_modules STATIC ${turbine_sources})
target_include_directories(turbine_modules PUBLIC "${PROJECT_SOURCE_DIR}/include")
find_package(Threads REQUIRED)
target_link_libraries(turbine_modules PUBLIC aeroacoustics Threads::Threads)
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${turbine_sources})
