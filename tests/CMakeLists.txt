find_package(Python3 3.11 REQUIRED COMPONENTS Interpreter)
option(AEROACOUSTICS_TEST_FORTRAN "Compare acoustic kernels with GFortran at test time" OFF)

function(add_probe name library)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE ${library})
endfunction()

# @WORK@ is a fresh directory per invocation, including repeated/parallel runs.
function(check name labels)
  add_test(NAME ${name} COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/run_test.py ${PROJECT_BINARY_DIR}/test-results ${name} ${ARGN})
  set_tests_properties(${name} PROPERTIES LABELS "${labels}" TIMEOUT 600)
  if(AEROACOUSTICS_USE_MKL)
    set_tests_properties(${name} PROPERTIES ENVIRONMENT "MKL_THREADING_LAYER=SEQUENTIAL")
  endif()
endfunction()

foreach(name auxiliary_probe mechanism_probe acoustic_workspace_probe)
  add_probe(${name} aeroacoustics)
endforeach()
check(mechanisms unit $<TARGET_FILE:mechanism_probe>)
check(acoustic_workspace unit $<TARGET_FILE:acoustic_workspace_probe> verify @WORK@/workspace.csv)

if(AEROACOUSTICS_BUILD_EXAMPLES)
  check(section_example example $<TARGET_FILE:aeroacoustics_example> @WORK@/spectrum.csv)
  check(bpm1989 example $<TARGET_FILE:aeroacoustics_bpm1989> @WORK@/bpm1989.csv --check)
  check(bpm1989_analysis example ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_bpm1989.py
    $<TARGET_FILE:aeroacoustics_bpm1989> @WORK@)
endif()

if(AEROACOUSTICS_TEST_FORTRAN)
  find_program(AEROACOUSTICS_GFORTRAN NAMES gfortran REQUIRED)
  check(fortran_kernels reference ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_fortran.py
    $<TARGET_FILE:aeroacoustics_shared> --fortran-compiler ${AEROACOUSTICS_GFORTRAN} --report @WORK@/report.json)
endif()

if(AEROACOUSTICS_BUILD_TURBINE)
  foreach(name optimization_probe optimization_benchmark trailing_edge_probe farm_probe validation_probe
      metrics_probe engineering_probe semantics_modules_probe performance_api_probe simulation_api_probe
      turbine_workspace_probe reliability_probe turbine_module_probe turbine_structure_probe
      turbine_mesh_probe turbine_solver_probe turbine_bem_probe)
    add_probe(${name} turbine_modules)
  endforeach()
  set(case_dir ${PROJECT_SOURCE_DIR}/examples/IEA_LB_RWT-AeroAcoustics)
  set(case ${case_dir}/IEA_LB_RWT-AeroAcoustics.fst)
  set(engineering ${PROJECT_SOURCE_DIR}/examples/engineering)
  check(allocations unit $<TARGET_FILE:optimization_probe>)
  check(validation unit $<TARGET_FILE:validation_probe>)
  check(trailing_edge "unit;concurrency" $<TARGET_FILE:trailing_edge_probe>)
  check(farm "unit;concurrency" $<TARGET_FILE:farm_probe> ${case})
  check(semantics unit $<TARGET_FILE:semantics_modules_probe> ${case})
  check(engineering "unit;concurrency" $<TARGET_FILE:engineering_probe> ${case} ${engineering}/controller.dat)
  check(metrics "unit;concurrency" $<TARGET_FILE:metrics_probe> ${case} ${engineering}/surfaces.dat @WORK@/metrics)
  check(performance_api "unit;concurrency" $<TARGET_FILE:performance_api_probe> ${case} @WORK@/performance)
  check(simulation_api unit $<TARGET_FILE:simulation_api_probe> ${case} @WORK@/simulation)
  check(turbine_workspace unit $<TARGET_FILE:turbine_workspace_probe> ${case})
  check(reliability unit $<TARGET_FILE:reliability_probe>
    ${case_dir}/Airfoils/RotorSE_FAST_IEA_landBased_RWT_AeroDyn_Polar_10.dat @WORK@/reliability)
  foreach(name reliability semantics engineering metrics)
    check(cli_${name} integration ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_${name}.py
      $<TARGET_FILE:aeroacoustics_turbine> @WORK@/output)
  endforeach()
  foreach(name performance validation farm trailing_edge)
    check(cli_${name} integration ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_${name}.py
      $<TARGET_FILE_DIR:aeroacoustics_turbine> @WORK@/output)
  endforeach()
  foreach(wind 8 9)
    check(reference_wind${wind} reference ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_reference.py
      $<TARGET_FILE:aeroacoustics_turbine> @WORK@/case --wind-speed ${wind} --report @WORK@/report.json)
  endforeach()
  if(AEROACOUSTICS_TEST_FORTRAN)
    check(fortran_trailing_edge reference ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_trailing_edge.py
      $<TARGET_FILE_DIR:trailing_edge_probe> @WORK@/output --fortran-compiler ${AEROACOUSTICS_GFORTRAN})
  endif()
endif()
