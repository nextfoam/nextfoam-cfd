##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

set(unit_tests
  UnitTestOpenMPArrayHandle.cxx
  UnitTestOpenMPArrayHandleFancy.cxx
  UnitTestOpenMPArrayHandleMultiplexer.cxx
  UnitTestOpenMPBitField.cxx
  UnitTestOpenMPCellLocatorRectilinearGrid.cxx
  UnitTestOpenMPCellLocatorTwoLevel.cxx
  UnitTestOpenMPCellLocatorUniformGrid.cxx
  UnitTestOpenMPColorTable.cxx
  UnitTestOpenMPComputeRange.cxx
  UnitTestOpenMPDataSetExplicit.cxx
  UnitTestOpenMPDataSetSingleType.cxx
  UnitTestOpenMPDeviceAdapter.cxx
  UnitTestOpenMPImplicitFunction.cxx
  UnitTestOpenMPPointLocatorSparseGrid.cxx
  UnitTestOpenMPRuntimeDeviceConfiguration.cxx
  )

if (NOT VTKm_NO_DEPRECATED_VIRTUAL)
  set(unit_tests ${unit_tests}
    UnitTestOpenMPVirtualObjectHandle.cxx
    )
endif()

vtkm_unit_tests(SOURCES ${unit_tests}
                LABEL "OPENMP"
                DEFINES VTKM_NO_ERROR_ON_MIXED_CUDA_CXX_TAG
                LIBRARIES vtkm_worklet)

#We need to have all OpenMP tests run serially as they
#will uses all the system cores, and we will cause a N*N thread
#explosion which causes the tests to run slower than when run
#serially
foreach (test ${unit_tests})
  if (NOT TEST "${test}")
    continue ()
  endif ()
  get_filename_component(tname ${test} NAME_WE)
  set_tests_properties(${tname} PROPERTIES
      RUN_SERIAL True
  )
endforeach()
