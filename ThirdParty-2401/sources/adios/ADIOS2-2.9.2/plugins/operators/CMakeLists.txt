#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_Sodium)
  add_library(EncryptionOperator MODULE
    EncryptionOperator.cpp
  )
  target_link_libraries(EncryptionOperator adios2_core sodium)
  install(TARGETS EncryptionOperator  EXPORT adios2Exports
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_core-runtime
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_core-libraries NAMELINK_COMPONENT adios2_core-development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT adios2_core-development
)
endif()
