FILE(GLOB TComb_Sources RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "*.c"
  "*.cpp"
  "*.hpp"
  "*.h"

  "avs/*.h"
)

IF( MSVC OR MINGW )
    # Export definitions in general are not needed on x64 and only cause warnings,
    # unfortunately we still must need a .def file for some COM functions.
    # NO C interface for this plugin
    # if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    #  LIST(APPEND TComb_Sources "TComb64.def")
    # else()
    #  LIST(APPEND TComb_Sources "TComb.def")
    # endif() 
ENDIF()

IF( MSVC_IDE )
    # Ninja, unfortunately, seems to have some issues with using rc.exe
    LIST(APPEND MedianBlur2_Sources "TComb.rc")
ENDIF()
