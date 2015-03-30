project(SPHARM_Mech)

cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()

find_package(VTK COMPONENTS
  vtkCommonCore
  vtkFiltersCore
  vtkInfovisCore
  vtkInteractionStyle
  vtkRenderingFreeTypeOpenGL
  vtkViewsQt
)
include(${VTK_USE_FILE})

# Use the include path and library for Qt that is used by VTK.
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  /Users/khairyk/VPCZ1_kk_share/c_code/shape_tools
  /Users/khairyk/VPCZ1_kk_share/c_code/buildnlopt/include
)

# Set your files and resources here
set( Srcs main.cpp simgast.cpp )

set( Hdrs simgast.h)

set( MOC_Hdrs simgast.h )

set( UIs simgast.ui )


# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# We have ui files, this will bring in the macro: qt5_wrap_ui
  find_package(Qt5Widgets REQUIRED QUIET)
  qt5_wrap_ui(UI_Srcs ${UIs})
  qt5_add_resources(QRC_Srcs ${QRCs} )

  source_group("Resources" FILES
    ${UIs}
    ${QRCs}
    ${EXE_ICON} # Not present
  )

  source_group("Generated" FILES
    ${UI_Srcs}
    ${MOC_Srcs}
    ${QRC_Srcs}
    ${QRC_Srcs}
  )

  add_executable(spharm_mech MACOSX_BUNDLE
    ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${QRC_Srcs})
  qt5_use_modules(spharm_mech Core Gui Widgets)
  target_link_libraries(spharm_mech ${VTK_LIBRARIES} /Users/khairyk/VPCZ1_kk_share/c_code/buildnlopt/lib/libnlopt.a)

