cmake_minimum_required(VERSION 3.13.4)

project(SlicerScope)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://slicerscope.github.io")
set(EXTENSION_CATEGORY "SlicerScope")
set(EXTENSION_CONTRIBUTORS "Yi Gao (Shenzhen University)")
set(EXTENSION_DESCRIPTION "SlicerScope is the extension hosting 3D Slicer based GB-level microscopic image viewing and computing.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/SlicerScope/SlicerScope/main/SlicerScope.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/SlicerScope/SlicerScope/main/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SlicerScopeViewer)
add_subdirectory(ColorDecomposition)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
