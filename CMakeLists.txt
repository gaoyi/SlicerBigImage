cmake_minimum_required(VERSION 3.13.4)

project(SlicerScope)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/SlicerScope")
set(EXTENSION_CATEGORY "SlicerScope")
set(EXTENSION_CONTRIBUTORS "Yi Gao (Shenzhen University)")
set(EXTENSION_DESCRIPTION "SlicerScope is the extension hosting 3D Slicer based GB-level microscopic image viewing and computing.")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/SlicerScope.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/SlicerScope/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SlicerScopeViewer)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
