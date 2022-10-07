cmake_minimum_required(VERSION 3.13.4)

project(BigImage)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://slicerscope.github.io")
set(EXTENSION_CATEGORY "BigImage")
set(EXTENSION_CONTRIBUTORS "Yi Gao (Shenzhen University), Andras Lasso (PerkLab, Queen's University)")
set(EXTENSION_DESCRIPTION "BigImage is the extension for viewing and analyzing big (>1GB) microscopy and other images.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/gaoyi/SlicerBigImage/main/SlicerBigImage.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/gaoyi/SlicerBigImage/main/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BigImageViewer)
add_subdirectory(ColorDecomposition)
add_subdirectory(NgffImageIO)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
