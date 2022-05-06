# SlicerScope

## Introduction

SlicerScope is an extension of 3D Slicer. It is used for large scale
microscopic image viewing and computing.

Computational histopathology is a fast emerging field which converts
the traditional glass slide based department to a new examination
platform. Such a paradigm shift also brings the in silico computation
to the field. Much research have been presented in the past decades on
the algorithm development for pathology image analysis. On the other
hand, a comprehensive software platform with advanced visualization
and computation capability, large developer community, flexible plugin
mechanism, and friendly transnational license, would be extremely
beneficial for the entire community.

SlicerScope is an open platform for whole slide histopathology image
computing based on the highly successful 3D Slicer.

In addition to viewing the giga-pixel whole slide image viewing,
currently, the extension also offers several specific analytical
modules for qualitative presentation, nucleus level analysis, tissue
scale computation, and 3D pathology. Thanks to the openess of Slicer,
the SlicerScope extension could also be further extended by you.

## Build
The module also needs the openslide library to be installed on your
computer. Under linux (apt system), this could be done by simply running:

sudo apt install libopenslide

The module in this extension will install the python wrapper of
openslide through pip. This is done automatically the first time the
module is run. But you may need to restart Slicer to use it.


## Citation

If you find this extension helpful please cite this paper:

Xiaxia Yu, Bingshuai Zhao, Haofan Huang, Mu Tian, Sai Zhang, Hongping Song, Zengshan Li, Kun Huang, Yi Gao, "An Open Source Platform for Computational Histopathology," in IEEE Access, vol. 9, pp. 73651-73661, 2021, doi: 10.1109/ACCESS.2021.3080429.


