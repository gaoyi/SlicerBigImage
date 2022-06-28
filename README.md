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

## Build and Installation

### Linux

You will need to **manully** install the openslide library to your OS. This can be done with, e.g., 
```
sudo apt install libopenslide0
```

This extension also need the following python modules (but they will be installed automatically):
* scikit-image 
* openslide-python

These modules, if not installed in your Slicer already, will be automatically installed for the first time the extension is run.

You may need to restart Slicer after the installation.

### Windows

You will need to **manully** install the openslide library to your OS. On Windows, download the pre-built binary library from [OpenSlide download site](https://openslide.org/download/) under the "Windows Binary" section. Pick the latest release (the 2017-11-22 one as of 20220628).

Note that OpenSlide winows binary library is built using MingW. And it's okay if you are building Slicer using other compiler (visual studio for example).

You can un-zip the openslide library to your system library directory. Or you can un-zip it to anywhere. In the later case, you will need to add the following lines in the SlicerScopeViewer.py, after the `import skimage.transform` on line 33:

```
import os
os.add_dll_directory("path-to-openslide/openslide-win64-20171122/bin")
```

This extension also need the following python modules (but they will be installed automatically):
* scikit-image 
* openslide-python

These modules, if not installed in your Slicer already, will be automatically installed for the first time the extension is run.

You may need to restart Slicer after the installation.

## Modules
This extension currently has two modules:

* SlicerScopeViewer: This module loads and views the WSI images.
* ColorDecomposition: This module performs the color/staining decomposition of the image.

Their detailed usages are listed below. More modules including nucleus segmentaiton and others will be uploaded soon.

## Usage

### Example data
Example large scale wholse slide image can be downloaded at, e.g., the
[OpenSlide
website](https://openslide.cs.cmu.edu/download/openslide-testdata/Aperio/CMU-1.svs
"Brightfield WSI")

### Large whole slide image viewing
Swith to the "SlicerScopeViewer" module in the SlicerScope category. As shown in the module panel, select the WSI file to view in the "Select WSI" box. Then click "Load WSI" below.

![image](https://user-images.githubusercontent.com/920557/174559913-77ccaee3-5063-4fa5-b562-dd1ad3b24236.png)

One can then view different region/scale of the image, using mouse dragging and mosue wheeling, as shown below:

![image](https://user-images.githubusercontent.com/89077084/174545844-83a5f601-32ca-4d88-b328-b3a0cba0e922.png)
![image](https://user-images.githubusercontent.com/89077084/174545870-063ae0a8-2e3d-49bd-8d61-08ca19c5dbb6.png)

### Staining decomposition
Histopathology images are often stained using different dyes. When a WSI is stained using multiple dyes, the different stains can be computationally separated using the module "ColorDecomposition", whose panel is shown below.

![image](https://user-images.githubusercontent.com/920557/174555656-1e227e15-2110-4bf1-8b9d-74b1fdddc823.png)

There are various options for color decomposition. This particular WSI is stained with H-E. Its original appearance is:
![image](https://user-images.githubusercontent.com/920557/174556082-81738b77-87f5-4111-bf31-bbca18501501.png)

After decomposition, the hematoxylin content is shown in gray-scale as:
![image](https://user-images.githubusercontent.com/920557/174556323-eb064126-c40b-48a2-95bd-f4a0c77d60b3.png)
where the dark regions corresponding to the high hematoxylin content.

If the eosin chanel is wanted, one can switch the output chanel in the module panel to the 2nd chanel, and the result will be like:
![image](https://user-images.githubusercontent.com/920557/174556464-e4e1d6d0-f1c3-4222-ad68-f490a520ae98.png)


## Citation

If you find this extension helpful please cite this paper:

Xiaxia Yu, Bingshuai Zhao, Haofan Huang, Mu Tian, Sai Zhang, Hongping Song, Zengshan Li, Kun Huang, Yi Gao, "An Open Source Platform for Computational Histopathology," in IEEE Access, vol. 9, pp. 73651-73661, 2021, doi: 10.1109/ACCESS.2021.3080429.


