import logging
import os

import vtk

import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin


#
# NgffImageIO
#

class NgffImageIO(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "NgffImageIO"
        self.parent.categories = ["Utilities"]
        # don't show this module - it is only for registering a reader
        self.parent.hidden = True
        self.parent.dependencies = []
        self.parent.contributors = ["Andras Lasso (PerkLab, Queen's University)"]
        self.parent.helpText = """
Read and write images in Zarr file format.
Array data is stored in the `image` array in the `xyz` axis-aligned physical coordinate system.
Transformation between this coordinate system (xyz) and the anatomical coordinate system (lps)
is stored in the `xyzToPhysicalTransform` attribute.
"""
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = """
"""

#
# NgffImageIOLogic
#

class NgffImageIOLogic(ScriptedLoadableModuleLogic):
    """
    This class contains utility functions for reading and writing volume nodes in zarr file format.
    """

    def __init__(self):
        """
        Called when the logic class is instantiated. Can be used for initializing member variables.
        """
        ScriptedLoadableModuleLogic.__init__(self)

    @staticmethod
    def installRequiredPythonPackages():
        try:
            import zarr
        except ImportError:
            slicer.util.pip_install('zarr')
        try:
            import dask
        except ImportError:
            slicer.util.pip_install('dask')
        try:
            import xarray
        except ImportError:
            slicer.util.pip_install('xarray')

    @staticmethod
    def xarrayFromVolume(volumeNode) -> "xr.DataArray":
        """Convert a volume node to an xarray.DataArray.
        Since image axes may be rotated in physical space but xarray accessors do not
        support rotated axes, the image in the xarray is defined in the "xyz" space,
        which is voxel space scaled with the image spacing.
        `xyzToPhysicalTransform` attribute stores a 4x4 homogeneous transformation matrix
        to transform between xyz to physical (LPS) coordinate systems.
        Spacing metadata is preserved in the xarray's coords.
        Dims are labeled as `x`, `y`, `z`, `t`, and `c`.
        This interface is and behavior is experimental and is subject to possible
        future changes."""

        NgffImageIOLogic.installRequiredPythonPackages()

        import xarray as xr
        import numpy as np
        array_view = slicer.util.arrayFromVolume(volumeNode)
        spacing = volumeNode.GetSpacing()
        origin_ras = volumeNode.GetOrigin()
        origin_lps = [-origin_ras[0], -origin_ras[1], origin_ras[2]]
        size = volumeNode.GetImageData().GetDimensions()
        image_dimension = 3
        image_dims = ("x", "y", "z", "t")
        coords = {}
        # When we export an image, xyz origin is always set to (0,0,0), but after processing
        # (such as resampling or extracting a subset of the data), the origin may change.
        origin_xyz = [0.0, 0.0, 0.0]
        for index, dim in enumerate(image_dims[:image_dimension]):
            coords[dim] = np.linspace(
                origin_xyz[index],
                origin_xyz[index] + (size[index] - 1) * spacing[index],
                size[index],
                dtype=np.float64,
            )
        dims = list(reversed(image_dims[:image_dimension]))
        components = volumeNode.GetImageData().GetNumberOfScalarComponents()
        if components > 1:
            dims.append("c")
            coords["c"] = np.arange(components, dtype=np.uint32)
        ijkToRasMatrixVtk = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(ijkToRasMatrixVtk)
        ijkToRasMatrix = slicer.util.arrayFromVTKMatrix(ijkToRasMatrixVtk)
        ijkToLpsMatrix = np.dot(ijkToRasMatrix, np.diag([-1.0, -1.0, 1.0, 1.0]))
        xyzToIjkMatrix = np.diag([1.0/spacing[0], 1.0/spacing[1], 1.0/spacing[2], 1.0])
        xyzToLpsMatrix = np.dot(ijkToLpsMatrix, xyzToIjkMatrix)
        print(f"xyzToPhysical={xyzToLpsMatrix}")
        attrs = {"xyzToPhysicalTransform": np.flip(xyzToLpsMatrix)}
        for attributeName in volumeNode.GetAttributeNames():
            attrs[key] = volumeNode.GetAttribute(attributeName)
        data_array = xr.DataArray(array_view, dims=dims, coords=coords, attrs=attrs)
        return data_array

    @staticmethod
    def volumeFromXarray(data_array: "xr.DataArray"):
        """Convert an xarray.DataArray to a MRML volume node.
        """

        NgffImageIOLogic.installRequiredPythonPackages()

        import numpy as np
        import builtins
        if not {"t", "z", "y", "x", "c"}.issuperset(data_array.dims):
            raise ValueError('Unsupported dims, supported dims: "t", "z", "y", "x", "c".')
        image_dims = list({"t", "z", "y", "x"}.intersection(set(data_array.dims)))
        image_dims.sort(reverse=True)
        image_dimension = len(image_dims)
        ordered_dims = ("t", "z", "y", "x")[-image_dimension:]
        is_vector = "c" in data_array.dims
        if is_vector:
            ordered_dims = ordered_dims + ("c",)
        values = data_array.values
        if ordered_dims != data_array.dims:
            dest = list(builtins.range(len(ordered_dims)))
            source = dest.copy()
            for ii in builtins.range(len(ordered_dims)):
                source[ii] = data_array.dims.index(ordered_dims[ii])
            values = np.moveaxis(values, source, dest).copy()
        volumeNode = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLScalarVolumeNode' if not is_vector else 'vtkMRMLVectorVolumeNode')
        slicer.util.updateVolumeFromArray(volumeNode, values) # is_vector)
        origin_xyz = [0.0] * image_dimension
        spacing = [1.0] * image_dimension
        for index, dim in enumerate(image_dims):
            coords = data_array.coords[dim]
            if coords.shape[0] > 1:
                origin_xyz[index] = float(coords[0])
                print(f'origin[{dim}] = {origin_xyz[index]}')
                spacing[index] = float(coords[-1] - coords[0]) / float(len(coords)-1)
                print(f'coords[{dim}] spacing = ({coords[-1]} - {coords[0]}) / {len(coords)-1} = {spacing[index]}')
        spacing.reverse()
        origin_xyz.reverse()
        ijkToXyz = np.diag([spacing[0], spacing[1], spacing[2], 1.0])
        ijkToXyz[:,3] = [-origin_xyz[0], -origin_xyz[1], origin_xyz[2], 1.0]  # TODO: it is not clear why first two components need sign inversion
        print(f"ijkToXyz={ijkToXyz}")
        if "xyzToPhysicalTransform" in data_array.attrs:
            xyzToPhysical = np.flip(data_array.attrs["xyzToPhysicalTransform"])
        else:
            xyzToPhysical = np.identity(4)
        ijkToLps = np.dot(xyzToPhysical, ijkToXyz)
        ijkToRas = np.dot(ijkToLps, np.diag([-1.0, -1.0, 1.0, 1.0]))
        print(f"xyzToPhysical={xyzToPhysical}")
        print(f"ijkToRas={ijkToRas}")
        volumeNode.SetIJKToRASMatrix(slicer.util.vtkMatrixFromArray(ijkToRas))
        ignore_keys = set(["xyzToPhysicalTransform"])
        for key in data_array.attrs:
            if not key in ignore_keys:
                volumeNode.SetAttribute(key, data_array.attrs[key])
        return volumeNode

#
# Reader plugin
#

class NgffImageIOFileReader:

    def __init__(self, parent):
        self.parent = parent

    def description(self):
        return 'Zarr image'

    def fileType(self):
        return 'ZarrImage'

    def extensions(self):
        return ['Zarr image (*.zarr)']

    def canLoadFile(self, filePath):
        return True

    def load(self, properties):
        try:
            filePath = properties['fileName']

            # Get node base name from filename
            if 'name' in properties.keys():
                baseName = properties['name']
            else:
                baseName = os.path.splitext(os.path.basename(filePath))[0]
                baseName = slicer.mrmlScene.GenerateUniqueName(baseName)

            # Read zarr image file
            import xarray as xr
            import zarr
            ds = xr.open_zarr(zarr.ZipStore(filePath))

            # Convert to volume
            volumeNode = NgffImageIOLogic.volumeFromXarray(ds.image)
            volumeNode.SetName(baseName)

        except Exception as e:
            logging.error('Failed to load file: ' + str(e))
            import traceback
            traceback.print_exc()
            return False

        # Show volume
        selectionNode = slicer.app.applicationLogic().GetSelectionNode()
        selectionNode.SetActiveVolumeID(volumeNode.GetID())
        slicer.app.applicationLogic().PropagateVolumeSelection()

        self.parent.loadedNodes = [volumeNode.GetID()]
        return True

#
# Writer plugin
#

class NgffImageIOFileWriter:

    def __init__(self, parent):
        self.parent = parent

    def description(self):
        return 'Zarr image'

    def fileType(self):
        return 'ZarrImage'

    def extensions(self, obj):
        return ['Zarr image (*.zarr)']

    def canWriteObject(self, obj):
        return bool(obj.IsA("vtkMRMLVolumeNode"))

    def write(self, properties):
        try:

            # Get node
            node = slicer.mrmlScene.GetNodeByID(properties["nodeID"])

            # Write node content to file
            filePath = properties['fileName']

            # Convert to xarray
            da = NgffImageIOLogic.xarrayFromVolume(node)
            ds = da.to_dataset(name='image')

            # Save as zarr
            import zarr
            with zarr.ZipStore(filePath) as store:
                ds.to_zarr(store)

        except Exception as e:
            logging.error('Failed to write file: ' + str(e))
            import traceback
            traceback.print_exc()
            return False

        self.parent.writtenNodes = [node.GetID()]
        return True


#
# NgffImageIOTest
#

class NgffImageIOTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """ Do whatever is needed to reset the state - typically a scene clear will be enough.
        """
        slicer.mrmlScene.Clear()

    def runTest(self):
        """Run as few or as many tests as needed here.
        """
        self.setUp()
        self.test_NgffImageIO1()

    def test_NgffImageIO1(self):
        """ Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("Starting the test")

        import SampleData
        inputVolume = SampleData.downloadSample('MRHead')

        self.delayDisplay('Save volume in zarr format')

        # Convert to xarray
        da = NgffImageIOLogic.xarrayFromVolume(inputVolume)
        print(da)
        self.assertEqual(len(da.x), 256)
        self.assertEqual(len(da.y), 256)
        self.assertEqual(len(da.z), 130)
        self.assertTrue(hasattr(da, 'xyzToPhysicalTransform'))
        # Save as zarr
        tempDir = slicer.util.tempDirectory()
        print(f"Write image in zarr format to {tempDir}")
        ds = da.to_dataset(name='image').chunk({'x': 20, 'y': -1})
        zs = ds.to_zarr(tempDir, mode='w')

        self.delayDisplay('Load volume from zarr')

        # Load zarr into xarray data set
        import xarray as xr
        ds = xr.open_dataset(tempDir, engine='zarr')
        self.assertEqual(len(ds.image.x), 256)
        self.assertEqual(len(ds.image.y), 256)
        self.assertEqual(len(ds.image.z), 130)
        self.assertTrue(hasattr(ds.image, 'xyzToPhysicalTransform'))
        # Convert to volume
        volumeNode = NgffImageIOLogic.volumeFromXarray(ds.image)
        # Display volume
        slicer.util.setSliceViewerLayers(volumeNode)

        self.delayDisplay('Load volume region from zarr')
        # This tests that xarray can correctly retrieve a region of the data at the specified resolution.
        # Efficient retrieval of a region of an image is a major feature of xarray/zarr, which allows
        # handling very large data sets.

        # Load zarr into xarray data set
        ds = xr.open_dataset(tempDir, engine='zarr', chunks={})
        import numpy as np
        dsPart = ds.sel(x=np.arange(50, 180, 2), y=np.arange(40, 160, 2), z=np.arange(50, 60, 2), method='nearest').image
        # Convert to volume
        volumePart = NgffImageIOLogic.volumeFromXarray(dsPart)
        # Display volume
        volumePart.CreateDefaultDisplayNodes()
        volumePart.GetScalarVolumeDisplayNode().SetAndObserveColorNodeID('vtkMRMLColorTableNodeRed')
        slicer.util.setSliceViewerLayers(foreground=volumePart, foregroundOpacity=0.5)

        self.delayDisplay('Clean up')

        # Remove temporary folder
        import shutil
        shutil.rmtree(tempDir, True)
