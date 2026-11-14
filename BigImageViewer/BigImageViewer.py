import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin

import logging


import SimpleITK as sitk
import sys
import numpy

import math

try:
    from tiffslide import TiffSlide
except ImportError:
    slicer.util.pip_install("tiffslide")
    from tiffslide import TiffSlide

#
# BigImageViewer
#

class BigImageViewer(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "BigImageViewer"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["BigImage"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["Yi Gao (Shenzhen Univ.)"]  # TODO: replace with "Firstname Lastname (Organization)"
    # TODO: update with short description of the module and a link to online module documentation
    self.parent.helpText = """
    BigImageViewer for viewing large image whose whole content is not able to be loaded into memory.
    """
    self.parent.acknowledgementText = """
    This file was developed by Yi Gao.
""" # replace with organization, grant and thanks.

#
# Register sample data sets in Sample Data module
#

# def registerSampleData():
#   """
#   Add data sets to Sample Data module.
#   """
#   # It is always recommended to provide sample data for users to make it easy to try the module,
#   # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

#   import SampleData
#   iconsPath = os.path.join(os.path.dirname(__file__), 'Resources/Icons')

#   # To ensure that the source code repository remains small (can be downloaded and installed quickly)
#   # it is recommended to store data sets that are larger than a few MB in a Github release.

#   # BigImageViewer1
#   SampleData.SampleDataLogic.registerCustomSampleDataSource(
#     # Category and sample name displayed in Sample Data module
#     category='BigImageViewer',
#     sampleName='BigImageViewer1',
#     # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
#     # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
#     thumbnailFileName=os.path.join(iconsPath, 'BigImageViewer1.png'),
#     # Download URL and target file name
#     uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
#     fileNames='BigImageViewer1.nrrd',
#     # Checksum to ensure file integrity. Can be computed by this command:
#     #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
#     checksums = 'SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95',
#     # This node name will be used when the data set is loaded
#     nodeNames='BigImageViewer1'
#   )

#   # BigImageViewer2
#   SampleData.SampleDataLogic.registerCustomSampleDataSource(
#     # Category and sample name displayed in Sample Data module
#     category='BigImageViewer',
#     sampleName='BigImageViewer2',
#     thumbnailFileName=os.path.join(iconsPath, 'BigImageViewer2.png'),
#     # Download URL and target file name
#     uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
#     fileNames='BigImageViewer2.nrrd',
#     checksums = 'SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97',
#     # This node name will be used when the data set is loaded
#     nodeNames='BigImageViewer2'
#   )

#
# BigImageViewerWidget
#

class BigImageViewerWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)  # needed for parameter node observation
    self.logic = None
    self._parameterNode = None
    self._updatingGUIFromParameterNode = False


    # self.kerasModelSegmentNuclei = keras.models.load_model("/tmp/unet_nucleus_20200120_at20200126.hdf5")
    # self.kerasModelSegmentColonGland = keras.models.load_model("/tmp/unet_glas_20200713.hdf5")


    # set to Red view only
    lm = slicer.app.layoutManager()
    lm.setLayout(6)

    #--------------------------------------------------------------------------------
    # Viewer variables
    self.topLeftX = 0
    self.topLeftY = 0

    self.patchSizeX = 0
    self.patchSizeY = 0

    self.MPP = 0.25 # um/pixel
    self.SliceThickness = 0.005 #mm

    self.BigRGBAImagePathname = ""

    self.BigRGBAImageNumberOfLevels = 0
    self.BigRGBAImageLevelToLoad = 0

    self.leftMouseButtonPos = (0, 0)

    self.WSISizesXAtAllLevels = []
    self.WSISizesYAtAllLevels = []

    self.H5FilePathname = ""
    self.h5FileLoaded = False

    self._zoomAnchor = None

    self.patchVolumeNode = None
    self.patchGrayVolumeNode = None
    self.patchLabelVolumeNode = None

    self.patchWsiOriginX = 0.0
    self.patchWsiOriginY = 0.0
    self.patchWsiScaleX = 1.0
    self.patchWsiScaleY = 1.0
    self.showWholeSlideOnLoad = False
    self.initialMagnification = 1.0


    # Interator
    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()

    # disable scroller in the Red View
    # redWidget.sliceController().setVisible(False)

    interactor = redView.interactorStyle().GetInteractor()
    interactor.AddObserver(vtk.vtkCommand.RightButtonPressEvent, self.onRightButtonPressed)

    interactor.AddObserver(vtk.vtkCommand.LeftButtonPressEvent, self.onLeftButtonPressed)
    interactor.AddObserver(vtk.vtkCommand.LeftButtonReleaseEvent, self.onLeftButtonReleased)

    interactor.AddObserver(vtk.vtkCommand.MouseWheelForwardEvent, self.onMouseWheelForwardEvent, 1.0)
    interactor.AddObserver(vtk.vtkCommand.MouseWheelBackwardEvent, self.onMouseWheelBackwardEvent, 1.0)

  #
  # Customized mouse right button pressed event
  #
  def onRightButtonPressed(self, obj, event=None):
    print ('onRightButtonPressed............................', event)

    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    style = redView.interactorStyle()
    interactor = style.GetInteractor()
    #print(interactor.GetEventPosition())
    # Do something here



  def onLeftButtonPressed(self, obj, event=None):
    print ('onLeftButtonPressed............................', event)

    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    style = redView.interactorStyle()
    interactor = style.GetInteractor()

    self.leftMouseButtonPos = interactor.GetEventPosition()
    #print(self.leftMouseButtonPos, type(self.leftMouseButtonPos), type(self.leftMouseButtonPos[0]))
    # Do something here


  def onLeftButtonReleased(self, obj, event=None):
    print ('onLeftButtonReleased............................', event)

    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    style = redView.interactorStyle()
    interactor = style.GetInteractor()
    #print(interactor.GetEventPosition())

    newPos = interactor.GetEventPosition()

    self.ui.topLeftXSliderWidget.setValue(self.ui.topLeftXSliderWidget.value + 10.0*(self.leftMouseButtonPos[0] - newPos[0]))
    self.ui.topLeftYSliderWidget.setValue(self.ui.topLeftYSliderWidget.value + 10.0*(newPos[1] - self.leftMouseButtonPos[1])) # x and y are different so the orders are different

    #self.ui.topLeftYSliderWidget.update()



  def onMouseWheelForwardEvent(self, obj, event=None):
    if self.ui.ObjectiveMagnificationSlicerWidget.isMaximized():
      return
    if not self.prepareZoomAtMouse():
      return

    self.ui.ObjectiveMagnificationSlicerWidget.setValue(self.ui.ObjectiveMagnificationSlicerWidget.value + 0.5*self.ui.ObjectiveMagnificationSlicerWidget.value)



  def onMouseWheelBackwardEvent(self, obj, event=None):
    if self.ui.ObjectiveMagnificationSlicerWidget.isMinimized():
      return
    if not self.prepareZoomAtMouse():
      return

    self.ui.ObjectiveMagnificationSlicerWidget.setValue(self.ui.ObjectiveMagnificationSlicerWidget.value - 0.5*self.ui.ObjectiveMagnificationSlicerWidget.value)

  def prepareZoomAtMouse(self):
    if not self.BigRGBAImagePathname or self.patchSizeX <= 0 or self.patchSizeY <= 0:
      return False
    patchVolumeNode = getattr(self, 'patchVolumeNode', None)
    if patchVolumeNode is None or patchVolumeNode.GetImageData() is None:
      return False
    objectiveMagnificationMax = getattr(self, 'objectiveMagnificationMax', 0)
    if objectiveMagnificationMax <= 0:
      return False

    currentMagnification = float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    if currentMagnification <= 0:
      return False

    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    interactor = redView.interactorStyle().GetInteractor()
    displayPosition = interactor.GetEventPosition()
    xyz = redView.convertDeviceToXYZ([int(displayPosition[0]), int(displayPosition[1])])
    ras = redView.convertXYZToRAS(xyz)

    ijkToRas = vtk.vtkMatrix4x4()
    patchVolumeNode.GetIJKToRASMatrix(ijkToRas)
    rasToIjk = vtk.vtkMatrix4x4()
    vtk.vtkMatrix4x4.Invert(ijkToRas, rasToIjk)
    ijk = rasToIjk.MultiplyPoint([ras[0], ras[1], ras[2], 1.0])

    dimensions = patchVolumeNode.GetImageData().GetDimensions()
    u = ijk[0] / float(dimensions[0])
    v = ijk[1] / float(dimensions[1])
    if not (0.0 <= u <= 1.0 and 0.0 <= v <= 1.0):
      u = 0.5
      v = 0.5

    wsiWidth = float(self.WSISizesXAtAllLevels[0])
    wsiHeight = float(self.WSISizesYAtAllLevels[0])
    anchorX = self.patchWsiOriginX + ijk[0] * self.patchWsiScaleX
    anchorY = self.patchWsiOriginY + ijk[1] * self.patchWsiScaleY

    if not (0.0 <= anchorX <= wsiWidth and 0.0 <= anchorY <= wsiHeight):
      u = 0.5
      v = 0.5
      anchorX = wsiWidth * 0.5
      anchorY = wsiHeight * 0.5

    self._zoomAnchor = (u, v, anchorX, anchorY)
    return True

  def applyZoomAtMouse(self):
    if self._zoomAnchor is None:
      return
    u, v, anchorX, anchorY = self._zoomAnchor
    self._zoomAnchor = None

    magnification = float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    objectiveMagnificationMax = getattr(self, 'objectiveMagnificationMax', 0)
    if magnification <= 0 or objectiveMagnificationMax <= 0:
      return

    level = self.BigRGBAImageLevelToLoad
    downsample = float(self.level_downsamples[level])
    levelMagnification = objectiveMagnificationMax / downsample
    ratio = levelMagnification / magnification
    spanX = round(self.patchSizeX * ratio) * downsample
    spanY = round(self.patchSizeY * ratio) * downsample

    topLeftX = anchorX - u * spanX
    topLeftY = anchorY - v * spanY
    topLeftX = round(topLeftX / downsample) * downsample
    topLeftY = round(topLeftY / downsample) * downsample

    xSlider = self.ui.topLeftXSliderWidget
    ySlider = self.ui.topLeftYSliderWidget
    topLeftX = min(max(topLeftX, xSlider.minimum), xSlider.maximum)
    topLeftY = min(max(topLeftY, ySlider.minimum), ySlider.maximum)

    xSlider.blockSignals(True)
    ySlider.blockSignals(True)
    xSlider.value = topLeftX
    ySlider.value = topLeftY
    xSlider.blockSignals(False)
    ySlider.blockSignals(False)

  def getPatchNodeName(self, suffix):
    base = os.path.splitext(os.path.basename(self.BigRGBAImagePathname))[0]
    if not base:
      base = "BigImage"
    return "%s_%s" % (base, suffix)

  def resetPatchNodeReferences(self):
    self.patchVolumeNode = None
    self.patchGrayVolumeNode = None
    self.patchLabelVolumeNode = None
    self._zoomAnchor = None
    self.patchWsiOriginX = 0.0
    self.patchWsiOriginY = 0.0
    self.patchWsiScaleX = 1.0
    self.patchWsiScaleY = 1.0
    self.showWholeSlideOnLoad = False
    hdf5File = getattr(self, 'hdf5File', None)
    if hdf5File is not None:
      try:
        hdf5File.close()
      except Exception:
        pass
      self.hdf5File = None
    self.h5FileLoaded = False

  def removePatchNodes(self):
    sh = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
    for attr in ("patchVolumeNode", "patchGrayVolumeNode", "patchLabelVolumeNode"):
      node = getattr(self, attr, None)
      if node is None:
        continue
      if sh is not None:
        itemID = sh.GetItemByDataNode(node)
        if itemID:
          sh.RemoveItem(itemID)
      displayNode = node.GetDisplayNode()
      displayNodeID = displayNode.GetID() if displayNode else None
      slicer.mrmlScene.RemoveNode(node)
      if displayNodeID:
        displayNode = slicer.mrmlScene.GetNodeByID(displayNodeID)
        if displayNode is not None:
          slicer.mrmlScene.RemoveNode(displayNode)
    self.resetPatchNodeReferences()





  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer).
    # Additional widgets can be instantiated manually and added to self.layout.
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/BigImageViewer.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create logic class. Logic implements all computations that should be possible to run
    # in batch mode, without a graphical user interface.
    self.logic = BigImageViewerLogic()

    self.wsiLevelToLoad = 0

    #--------------------------------------------------------------------------------
    # connections
    self.ui.loadWSIMetaInfoButton.connect('clicked(bool)', self.onLoadWSIMetaInfoButton)

    #self.loadBigRGBAImageButton.connect('clicked(bool)', self.onLoadBigRGBAImageButton)
    self.ui.loadH5FileButton.connect('clicked(bool)', self.onLoadH5FileButton)

    self.ui.topLeftXSliderWidget.connect("valueChanged(double)", self.loadPatchFromBigRGBAImage)
    self.ui.topLeftYSliderWidget.connect("valueChanged(double)", self.loadPatchFromBigRGBAImage)

    self.ui.ObjectiveMagnificationSlicerWidget.connect("valueChanged(double)", self.onWSILevelChanged)

    self.ui.segmentNucleiButton.connect('clicked(bool)', self.onSegmentNucleiButton)
    self.ui.decomposeStainButton.connect('clicked(bool)', self.onDecomposeStainButton)
    self.ui.detectGlandButton.connect('clicked(bool)', self.onDetectGlandButton)
    # connections
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Connections

    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

    # # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # # (in the selected parameter node).
    # self.ui.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    # self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    # self.ui.imageThresholdSliderWidget.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    # self.ui.invertOutputCheckBox.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    # self.ui.invertedOutputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)

    # # Buttons
    # self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)

    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()




    # Add vertical spacer
    self.layout.addStretch(1)

  def setupForScalarPatch(self):
    imageSize = [self.patchSizeX, self.patchSizeY, 1]
    imageSpacing = [self.MPP/1000., self.MPP/1000., self.SliceThickness]
    voxelType = vtk.VTK_UNSIGNED_CHAR

    #--------------------------------------------------------------------------------
    # Allocate space and create node to store the RGB patch
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType, 1)

    # Create volume node
    IJKToRASDirectionMatrix = vtk.vtkMatrix4x4()
    IJKToRASDirectionMatrix.SetElement(0, 0, -1)
    IJKToRASDirectionMatrix.SetElement(1, 1, -1)
    # now input patch is using LPS (which does not make sense for
    # 2D image....... but in CLI for color decomposition, the ITK
    # image is by default LPS. so the output image is not aligned
    # with the color image (ras)

    volumeNode = slicer.vtkMRMLScalarVolumeNode()
    volumeNode.SetName(self.getPatchNodeName("patch_gray"))
    volumeNode.SetSaveWithScene(False)
    volumeNode.SetHideFromEditors(True)
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
    volumeNode.SetAndObserveImageData(imageData)
    slicer.mrmlScene.AddNode(volumeNode)

    # Add volume to scene
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    displayNode.SetSaveWithScene(False)

    slicer.mrmlScene.AddNode(displayNode)
    #colorNode = slicer.util.getNode('Green')
    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())

    # store the loaded patch to widget member
    self.patchGrayVolumeNode = volumeNode



  def setupVolumeNodeToStoreRGBPatch(self):
    imageSize = [self.patchSizeX, self.patchSizeY, 1]
    imageSpacing = [self.MPP/1000., self.MPP/1000., self.SliceThickness]
    voxelType = vtk.VTK_UNSIGNED_CHAR

    #--------------------------------------------------------------------------------
    # Allocate space and create node to store the RGB patch
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType, 3)

    # Create volume node
    IJKToRASDirectionMatrix = vtk.vtkMatrix4x4()
    IJKToRASDirectionMatrix.SetElement(0, 0, -1)
    IJKToRASDirectionMatrix.SetElement(1, 1, -1)
    # now input patch is using LPS (which does not make sense for
    # 2D image....... but in CLI for color decomposition, the ITK
    # image is by default LPS. so the output image is not aligned
    # with the color image (ras)

    volumeNode = slicer.vtkMRMLVectorVolumeNode()
    volumeNode.SetName(self.getPatchNodeName("patch_rgb"))
    volumeNode.SetSaveWithScene(False)
    volumeNode.SetHideFromEditors(True)
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
    volumeNode.SetAndObserveImageData(imageData)
    slicer.mrmlScene.AddNode(volumeNode)

    # Add volume to scene
    displayNode = slicer.vtkMRMLVectorVolumeDisplayNode()
    displayNode.SetSaveWithScene(False)
    slicer.mrmlScene.AddNode(displayNode)

    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())

    # store the loaded patch to widget member
    self.patchVolumeNode = volumeNode

    #--------------------------------------------------------------------------------
    # Set the patch volume to BigRGBAImage loader so it knows where to store
    # patch data

    selectionNode = slicer.app.applicationLogic().GetSelectionNode()

    # This will select the RGB patch to the bacgkround of the Red view
    selectionNode.SetReferenceActiveVolumeID(self.patchVolumeNode.GetID())

    # This will select the label mask to the label of the Red
    # view. Without this line, the label node will be in slicer but
    # not selected. You need to manually pick everytime
    #selectionNode.SetReferenceActiveLabelVolumeID(self.patchLabelVolumeNode.GetID())
    slicer.app.applicationLogic().PropagateVolumeSelection(0)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  def setupForLoadingPatchFromH5File(self):

    if self.ui.H5FileFileNameEditor.currentPath and os.path.isfile(self.ui.H5FileFileNameEditor.currentPath):
      self.ui.H5FilePathname = self.ui.H5FileFileNameEditor.currentPath

      #--------------------------------------------------------------------------------
      # Read some meta info from H5File to populate slicer UI

      # self.H5FileLoader.SetH5FileFileName(self.H5FilePathname)
      # self.H5FileLoader.Initialization()

      # self.H5FileNumberOfLevels = self.H5FileLoader.GetH5FileLevel()
      # self.H5FileSizeX0 = self.H5FileLoader.GetSizeX0()
      # self.H5FileSizeY0 = self.H5FileLoader.GetSizeY0()
      # self.MPP = self.H5FileLoader.GetMPP()

      self.hdf5File = h5py.File(self.H5FilePathname, 'r')

      # #--------------------------------------------------------------------------------
      # # Read UI info to determine the patch size to be extracted
      # lm = slicer.app.layoutManager()
      # redWidget = lm.sliceWidget('Red')
      # redView = redWidget.sliceView()

      # self.patchSizeX = redView.width
      # self.patchSizeY = redView.height

      imageSize = [self.patchSizeX, self.patchSizeY, 1]
      imageSpacing = [self.MPP/1000., self.MPP/1000., self.SliceThickness]
      voxelType = vtk.VTK_UNSIGNED_CHAR

      #--------------------------------------------------------------------------------
      # Allocate space and create node to store the label patch
      labelImageData = vtk.vtkImageData()
      labelImageData.SetDimensions(imageSize)
      labelImageData.AllocateScalars(voxelType, 1)

      IJKToRASDirectionMatrix = vtk.vtkMatrix4x4()
      IJKToRASDirectionMatrix.SetElement(0, 0, -1)
      IJKToRASDirectionMatrix.SetElement(1, 1, -1)
      # now input patch is using LPS (which does not make sense for
      # 2D image....... but in CLI for color decomposition, the ITK
      # image is by default LPS. so the output image is not aligned
      # with the color image (ras)

      labelNode = slicer.vtkMRMLLabelMapVolumeNode()
      labelNode.SetName("currentPatchFromH5File-label")
      labelNode.SetSpacing(imageSpacing)
      labelNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
      labelNode.SetAndObserveImageData(labelImageData)
      slicer.mrmlScene.AddNode(labelNode)

      a = slicer.util.array("currentPatchFromH5File-label")
      a[:] = 0

      # Add volume to scene
      labelDisplayNode = slicer.vtkMRMLLabelMapVolumeDisplayNode()
      slicer.mrmlScene.AddNode(labelDisplayNode)
      #labelColorNode = slicer.util.getNode('Random')
      labelColorNode = slicer.util.getNode('GenericColors')
      labelDisplayNode.SetAndObserveColorNodeID(labelColorNode.GetID())
      labelNode.SetAndObserveDisplayNodeID(labelDisplayNode.GetID())
      labelNode.CreateDefaultStorageNode()

      red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
      red_cn = red_logic.GetSliceCompositeNode()
      red_cn.SetLabelVolumeID(labelNode.GetID())

      # store the loaded label patch to widget member
      self.patchLabelVolumeNode = labelNode


      #--------------------------------------------------------------------------------
      # Set the patch volume to H5File loader so it knows where to store
      # patch data
      #self.H5FileLoader.SetOutputPatchImage(self.patchLabelVolumeNode.GetImageData())

      selectionNode = slicer.app.applicationLogic().GetSelectionNode()

      # This will select the RGB patch to the bacgkround of the Red view
      #selectionNode.SetReferenceActiveVolumeID(self.patchVolumeNode.GetID())

      # This will select the label mask to the label of the Red
      # view. Without this line, the label node will be in slicer but
      # not selected. You need to manually pick everytime
      selectionNode.SetReferenceActiveLabelVolumeID(self.patchLabelVolumeNode.GetID())
      slicer.app.applicationLogic().PropagateVolumeSelection(0)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.removeObservers()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()






  def onWSILevelChanged(self):
    # Change from: import openslide
    import tiffslide

    ds = (self.objectiveMagnificationMax)/(self.ui.ObjectiveMagnificationSlicerWidget.value)
    slide = tiffslide.TiffSlide(self.BigRGBAImagePathname) # Use TiffSlide instead of OpenSlide
    self.BigRGBAImageLevelToLoad = slide.get_best_level_for_downsample(ds)
    slide.close()

    #print(self.BigRGBAImageLevelToLoad)

    self.updateUIWidget()

    self.applyZoomAtMouse()

    self.loadPatchFromBigRGBAImage()

  def enableAndInitUIWidget(self):
    #self.loadBigRGBAImageButton.enabled = True

    #--------------------------------------------------------------------------------
    # Set the max of slider for top-left corner
    for slider in (self.ui.topLeftXSliderWidget, self.ui.topLeftYSliderWidget):
      slider.blockSignals(True)
      slider.value = 0
      slider.maximum = 0
      slider.blockSignals(False)
      slider.enabled = True

    #--------------------------------------------------------------------------------
    # Set the max of number of levels
    magnificationSlider = self.ui.ObjectiveMagnificationSlicerWidget
    magnificationSlider.maximum = self.objectiveMagnificationMax
    magnificationSlider.minimum = self.objectiveMagnificationMin
    magnificationSlider.blockSignals(True)
    magnificationSlider.value = self.initialMagnification
    magnificationSlider.blockSignals(False)
    magnificationSlider.enabled = True

    self.showWholeSlideOnLoad = True
    self.onWSILevelChanged()


  def updateUIWidget(self):
    # self.ui.topLeftXSliderWidget.value = self.topLeftX
    # self.ui.topLeftYSliderWidget.value = self.topLeftY

    magnification = float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    if magnification > 0:
      downsample = self.objectiveMagnificationMax / magnification
      self.ui.topLeftXSliderWidget.maximum = max(0, self.WSISizesXAtAllLevels[0] - self.patchSizeX*downsample - 1)
      self.ui.topLeftYSliderWidget.maximum = max(0, self.WSISizesYAtAllLevels[0] - self.patchSizeY*downsample - 1)

    #self.ui.ObjectiveMagnificationSlicerWidget.maximum = self.BigRGBAImageNumberOfLevels - 1


  def getRegionFromFileAsRGBNumpyArray(self, WSIName, level, topX0, topY0, width, height):
    # Change from: import openslide
    import tiffslide

    slide = tiffslide.TiffSlide(WSIName) # Use TiffSlide instead of OpenSlide
    thisTilePilIm = slide.read_region((topX0, topY0), level, (width, height))
    if thisTilePilIm.mode != "RGB":
        thisTilePilIm = thisTilePilIm.convert("RGB")
    imRGBNdarray = numpy.asarray(thisTilePilIm)
    slide.close()

    return imRGBNdarray

  def getWholeSlidePatch(self):
    from PIL import Image

    level = self.BigRGBAImageLevelToLoad
    downsample = float(self.level_downsamples[level])
    levelWidth = int(self.WSISizesXAtAllLevels[level])
    levelHeight = int(self.WSISizesYAtAllLevels[level])
    q = self.getRegionFromFileAsRGBNumpyArray(self.BigRGBAImagePathname, level, 0, 0, levelWidth, levelHeight)

    wsiWidth = float(self.WSISizesXAtAllLevels[0])
    wsiHeight = float(self.WSISizesYAtAllLevels[0])
    scale = min(float(self.patchSizeX)/wsiWidth, float(self.patchSizeY)/wsiHeight)
    displayWidth = max(1, int(round(wsiWidth*scale)))
    displayHeight = max(1, int(round(wsiHeight*scale)))

    resample = getattr(Image, 'Resampling', Image).LANCZOS
    resized = numpy.asarray(Image.fromarray(q).resize((displayWidth, displayHeight), resample))
    p = numpy.full((self.patchSizeY, self.patchSizeX, 3), 255, dtype=numpy.uint8)
    offsetX = int((self.patchSizeX - displayWidth)/2)
    offsetY = int((self.patchSizeY - displayHeight)/2)
    p[offsetY:offsetY+displayHeight, offsetX:offsetX+displayWidth, :] = resized

    self.patchWsiScaleX = (levelWidth*downsample) / float(displayWidth)
    self.patchWsiScaleY = (levelHeight*downsample) / float(displayHeight)
    self.patchWsiOriginX = -offsetX*self.patchWsiScaleX
    self.patchWsiOriginY = -offsetY*self.patchWsiScaleY

    self.topLeftX = 0
    self.topLeftY = 0
    return p

  def loadPatchFromBigRGBAImage(self):
    import skimage

    #--------------------------------------------------------------------------------
    # Load RGB image
    #
    # The image to view is of the size (self.patchSizeX,
    # self.patchSizeY) on the screen. However, i need to compute how
    # at the current zooming level (self.BigRGBAImageLevelToLoad), how
    # big a region i need to load, and then zoom to (self.patchSizeX,
    # self.patchSizeY)
    #
    #

    #print(self.BigRGBAImageLevelToLoad)
    # print(int(self.ui.topLeftXSliderWidget.value))
    # print(int(self.ui.topLeftYSliderWidget.value))
    # print(int(self.patchSizeX))
    # print(int(self.patchSizeY))


    magnification = float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    atWholeSlideView = (self.initialMagnification > 0
                        and magnification <= self.initialMagnification*(1.0 + 1e-6))

    if self.showWholeSlideOnLoad or atWholeSlideView:
      self.showWholeSlideOnLoad = False
      for slider in (self.ui.topLeftXSliderWidget, self.ui.topLeftYSliderWidget):
        slider.blockSignals(True)
        slider.value = 0
        slider.blockSignals(False)
      p = self.getWholeSlidePatch()
      ratio = 1.0
    else:
      magAtThisLevel = float(self.slideInfo["objectiveMagnification"])/float(self.level_downsamples[self.BigRGBAImageLevelToLoad])
      ratio = magAtThisLevel/float(self.ui.ObjectiveMagnificationSlicerWidget.value)

      sizeToLoadX = round(self.patchSizeX*ratio)
      sizeToLoadY = round(self.patchSizeY*ratio)

      self.topLeftX = int(self.ui.topLeftXSliderWidget.value)
      self.topLeftY = int(self.ui.topLeftYSliderWidget.value)

      q = self.getRegionFromFileAsRGBNumpyArray(self.BigRGBAImagePathname, self.BigRGBAImageLevelToLoad, int(self.ui.topLeftXSliderWidget.value), int(self.ui.topLeftYSliderWidget.value), int(sizeToLoadX), int(sizeToLoadY))

      p = skimage.transform.resize(q, (int(self.patchSizeY), int(self.patchSizeX)), preserve_range=True, mode='wrap', order=0, anti_aliasing=False)

      downsample = float(self.level_downsamples[self.BigRGBAImageLevelToLoad])
      self.patchWsiScaleX = (sizeToLoadX*downsample) / float(self.patchSizeX)
      self.patchWsiScaleY = (sizeToLoadY*downsample) / float(self.patchSizeY)
      self.patchWsiOriginX = int(self.topLeftX / downsample) * downsample
      self.patchWsiOriginY = int(self.topLeftY / downsample) * downsample

    #a = slicer.util.array('currentPatchFromBigRGBAImage')
    volumeNode = self.patchVolumeNode
    a = slicer.util.arrayFromVolume(volumeNode)

    # resizing is needed for python (not for CLI version) coz we have
    # to fill the numpy array here
    a[:] = p

    resolutionNow = self.MPP/1000.0*float(self.slideInfo["objectiveMagnification"])/float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    imageSpacing = [resolutionNow, resolutionNow, self.SliceThickness]

    n = self.patchVolumeNode
    n.SetSpacing(imageSpacing)

    n.GetImageData().Modified()

    displayNode = n.GetDisplayNode()
    if displayNode is not None:
      displayNode.Modified()

    #--------------------------------------------------------------------------------
    # Load label image
    if self.h5FileLoaded:
      group = self.hdf5File[self.ui.h5DatasetOption.currentText]

      dsetName = "data" + str(self.BigRGBAImageLevelToLoad)
      dset = group[dsetName]

      start0 = int(round(self.ui.topLeftYSliderWidget.value/self.level_downsamples[self.BigRGBAImageLevelToLoad]))
      start1 = int(round(self.ui.topLeftXSliderWidget.value/self.level_downsamples[self.BigRGBAImageLevelToLoad]))
      size0 = int(round(self.patchSizeY*ratio))
      size1 = int(round(self.patchSizeX*ratio))

      # print(dsetName)
      # print((start0, start1))
      #print(size1)

      #print("requested size = ", (size0, size1))

      labelArrayInLevel = dset[start0:(start0 + size0), start1:(start1 + size1)]
      #print("max label labelArrayInLevel = " + str(labelArrayInLevel.max()))

      #---------------------------------------------------------------
      # When at the borader, when requested region is to the right or
      # below the read image, the openslide returns an array the same
      # size as the requested size. But h5 only return the intersected
      # region. So the return image size is smaller than the requested
      # size. If not padded first, the subsequent resize will stretch
      # the label image
      if labelArrayInLevel.shape != (size0, size1):
        tmp = numpy.zeros((size0, size1))
        tmp[:labelArrayInLevel.shape[0], :labelArrayInLevel.shape[1]] = labelArrayInLevel
        labelArrayInLevel = tmp
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      labelArray = skimage.transform.resize(labelArrayInLevel, (int(self.patchSizeY), int(self.patchSizeX)), preserve_range=True, mode='wrap', order=0, anti_aliasing=False)
      #print("got size = ", labelArrayInLevel.shape, "resize to = ", labelArray.shape, "ratio = ", labelArrayInLevel.shape[0]/labelArray.shape[0], labelArrayInLevel.shape[1]/labelArray.shape[1])
      #print("max label labelArray = " + str(labelArray.max()))

      #labelArray = dset[int(self.ui.topLeftYSliderWidget.value):int(self.ui.topLeftYSliderWidget.value + self.patchSizeY), int(self.ui.topLeftXSliderWidget.value):int(self.ui.topLeftXSliderWidget.value + self.patchSizeX)]

      # print(self.ui.topLeftXSliderWidget.value, self.ui.topLeftXSliderWidget.value, self.patchSizeX, self.patchSizeY)
      # print(labelArray.max())

      # labelArray = numpy.zeros((self.patchSizeY, self.patchSizeX), dtype='uint8')
      # labelArray[int(self.patchSizeY/3):int(2*self.patchSizeY/3), int(self.patchSizeX/3):int(2*self.patchSizeX/3)] = 1

      n = slicer.util.getNode('currentPatchFromH5File-label')
      a = slicer.util.array('currentPatchFromH5File-label')
      a[:] = labelArray

      n.SetSpacing(imageSpacing)
      n.GetImageData().Modified()

      b = slicer.util.array('currentPatchFromH5File-label')

      self.patchLabelVolumeNode.GetImageData().Modified()

    #--------------------------------------------------------------------------------
    # If need
    if self.ui.extractHematoxylinOnFlyCheckBox.checked:
      parameters = {}
      parameters['inputVolume'] = self.patchVolumeNode.GetID()
      parameters['outputVolume'] = self.patchGrayVolumeNode.GetID()
      slicer.cli.run( slicer.modules.colordecomposition, None, parameters, wait_for_completion=True )



    #--------------------------------------------------------------------------------
    # magnitude = vtk.vtkImageMagnitude()
    # magnitude.SetInputData(self.patchVolumeNode.GetImageData())
    # magnitude.Update()

    # grayNode = slicer.vtkMRMLScalarVolumeNode()
    # grayNode.SetImageDataConnection(magnitude.GetOutputPort())
    # grayNode.SetName("tempGrayNode")

    # labelNode = slicer.vtkMRMLLabelMapVolumeNode()
    # logic = BigViewerModuleLogic()
    # logic.run(grayNode, labelNode, 100)

    # slicer.mrmlScene.AddNode(labelNode)

    # # Add volume to scene
    # displayNode = slicer.vtkMRMLVectorVolumeDisplayNode()

    # slicer.mrmlScene.AddNode(displayNode)
    # colorNode = slicer.util.getNode('FreeSurferLabels')
    # displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    # labelNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    # labelNode.CreateDefaultStorageNode()

    # red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
    # red_cn = red_logic.GetSliceCompositeNode()
    # red_cn.SetLabelVolumeID(labelNode.GetID())


    # magnitude = vtk.vtkImageMagnitude()
    # magnitude.SetInputData(self.patchVolumeNode.GetImageData())
    # magnitude.Update()

    # grayNode = slicer.vtkMRMLScalarVolumeNode()
    # grayNode.SetImageDataConnection(magnitude.GetOutputPort())
    # grayNode.SetName("tempGrayNode")

    # slicer.mrmlScene.AddNode(grayNode)

    # displayNode = slicer.vtkMRMLVectorVolumeDisplayNode()

    # slicer.mrmlScene.AddNode(displayNode)
    # colorNode = slicer.util.getNode('FreeSurferLabels')
    # displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    # grayNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    # grayNode.CreateDefaultStorageNode()

    # red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
    # red_cn = red_logic.GetSliceCompositeNode()
    # red_cn.SetForegroundVolumeID(grayNode.GetID())


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # #bgrdNode.SetIJKToRASDirectionMatrix(fMat)
    # slicer.mrmlScene.AddNode(bgrdNode)
    # bgrdVolID = bgrdNode.GetID()
    # red_cn.SetForegroundVolumeID(fgrdVolID)
    # red_cn.SetBackgroundVolumeID(bgrdVolID)
    # red_cn.SetForegroundOpacity(1)



    # bgrdNode.SetImageDataConnection(magnitude.GetOutputPort())
    # bgrdNode.SetName(bgrdName)
    # #bgrdNode.SetIJKToRASDirectionMatrix(fMat)
    # slicer.mrmlScene.AddNode(bgrdNode)
    # bgrdVolID = bgrdNode.GetID()
    # red_cn.SetForegroundVolumeID(fgrdVolID)
    # red_cn.SetBackgroundVolumeID(bgrdVolID)
    # red_cn.SetForegroundOpacity(1)


#    enableScreenshotsFlag = self.ui.extractHematoxylinOnFlyCheckBox.checked



    # redViewWidget = slicer.app.layoutManager().sliceWidget("Red")
    # redView = redViewWidget.sliceView()
    # redView.

    lm = slicer.app.layoutManager()
    redWidget = lm.sliceWidget('Red')
    redView = redWidget.sliceView()

    redController = redWidget.sliceController()
    redController.fitSliceToBackground()
    #redView.resetFocalPoint()
    redView.forceRender()

    #    r = slicer.app.layoutManager().sliceWidget("Red").sliceController()
#    r.fitSliceToBackground()


    # wti = vtk.vtkWindowToImageFilter()
    # wti.SetInput(redView.renderWindow())
    # wti.Update()
    # v = vtk.vtkImageViewer()
    # v.SetColorWindow(255)
    # v.SetColorLevel(128)
    # v.SetInputConnection(wti.GetOutputPort())
    # v.Render()

    #slicer.app.processEvents()
    #qt.QApplication.processEvents()

  def determinPatchSizeByViewerSize(self):
    #--------------------------------------------------------------------------------
    # Read UI info to determine the patch size to be extracted
    lm = slicer.app.layoutManager()
    redWidget = lm.sliceWidget('Red')
    redView = redWidget.sliceView()

    self.patchSizeX = max(1, int(redView.width/2))
    self.patchSizeY = max(1, int(redView.height/2))

    if not self.WSISizesXAtAllLevels or not self.WSISizesYAtAllLevels:
      return

    wsiWidth = float(self.WSISizesXAtAllLevels[0])
    wsiHeight = float(self.WSISizesYAtAllLevels[0])
    scale = min(self.patchSizeX / wsiWidth, self.patchSizeY / wsiHeight)
    wholeSlideDownsample = 1.0 / scale

    self.initialMagnification = min(self.objectiveMagnificationMax,
                                    self.objectiveMagnificationMax / wholeSlideDownsample)
    self.objectiveMagnificationMin = min(self.objectiveMagnificationMin,
                                         self.initialMagnification)

    return

  def getMetaInfoFromWSI(self):

    self.logic.installPreRequisites()
    self.slideInfo = self.logic.getSlideInfo(self.BigRGBAImagePathname)

    print("returned from logic.getSlideInfo")

    self.BigRGBAImageNumberOfLevels = self.slideInfo["level_count"]
    self.level_downsamples = self.slideInfo["level_downsamples"]
    self.MPP = self.slideInfo["mpp"][0]

    self.WSISizesXAtAllLevels = []
    self.WSISizesYAtAllLevels = []
    for it in range(self.BigRGBAImageNumberOfLevels):
      self.WSISizesXAtAllLevels.append(self.slideInfo["level_dimensions"][it][0])
      self.WSISizesYAtAllLevels.append(self.slideInfo["level_dimensions"][it][1])

    self.objectiveMagnificationMax = self.slideInfo["objectiveMagnification"]
    self.objectiveMagnificationMin = self.slideInfo["objectiveMagnification"]/self.level_downsamples[-1]/3.0 # so to include more

    # print("======== self.objectiveMagnificationMax = ", self.objectiveMagnificationMax)
    # print("======== self.MPP = ", self.MPP)

    # This will trigger the onWSILevelChanged in which the
    # loadPatchFromBigRGBAImage will be called, where the
    # currentPatchGrayChannel node will be filled--- but it's not
    # allocated yet. So this should not be set here.
    #
    # self.ui.ObjectiveMagnificationSlicerWidget.value = self.objectiveMagnificationMax

    # print(self.objectiveMagnificationMax)
    # print(self.objectiveMagnificationMin)

    return

  def onLoadWSIMetaInfoButton(self):

    #--------------------------------------------------------------------------------
    # Read meta infrom from WSI image
    #if self.BigRGBAImageFileNameEditor.currentPath and os.path.isfile(self.ui.BigRGBAImageFileNameEditor.currentPath):
    if not os.path.isfile(self.ui.BigRGBAImageFileNameEditor.currentPath):
      print("File not exist")
      return

    self.logic.installPreRequisites()

    self.BigRGBAImagePathname = self.ui.BigRGBAImageFileNameEditor.currentPath
    self.removePatchNodes()
    self.getMetaInfoFromWSI()
    self.determinPatchSizeByViewerSize()

    self.setupVolumeNodeToStoreRGBPatch()
    self.setupForScalarPatch()

    self.enableAndInitUIWidget()

  # def onLoadBigRGBAImageButton(self):
  #   self.setupVolumeNodeToStoreRGBPatch()
  #   self.loadPatchFromBigRGBAImage()

  #   self.enableAndInitUIWidget()



  def onDecomposeStainButton(self):
    parameters = {}
    parameters['inputVolume'] = self.patchVolumeNode.GetID()
    parameters['outputVolume'] = self.patchGrayVolumeNode.GetID()
    slicer.cli.run( slicer.modules.colordecomposition, None, parameters, wait_for_completion=True )

    volumeNode = self.patchGrayVolumeNode
    slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = self.patchVolumeNode
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)


  def onDetectMitosisButton(self):
    pass

  def onDetectGlandButton(self):
    self.logic.installPreRequisites()
    self.logic.processSegmentGland(self.patchVolumeNode, self.patchGrayVolumeNode, self.kerasModelSegmentColonGland)

    volumeNode = self.patchGrayVolumeNode
    slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = self.patchVolumeNode
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)

    return

  def onSegmentNucleiButton(self):
    #self.logic.installPreRequisites()

    nucleusLabelImageNode = slicer.vtkMRMLLabelMapVolumeNode()
    nucleusLabelImageNode.SetName("nucleusLabelImage")
    slicer.mrmlScene.AddNode( nucleusLabelImageNode )
    #parameters['outputAxisMaskVolume'] =

    volumeNode = self.patchGrayVolumeNode
    #slicer.util.setSliceViewerLayers(foreground=volumeNode)

    parameters = {}
    parameters['inputfileName'] = self.patchVolumeNode.GetID()
    parameters['outputEImage'] = volumeNode.GetID()
    parameters['outputHImage'] = volumeNode.GetID()
    parameters['outputImageName'] = nucleusLabelImageNode.GetID()
    parameters['mpp'] = self.MPP/self.objectiveMagnificationMax*float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    slicer.cli.run( slicer.modules.nucleussegmentation, None, parameters, wait_for_completion=True )


    # volumeNode = slicer.util.getNode('currentPatchGrayChannel') #self.patchGrayVolumeNode
    # slicer.util.setSliceViewerLayers(foreground=volumeNode)
    # volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage') #self.patchVolumeNode
    # slicer.util.setSliceViewerLayers(background = volumeNode1)

    # slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)


    # self.logic.processSegmentNuclei(self.patchVolumeNode, nucleusLabelImageNode.GetID())


    # volumeNode = slicer.util.getNode('currentPatchGrayChannel')
    # slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = self.patchVolumeNode
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.getNode('vtkMRMLSliceNodeRed').SetUseLabelOutline(1)

    # slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)

  def onLoadH5FileButton(self):
    if not self.patchVolumeNode:
      return

    self.setupForLoadingPatchFromH5File()

    self.h5FileLoaded = True

    f = h5py.File(self.H5FilePathname, 'r')

    self.groupNamesInH5 = []

    for key in f.keys():
        if isinstance(f[key], h5py.Group):
          # node is a dataset
          print("Group ", key)
          self.groupNamesInH5.append(key)
          self.h5DatasetOption.addItem(key)
        elif isinstance(f[key], h5py.Dataset):
          #self.groupNamesInH5.append(key)
          print("Dataset ", key)

    f.close()

    self.h5DatasetOption.enabled = True





  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)
    self.resetPatchNodeReferences()

  def onSceneEndClose(self, caller, event):
    """
    Called just after the scene is closed.
    """
    # If this module is shown while the scene is closed then recreate a new parameter node immediately
    if self.parent.isEntered:
      self.initializeParameterNode()

  def initializeParameterNode(self):
    """
    Ensure parameter node exists and observed.
    """
    # Parameter node stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.

    self.setParameterNode(self.logic.getParameterNode())

    # Select default input nodes if nothing is selected yet to save a few clicks for the user
    if not self._parameterNode.GetNodeReference("InputVolume"):
      firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
      if firstVolumeNode:
        self._parameterNode.SetNodeReferenceID("InputVolume", firstVolumeNode.GetID())

  def setParameterNode(self, inputParameterNode):
    """
    Set and observe parameter node.
    Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
    """

    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    # # Update node selectors and sliders
    # self.ui.inputSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputVolume"))
    # self.ui.outputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputVolume"))
    # self.ui.invertedOutputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputVolumeInverse"))
    # self.ui.imageThresholdSliderWidget.value = float(self._parameterNode.GetParameter("Threshold"))
    # self.ui.invertOutputCheckBox.checked = (self._parameterNode.GetParameter("Invert") == "true")

    # # Update buttons states and tooltips
    # if self._parameterNode.GetNodeReference("InputVolume") and self._parameterNode.GetNodeReference("OutputVolume"):
    #   self.ui.applyButton.toolTip = "Compute output volume"
    #   self.ui.applyButton.enabled = True
    # else:
    #   self.ui.applyButton.toolTip = "Select input and output volume nodes"
    #   self.ui.applyButton.enabled = False

    # All the GUI updates are done
    self._updatingGUIFromParameterNode = False

  def updateParameterNodeFromGUI(self, caller=None, event=None):
    """
    This method is called when the user makes any change in the GUI.
    The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    # self._parameterNode.SetNodeReferenceID("InputVolume", self.ui.inputSelector.currentNodeID)
    # self._parameterNode.SetNodeReferenceID("OutputVolume", self.ui.outputSelector.currentNodeID)
    # self._parameterNode.SetParameter("Threshold", str(self.ui.imageThresholdSliderWidget.value))
    # self._parameterNode.SetParameter("Invert", "true" if self.ui.invertOutputCheckBox.checked else "false")
    # self._parameterNode.SetNodeReferenceID("OutputVolumeInverse", self.ui.invertedOutputSelector.currentNodeID)

    self._parameterNode.EndModify(wasModified)

  # def onApplyButton(self):
  #   """
  #   Run processing when user clicks "Apply" button.
  #   """
  #   try:

  #     # Compute output
  #     self.logic.process(self.ui.inputSelector.currentNode(), self.ui.outputSelector.currentNode(),
  #       self.ui.imageThresholdSliderWidget.value, self.ui.invertOutputCheckBox.checked)

  #     # Compute inverted output (if needed)
  #     if self.ui.invertedOutputSelector.currentNode():
  #       # If additional output volume is selected then result with inverted threshold is written there
  #       self.logic.process(self.ui.inputSelector.currentNode(), self.ui.invertedOutputSelector.currentNode(),
  #         self.ui.imageThresholdSliderWidget.value, not self.ui.invertOutputCheckBox.checked, showResult=False)

  #   except Exception as e:
  #     slicer.util.errorDisplay("Failed to compute results: "+str(e))
  #     import traceback
  #     traceback.print_exc()


#
# BigImageViewerLogic
#

class BigImageViewerLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    """
    Called when the logic class is instantiated. Can be used for initializing member variables.
    """
    ScriptedLoadableModuleLogic.__init__(self)

  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    if not parameterNode.GetParameter("Threshold"):
      parameterNode.SetParameter("Threshold", "100.0")
    if not parameterNode.GetParameter("Invert"):
      parameterNode.SetParameter("Invert", "false")

  def installPreRequisites(self):
    import os.path
    from os.path import expanduser
    homeDir = expanduser("~")

    try:
        import skimage
    except ImportError:
        slicer.util.pip_install("scikit-image")

    import skimage.transform

    # try:
    #     import h5py
    # except ImportError:
    #     slicer.util.pip_install("h5py")

    # try:
    #     import tensorflow
    # except ImportError:
    #     slicer.util.pip_install("tensorflow")


    # try:
    #     import keras.models
    # except ImportError:
    #     slicer.util.pip_install("keras")

    os.environ["CUDA_VISIBLE_DEVICES"] = "0"



  def getSlideInfo(self, slidePath):
    slide = TiffSlide(slidePath)
    info = {}

    # 基础金字塔信息
    info["level_count"] = slide.level_count
    info["level_dimensions"] = slide.level_dimensions
    info["level_downsamples"] = slide.level_downsamples

    # --------------------------
    # 修复：正确获取 MPP
    # --------------------------
    try:
        # 方式 1：新版 tiffslide
        mpp_x = slide.properties.get("tiffslide.mpp-x", None)
        mpp_y = slide.properties.get("tiffslide.mpp-y", None)
        if mpp_x is None:
            # 方式 2：从 OpenSlide 兼容属性读取
            mpp_x = slide.properties.get("openslide.mpp-x", 0.5)
            mpp_y = slide.properties.get("openslide.mpp-y", 0.5)
    except:
        mpp_x = 0.5
        mpp_y = 0.5

    info["mpp"] = (float(mpp_x), float(mpp_y))

    # 放大倍数
    info["objectiveMagnification"] = 20.0
    try:
        obj = slide.properties.get("tiffslide.objective-power", None)
        if obj is None:
            obj = slide.properties.get("openslide.objective-power", 20.0)
        info["objectiveMagnification"] = float(obj)
    except Exception:
        pass

    slide.close()
    return info
    # import tiffslide as openslide
    # slide = openslide.OpenSlide(slidePath)
    
    # # 1. Safely extract metadata with defaults (The previous fix)
    # props = slide.properties
    # mpp_x = float(props.get('openslide.mpp-x', 1.0))
    # mpp_y = float(props.get('openslide.mpp-y', 1.0))
    # obj_power = float(props.get('openslide.objective-power', 40.0))
    # vendor = props.get('openslide.vendor', 'Unknown')
    
    # # 2. Extract intrinsic slide object attributes (The missing piece)
    # width, height = slide.dimensions
    
    # slideInfo = {
    #     'width': width,
    #     'height': height,
    #     'mpp': [mpp_x, mpp_y],
    #     'objectiveMagnification': obj_power,
    #     'vendor': vendor,
        
    #     # Pyramid structural data required by BigImageViewer UI
    #     'level_count': slide.level_count,
    #     'level_dimensions': slide.level_dimensions,
    #     'level_downsamples': slide.level_downsamples
    # }
    
    # slide.close()
    # return slideInfo



  def segment2DRGBPatchBatch(self, model, inputImagePatchBatchArray):
    # This segments a list of RGB 2D patches. The input
    # inputImagePatchArray is a (#batches, patchSideLen, patchSideLen,
    # numChannel) numpy array. Next this needs to be reshaped to
    # (#batches, patchSideLen, patchSideLen, #channels)
    # where and #channels=1 in this case

    sz = inputImagePatchBatchArray.shape

    inputImagePatchBatchArray /= 255.0

    testBatchSize = sz[0]
    results = model.predict(inputImagePatchBatchArray, testBatchSize, verbose=1)

    outputSegBatchArray = results[:, :, :, :]

    return outputSegBatchArray


  def processSegmentGland(self, inputVolume, outputVolume, kerasModel):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputVolume: volume to be thresholded
    :param outputVolume: thresholding result
    :param imageThreshold: values above/below this threshold will be set to 0
    :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
    :param showResult: show output volume in slice viewers
    """

    if not inputVolume or not outputVolume:
      raise ValueError("Input or output volume is invalid")

    import time
    startTime = time.time()
    logging.info('Processing started')

    ################################################################################
    # keras seg
    imgArray = slicer.util.arrayFromVolume(inputVolume)
    imgArray = imgArray[0, :, :, :]

    testImgSeg = self.segment2DRGBImageRandomSampleDividePrior(model = kerasModel, imageArray = imgArray.astype(numpy.float), patchSideLen = 64, numPatchSampleFactor = 5, batch_size = 10)

    testImgSeg = testImgSeg[numpy.newaxis, :, :]

    # volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    # outputVolume.CreateDefaultDisplayNodes()
    slicer.util.updateVolumeFromArray(outputVolume, testImgSeg)
    #setSliceViewerLayers(background=volumeNode)
    slicer.util.setSliceViewerLayers(foreground=outputVolume)


    # keras seg, end
    ################################################################################

    stopTime = time.time()
    logging.info('Processing completed in {0:.2f} seconds'.format(stopTime-startTime))



  def processSegmentNuclei(self, inputVolume, outputVolume):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputVolume: volume to be thresholded
    :param outputVolume: thresholding result
    :param imageThreshold: values above/below this threshold will be set to 0
    :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
    :param showResult: show output volume in slice viewers
    """

    if not inputVolume or not outputVolume:
      raise ValueError("Input or output volume is invalid")

    import time
    startTime = time.time()
    logging.info('Processing started')

    # # Compute the thresholded output volume using the "Threshold Scalar Volume" CLI module
    # cliParams = {
    #   'InputVolume': inputVolume.GetID(),
    #   'OutputVolume': outputVolume.GetID(),
    #   'ThresholdValue' : imageThreshold,
    #   'ThresholdType' : 'Above' if invert else 'Below'
    #   }
    # cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True, update_display=showResult)
    # # We don't need the CLI module node anymore, remove it to not clutter the scene with it
    # slicer.mrmlScene.RemoveNode(cliNode)


    ################################################################################
    # keras seg
    imgArray = slicer.util.arrayFromVolume(inputVolume)
    imgArray = imgArray[0, :, :, :]

    testImgSeg = self.segment2DRGBImageRandomSampleDividePrior(model = kerasModel, imageArray = imgArray.astype(numpy.float), patchSideLen = 64, numPatchSampleFactor = 5, batch_size = 10)

    testImgSeg = testImgSeg[numpy.newaxis, :, :]

    # volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    # outputVolume.CreateDefaultDisplayNodes()
    slicer.util.updateVolumeFromArray(outputVolume, testImgSeg)
    #setSliceViewerLayers(background=volumeNode)
    slicer.util.setSliceViewerLayers(foreground=outputVolume)


    # keras seg, end
    ################################################################################

    stopTime = time.time()
    logging.info('Processing completed in {0:.2f} seconds'.format(stopTime-startTime))




  def process(self, inputVolume, outputVolume, imageThreshold, invert=False, showResult=True):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputVolume: volume to be thresholded
    :param outputVolume: thresholding result
    :param imageThreshold: values above/below this threshold will be set to 0
    :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
    :param showResult: show output volume in slice viewers
    """

    if not inputVolume or not outputVolume:
      raise ValueError("Input or output volume is invalid")

    import time
    startTime = time.time()
    logging.info('Processing started')

    # Compute the thresholded output volume using the "Threshold Scalar Volume" CLI module
    cliParams = {
      'InputVolume': inputVolume.GetID(),
      'OutputVolume': outputVolume.GetID(),
      'ThresholdValue' : imageThreshold,
      'ThresholdType' : 'Above' if invert else 'Below'
      }
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True, update_display=showResult)
    # We don't need the CLI module node anymore, remove it to not clutter the scene with it
    slicer.mrmlScene.RemoveNode(cliNode)

    stopTime = time.time()
    logging.info(f'Processing completed in {stopTime-startTime:.2f} seconds')


#
# BigImageViewerTest
#

class BigImageViewerTest(ScriptedLoadableModuleTest):
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
    self.test_BigImageViewer1()

  def test_BigImageViewer1(self):
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

    # Get/create input data

    import SampleData
    #registerSampleData()
    inputVolume = SampleData.downloadSample('BigImageViewer1')
    self.delayDisplay('Loaded test data set')

    inputScalarRange = inputVolume.GetImageData().GetScalarRange()
    self.assertEqual(inputScalarRange[0], 0)
    self.assertEqual(inputScalarRange[1], 695)

    outputVolume = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    threshold = 100

    # Test the module logic

    logic = BigImageViewerLogic()

    # Test algorithm with non-inverted threshold
    logic.process(inputVolume, outputVolume, threshold, True)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.process(inputVolume, outputVolume, threshold, False)
    outputScalarRange = outputVolume.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')
