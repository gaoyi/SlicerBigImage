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

    interactor.AddObserver(vtk.vtkCommand.MouseWheelForwardEvent, self.onMouseWheelForwardEvent)
    interactor.AddObserver(vtk.vtkCommand.MouseWheelBackwardEvent, self.onMouseWheelBackwardEvent)

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
    print ('onMouseWheelForwardEvent............................')
    if not self.ui.ObjectiveMagnificationSlicerWidget.isMaximized():
      #self.ui.ObjectiveMagnificationSlicerWidget.setValue(self.ui.ObjectiveMagnificationSlicerWidget.value + 2.0*self.ui.ObjectiveMagnificationSlicerWidget.singleStep)

      # set the increment to be 1/2 of the current values accelerates
      # the zooming well when the mag is large. It also make the
      # zooming stable when the value is small. This is better than
      # setting the increament to a fixed value
      self.ui.ObjectiveMagnificationSlicerWidget.setValue(self.ui.ObjectiveMagnificationSlicerWidget.value + 0.5*self.ui.ObjectiveMagnificationSlicerWidget.value)

      # If send this signal every time the slider changes, will trigger loading patch. This will make the zooming very slow
      #self.ui.ObjectiveMagnificationSlicerWidget.valueChanged(self.ui.ObjectiveMagnificationSlicerWidget.value)


  def onMouseWheelBackwardEvent(self, obj, event=None):
    print ('onMouseWheelBackwardEvent............................')
    if not self.ui.ObjectiveMagnificationSlicerWidget.isMinimized():
      #self.ui.ObjectiveMagnificationSlicerWidget.setValue(self.ui.ObjectiveMagnificationSlicerWidget.value - 2.0*self.ui.ObjectiveMagnificationSlicerWidget.singleStep)

      # set the increment to be 1/2 of the current values accelerates
      # the zooming well when the mag is large. It also make the
      # zooming stable when the value is small. This is better than
      # setting the increament to a fixed value
      self.ui.ObjectiveMagnificationSlicerWidget.setValue(self.ui.ObjectiveMagnificationSlicerWidget.value - 0.5*self.ui.ObjectiveMagnificationSlicerWidget.value)

      # If send this signal every time the slider changes, will trigger loading patch. This will make the zooming very slow
      #self.ui.ObjectiveMagnificationSlicerWidget.valueChanged(self.ui.ObjectiveMagnificationSlicerWidget.value)





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

    if os.name == 'nt':
      openSlidePath = self.logic.getOpenSlidePath()
      self.ui.OpenSlidePathLineEdit.currentPath = openSlidePath
    else:
      self.ui.OpenSlidePathLineEdit.hide()
      # There are no other advanced settings for now, so hide the whole section
      self.ui.advancedCollapsibleButton.hide()

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
    volumeNode.SetName("currentPatchGrayChannel")
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
    volumeNode.SetAndObserveImageData(imageData)
    slicer.mrmlScene.AddNode(volumeNode)

    # Add volume to scene
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()

    slicer.mrmlScene.AddNode(displayNode)
    #colorNode = slicer.util.getNode('Green')
    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    volumeNode.CreateDefaultStorageNode()

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
    volumeNode.SetName("currentPatchFromBigRGBAImage")
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
    volumeNode.SetAndObserveImageData(imageData)
    slicer.mrmlScene.AddNode(volumeNode)

    # Add volume to scene
    displayNode = slicer.vtkMRMLVectorVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)

    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    volumeNode.CreateDefaultStorageNode()

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
    import openslide
    ds = (self.objectiveMagnificationMax)/(self.ui.ObjectiveMagnificationSlicerWidget.value)

    slide = openslide.OpenSlide(self.BigRGBAImagePathname)
    self.BigRGBAImageLevelToLoad = slide.get_best_level_for_downsample(ds)
    slide.close()

    #print(self.BigRGBAImageLevelToLoad)

    self.updateUIWidget()

    self.loadPatchFromBigRGBAImage()

  def enableAndInitUIWidget(self):
    #self.loadBigRGBAImageButton.enabled = True

    #--------------------------------------------------------------------------------
    # Set the max of slider for top-left corner
    # self.ui.topLeftXSliderWidget.value = self.topLeftX
    # self.ui.topLeftYSliderWidget.value = self.topLeftY

    self.ui.topLeftXSliderWidget.value = 0
    self.ui.topLeftYSliderWidget.value = 0

    self.ui.topLeftXSliderWidget.maximum = max(0, self.WSISizesXAtAllLevels[0] - self.patchSizeX*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1)
    self.ui.topLeftYSliderWidget.maximum = max(0, self.WSISizesYAtAllLevels[0] - self.patchSizeY*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1)

    self.ui.topLeftXSliderWidget.enabled = True
    self.ui.topLeftYSliderWidget.enabled = True

    #--------------------------------------------------------------------------------
    # Set the max of number of levels
    self.ui.ObjectiveMagnificationSlicerWidget.maximum = self.objectiveMagnificationMax
    self.ui.ObjectiveMagnificationSlicerWidget.minimum = self.objectiveMagnificationMin
    self.ui.ObjectiveMagnificationSlicerWidget.value = 1.0
    self.ui.ObjectiveMagnificationSlicerWidget.enabled = True


  def updateUIWidget(self):
    # self.ui.topLeftXSliderWidget.value = self.topLeftX
    # self.ui.topLeftYSliderWidget.value = self.topLeftY

    self.ui.topLeftXSliderWidget.maximum = max(0, self.WSISizesXAtAllLevels[0] - self.patchSizeX*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1)
    self.ui.topLeftYSliderWidget.maximum = max(0, self.WSISizesYAtAllLevels[0] - self.patchSizeY*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1)

    #self.ui.ObjectiveMagnificationSlicerWidget.maximum = self.BigRGBAImageNumberOfLevels - 1


  def getRegionFromFileAsRGBNumpyArray(self, WSIName, level, topX0, topY0, width, height):
    import openslide

    # I currently do not have range check. May add later
    slide = openslide.OpenSlide(WSIName)

    thisTilePilIm = slide.read_region((topX0, topY0), level, (width, height))
    if thisTilePilIm.mode != "RGB":
        thisTilePilIm = thisTilePilIm.convert("RGB")

    imRGBNdarray = numpy.asarray(thisTilePilIm)

    slide.close()

    return imRGBNdarray

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


    magAtThisLevel = float(self.slideInfo["objectiveMagnification"])/float(self.level_downsamples[self.BigRGBAImageLevelToLoad])

    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print(self.ui.ObjectiveMagnificationSlicerWidget.value)
    ratio = magAtThisLevel/float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    # print(ratio)
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    sizeToLoadX = round(self.patchSizeX*ratio)
    sizeToLoadY = round(self.patchSizeY*ratio)

    self.topLeftX = int(self.ui.topLeftXSliderWidget.value)

    #print("RGB requested size = ", (sizeToLoadY, sizeToLoadX))
    q = self.getRegionFromFileAsRGBNumpyArray(self.BigRGBAImagePathname, self.BigRGBAImageLevelToLoad, int(self.ui.topLeftXSliderWidget.value), int(self.ui.topLeftYSliderWidget.value), int(sizeToLoadX), int(sizeToLoadY))

    # print(q.max())

    p = skimage.transform.resize(q, (int(self.patchSizeY), int(self.patchSizeX)), preserve_range=True, mode='wrap', order=0, anti_aliasing=False)
    #print("RGB Got size = ", q.shape, "resize to ", p.shape, q.shape[0]/p.shape[0], q.shape[1]/p.shape[1])
    # print(p.max())

    #a = slicer.util.array('currentPatchFromBigRGBAImage')
    volumeNode = slicer.util.getNode('currentPatchFromBigRGBAImage')
    a = slicer.util.arrayFromVolume(volumeNode)

    # resizing is needed for python (not for CLI version) coz we have
    # to fill the numpy array here
    a[:] = p

    resolutionNow = self.MPP/1000.0*float(self.slideInfo["objectiveMagnification"])/float(self.ui.ObjectiveMagnificationSlicerWidget.value)
    imageSpacing = [resolutionNow, resolutionNow, self.SliceThickness]

    n = slicer.util.getNode('currentPatchFromBigRGBAImage')
    n.SetSpacing(imageSpacing)

    n.GetImageData().Modified()

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

    self.patchSizeX = int(redView.width/2)
    self.patchSizeY = int(redView.height/2)

    return

  def getMetaInfoFromWSI(self):

    self.logic.installPreRequisites()
    self.slideInfo = self.logic.getSlideInfo(self.BigRGBAImagePathname)

    print("returned from logic.getSlideInfo")

    self.BigRGBAImageNumberOfLevels = self.slideInfo["level_count"]
    self.level_downsamples = self.slideInfo["level_downsamples"]
    self.MPP = self.slideInfo["mpp"][0]

    for it in range(self.BigRGBAImageNumberOfLevels):
      self.WSISizesXAtAllLevels.append(self.slideInfo["level_dimensions"][it][0])
      self.WSISizesYAtAllLevels.append(self.slideInfo["level_dimensions"][it][1])

    self.objectiveMagnificationMax = self.slideInfo["objectiveMagnification"]
    self.objectiveMagnificationMin = self.slideInfo["objectiveMagnification"]/self.level_downsamples[-1]/3.0 # so to include more

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

    # For now OpenSlide is only installed on Windows
    if os.name == 'nt':
      self.logic.setOpenSlidePath(self.ui.OpenSlidePathLineEdit.currentPath)
      if not self.logic.isOpenSlidePathValid():
        self.logic.findOpenSlide()
        if not self.logic.isOpenSlidePathValid() and os.name == 'nt':  # TODO: implement download for Linux/MacOS?
          # OpenSlide not found, offer downloading it
          if slicer.util.confirmOkCancelDisplay(
              'OpenSlide not detected on your system. '
              'Download OpenSlide library?',
                  windowTitle='Download confirmation'):
              if not self.logic.openSlideDownload():
                  slicer.util.errorDisplay("OpenSlide download failed")
      if not self.logic.isOpenSlidePathValid():
          # still not found, user has to specify path manually
          self.advancedCollapsibleButton.collapsed = False
          return
      self.ui.OpenSlidePathLineEdit.currentPath = self.logic.getOpenSlidePath()

    self.determinPatchSizeByViewerSize()

    #--------------------------------------------------------------------------------
    # Read meta infrom from WSI image
    #if self.BigRGBAImageFileNameEditor.currentPath and os.path.isfile(self.ui.BigRGBAImageFileNameEditor.currentPath):
    if not os.path.isfile(self.ui.BigRGBAImageFileNameEditor.currentPath):
      print("File not exist")
      return

    self.logic.installPreRequisites()

    self.BigRGBAImagePathname = self.ui.BigRGBAImageFileNameEditor.currentPath
    self.getMetaInfoFromWSI()

    self.setupVolumeNodeToStoreRGBPatch()
    self.setupForScalarPatch()

    # Setting the magnification slide bar will automatically trigger
    # the onWSILevelChanged in which will call the
    # loadPatchFromBigRGBAImage
    #self.ui.ObjectiveMagnificationSlicerWidget.value = self.objectiveMagnificationMax
    self.ui.ObjectiveMagnificationSlicerWidget.value = self.objectiveMagnificationMin

    # Setting the magnification slide bar will automatically trigger
    # the onWSILevelChanged in which will call the
    # loadPatchFromBigRGBAImage. So no need to do again
    # self.loadPatchFromBigRGBAImage()

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

    volumeNode = slicer.util.getNode('currentPatchGrayChannel') #self.patchGrayVolumeNode
    slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage') #self.patchVolumeNode
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)


  def onDetectMitosisButton(self):
    pass

  def onDetectGlandButton(self):
    self.logic.installPreRequisites()
    self.logic.processSegmentGland(self.patchVolumeNode, self.patchGrayVolumeNode, self.kerasModelSegmentColonGland)

    volumeNode = slicer.util.getNode('currentPatchGrayChannel')
    slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage')
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)

    return

  def onSegmentNucleiButton(self):
    #self.logic.installPreRequisites()

    nucleusLabelImageNode = slicer.vtkMRMLLabelMapVolumeNode()
    nucleusLabelImageNode.SetName("nucleusLabelImage")
    slicer.mrmlScene.AddNode( nucleusLabelImageNode )
    #parameters['outputAxisMaskVolume'] =

    volumeNode = slicer.util.getNode('currentPatchGrayChannel')
    #slicer.util.setSliceViewerLayers(foreground=volumeNode)

    parameters = {}
    parameters['inputfileName'] = self.patchVolumeNode.GetID()
    parameters['outputEImage'] = volumeNode.GetID()
    parameters['outputHImage'] = volumeNode.GetID()
    parameters['outputImageName'] = nucleusLabelImageNode.GetID()
    parameters['mpp'] = 0.5/20.0*float(self.ui.ObjectiveMagnificationSlicerWidget.value) # hack: assume mpp=0.5 at 20x
    slicer.cli.run( slicer.modules.nucleussegmentation, None, parameters, wait_for_completion=True )

    # volumeNode = slicer.util.getNode('currentPatchGrayChannel') #self.patchGrayVolumeNode
    # slicer.util.setSliceViewerLayers(foreground=volumeNode)
    # volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage') #self.patchVolumeNode
    # slicer.util.setSliceViewerLayers(background = volumeNode1)

    # slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)


    # self.logic.processSegmentNuclei(self.patchVolumeNode, nucleusLabelImageNode.GetID())


    # volumeNode = slicer.util.getNode('currentPatchGrayChannel')
    # slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage')
    slicer.util.setSliceViewerLayers(background = volumeNode1)
    #nucleusLabelImageNode

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
        import openslide
    except ImportError:
        slicer.util.pip_install("openslide-python")

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

    if os.name == 'nt':
      import os
      os.add_dll_directory(os.path.dirname(self.getOpenSlidePath()))



  def getSlideInfo(self, svsPathname):
    import openslide

    slide = openslide.OpenSlide(svsPathname)

    slideInfo = {'slidePathname':svsPathname,
                 'level_count':slide.level_count,
                 'level_dimensions': slide.level_dimensions,
                 'level_downsamples': slide.level_downsamples,
                 'mpp': [float(slide.properties['openslide.mpp-x']), float(slide.properties['openslide.mpp-y'])],
                 'objectiveMagnification': float(slide.properties['openslide.objective-power'])}
    slide.close()

    return slideInfo



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
    imgArray = slicer.util.array("currentPatchFromBigRGBAImage")
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
    imgArray = slicer.util.array("currentPatchFromBigRGBAImage")
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


  def isOpenSlidePathValid(self):
      import os
      openSlidePath = self.getOpenSlidePath()
      return os.path.isfile(openSlidePath)

  def getDownloadedOpenSlideDirectory(self):
      return os.path.dirname(slicer.app.slicerUserSettingsFilePath) + '/OpenSlide'

  def getOpenSlideExecutableFilename(self):
      return 'openslide.jar'

  def findOpenSlide(self):
      # Try to find the executable at specific paths
      commonOpenSlidePaths = [
          '/usr/local/bin/openslide',
          '/usr/bin/openslide'
      ]
      for openSlidePath in commonOpenSlidePaths:
          if os.path.isfile(openSlidePath):
              # found one
              self.setOpenSlidePath(openSlidePath)
              return True
      # Search for the executable in directories
      commonOpenSlideDirs = [
          self.getDownloadedOpenSlideDirectory()
      ]
      for openSlideDir in commonOpenSlideDirs:
          if self.findOpenSlideInDirectory(openSlideDir):
              # found it
              return True
      # Not found
      return False

  def findOpenSlideInDirectory(self, openSlideDir):
      openSlideExecutableFilename = self.getOpenSlideExecutableFilename()
      for dirpath, dirnames, files in os.walk(openSlideDir):
          for name in files:
              if name == openSlideExecutableFilename:
                  openSlideExecutablePath = (dirpath + '/' + name).replace('\\', '/')
                  self.setOpenSlidePath(openSlideExecutablePath)
                  return True
      return False

  def unzipOpenSlide(self, filePath, openSlideTargetDirectory):
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
          logging.info('openSlide package is not found at ' + filePath)
          return False

      logging.info('Unzipping OpenSlide package ' + filePath)
      qt.QDir().mkpath(openSlideTargetDirectory)
      slicer.app.applicationLogic().Unzip(filePath, openSlideTargetDirectory)
      success = self.findOpenSlideInDirectory(openSlideTargetDirectory)
      return success

  def openSlideDownload(self):
      openSlideTargetDirectory = self.getDownloadedOpenSlideDirectory()
      # The number in the filePath can be incremented each time a significantly different OpenSlide version
      # is to be introduced (it prevents reusing a previously downloaded package).
      filePath = slicer.app.temporaryPath + '/OpenSlide-package-slicer-01.zip'
      success = self.unzipOpenSlide(filePath, openSlideTargetDirectory)
      if success:
          # there was a valid downloaded package already
          return True

      # List of mirror sites to attempt download OpenSlide pre-built binaries from
      urls = []
      if os.name == 'nt':
          urls.append('https://github.com/openslide/openslide-winbuild/releases/download/v20160612/openslide-win64-20160612.zip')
      else:
          # TODO: implement downloading for Linux/MacOS?
          pass

      success = False
      qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

      for url in urls:

          success = True
          try:
              logging.info('Requesting download OpenSlide from %s...' % url)
              import urllib.request, urllib.error, urllib.parse
              req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
              data = urllib.request.urlopen(req).read()
              with open(filePath, "wb") as f:
                  f.write(data)

              success = self.unzipOpenSlide(filePath, openSlideTargetDirectory)
          except:
              success = False

          if success:
              break

      qt.QApplication.restoreOverrideCursor()
      return success

  def getOpenSlidePath(self):
      settings = qt.QSettings()
      if settings.contains('BigImage/OpenSlidePath'):
          return slicer.app.toSlicerHomeAbsolutePath(settings.value('BigImage/OpenSlidePath'))
      return ''

  def setOpenSlidePath(self, openSlidePath):
      # don't save it if already saved
      settings = qt.QSettings()
      if settings.contains('BigImage/OpenSlidePath'):
          if openSlidePath == slicer.app.toSlicerHomeAbsolutePath(settings.value('BigImage/OpenSlidePath')):
              return
      settings.setValue('BigImage/OpenSlidePath', slicer.app.toSlicerHomeRelativePath(openSlidePath))

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
