#include "vtkMRMLBigImageNode.h"

#include "vtkMRMLNodePropertyMacros.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkMRMLBigImageNode);

//----------------------------------------------------------------------------
vtkMRMLBigImageNode::vtkMRMLBigImageNode() = default;

//----------------------------------------------------------------------------
vtkMRMLBigImageNode::~vtkMRMLBigImageNode()
{
  this->SetSourceURI(nullptr);
  this->SetImageType(nullptr);
  this->SetAccessMode(nullptr);
  this->SetPixelType(nullptr);
  this->SetFrameUID(nullptr);
}

//----------------------------------------------------------------------------
vtkMRMLNode* vtkMRMLBigImageNode::CreateNodeInstance()
{
  return vtkMRMLBigImageNode::New();
}

//----------------------------------------------------------------------------
void vtkMRMLBigImageNode::ReadXMLAttributes(const char** atts)
{
  int disabledModify = this->StartModify();

  Superclass::ReadXMLAttributes(atts);

  vtkMRMLReadXMLBeginMacro(atts);
  vtkMRMLReadXMLStringMacro(sourceURI, SourceURI);
  vtkMRMLReadXMLStringMacro(imageType, ImageType);
  vtkMRMLReadXMLStringMacro(accessMode, AccessMode);
  vtkMRMLReadXMLStringMacro(pixelType, PixelType);
  vtkMRMLReadXMLStringMacro(frameUID, FrameUID);
  vtkMRMLReadXMLVectorMacro(dimensions, Dimensions, int, 3);
  vtkMRMLReadXMLVectorMacro(spacing, Spacing, double, 3);
  vtkMRMLReadXMLVectorMacro(mpp, MPP, double, 2);
  vtkMRMLReadXMLIntMacro(levelCount, LevelCount);
  vtkMRMLReadXMLIntMacro(channelCount, ChannelCount);
  vtkMRMLReadXMLBooleanMacro(pyramid, Pyramid);
  vtkMRMLReadXMLBooleanMacro(mppEstimated, MPPEstimated);
  vtkMRMLReadXMLEndMacro();

  this->EndModify(disabledModify);
}

//----------------------------------------------------------------------------
void vtkMRMLBigImageNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of, nIndent);

  vtkMRMLWriteXMLBeginMacro(of);
  vtkMRMLWriteXMLStringMacro(sourceURI, SourceURI);
  vtkMRMLWriteXMLStringMacro(imageType, ImageType);
  vtkMRMLWriteXMLStringMacro(accessMode, AccessMode);
  vtkMRMLWriteXMLStringMacro(pixelType, PixelType);
  vtkMRMLWriteXMLStringMacro(frameUID, FrameUID);
  vtkMRMLWriteXMLVectorMacro(dimensions, Dimensions, int, 3);
  vtkMRMLWriteXMLVectorMacro(spacing, Spacing, double, 3);
  vtkMRMLWriteXMLVectorMacro(mpp, MPP, double, 2);
  vtkMRMLWriteXMLIntMacro(levelCount, LevelCount);
  vtkMRMLWriteXMLIntMacro(channelCount, ChannelCount);
  vtkMRMLWriteXMLBooleanMacro(pyramid, Pyramid);
  vtkMRMLWriteXMLBooleanMacro(mppEstimated, MPPEstimated);
  vtkMRMLWriteXMLEndMacro();
}

//----------------------------------------------------------------------------
void vtkMRMLBigImageNode::Copy(vtkMRMLNode* node)
{
  if (this == node)
  {
    return;
  }

  Superclass::Copy(node);

  vtkMRMLBigImageNode* sourceNode = vtkMRMLBigImageNode::SafeDownCast(node);
  if (!sourceNode)
  {
    return;
  }

  this->SetSourceURI(sourceNode->GetSourceURI());
  this->SetImageType(sourceNode->GetImageType());
  this->SetAccessMode(sourceNode->GetAccessMode());
  this->SetPixelType(sourceNode->GetPixelType());
  this->SetFrameUID(sourceNode->GetFrameUID());
  this->SetDimensions(sourceNode->GetDimensions());
  this->SetSpacing(sourceNode->GetSpacing());
  this->SetMPP(sourceNode->GetMPP());
  this->SetLevelCount(sourceNode->GetLevelCount());
  this->SetChannelCount(sourceNode->GetChannelCount());
  this->SetPyramid(sourceNode->GetPyramid());
  this->SetMPPEstimated(sourceNode->GetMPPEstimated());
}

//----------------------------------------------------------------------------
void vtkMRMLBigImageNode::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os, indent);

  vtkMRMLPrintBeginMacro(os, indent);
  vtkMRMLPrintStringMacro(SourceURI);
  vtkMRMLPrintStringMacro(ImageType);
  vtkMRMLPrintStringMacro(AccessMode);
  vtkMRMLPrintStringMacro(PixelType);
  vtkMRMLPrintStringMacro(FrameUID);
  vtkMRMLPrintVectorMacro(Dimensions, int, 3);
  vtkMRMLPrintVectorMacro(Spacing, double, 3);
  vtkMRMLPrintVectorMacro(MPP, double, 2);
  vtkMRMLPrintIntMacro(LevelCount);
  vtkMRMLPrintIntMacro(ChannelCount);
  vtkMRMLPrintBooleanMacro(Pyramid);
  vtkMRMLPrintBooleanMacro(MPPEstimated);
  vtkMRMLPrintEndMacro();
}
