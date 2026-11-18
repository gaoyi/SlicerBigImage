#ifndef __vtkMRMLBigImageNode_h
#define __vtkMRMLBigImageNode_h

#include "vtkSlicerBigImageCoreModuleMRMLExport.h"
#include "vtkMRMLTransformableNode.h"

/// MRML node representing a two-dimensional big image (for example a WSI).
///
/// Pixel data is not stored in this node. The node stores the source URI,
/// image geometry and acquisition metadata, and references a pose plane and
/// transform. Pixels are read on demand by the viewer/reader backend.
class VTK_SLICER_BIGIMAGECORE_MODULE_MRML_EXPORT vtkMRMLBigImageNode : public vtkMRMLTransformableNode
{
public:
  static vtkMRMLBigImageNode* New();
  vtkTypeMacro(vtkMRMLBigImageNode, vtkMRMLTransformableNode);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkMRMLNode* CreateNodeInstance() override;
  const char* GetNodeTagName() override { return "BigImage"; }

  vtkSetStringMacro(SourceURI);
  vtkGetStringMacro(SourceURI);

  vtkSetStringMacro(ImageType);
  vtkGetStringMacro(ImageType);

  vtkSetStringMacro(AccessMode);
  vtkGetStringMacro(AccessMode);

  vtkSetStringMacro(PixelType);
  vtkGetStringMacro(PixelType);

  vtkSetStringMacro(FrameUID);
  vtkGetStringMacro(FrameUID);

  vtkSetVector3Macro(Dimensions, int);
  vtkGetVector3Macro(Dimensions, int);

  vtkSetVector3Macro(Spacing, double);
  vtkGetVector3Macro(Spacing, double);

  vtkSetVector2Macro(MPP, double);
  vtkGetVector2Macro(MPP, double);

  vtkSetMacro(LevelCount, int);
  vtkGetMacro(LevelCount, int);

  vtkSetMacro(ChannelCount, int);
  vtkGetMacro(ChannelCount, int);

  vtkSetMacro(Pyramid, bool);
  vtkGetMacro(Pyramid, bool);
  vtkBooleanMacro(Pyramid, bool);

  vtkSetMacro(MPPEstimated, bool);
  vtkGetMacro(MPPEstimated, bool);
  vtkBooleanMacro(MPPEstimated, bool);

  void SetPosePlaneNodeID(const char* nodeID) { this->SetNodeReferenceID("PosePlane", nodeID); }
  const char* GetPosePlaneNodeID() { return this->GetNodeReferenceID("PosePlane"); }
  vtkMRMLNode* GetPosePlaneNode() { return this->GetNodeReference("PosePlane"); }

  void ReadXMLAttributes(const char** atts) override;
  void WriteXML(ostream& of, int indent) override;
  void Copy(vtkMRMLNode* node) override;

protected:
  vtkMRMLBigImageNode();
  ~vtkMRMLBigImageNode() override;

  char* SourceURI{ nullptr };
  char* ImageType{ nullptr };
  char* AccessMode{ nullptr };
  char* PixelType{ nullptr };
  char* FrameUID{ nullptr };
  int Dimensions[3]{ 0, 0, 1 };
  double Spacing[3]{ 1.0, 1.0, 1.0 };
  double MPP[2]{ 1.0, 1.0 };
  int LevelCount{ 0 };
  int ChannelCount{ 1 };
  bool Pyramid{ false };
  bool MPPEstimated{ false };

private:
  vtkMRMLBigImageNode(const vtkMRMLBigImageNode&) = delete;
  void operator=(const vtkMRMLBigImageNode&) = delete;
};

#endif
