#include "vtkSlicerBigImageCoreModuleLogic.h"

#include "vtkMRMLBigImageNode.h"

#include "vtkCallbackCommand.h"
#include "vtkCollection.h"
#include "vtkMatrix4x4.h"
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLMarkupsPlaneDisplayNode.h"
#include "vtkMRMLMarkupsPlaneNode.h"
#include "vtkMRMLScene.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"

vtkStandardNewMacro(vtkSlicerBigImageCoreModuleLogic);

//----------------------------------------------------------------------------
vtkSlicerBigImageCoreModuleLogic::vtkSlicerBigImageCoreModuleLogic() = default;

//----------------------------------------------------------------------------
vtkSlicerBigImageCoreModuleLogic::~vtkSlicerBigImageCoreModuleLogic()
{
  if (this->SceneEventsCallback)
  {
    if (this->GetMRMLScene())
    {
      this->GetMRMLScene()->RemoveObserver(this->SceneEventsCallback);
    }
    this->SceneEventsCallback->Delete();
    this->SceneEventsCallback = nullptr;
  }
}

//----------------------------------------------------------------------------
void vtkSlicerBigImageCoreModuleLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkSlicerBigImageCoreModuleLogic::RegisterNodes()
{
  vtkMRMLScene* scene = this->GetMRMLScene();
  if (!scene)
  {
    return;
  }

  vtkSmartPointer<vtkMRMLBigImageNode> bigImageNode = vtkSmartPointer<vtkMRMLBigImageNode>::New();
  scene->RegisterNodeClass(bigImageNode);

  if (!this->SceneEventsCallback)
  {
    this->SceneEventsCallback = vtkCallbackCommand::New();
    this->SceneEventsCallback->SetCallback(&vtkSlicerBigImageCoreModuleLogic::OnSceneEvents);
    this->SceneEventsCallback->SetClientData(this);
  }

  scene->AddObserver(vtkMRMLScene::EndImportEvent, this->SceneEventsCallback);
  scene->AddObserver(vtkMRMLScene::EndRestoreEvent, this->SceneEventsCallback);
}

//----------------------------------------------------------------------------
void vtkSlicerBigImageCoreModuleLogic::OnSceneEvents(
  vtkObject* vtkNotUsed(caller), unsigned long event, void* clientData, void* vtkNotUsed(callData))
{
  vtkSlicerBigImageCoreModuleLogic* self =
    static_cast<vtkSlicerBigImageCoreModuleLogic*>(clientData);
  if (!self)
  {
    return;
  }

  if (event == vtkMRMLScene::EndImportEvent || event == vtkMRMLScene::EndRestoreEvent)
  {
    self->UpdateBigImagePlanes();
  }
}

//----------------------------------------------------------------------------
void vtkSlicerBigImageCoreModuleLogic::UpdateBigImagePlanes()
{
  vtkMRMLScene* scene = this->GetMRMLScene();
  if (!scene)
  {
    return;
  }

  vtkSmartPointer<vtkCollection> nodes = vtkSmartPointer<vtkCollection>::Take(
    scene->GetNodesByClass("vtkMRMLBigImageNode"));
  if (!nodes)
  {
    return;
  }

  vtkCollectionSimpleIterator it;
  nodes->InitTraversal(it);
  while (vtkMRMLNode* node = vtkMRMLNode::SafeDownCast(nodes->GetNextItemAsObject(it)))
  {
    vtkMRMLBigImageNode* bigImageNode = vtkMRMLBigImageNode::SafeDownCast(node);
    if (!bigImageNode)
    {
      continue;
    }

    vtkMRMLMarkupsPlaneNode* planeNode = vtkMRMLMarkupsPlaneNode::SafeDownCast(
      scene->GetNodeByID(bigImageNode->GetPosePlaneNodeID()));
    if (!planeNode)
    {
      continue;
    }

    int* dimensions = bigImageNode->GetDimensions();
    double* spacing = bigImageNode->GetSpacing();
    double width = dimensions[0] * spacing[0];
    double height = dimensions[1] * spacing[1];
    if (width <= 0.0 || height <= 0.0)
    {
      continue;
    }

    planeNode->SetSizeMode(vtkMRMLMarkupsPlaneNode::SizeModeAbsolute);
    planeNode->SetSize(width, height);

    vtkMRMLMarkupsPlaneDisplayNode* planeDisplayNode =
      vtkMRMLMarkupsPlaneDisplayNode::SafeDownCast(planeNode->GetDisplayNode());
    if (planeDisplayNode)
    {
      planeDisplayNode->SetScaleHandleVisibility(false);
      planeDisplayNode->SetVisibility2D(false);
    }

    vtkMRMLLinearTransformNode* transformNode = vtkMRMLLinearTransformNode::SafeDownCast(
      scene->GetNodeByID(bigImageNode->GetTransformNodeID()));
    if (!transformNode)
    {
      continue;
    }

    vtkSmartPointer<vtkMatrix4x4> transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
    transformNode->GetMatrixTransformToWorld(transformMatrix);

    double xAxis[3] = { transformMatrix->GetElement(0, 0),
                        transformMatrix->GetElement(1, 0),
                        transformMatrix->GetElement(2, 0) };
    double yAxis[3] = { transformMatrix->GetElement(0, 1),
                        transformMatrix->GetElement(1, 1),
                        transformMatrix->GetElement(2, 1) };
    double zAxis[3] = { transformMatrix->GetElement(0, 2),
                        transformMatrix->GetElement(1, 2),
                        transformMatrix->GetElement(2, 2) };
    if (vtkMath::Normalize(xAxis) == 0.0 ||
        vtkMath::Normalize(yAxis) == 0.0 ||
        vtkMath::Normalize(zAxis) == 0.0)
    {
      continue;
    }

    double origin[3] = { transformMatrix->GetElement(0, 3),
                         transformMatrix->GetElement(1, 3),
                         transformMatrix->GetElement(2, 3) };
    double center[3] = {
      origin[0] + 0.5 * width * xAxis[0] + 0.5 * height * yAxis[0],
      origin[1] + 0.5 * width * xAxis[1] + 0.5 * height * yAxis[1],
      origin[2] + 0.5 * width * xAxis[2] + 0.5 * height * yAxis[2]
    };

    planeNode->SetAxesWorld(xAxis, yAxis, zAxis);
    planeNode->SetOriginWorld(center);
  }

}
