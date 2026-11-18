#ifndef __vtkSlicerBigImageCoreModuleLogic_h
#define __vtkSlicerBigImageCoreModuleLogic_h

#include "vtkSlicerBigImageCoreModuleLogicExport.h"
#include "vtkSlicerModuleLogic.h"

class vtkCallbackCommand;

class VTK_SLICER_BIGIMAGECORE_MODULE_LOGIC_EXPORT vtkSlicerBigImageCoreModuleLogic : public vtkSlicerModuleLogic
{
public:
  static vtkSlicerBigImageCoreModuleLogic* New();
  vtkTypeMacro(vtkSlicerBigImageCoreModuleLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  static void OnSceneEvents(vtkObject* caller, unsigned long event, void* clientData, void* callData);

protected:
  vtkSlicerBigImageCoreModuleLogic();
  ~vtkSlicerBigImageCoreModuleLogic() override;

  void RegisterNodes() override;

  void UpdateBigImagePlanes();

  vtkCallbackCommand* SceneEventsCallback{ nullptr };

private:
  vtkSlicerBigImageCoreModuleLogic(const vtkSlicerBigImageCoreModuleLogic&) = delete;
  void operator=(const vtkSlicerBigImageCoreModuleLogic&) = delete;
};

#endif
