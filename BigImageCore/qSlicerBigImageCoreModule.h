#ifndef __qSlicerBigImageCoreModule_h
#define __qSlicerBigImageCoreModule_h

#include "qSlicerLoadableModule.h"
#include "qSlicerBigImageCoreModuleExport.h"

class Q_SLICER_QTMODULES_BIGIMAGECORE_EXPORT qSlicerBigImageCoreModule : public qSlicerLoadableModule
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "org.slicer.modules.loadable.qSlicerLoadableModule/1.0");
  Q_INTERFACES(qSlicerLoadableModule);

public:
  typedef qSlicerLoadableModule Superclass;
  explicit qSlicerBigImageCoreModule(QObject* parent = nullptr);
  ~qSlicerBigImageCoreModule() override;

  qSlicerGetTitleMacro(tr("BigImage Core"));

  QString helpText() const override;
  QString acknowledgementText() const override;
  QStringList contributors() const override;
  QIcon icon() const override;
  QStringList categories() const override;
  QStringList dependencies() const override;
  bool isHidden() const override;

protected:
  qSlicerAbstractModuleRepresentation* createWidgetRepresentation() override;
  vtkMRMLAbstractLogic* createLogic() override;
};

#endif
