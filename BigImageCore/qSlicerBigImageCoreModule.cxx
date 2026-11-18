#include "qSlicerBigImageCoreModule.h"

#include "vtkSlicerBigImageCoreModuleLogic.h"

#include <vtkMRMLAbstractLogic.h>

//----------------------------------------------------------------------------
qSlicerBigImageCoreModule::qSlicerBigImageCoreModule(QObject* parent)
  : Superclass(parent)
{
}

//----------------------------------------------------------------------------
qSlicerBigImageCoreModule::~qSlicerBigImageCoreModule() = default;

//----------------------------------------------------------------------------
QString qSlicerBigImageCoreModule::helpText() const
{
  return tr("Core MRML node types for SlicerBigImage.");
}

//----------------------------------------------------------------------------
QString qSlicerBigImageCoreModule::acknowledgementText() const
{
  return QString();
}

//----------------------------------------------------------------------------
QStringList qSlicerBigImageCoreModule::contributors() const
{
  return QStringList() << QString("Yi Gao (Shenzhen University)");
}

//----------------------------------------------------------------------------
QIcon qSlicerBigImageCoreModule::icon() const
{
  return QIcon();
}

//----------------------------------------------------------------------------
QStringList qSlicerBigImageCoreModule::categories() const
{
  return QStringList() << QString("BigImage");
}

//----------------------------------------------------------------------------
QStringList qSlicerBigImageCoreModule::dependencies() const
{
  return QStringList();
}

//----------------------------------------------------------------------------
bool qSlicerBigImageCoreModule::isHidden() const
{
  return true;
}

//----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerBigImageCoreModule::createWidgetRepresentation()
{
  return nullptr;
}

//----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerBigImageCoreModule::createLogic()
{
  return vtkSlicerBigImageCoreModuleLogic::New();
}
