#ifndef utilitiesImage_h_
#define utilitiesImage_h_


#include <string>

// itk
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"

// local
#include "itkTypedefs.h"

namespace gth818n
{
  void computeNucleiSizes(itkUIntImage2DType::Pointer itkLabelImage, itkUIntImage2DType::PixelType totalNumberOfObjects, std::vector<long>& sizeOfAllObjects, int magnification);
  void computeNucleiFractionalAnisotrpy(itkUIntImage2DType::Pointer itkLabelImage, itkUIntImage2DType::PixelType totalNumberOfObjects, std::vector<float>& fractionalAnisotrpyOfAllObjects, int magnification);

  itkFloatImage2DType::Pointer bw2dist(itkShortImage2DType::Pointer bwImage);

  template<typename InputImageType, typename OutputImageType>
  typename OutputImageType::Pointer
  otsuThresholdImageT(const InputImageType* image, typename OutputImageType::PixelType maskValue = 1)
  {
    typedef itk::OtsuThresholdImageFilter<InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer otsuFilter = FilterType::New();
    otsuFilter->SetInput(image);
    otsuFilter->SetInsideValue(maskValue);
    otsuFilter->Update();

    return otsuFilter->GetOutput();
  }

  itkUCharImage2DType::Pointer otsuThresholdImage(itkFloatImage2DType::Pointer image, itkUCharImage2DType::PixelType maskValue = 1, float ratio = 1.0);

  itkShortImage2DType::Pointer thresholdHematoxylinImage(itkFloatImage2DType::Pointer image, itkShortImage2DType::PixelType maskValue = 1);

  itkShortImage2DType::Pointer smoothContour(itkShortImage2DType::Pointer bwImage);

  itkUIntImage2DType::Pointer binaryImageToConnectedComponentLabelImage(itkShortImage2DType::Pointer bwImage);

  itkUIntImage2DType::Pointer labelImageToSizeLabeledImage(itkUIntImage2DType::Pointer labelImage, itkUIntImage2DType::PixelType lowerThreshold = 3); ///< if an area is smaller than lowerThreshold pixels, discard

  itkFloatImage2DType::Pointer PeronaMalikDiffusion(itkFloatImage2DType::Pointer img, int numberOfIteration = 10, double conductanceParameter = 0.5);
  itkFloatImage2DType::Pointer CurvatureFlowSmoothing(itkFloatImage2DType::Pointer img, int numberOfIteration = 10);
  itkFloatImage2DType::Pointer GaussianSmoothing(itkFloatImage2DType::Pointer img, float variance = 3);

  itkUCharImage2DType::Pointer edgesOfDifferentLabelRegion(itkUIntImage2DType::Pointer labelImage);

  itkUCharImage2DType::Pointer fillHole(const itkUCharImage2DType* maskWithHoles);

  template< typename InputImageType, typename OutputImageType >
  typename OutputImageType::Pointer
  castItkImage( const InputImageType* inputImage )
  {
    typedef itk::CastImageFilter< InputImageType, OutputImageType > itkCastFilter_t;

    typename itkCastFilter_t::Pointer caster = itkCastFilter_t::New();
    caster->SetInput( inputImage );
    caster->Update();

    return caster->GetOutput();
  }

}// namespace gth818n


#endif
