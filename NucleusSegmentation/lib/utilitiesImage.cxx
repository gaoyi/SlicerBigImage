#include <csignal>
#include <string>

//#include <omp.h>

// itk
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkThresholdSegmentationLevelSetImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkValuedRegionalMinimaImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkZeroCrossingBasedEdgeDetectionImageFilter.h"
#include "itkZeroCrossingImageFilter.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"


#include "utilitiesImage.h"


namespace gth818n
{
  itkShortImage2DType::Pointer smoothContour(itkShortImage2DType::Pointer bwImage)
  {
    typedef itk::CastImageFilter< itkShortImage2DType, itkFloatImage2DType > CastFilterType;
    CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(bwImage);
    castFilter->Update();

    itkFloatImage2DType::Pointer itkImage = castFilter->GetOutput();
    /// Convert cv::Mat image to itk image
    ////////////////////////////////////////////////////////////////////////////////

    typedef itk::SignedDanielssonDistanceMapImageFilter< itkShortImage2DType, itkFloatImage2DType >  SignedDanielssonDistanceMapImageFilterType;
    SignedDanielssonDistanceMapImageFilterType::Pointer sdfFilter = SignedDanielssonDistanceMapImageFilterType::New();
    sdfFilter->SetInput( bwImage );
    sdfFilter->UseImageSpacingOn();
    sdfFilter->Update();

    typedef itk::ThresholdSegmentationLevelSetImageFilter< itkFloatImage2DType, itkFloatImage2DType > ThresholdSegmentationLevelSetImageFilterType;
    ThresholdSegmentationLevelSetImageFilterType::Pointer thresholdSegmentation = ThresholdSegmentationLevelSetImageFilterType::New();
    thresholdSegmentation->SetCurvatureScaling( 50.0 );
    thresholdSegmentation->SetPropagationScaling( 0.1 ); ///< this param and the SetCurvatureScaling one, as long as having the same ratio, (5, 1) and (50, 10) seem to be the same.
    thresholdSegmentation->SetMaximumRMSError( 1e-6 );
    thresholdSegmentation->SetNumberOfIterations( 100 ); ///< 1000 may be good
    thresholdSegmentation->SetUpperThreshold( 1 ); ///< change to 75 will significantly shrink to darker area
    thresholdSegmentation->SetLowerThreshold( 0 ); ///< change this to -10 does not affect result
    thresholdSegmentation->SetIsoSurfaceValue(0.5); ///< the output of sdfFilter has 0.5-isocontour the same as sdf's input binary
    thresholdSegmentation->SetInput( sdfFilter->GetOutput() );
    thresholdSegmentation->SetFeatureImage( itkImage );
    thresholdSegmentation->Update();

    typedef itk::BinaryThresholdImageFilter<itkFloatImage2DType, itkShortImage2DType> ThresholdingFilterType;
    ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
    thresholder->SetLowerThreshold( -1000.0 );
    thresholder->SetUpperThreshold(     0.0 );

    thresholder->SetOutsideValue(  0  );
    thresholder->SetInsideValue(  1 );
    thresholder->SetInput( thresholdSegmentation->GetOutput() );
    thresholder->Update();

    return thresholder->GetOutput();
  }


  void computeNucleiSizes(itkUIntImage2DType::Pointer itkLabelImage, itkUIntImage2DType::PixelType totalNumberOfObjects, std::vector<long>& sizeOfAllObjects, int magnification)
  {
    sizeOfAllObjects.resize(totalNumberOfObjects+1);

    const itkUIntImage2DType::PixelType* itkLabelImagePtr = static_cast<itkUIntImage2DType::PixelType*>(itkLabelImage->GetBufferPointer());

    unsigned long numberOfPixels = itkLabelImage->GetLargestPossibleRegion().GetNumberOfPixels();

    float magnificationAdjustRatio = static_cast<float>(magnification)/20.0; ///< So that the size will the equivalent to under the magnification of 20x

    for (unsigned long it = 0; it < numberOfPixels; ++it)
      {
        ++sizeOfAllObjects[itkLabelImagePtr[it]];
      }

    for (unsigned long it = 0; it < sizeOfAllObjects.size(); ++it)
      {
        sizeOfAllObjects[it] /= magnificationAdjustRatio;
      }

    return;
  }

  void computeNucleiFractionalAnisotrpy(itkUIntImage2DType::Pointer itkLabelImage, itkUIntImage2DType::PixelType totalNumberOfObjects, std::vector<float>& fractionalAnisotrpyOfAllObjects, int magnification)
  {
    /// go thru each pixel, put its coordinate to the corresponding list
    std::vector< std::vector<itk2DIndexType> > listOfCoordinatesOfEachObject(totalNumberOfObjects+1);

    const itkUIntImage2DType::PixelType* itkLabelImagePtr = static_cast<itkUIntImage2DType::PixelType*>(itkLabelImage->GetBufferPointer());

    unsigned long nx = itkLabelImage->GetLargestPossibleRegion().GetSize()[0];
    unsigned long ny = itkLabelImage->GetLargestPossibleRegion().GetSize()[1];

    {
      unsigned long it = 0;
      itk2DIndexType idx;

      for (unsigned long iy = 0; iy < ny; ++iy)
        {
          idx[1] = iy;
          for (unsigned long ix = 0; ix < nx; ++ix)
            {
              idx[0] = ix;
              itkUIntImage2DType::PixelType objectId = itkLabelImagePtr[it++];
              listOfCoordinatesOfEachObject[objectId].push_back(idx);
            }
        }
    }

    /// For each object, if it contains more than pixels in [10, 200],
    /// compute the fractional anisotrpoy. If we have 3 pixels, we can
    /// compute the co-variance matrix and then the eigen
    /// decomposition. But only 10 to 200 pixel area may be
    /// meaningful.
    std::vector<long> sizeOfAllObjects;
    computeNucleiSizes(itkLabelImage, totalNumberOfObjects, sizeOfAllObjects, magnification);

    fractionalAnisotrpyOfAllObjects.resize(totalNumberOfObjects+1);

    for (itkUIntImage2DType::PixelType itObject = 0; itObject < totalNumberOfObjects; ++itObject)
      {
        if (10 > sizeOfAllObjects[itObject] || 200 < sizeOfAllObjects[itObject])
          {
            fractionalAnisotrpyOfAllObjects[itObject] = -1;
            continue;
          }

        const std::vector<itk2DIndexType>& thePointsInThisObject = listOfCoordinatesOfEachObject[itObject];
        int numOfPointsInThisObject = thePointsInThisObject.size();

        float meanX = 0;
        float meanY = 0;

        for (int itPoint = 0; itPoint < numOfPointsInThisObject; ++itPoint)
          {
            meanX += thePointsInThisObject[itPoint][0];
            meanY += thePointsInThisObject[itPoint][1];
          }

        meanX /= static_cast<float>(numOfPointsInThisObject);
        meanY /= static_cast<float>(numOfPointsInThisObject);

        vnlFloatMatrixType M(2, numOfPointsInThisObject);
        for (int itPoint = 0; itPoint < numOfPointsInThisObject; ++itPoint)
          {
            M(0, itPoint) = thePointsInThisObject[itPoint][0] - meanX;
            M(1, itPoint) = thePointsInThisObject[itPoint][1] - meanY;
          }

        vnlFloatMatrixType covMatrix = M*( M.transpose() );

        vnl_symmetric_eigensystem<float> eig(covMatrix);
        float ev0 = eig.get_eigenvalue(0);
        float ev1 = eig.get_eigenvalue(1);

        float fa = fabs(ev0 - ev1)/sqrt(ev0*ev0 + ev1*ev1);

        fractionalAnisotrpyOfAllObjects[itObject] = fa;
      }

    return;
  }


  itkFloatImage2DType::Pointer bw2dist(itkShortImage2DType::Pointer bwImage)
  {
    typedef itk::SignedDanielssonDistanceMapImageFilter< itkShortImage2DType, itkFloatImage2DType >  SignedDanielssonDistanceMapImageFilterType;
    SignedDanielssonDistanceMapImageFilterType::Pointer sdfFilter = SignedDanielssonDistanceMapImageFilterType::New();
    sdfFilter->SetInput( bwImage );
    sdfFilter->UseImageSpacingOn();
    sdfFilter->Update();

    return sdfFilter->GetOutput();
  }

  itkUCharImage2DType::Pointer otsuThresholdImage(itkFloatImage2DType::Pointer image, itkUCharImage2DType::PixelType maskValue, float ratio)
  {
    typedef itk::OtsuThresholdImageFilter<itkFloatImage2DType, itkUCharImage2DType> FilterType;
    FilterType::Pointer otsuFilter = FilterType::New();
    otsuFilter->SetInput(image);
    otsuFilter->SetInsideValue(maskValue);
    otsuFilter->Update();

    itkUCharImage2DType::Pointer mask = otsuFilter->GetOutput();

    itkFloatImage2DType::PixelType otsuThld = otsuFilter->GetThreshold();

    otsuThld *= ratio;

    const itkFloatImage2DType::PixelType* imageBufferPointer = image->GetBufferPointer();
    itkUCharImage2DType::PixelType* maskBufferPointer = mask->GetBufferPointer();
    long numPixels = mask->GetLargestPossibleRegion().GetNumberOfPixels();
    for (long it = 0; it < numPixels; ++it)
      {
        maskBufferPointer[it] = imageBufferPointer[it]<=otsuThld?maskValue:0;
      }

    return mask;
  }


  itkShortImage2DType::Pointer thresholdHematoxylinImage(itkFloatImage2DType::Pointer image, itkShortImage2DType::PixelType maskValue)
  {
    typedef itk::OtsuThresholdImageFilter<itkFloatImage2DType, itkShortImage2DType> FilterType;
    FilterType::Pointer otsuFilter = FilterType::New();
    otsuFilter->SetInput(image);
    otsuFilter->SetInsideValue(maskValue);
    otsuFilter->Update();

    itkFloatImage2DType::PixelType otsuThreshold = otsuFilter->GetThreshold();
    itkFloatImage2DType::PixelType modifiedThreshold = otsuThreshold*0.8; ///< Otsu seems to be too high so include too much

    itkShortImage2DType::Pointer mask = otsuFilter->GetOutput();
    itkShortImage2DType::PixelType* maskBufferPointer = mask->GetBufferPointer();

    const itkFloatImage2DType::PixelType* imageBufferPointer = image->GetBufferPointer();

    itkShortImage2DType::SizeValueType numPixels = mask->GetLargestPossibleRegion().GetNumberOfPixels();
    for (itkShortImage2DType::SizeValueType it = 0; it < numPixels; ++it)
      {
        maskBufferPointer[it] = imageBufferPointer[it] < modifiedThreshold?1:0;
      }

    return mask;
  }

  itkFloatImage2DType::Pointer PeronaMalikDiffusion(itkFloatImage2DType::Pointer img, int numberOfIteration, double conductanceParameter)
  {
    typedef itk::GradientAnisotropicDiffusionImageFilter< itkFloatImage2DType, itkFloatImage2DType > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( img );
    filter->SetNumberOfIterations( numberOfIteration );
    filter->SetTimeStep( 0.125 );
    filter->SetConductanceParameter( conductanceParameter );
    filter->Update();

    return filter->GetOutput();
  }


  itkFloatImage2DType::Pointer CurvatureFlowSmoothing(itkFloatImage2DType::Pointer img, int numberOfIteration)
  {
    typedef itk::CurvatureFlowImageFilter< itkFloatImage2DType,itkFloatImage2DType >CurvatureFlowImageFilterType;
    CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();
    smoothing->SetInput( img );
    smoothing->SetNumberOfIterations( numberOfIteration );
    smoothing->SetTimeStep( 0.125 );
    smoothing->Update();

    return smoothing->GetOutput();
  }

  itkUIntImage2DType::Pointer binaryImageToConnectedComponentLabelImage(itkShortImage2DType::Pointer bwImage)
  {
    typedef itk::ConnectedComponentImageFilter <itkShortImage2DType, itkUIntImage2DType > ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
    connected->SetInput(bwImage);
    connected->Update();
    /// Connected component to separate each nucleus
    ////////////////////////////////////////////////////////////////////////////////

    //m_totalNumberOfConnectedComponents = connected->GetObjectCount();

    return connected->GetOutput();
  }

  itkUIntImage2DType::Pointer labelImageToSizeLabeledImage(itkUIntImage2DType::Pointer labelImage, itkUIntImage2DType::PixelType lowerThreshold)
  {
    typedef itk::MinimumMaximumImageCalculator<itkUIntImage2DType> ImageCalculatorFilterType;
    ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
    imageCalculatorFilter->SetImage(labelImage);
    imageCalculatorFilter->Compute();
    //std::cout<<imageCalculatorFilter->GetMinimum()<<"\t"<<imageCalculatorFilter->GetMaximum()<<std::endl;

    itkUIntImage2DType::PixelType totalNumberOfConnectedComponents = imageCalculatorFilter->GetMaximum();

    std::vector<itkUIntImage2DType::PixelType> sizesOfAllNuclei(totalNumberOfConnectedComponents + 1);



    itkUIntImage2DType::Pointer sizeValuedLabelImage = itkUIntImage2DType::New();
    sizeValuedLabelImage->SetRegions( labelImage->GetLargestPossibleRegion() );
    sizeValuedLabelImage->Allocate();
    sizeValuedLabelImage->FillBuffer(0);

    const itkUIntImage2DType::PixelType* labelImageOfCurrentTilePtr = static_cast<itkUIntImage2DType::PixelType*>(labelImage->GetBufferPointer());

    unsigned long numberOfPixels = sizeValuedLabelImage->GetLargestPossibleRegion().GetNumberOfPixels();

    for (unsigned long it = 0; it < numberOfPixels; ++it)
      {
        ++sizesOfAllNuclei[labelImageOfCurrentTilePtr[it]];
      }


    itkUIntImage2DType::PixelType* sizeValuedLabelImagePtr = static_cast<itkUIntImage2DType::PixelType*>(sizeValuedLabelImage->GetBufferPointer());
    for (unsigned long it = 0; it < numberOfPixels; ++it)
      {
        if (labelImageOfCurrentTilePtr[it])
          {
            itkUIntImage2DType::PixelType s = sizesOfAllNuclei[ labelImageOfCurrentTilePtr[it] ];

            if (s >= lowerThreshold)
              {
                sizeValuedLabelImagePtr[it] = s;
              }
            else
              {
                sizeValuedLabelImagePtr[it] = 0;
              }
          }
      }

    return sizeValuedLabelImage;
  }


  itkUCharImage2DType::Pointer edgesOfDifferentLabelRegion(itkUIntImage2DType::Pointer labelImage)
  {
    typedef itk::GradientMagnitudeImageFilter<itkUIntImage2DType, itkFloatImage2DType >  GradientMagnitudeImageFilterType;
    GradientMagnitudeImageFilterType::Pointer gradientMagnitudeImageFilter = GradientMagnitudeImageFilterType::New();
    gradientMagnitudeImageFilter->SetInput(labelImage);
    gradientMagnitudeImageFilter->Update();

    itkFloatImage2DType::Pointer gradImage = gradientMagnitudeImageFilter->GetOutput();
    itkFloatImage2DType::PixelType* gradImageBufferPointer = gradImage->GetBufferPointer();

    itkUCharImage2DType::Pointer mask = itkUCharImage2DType::New();
    mask->SetRegions( gradImage->GetLargestPossibleRegion() );
    mask->Allocate();
    mask->FillBuffer(0);

    itkUCharImage2DType::PixelType* maskBufferPointer = mask->GetBufferPointer();

    long numPixels = mask->GetLargestPossibleRegion().GetNumberOfPixels();

    for (long it = 0; it < numPixels; ++it)
      {
        if (gradImageBufferPointer[it] >= 0.5)
          {
            maskBufferPointer[it] = 1;
          }
      }


    return mask;
  }

  itkFloatImage2DType::Pointer GaussianSmoothing(itkFloatImage2DType::Pointer img, float variance)
  {
    typedef itk::DiscreteGaussianImageFilter<itkFloatImage2DType, itkFloatImage2DType >  filterType;
    filterType::Pointer gaussianFilter = filterType::New();
    gaussianFilter->SetInput( img );
    gaussianFilter->SetVariance(variance);
    gaussianFilter->Update();

    return gaussianFilter->GetOutput();
  }

  itkUCharImage2DType::Pointer fillHole(const itkUCharImage2DType* maskWithHoles)
  {
    typedef itk::BinaryFillholeImageFilter< itkUCharImage2DType > FilterType;
    FilterType::Pointer filler = FilterType::New();
    filler->SetFullyConnected( 1 );
    filler->SetForegroundValue( 1 );
    filler->SetInput(maskWithHoles);
    filler->Update();

    return filler->GetOutput();
  }

}// namespace
