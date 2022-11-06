// itk
#include "itkImage.h"
#include "itkWatershedImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"

// local
// #include "itkTypedefs.h"
#include "utilitiesImage.h"

#include "HematoxylinImageSegmentationFilter.h"

#include "../include/SFLSLocalChanVeseSegmentor2D.h"

#include "itkMinimumMaximumImageCalculator.h"

#include "BinaryMaskAnalysisFilter.h"


namespace gth818n
{
  HematoxylinImageSegmentationFilter::HematoxylinImageSegmentationFilter()
  {
    m_nucleusSizeThreshold = 0;

    m_mpp = 0.25;

    m_allDone = false;

    return;
  }


  void HematoxylinImageSegmentationFilter::update()
  {
    _segmentHematoxylinImage();

    m_allDone = true;

    return;
  }

  void HematoxylinImageSegmentationFilter::_segmentHematoxylinImage()
  {

    //_segmentHematoxylinImage_otsu();

    _segmentHematoxylinImage_otsu_ChanVese_removeSmallIsland();

    return;
  }

  void HematoxylinImageSegmentationFilter::setMPP(float mpp)
  {
    if (mpp > 0)
      {
        m_mpp = mpp;
      }
    else
      {
        std::cerr<<"Error: mpp should be > 0. But got "<<mpp<<std::endl;
      }

    return;
  }

  void HematoxylinImageSegmentationFilter::_segmentHematoxylinImage_otsu()
  {

    // int numberOfIteration = 10;
    // double conductanceParameter = 0.5;
    // itkFloatImage2DType::Pointer image = PeronaMalikDiffusion(castItkImage<itkUCharImage2DType, itkFloatImage2DType>(m_hematoxylinImage), numberOfIteration, conductanceParameter);

    // // typedef itk::GradientRecursiveGaussianImageFilter<itkFloatImage2DType, itkFloatImage2DType> GradientRecursiveGaussianImageFilterType;
    // // GradientRecursiveGaussianImageFilterType::Pointer gradfilter = GradientRecursiveGaussianImageFilterType::New();
    // // //sigma is specified in millimeters
    // // gradfilter->SetSigma( 1.5 );

    // {
    //   typedef itk::GradientMagnitudeImageFilter<itkFloatImage2DType, itkFloatImage2DType >  GradientMagnitudeImageFilterType;
    //   GradientMagnitudeImageFilterType::Pointer gradientMagnitudeImageFilter = GradientMagnitudeImageFilterType::New();
    //   gradientMagnitudeImageFilter->SetInput(image);
    //   gradientMagnitudeImageFilter->Update();
    //   image = gradientMagnitudeImageFilter->GetOutput();
    // }

    // double threshold = 0.00;
    // double level = 0.2;

    // itkUIntImage2DType::Pointer waterShedResults;
    // {
    //   typedef itk::WatershedImageFilter<itkFloatImage2DType> WatershedFilterType;
    //   WatershedFilterType::Pointer watershed = WatershedFilterType::New();
    //   watershed->SetThreshold(threshold);
    //   watershed->SetLevel(level);
    //   watershed->SetInput(image);
    //   watershed->Update();
    //   waterShedResults = castItkImage<itkULongImageType, itkUIntImage2DType>(watershed->GetOutput());
    // }

    // m_nucleusBinaryMask = gth818n::edgesOfDifferentLabelRegion(waterShedResults);
    // itkUCharImage2DType::PixelType* m_nucleusBinaryMaskBufferPointer = m_nucleusBinaryMask->GetBufferPointer();

    // itkShortImage2DType::PixelType* maskBufferPointer = mask->GetBufferPointer();

    // for (long it = 0; it < mask->GetLargestPossibleRegion().GetNumberOfPixels(); ++it )
    //   {
    //     if (maskBufferPointer[it] == 1 && m_nucleusBinaryMaskBufferPointer[it] != 0)
    //       {
    //         m_nucleusBinaryMaskBufferPointer[it] = 1;
    //       }
    //     else
    //       {
    //         m_nucleusBinaryMaskBufferPointer[it] = 0;
    //       }
    //   }


    m_nucleusBinaryMask = otsuThresholdImageT<itkUCharImage2DType, itkUCharImage2DType>(m_hematoxylinImage);
    m_nucleusBinaryMask = fillHole(m_nucleusBinaryMask);

    return;
  }

  void HematoxylinImageSegmentationFilter::_segmentHematoxylinImage_otsu_ChanVese_removeSmallIsland()
  {
    short maskValue = 1;
    float otsuRatio = 0.7;
    itkFloatImage2DType::Pointer hemaFlow = castItkImage<itkUCharImage2DType, itkFloatImage2DType>(m_hematoxylinImage);

    m_nucleusBinaryMask = otsuThresholdImage(hemaFlow, maskValue, otsuRatio);
    long numPixels = m_nucleusBinaryMask->GetLargestPossibleRegion().GetNumberOfPixels();

    /// Check if the current mask is all zero, if so return coz that will cause problem in CV (plus there is no point to run Chan Vese
    {
      bool imageIsAllZero = true;
      const itkBinaryMaskImage2DType::PixelType* nucleusBinaryMaskBufferPointer = m_nucleusBinaryMask->GetBufferPointer();
      for (long it = 0; it < numPixels; ++it)
        {
          if (nucleusBinaryMaskBufferPointer[it] >= 1)
            {
              imageIsAllZero = false;
              break;
            }
        }

      if (imageIsAllZero)
        {
          return;
        }
    }

    m_nucleusBinaryMask = gth818n::fillHole(m_nucleusBinaryMask);

    int numiter = 100;
    CSFLSLocalChanVeseSegmentor2D< itkFloatImage2DType::PixelType > cv;
    cv.setImage(hemaFlow);
    cv.setMask( m_nucleusBinaryMask );
    cv.setNumIter(numiter);
    float lambda = 0.8;
    cv.setCurvatureWeight(lambda);
    cv.doSegmenation();

    CSFLSLocalChanVeseSegmentor2D< itkFloatImage2DType::PixelType >::LSImageType::Pointer phi = cv.mp_phi;

    itkUCharImage2DType::PixelType* nucleusBinaryMaskBufferPointer = m_nucleusBinaryMask->GetBufferPointer();
    CSFLSLocalChanVeseSegmentor2D< itkFloatImage2DType::PixelType >::LSImageType::PixelType* phiBufferPointer = phi->GetBufferPointer();

    for (long it = 0; it < numPixels; ++it)
      {
        nucleusBinaryMaskBufferPointer[it] = phiBufferPointer[it]<=1.0?1:0;
      }


    gth818n::BinaryMaskAnalysisFilter binaryMaskAnalyzer;
    binaryMaskAnalyzer.setMaskImage( m_nucleusBinaryMask );
    binaryMaskAnalyzer.setMPP(m_mpp);
    binaryMaskAnalyzer.update();

    // //gth818n::writeImage<gth818n::BinaryMaskAnalysisFilter::itkFloatImage2DType>(binaryMaskAnalyzer.getSizeColoredImage(), "color.nrrd" );
    // gth818n::writeImage<gth818n::BinaryMaskAnalysisFilter::itkFloatImage2DType>(binaryMaskAnalyzer.getFeatureColoredImage(2), "ratio1.nrrd" );
    // //gth818n::writeImage<gth818n::BinaryMaskAnalysisFilter::itkFloatImage2DType>(binaryMaskAnalyzer.getFeatureColoredImage(1), "color1.nrrd" );
    //gth818n::writeImage<gth818n::BinaryMaskAnalysisFilter::itkIntImage2DType>(binaryMaskAnalyzer.getConnectedComponentLabelImage(), outputImageName.c_str() );

    itkUIntImage2DType::Pointer sizeLabel = binaryMaskAnalyzer.getConnectedComponentLabelImage();
    itkUCharImage2DType::Pointer edgeBetweenLabelsMask = edgesOfDifferentLabelRegion(gth818n::castItkImage<itkUIntImage2DType, itkUIntImage2DType>(binaryMaskAnalyzer.getConnectedComponentLabelImage()));
    itkUCharImage2DType::PixelType* edgeBetweenLabelsMaskBufferPointer = edgeBetweenLabelsMask->GetBufferPointer();

    // itkUIntImage2DType::Pointer ccLabel = gth818n::binaryImageToConnectedComponentLabelImage( gth818n::castItkImage<itkUCharImage2DType, itkShortImage2DType>(m_nucleusBinaryMask) );
    // itkUIntImage2DType::Pointer sizeLabel = gth818n::labelImageToSizeLabeledImage(ccLabel, m_nucleusSizeThreshold); ///< if an area is smaller than lowerThreshold pixels, discard

    const itkUIntImage2DType::PixelType* sizeLabelBufferPointer = sizeLabel->GetBufferPointer();
    for (long it = 0; it < numPixels; ++it)
      {
        nucleusBinaryMaskBufferPointer[it] = sizeLabelBufferPointer[it] >= 1?1:0;
        nucleusBinaryMaskBufferPointer[it] *= (1 - edgeBetweenLabelsMaskBufferPointer[it]);
      }

    return;
  }


  HematoxylinImageSegmentationFilter::itkBinaryMaskImage2DType::Pointer
  HematoxylinImageSegmentationFilter::getNucleiBinaryMaskImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not done.\n";
        abort();
      }

    return m_nucleusBinaryMask;
  }



}// namespace gth818n
