// itk
#include "itkImage.h"

// local
#include "utilitiesImage.h"
#include "HematoxylinImageSegmentationFilter.h"
#include "SFLSLocalChanVeseSegmentor2D.h"
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

    itkUIntImage2DType::Pointer sizeLabel = binaryMaskAnalyzer.getConnectedComponentLabelImage();
    itkUCharImage2DType::Pointer edgeBetweenLabelsMask = edgesOfDifferentLabelRegion(gth818n::castItkImage<itkUIntImage2DType, itkUIntImage2DType>(binaryMaskAnalyzer.getConnectedComponentLabelImage()));
    itkUCharImage2DType::PixelType* edgeBetweenLabelsMaskBufferPointer = edgeBetweenLabelsMask->GetBufferPointer();

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
