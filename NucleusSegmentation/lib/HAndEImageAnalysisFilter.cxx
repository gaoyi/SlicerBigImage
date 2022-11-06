#include <vector>

// itk
#include "itkImage.h"
#include "itkConnectedComponentImageFilter.h"

// local
#include "utilitiesImage.h"
#include "HAndEImageAnalysisFilter.h"
#include "SFLSLocalChanVeseSegmentor2D.h"
#include "ColorDecomposition.h"
#include "HematoxylinImageSegmentationFilter.h"

namespace gth818n
{
  void HAndEImageAnalysisFilter::update()
  {
    _extractHematoxylinImage();

    _segmentHematoxylinImage();

    _computeNucleiLabelImage();

    m_allDone = true;

    return;
  }

  void HAndEImageAnalysisFilter::setMPP(float mpp)
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

  void HAndEImageAnalysisFilter::_segmentHematoxylinImage()
  {
    HematoxylinImageSegmentationFilter imageAnalyzer;
    imageAnalyzer.setHematoxylinImage( m_hematoxylinImage );
    imageAnalyzer.setNucleusSizeThreshold(m_nucleusSizeThreshold);
    imageAnalyzer.setMPP(m_mpp);
    imageAnalyzer.update();

    m_nucleusBinaryMask = imageAnalyzer.getNucleiBinaryMaskImage();

    return;
  }

  void HAndEImageAnalysisFilter::_computeNucleiLabelImage()
  {
    ////////////////////////////////////////////////////////////////////////////////
    /// Connected component to separate each nucleus
    ///
    /// The max of <unsigned int>, as tested in
    /// ../test/mainNumericTraits.cxx, is 4G. Should be large enough for
    /// a tile. Enough for a WSI?
    typedef itk::ConnectedComponentImageFilter <itkUCharImage2DType, itkUIntImage2DType > ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
    connected->SetInput(m_nucleusBinaryMask);
    connected->Update();
    /// Connected component to separate each nucleus
    ////////////////////////////////////////////////////////////////////////////////

    m_totalNumberOfConnectedComponents = connected->GetObjectCount();

    m_nucleiLabelImage = connected->GetOutput();

    return;
  }

  HAndEImageAnalysisFilter::HAndEImageAnalysisFilter()
  {
    m_nucleusSizeThreshold = 0;

    m_mpp = 0.25;

    m_allDone = false;

    m_threshold = 50; ///< 75 is small, miss a lot of nuclei. But of course the false positive is low. To make the false positive even lower, let me use 50
    //m_threshold = 150; ///< 150 till now over seg all luad tiles. But i think it safer to put to 175 if we want small false negative

    return;
  }


  HAndEImageAnalysisFilter::itkIntImage2DType::Pointer
  HAndEImageAnalysisFilter::getNucleiLabelImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_nucleiLabelImage;
  }

  HAndEImageAnalysisFilter::itkBinaryMaskImage2DType::Pointer
  HAndEImageAnalysisFilter::getNucleiBinaryMaskImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_nucleusBinaryMask;
  }

  std::vector<long>
  HAndEImageAnalysisFilter::getSizesOfAllNuclei()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_sizesOfAllNuclei;
  }


  std::vector<float>
  HAndEImageAnalysisFilter::getRoundnessOfAllNuclei()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_roundnessOfAllNuclei;
  }

  int64_t
  HAndEImageAnalysisFilter::getTotalNumberOfConnectedComponents()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_totalNumberOfConnectedComponents;
  }

  void
  HAndEImageAnalysisFilter::_extractHematoxylinImage()
  {
    m_hematoxylinImage = gth818n::ExtractHematoxylinChannel(m_HAndEImage);

    return;
  }

}// namespace gth818n
