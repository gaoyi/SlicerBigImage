#ifndef MeanshiftImageClusteringFilter_hxx_
#define MeanshiftImageClusteringFilter_hxx_

#include "MeanshiftImageClusteringFilter.h"

//dbg
#include "Image/gth818nImage.h"
//dbg, end

namespace gth818n
{

  template< typename TInputImage, typename TOutputImage >
  MeanshiftImageClusteringFilter< TInputImage, TOutputImage >::MeanshiftImageClusteringFilter()
  {
    m_inputImage = 0;
    m_inputLabelImage = 0;
    m_outputImage = 0;

    m_inputLabelOfInterest = 0;
    m_inputLabelOfInterestSet = false;


    m_allDone = false;
  }


  template< typename TInputImage, typename TOutputImage >
  void MeanshiftImageClusteringFilter< TInputImage, TOutputImage >::update()
  {
    _constructFeatureVectors();
    _doMeanShiftClustering();
    _meanshiftResultToLabelImage();


    m_allDone = true;

    return;
  }

  template< typename TInputImage, typename TOutputImage >
  void MeanshiftImageClusteringFilter< TInputImage, TOutputImage >::_constructFeatureVectors()
  {
    if (!m_inputLabelImage)
      {
        /// clustering on entire image
      }
    else if (!m_inputLabelOfInterestSet)
      {
        /// clustering on region with non-zero labels
      }
    else
      {
        /// clustering on region with label == m_inputLabelOfInterest
      }


  TODO: populate these two:
    typename SampleType::Pointer m_samples;
    std::vector<IndexType> m_indexOfSamplePoints;





    return;
  }

  template< typename TInputImage, typename TOutputImage >
  void MeanshiftImageClusteringFilter< TInputImage, TOutputImage >::_doMeanShiftClustering()

    return;
  }

  template< typename TInputImage, typename TOutputImage >
  void MeanshiftImageClusteringFilter< TInputImage, TOutputImage >::_meanshiftResultToLabelImage()

    return;
  }


  template< typename TInputImage, typename TOutputImage >
  typename MeanshiftImageClusteringFilter< TInputImage, TOutputImage >::OutputImageType::Pointer
  MeanshiftImageClusteringFilter< TInputImage, TOutputImage >::getOutputImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_outputImage;
  }





}// namespace gth818n


#endif
