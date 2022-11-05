#ifndef MeanshiftImageClusteringFilter_h_
#define MeanshiftImageClusteringFilter_h_


// itk
#include "itkVector.h"
#include "itkListSample.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistanceMetric.h"

// local
#include "HierarchicalMeanshiftClusteringFilter.h"

namespace gth818n
{
  template< typename TInputImage, typename TInputLabelImage, typename TOutputImage >
  class MeanshiftImageClusteringFilter
  {
  public:
    /*--------------------------------------------------------------------------------*/
    /// typedef

    /** Standard class typedefs. */
    typedef MeanshiftImageClusteringFilter Self;

    typedef TInputImage InputImageType;
    typedef TInputLabelImage InputLabelImageType;
    typedef TOutputImage OutputImageType;

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename InputLabelImageType::PixelType InputLabelPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    typedef typename InputImageType::IndexType IndexType;
    /// typedef, end
    /*--------------------------------------------------------------------------------*/


    /*--------------------------------------------------------------------------------*/
    /// ctor
    MeanshiftImageClusteringFilter();
    ~MeanshiftImageClusteringFilter() {}
    /// ctor, end
    /*--------------------------------------------------------------------------------*/


    /*--------------------------------------------------------------------------------*/
    /// public fn
    /** Get the dimension (size) of the point. */
    void setInputImage(typename InputImageType::Pointer img) {m_inputImage = img;}
    void setLabelImage(typename InputLabelImageType::Pointer l) {m_inputLabelImage = l;}

    void setLabelOfInterest(InputLabelPixelType l)
    {
      m_inputLabelOfInterest = l;
      m_inputLabelOfInterestSet = true;
    }

    typename OutputImageType::Pointer getOutputImage();

    void update();
    /// public fn, end
    /*--------------------------------------------------------------------------------*/

  private:
    static const unsigned int FeatureDimension = TInputImage::ImageDimension + 1; ///< 1 for scalar intensity

    typename InputImageType::Pointer m_inputImage;
    typename InputLabelImageType::Pointer m_inputLabelImage;
    typename OutputImageType::Pointer m_outputImage;

    typedef gth818n::HierarchicalMeanshiftClusteringFilter<float, FeatureDimension> MeanshiftClusteringFilterType;
    typedef typename MeanshiftClusteringFilterType::VectorType VectorType;
    typedef typename itk::Statistics::ListSample< VectorType > SampleType;


    /*--------------------------------------------------------------------------------*/
    /// private data
    InputLabelPixelType m_inputLabelOfInterest;
    bool m_inputLabelOfInterestSet;

    typename SampleType::Pointer m_samples;
    std::vector<IndexType> m_indexOfSamplePoints;


    bool m_allDone;
    /// private data, end
    /*--------------------------------------------------------------------------------*/


    /*--------------------------------------------------------------------------------*/
    /// private fn
    void _constructFeatureVectors();
    void _doMeanShiftClustering();
    void _meanshiftResultToLabelImage();
    /// private fn
    /*--------------------------------------------------------------------------------*/
  };

}// namespace gth818n

#include "MeanshiftImageClusteringFilter.hxx"

#endif // MeanshiftImageClusteringFilter_h_
