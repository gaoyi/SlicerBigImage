#ifndef BinaryMaskWatershedClusteringFilter_h_
#define BinaryMaskWatershedClusteringFilter_h_


// itk
#include "itkImage.h"
#include "itkWatershedImageFilter.h"


namespace gth818n
{
  template< typename TInputImage, typename TOutputImage >
  class BinaryMaskWatershedClusteringFilter
  {
  public:
    typedef BinaryMaskWatershedClusteringFilter Self;

    ////////////////////////////////////////////////////////////////////////////////
    /// ctor
    BinaryMaskWatershedClusteringFilter();
    ~BinaryMaskWatershedClusteringFilter() {}
    /// ctor, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// typedef
    typedef TInputImage InputImageType;
    typedef TOutputImage OutputImageType;

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    typedef itk::Image<unsigned char, InputImageType::ImageDimension> UCharImageType;
    typedef itk::Image<float, InputImageType::ImageDimension> FloatImageType;

    typedef typename InputImageType::IndexType IndexType;
    /// typedef, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// public fn
    void setInputImage(typename InputImageType::Pointer img) {m_inputImage = img;}
    typename OutputImageType::Pointer getOutputImage();

    void setThreshold(double thld) {m_threshold = thld;}
    void setLevel(double level) {m_level = level;}

    void update();
    /// public fn, end
    ////////////////////////////////////////////////////////////////////////////////


  private:

    typedef itk::WatershedImageFilter<FloatImageType> WatershedFilterType;
    typedef typename WatershedFilterType::OutputImageType WSOutputImageType;


    ////////////////////////////////////////////////////////////////////////////////
    /// private data
    typename InputImageType::Pointer m_inputImage;
    typename UCharImageType::Pointer m_binaryMaskOfInputImage;

    typename OutputImageType::Pointer m_outputImage;

    unsigned int m_numberOfObjects; ///< I will use "Object" as well as "Connected Component"

    double m_threshold;
    double m_level;

    bool m_allDone;
    /// private data, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// private fn
    void _computeBinaryMask();
    void _watershedOnSignedDistanceMap();
    /// private fn
    ////////////////////////////////////////////////////////////////////////////////
  };

}// namespace gth818n

#include "BinaryMaskWatershedClusteringFilter.hxx"

#endif // BinaryMaskWatershedClusteringFilter_h_
