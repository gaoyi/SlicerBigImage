#ifndef GrayImageWatershedSegmentationFilter_h_
#define GrayImageWatershedSegmentationFilter_h_


// itk
#include "itkImage.h"
#include "itkWatershedImageFilter.h"


namespace gth818n
{
  template< typename TInputImage, typename TOutputImage >
  class GrayImageWatershedSegmentationFilter
  {
  public:
    typedef GrayImageWatershedSegmentationFilter Self;

    ////////////////////////////////////////////////////////////////////////////////
    /// ctor
    GrayImageWatershedSegmentationFilter();
    ~GrayImageWatershedSegmentationFilter() {}
    /// ctor, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// typedef
    typedef TInputImage InputImageType;
    typedef TOutputImage OutputImageType;

    typedef typename InputImageType::PixelType InputPixelType;
    typedef typename OutputImageType::PixelType OutputPixelType;
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
    ////////////////////////////////////////////////////////////////////////////////
    /// private data
    typename InputImageType::Pointer m_inputImage;
    typename OutputImageType::Pointer m_outputImage;

    unsigned int m_numberOfObjects; ///< I will use "Object" as well as "Connected Component"

    double m_threshold;
    double m_level;

    bool m_allDone;
    /// private data, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// private fn
    void _watershed();
    /// private fn
    ////////////////////////////////////////////////////////////////////////////////
  };

}// namespace gth818n

#include "GrayImageWatershedSegmentationFilter.hxx"

#endif // GrayImageWatershedSegmentationFilter_h_
