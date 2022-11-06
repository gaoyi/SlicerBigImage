#ifndef HematoxylinImageSegmentationFilter_h_
#define HematoxylinImageSegmentationFilter_h_

// itk
#include "itkImage.h"
#include "itkRGBPixel.h"

namespace gth818n
{
  class HematoxylinImageSegmentationFilter
  {
  public:
    //typedef HematoxylinImageSegmentationFilter Self;

    ////////////////////////////////////////////////////////////////////////////////
    /// ctor
    HematoxylinImageSegmentationFilter();
    ~HematoxylinImageSegmentationFilter() {}
    /// ctor, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// typedef
    static const unsigned int ImageDimension = 2;
    typedef itk::Image<unsigned char, ImageDimension> itkUCharImage2DType;
    typedef itkUCharImage2DType itkBinaryMaskImage2DType;
    /// typedef, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// public fn
    void setHematoxylinImage(const itkUCharImage2DType* hematoxylinImage) {m_hematoxylinImage = hematoxylinImage;}

    void setNucleusSizeThreshold(float sizeThld) {m_nucleusSizeThreshold = sizeThld;}

    void setMPP(float mpp);

    void update();

    itkBinaryMaskImage2DType::Pointer getNucleiBinaryMaskImage();
    /// public fn, end
    ////////////////////////////////////////////////////////////////////////////////


  private:
    typedef itk::Image<float, ImageDimension> itkFloatImage2DType;
    typedef itk::Index<ImageDimension> itk2DIndexType;
    ////////////////////////////////////////////////////////////////////////////////
    /// The max of <unsigned int>, as tested in
    /// ../test/mainNumericTraits.cxx, is 4G. Should be large enough for
    /// representing the label of each nucleus in a tile. Enough for a
    /// WSI?
    ///
    /// typedef unsigned int LabelType;
    /// typedef itk::Image< LabelType, Dimension > LabelImageType;
    typedef itk::Image<unsigned long, ImageDimension> itkULongImageType;
    typedef itk::Image<unsigned int, ImageDimension> itkUIntImage2DType;
    typedef itk::Image<int, ImageDimension> itkIntImageType;
    typedef itkUIntImage2DType itkIntImage2DType;

    ////////////////////////////////////////////////////////////////////////////////
    /// private data
    const itkUCharImage2DType* m_hematoxylinImage;
    itkBinaryMaskImage2DType::Pointer m_nucleusBinaryMask;

    float m_nucleusSizeThreshold; ///< region smaller than this will be removed, unit: um^2
    float m_mpp; ///< Micron Per Pixel

    bool m_allDone;
    /// private data, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// private fn
    void _segmentHematoxylinImage();
    void _segmentHematoxylinImage_otsu_ChanVese_removeSmallIsland();
    void _computeNucleiLabelImage();
    /// private fn
    ////////////////////////////////////////////////////////////////////////////////
  };

}// namespace gth818n


#endif // HematoxylinImageSegmentationFilter_h_
