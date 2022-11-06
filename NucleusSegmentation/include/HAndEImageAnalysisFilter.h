#ifndef HAndEImageAnalysisFilter_h_
#define HAndEImageAnalysisFilter_h_

#include <vector>

// itk
#include "itkImage.h"
#include "itkRGBPixel.h"

// openCV
//#include <opencv2/opencv.hpp>


namespace gth818n
{
  class HAndEImageAnalysisFilter
  {
  public:
    typedef HAndEImageAnalysisFilter Self;

    ////////////////////////////////////////////////////////////////////////////////
    /// ctor
    HAndEImageAnalysisFilter();
    ~HAndEImageAnalysisFilter() {}
    /// ctor, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// typedef
    static const unsigned int ImageDimension = 2;
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
    typedef itk::Image<unsigned int, ImageDimension> itkUIntImage2DType;
    typedef itk::Image<int, ImageDimension> itkIntImageType;

    typedef itk::Image<unsigned char, ImageDimension> itkUCharImage2DType;
    typedef itk::Image<char, ImageDimension> itkCharImageType;

    typedef itk::Image<unsigned short, ImageDimension> itkUShortImageType;
    typedef itk::Image<short, ImageDimension> itkShortImage2DType;

    typedef itk::RGBPixel<unsigned char> RGBPixelType;
    typedef itk::Image<RGBPixelType, ImageDimension> itkRGBImage2DType;

    typedef itkUCharImage2DType itkBinaryMaskImage2DType;
    typedef itkUIntImage2DType itkIntImage2DType;
    /// typedef, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// public fn
    //void setInputHAndEImage(cv::Mat inputHAndEImage) {m_HAndEImage = inputHAndEImage;}
    void setInputHAndEImage(itkRGBImage2DType::Pointer inputHAndEImage) {m_HAndEImage = inputHAndEImage;}

    void setMPP(float mpp);
    void setMagnification(int magnification) {m_magnification = magnification;}
    void setThreshold(float threshold) { m_threshold = threshold;}

    void setNucleusSizeThreshold(float sizeThld) {m_nucleusSizeThreshold = sizeThld;}

    void update();

    //cv::Mat getNormalizedImage();
    itkIntImage2DType::Pointer getNucleiLabelImage();
    itkBinaryMaskImage2DType::Pointer getNucleiBinaryMaskImage();

    int64_t getTotalNumberOfConnectedComponents();

    std::vector<long> getSizesOfAllNuclei(); ///< For this version I will just return this which causes memory copy. See if this is okay. Right way may use shared_ptr
    std::vector<float> getRoundnessOfAllNuclei(); ///< For this version I will just return this which causes memory copy. See if this is okay. Right way may use shared_ptr
    /// public fn, end
    ////////////////////////////////////////////////////////////////////////////////


  private:
    ////////////////////////////////////////////////////////////////////////////////
    /// typedefs
    typedef vnl_matrix< float > vnlFloatMatrixType;
    typedef vnl_matrix< double > vnlDoubleMatrixType;
    /// typedefs, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// private data
    int m_magnification;

    float m_mpp;

    // cv::Mat m_HAndEImage;
    // cv::Mat m_HAndEImageNormalized;
    itkRGBImage2DType::Pointer m_HAndEImage;
    itkUCharImage2DType::Pointer m_hematoxylinImage;

    float m_threshold;

    itkBinaryMaskImage2DType::Pointer m_nucleusBinaryMask;
    itkIntImage2DType::Pointer m_nucleiLabelImage;

    void _extractHematoxylinImage();


    int64_t m_totalNumberOfConnectedComponents;

    std::vector<long> m_sizesOfAllNuclei;
    std::vector<float> m_roundnessOfAllNuclei;

    float m_nucleusSizeThreshold; ///< region smaller than this number of pixels will be removed

    bool m_allDone;
    /// private data, end
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    /// private fn
    void _normalizeHAndEImage();
    void _segment();
    void _segmentHematoxylinImage();
    //void _segmentHematoxylinImage_otsu_watershed();

    void _computeNucleiLabelImage();

    void _computeAllNucleiSizes();
    void _computeAllNucleiRoundness();
    /// private fn
    ////////////////////////////////////////////////////////////////////////////////
  };

}// namespace gth818n


#endif
