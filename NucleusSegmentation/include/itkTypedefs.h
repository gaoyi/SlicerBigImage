#ifndef itkTypedefs_h_
#define itkTypedefs_h_


// itk
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "vnl/vnl_matrix.h"


namespace gth818n
{
  const unsigned int ImageDimension = 2;
  typedef itk::Image<float, ImageDimension> itkFloatImage2DType;

  typedef itk::RGBPixel<unsigned char> RGBPixelType;
  typedef itk::Image<RGBPixelType, ImageDimension> itkRGBImage2DType;

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

  typedef itk::Vector< float, 2 > itkVectorType;


  ////////////////////////////////////////////////////////////////////////////////
  /// VNL types
  typedef vnl_matrix< float > vnlFloatMatrixType;
  typedef vnl_matrix< double > vnlDoubleMatrixType;


}// namespace gth818n



#endif
