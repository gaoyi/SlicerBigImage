/*******************************************************************************/
/*                                                                             */
/*  This program is distributed in the hope that it will be useful, but        */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY */
/*  or FITNESS FOR A PARTICULAR PURPOSE.                                       */
/*                                                                             */
/*  Please contact the author for reusing or redistributing this program.      */
/*                                                                             */
/*                                                  Copyright (c) 2010, Yi Gao */
/*                                                            gaoyi@gatech.edu */
/*                                                                             */
/*******************************************************************************/

#ifndef gth818nImageIO_h_
#define gth818nImageIO_h_


#include <string>
#include <vector>

// itk
#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkImageIOBase.h"



namespace gth818n
{
  //--------------------------------------------------------------------------------
  // IO, image
  template< typename TItkImage >
  typename TItkImage::Pointer
  readImage(const char *fileName);

  template< typename TItkImage >
  void
  writeImage(typename TItkImage::Pointer img, const char *fileName, bool compress = true);

  template< typename TItkImage >
  std::vector< typename TItkImage::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList );

  template< typename TItkVectorImage >
  void
  writeVectorImage(typename TItkVectorImage::Pointer img, const char *fileName, int component);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  //-----------------------------------------------------------------------------
  /// Get the PixelType and ComponentType from fileName
  void GetImageType (std::string fileName, itk::ImageIOBase::IOPixelType &pixelType, itk::ImageIOBase::IOComponentType &componentType);


  //--------------------------------------------------------------------------------
  // convert itk image with c-array. This is for cython
  template<typename RGBPixelType = unsigned char>
  typename itk::Image<itk::RGBPixel<RGBPixelType>, 2>::Pointer
  ucharArraryToItkRGBImage2D(unsigned char* rgbByteFlow, long nx, long ny, double spacing = 1.0); // assume isotropic spacing unit: um

  template<typename RGBPixelType = unsigned char>
  void
  itkRGBImageToUcharArrary2D(typename itk::Image<itk::RGBPixel<RGBPixelType>, 2>::Pointer img, unsigned char* rgbByteFlow);

  template<typename ScalarPixelType>
  void
  itkScalarImageToArray2D(typename itk::Image<ScalarPixelType, 2>::Pointer img, ScalarPixelType* array);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


}// gth818n


#include "gth818nImageIO.hxx"

#endif
