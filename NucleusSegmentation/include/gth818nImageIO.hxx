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

#ifndef gth818nImageIO_hxx_
#define gth818nImageIO_hxx_

#include <csignal>
#include <string>

// itk
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"

#include "itkImageRegionIterator.h"

#include "itkVector.h"

#include "itkTransformFactoryBase.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"



// vnl
#include "vnl/vnl_matrix.h"


// local
#include "gth818nImageIO.h"

namespace gth818n
{
  //--------------------------------------------------------------------------------
  // IO, image
  template< typename TItkImage >
  typename TItkImage::Pointer
  readImage(const char *fileName)
  {
    typedef itk::ImageFileReader< TItkImage > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename TItkImage::Pointer image;

    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl;
        std::cerr<< err << std::endl;
        raise(SIGABRT);
      }

    return image;
  }


  template< typename TItkImage >
  void
  writeImage(typename TItkImage::Pointer img, const char *fileName, bool compress)
  {
    typedef itk::ImageFileWriter< TItkImage > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(img);
    if (compress)
      {
        writer->UseCompressionOn();
      }
    else
      {
        writer->UseCompressionOff();
      }

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        raise(SIGABRT);
      }
  }

  template< typename TItkImage >
  std::vector< typename TItkImage::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList )
  {
    typedef typename TItkImage::Pointer itkImagePointer_t;
    typedef std::vector< itkImagePointer_t > itkImageList_t;
    typedef itk::ImageFileReader< TItkImage > itkImageReader_t;


    itkImageList_t imageSeries;

    long n = imageNameList.size();
    for (long i = 0; i < n; ++i)
      {
        std::string thisName = imageNameList[i];

        typename itkImageReader_t::Pointer reader = itkImageReader_t::New();
        reader->SetFileName(thisName);

        itkImagePointer_t img;

        try
          {
            reader->Update();
            img = reader->GetOutput();
          }
        catch ( itk::ExceptionObject &err)
          {
            std::cerr<< "ExceptionObject caught !" << std::endl;
            std::cerr<< err << std::endl;
            raise(SIGABRT);
          }


        imageSeries.push_back(img);
      }

    return imageSeries;
  }


  template< typename TItkVectorImage >
  void
  writeVectorImage(typename TItkVectorImage::Pointer img, const char *fileName, int component)
  {
    typedef itk::Image<double, TItkVectorImage::ImageDimension> ItkImageType;
    typename ItkImageType::Pointer componentImg = ItkImageType::New();
    componentImg->SetRegions(img->GetLargestPossibleRegion() );
    componentImg->Allocate();


    typedef itk::ImageRegionIterator<TItkVectorImage> VectorIteratorType;
    VectorIteratorType vIter(img, img->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<ItkImageType> IteratorType;
    IteratorType iter(componentImg, componentImg->GetLargestPossibleRegion());

    for (vIter.GoToBegin(), iter.GoToBegin(); !vIter.IsAtEnd(); ++iter, ++vIter)
      {
        iter.Set(vIter.Get()[component]);
      }

    typedef itk::ImageFileWriter< ItkImageType > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(componentImg);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        abort();
      }
  }




  //-----------------------------------------------------------------------------
  /// Get the PixelType and ComponentType from fileName
  void GetImageType(std::string fileName, itk::ImageIOBase::IOPixelType &pixelType, itk::ImageIOBase::IOComponentType &componentType)
  {
    typedef itk::Image<unsigned char, 3> ImageType;
    itk::ImageFileReader<ImageType>::Pointer imageReader = itk::ImageFileReader<ImageType>::New();
    imageReader->SetFileName(fileName.c_str());
    imageReader->UpdateOutputInformation();

    pixelType = imageReader->GetImageIO()->GetPixelType();
    componentType = imageReader->GetImageIO()->GetComponentType();

    return;
  }





  //--------------------------------------------------------------------------------
  // For cython
  template<typename ScalarPixelType>
  void
  itkScalarImageToArray2D(typename itk::Image<ScalarPixelType, 2>::Pointer img, ScalarPixelType* array)
  {
    // Assuming the array space has already been allocated and owned outside here
    memcpy(array, img->GetBufferPointer(), img->GetLargestPossibleRegion().GetNumberOfPixels()*sizeof(ScalarPixelType));

    return;
  }

  template<typename RGBPixelType>
  typename itk::Image<itk::RGBPixel<RGBPixelType>, 2>::Pointer
  ucharArraryToItkRGBImage2D(unsigned char* rgbByteFlow, long nx, long ny, double spacing)
  {
    typedef itk::Image<itk::RGBPixel<RGBPixelType>, 2> itkRGBImage2DType;

    typename itkRGBImage2DType::SizeType size;
    size[0]  = nx;  // size along X
    size[1]  = ny;  // size along Y

    typename itkRGBImage2DType::IndexType start;
    start.Fill( 0 );

    typename itkRGBImage2DType::RegionType region;
    region.SetIndex( start );
    region.SetSize(  size  );

    typename itkRGBImage2DType::Pointer img = itkRGBImage2DType::New();
    img->SetRegions(region);
    img->Allocate();

    img->SetSpacing(spacing);

    typename itkRGBImage2DType::PixelType* imgBufferPointer = img->GetBufferPointer();
    long it = 0;
    for (long iy = 0; iy < ny; ++iy)
      {
        for (long ix = 0; ix < nx; ++ix)
          {
            imgBufferPointer[it].SetRed(static_cast<RGBPixelType>(rgbByteFlow[3*it]));
            imgBufferPointer[it].SetGreen(static_cast<RGBPixelType>(rgbByteFlow[3*it + 1]));
            imgBufferPointer[it].SetBlue(static_cast<RGBPixelType>(rgbByteFlow[3*it + 2]));

            ++it;
          }
      }

    return img;
  }

  template<typename RGBPixelType>
  void
  itkRGBImageToUcharArrary2D(typename itk::Image<itk::RGBPixel<RGBPixelType>, 2>::Pointer img, unsigned char* rgbByteFlow)
  {
    typedef itk::Image<itk::RGBPixel<RGBPixelType>, 2> itkRGBImage2DType;

    // we assume the space (nx*ny*3) for the image has been allocated in rgbByteFlow
    typename itkRGBImage2DType::PixelType* imgBufferPointer = img->GetBufferPointer();

    long n = static_cast<long>(img->GetLargestPossibleRegion().GetNumberOfPixels());
    for (long it = 0; it < n; ++it)
      {
        rgbByteFlow[3*it] = static_cast<unsigned char>(imgBufferPointer[it].GetRed());
        rgbByteFlow[3*it + 1] = static_cast<unsigned char>(imgBufferPointer[it].GetGreen());
        rgbByteFlow[3*it + 2] = static_cast<unsigned char>(imgBufferPointer[it].GetBlue());
      }

    return;
  }
  // For cython
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


}// gth818n

#endif
