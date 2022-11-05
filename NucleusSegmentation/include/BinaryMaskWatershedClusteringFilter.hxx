#ifndef BinaryMaskWatershedClusteringFilter_hxx_
#define BinaryMaskWatershedClusteringFilter_hxx_

#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkWatershedImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRelabelComponentImageFilter.h"

#include "BinaryMaskWatershedClusteringFilter.h"

//dbg
#include "Image/gth818nImage.h"
//dbg, end

namespace gth818n
{

  template< typename TInputImage, typename TOutputImage >
  BinaryMaskWatershedClusteringFilter< TInputImage, TOutputImage >::BinaryMaskWatershedClusteringFilter()
  {
    m_threshold = 0;
    m_level = 0.02;

    m_allDone = false;
  }


  template< typename TInputImage, typename TOutputImage >
  void BinaryMaskWatershedClusteringFilter< TInputImage, TOutputImage >::update()
  {
    _computeBinaryMask();
    _watershedOnSignedDistanceMap();

    m_allDone = true;

    return;
  }

  template< typename TInputImage, typename TOutputImage >
  void BinaryMaskWatershedClusteringFilter< TInputImage, TOutputImage >::_computeBinaryMask()
  {
    typedef itk::BinaryThresholdImageFilter<InputImageType, UCharImageType> BinaryThresholdImageFilterType;
    typename BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
    thresholdFilter->SetInput(m_inputImage);
    thresholdFilter->SetLowerThreshold(1);
    thresholdFilter->SetUpperThreshold( itk::NumericTraits< InputPixelType >::max() );
    thresholdFilter->SetInsideValue(1);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();

    m_binaryMaskOfInputImage = thresholdFilter->GetOutput();


    //dbg
    IO::writeImage<UCharImageType>(m_binaryMaskOfInputImage, "bin.nrrd" );
    //dbg, end

    return;
  }

  template< typename TInputImage, typename TOutputImage >
  typename BinaryMaskWatershedClusteringFilter< TInputImage, TOutputImage >::OutputImageType::Pointer
  BinaryMaskWatershedClusteringFilter< TInputImage, TOutputImage >::getOutputImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_outputImage;
  }


  template< typename TInputImage, typename TOutputImage >
  void BinaryMaskWatershedClusteringFilter< TInputImage, TOutputImage >::_watershedOnSignedDistanceMap()
  {
    typedef itk::SignedDanielssonDistanceMapImageFilter< UCharImageType, FloatImageType >  SignedDanielssonDistanceMapImageFilter;
    typename SignedDanielssonDistanceMapImageFilter::Pointer sdfFilter = SignedDanielssonDistanceMapImageFilter::New();
    sdfFilter->SetInput( m_binaryMaskOfInputImage );
    sdfFilter->SetInsideIsPositive(false);
    sdfFilter->UseImageSpacingOn();
    sdfFilter->Update();

    //dbg
    IO::writeImage<FloatImageType>(sdfFilter->GetOutput(), "sdf.nrrd" );
    //dbg, end

    typedef itk::WatershedImageFilter<FloatImageType> WatershedFilterType;
    typename WatershedFilterType::Pointer watershed = WatershedFilterType::New();
    watershed->SetThreshold(m_threshold);
    watershed->SetLevel(m_level);
    watershed->SetInput( sdfFilter->GetOutput() );
    watershed->Update();

    typename WSOutputImageType::Pointer wsOutput = watershed->GetOutput();

    typedef itk::RelabelComponentImageFilter<WSOutputImageType, WSOutputImageType> FilterType;
    typename FilterType::Pointer relabelFilter = FilterType::New();
    relabelFilter->SetInput(wsOutput);
    //relabelFilter->SetMinimumObjectSize(static_cast<FilterType::ObjectSizeType>(m_objectSizeThreshold/m_mpp/m_mpp)); // This command takes number of pixels as input
    relabelFilter->Update();
    wsOutput = relabelFilter->GetOutput();

    typedef itk::MinimumMaximumImageCalculator<WSOutputImageType> ImageCalculatorFilterType;
    typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
    imageCalculatorFilter->SetImage(wsOutput);
    imageCalculatorFilter->Compute();

    if (static_cast<double>(imageCalculatorFilter->GetMaximum()) > static_cast<double>(itk::NumericTraits< OutputPixelType >::max()) )
      {
        std::cout<<"Warning: water shed output exceed output pixel holder max.\n";
        std::cout<<imageCalculatorFilter->GetMaximum()<<"\t"<<itk::NumericTraits< OutputPixelType >::max()<<std::endl;
      }

    if (static_cast<double>(imageCalculatorFilter->GetMinimum()) < static_cast<double>(itk::NumericTraits< OutputPixelType >::min()) )
      {
        /// if output pixel type is float or double, min() is 0. This may cause problem
        std::cout<<"Warning: water shed output exceed output pixel holder min.\n";
        std::cout<<imageCalculatorFilter->GetMinimum()<<"\t"<<itk::NumericTraits< OutputPixelType >::min()<<std::endl;
      }

    typedef itk::CastImageFilter< WSOutputImageType, OutputImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(wsOutput);
    castFilter->Update();


    m_outputImage = castFilter->GetOutput();

    return;
  }



}// namespace gth818n


#endif
