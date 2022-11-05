#ifndef GrayImageWatershedSegmentationFilter_hxx_
#define GrayImageWatershedSegmentationFilter_hxx_

#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkWatershedImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"

#include "GrayImageWatershedSegmentationFilter.h"

//dbg
#include "Image/gth818nImage.h"
//dbg, end

namespace gth818n
{

  template< typename TInputImage, typename TOutputImage >
  GrayImageWatershedSegmentationFilter< TInputImage, TOutputImage >::GrayImageWatershedSegmentationFilter()
  {
    m_threshold = 0.1;
    m_level = 0.33;

    m_allDone = false;
  }


  template< typename TInputImage, typename TOutputImage >
  void GrayImageWatershedSegmentationFilter< TInputImage, TOutputImage >::update()
  {
    _watershed();

    m_allDone = true;

    return;
  }


  template< typename TInputImage, typename TOutputImage >
  typename GrayImageWatershedSegmentationFilter< TInputImage, TOutputImage >::OutputImageType::Pointer
  GrayImageWatershedSegmentationFilter< TInputImage, TOutputImage >::getOutputImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_outputImage;
  }


  template< typename TInputImage, typename TOutputImage >
  void GrayImageWatershedSegmentationFilter< TInputImage, TOutputImage >::_watershed()
  {
    typedef itk::Image<float, TInputImage::ImageDimension> FloatImageType;
    typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType, FloatImageType > PMFilterType;
    typename PMFilterType::Pointer pmfilter = PMFilterType::New();
    pmfilter->SetInput( m_inputImage );
    pmfilter->SetNumberOfIterations( 50 );
    pmfilter->SetTimeStep( 0.125 );
    double conductanceParameter = 0.5;
    pmfilter->SetConductanceParameter( conductanceParameter );
    pmfilter->Update();

    typedef itk::GradientMagnitudeImageFilter<FloatImageType, FloatImageType >  GradientMagnitudeImageFilterType;
    typename GradientMagnitudeImageFilterType::Pointer gradientMagnitudeImageFilter = GradientMagnitudeImageFilterType::New();
    gradientMagnitudeImageFilter->SetInput(pmfilter->GetOutput());
    gradientMagnitudeImageFilter->Update();

    typedef itk::WatershedImageFilter<FloatImageType> WatershedFilterType;
    typename WatershedFilterType::Pointer watershed = WatershedFilterType::New();
    watershed->SetThreshold(m_threshold);
    watershed->SetLevel(m_level);
    watershed->SetInput( gradientMagnitudeImageFilter->GetOutput() );
    watershed->Update();

    // typedef itk::WatershedImageFilter<InputImageType> WatershedFilterType;
    // typename WatershedFilterType::Pointer watershed = WatershedFilterType::New();
    // watershed->SetThreshold(m_threshold);
    // watershed->SetLevel(m_level);
    // watershed->SetInput( m_inputImage );
    // watershed->Update();

    typedef typename WatershedFilterType::OutputImageType WSOutputImageType;
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
