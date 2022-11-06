#include <vector>

// itk
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkThresholdSegmentationLevelSetImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"


//#include "Image/gth818nImage.h"
#include "utilitiesImage.h"

#include "itkCastImageFilter.h"

#include "itkWatershedImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

//#include "itkGradientMagnitudeImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

#include "ColorDecomposition.h"

// itk opencv bridge
//#include "itkOpenCVImageBridge.h"

// openCV
//#include <opencv2/opencv.hpp>

// local
#include "HAndEImageAnalysisFilter.h"
#include "../include/SFLSLocalChanVeseSegmentor2D.h"

#include "FastGrowCutSegmenter2D.h"
#include "HematoxylinImageSegmentationFilter.h"

// for debug
//#include "Image/gth818nImage.h"

namespace gth818n
{
  void HAndEImageAnalysisFilter::update()
  {
    // _normalizeHAndEImage();
    // _segment();

    _extractHematoxylinImage();

    _segmentHematoxylinImage();
    //_segmentHematoxylinImage_otsu_watershed();

    _computeNucleiLabelImage();
    // _computeAllNucleiSizes(); 
    // _computeAllNucleiRoundness();

    m_allDone = true;

    return;
  }

  // void HAndEImageAnalysisFilter::_normalizeHAndEImage()
  // {
  //   float meanT[3] = {-0.176, -0.0033, 0.0143};
  //   float stdT[3] = {0.2317, 0.0491, 0.0156};

  //   m_HAndEImageNormalized =  nscale::Normalization::normalization(m_HAndEImage, meanT, stdT);

  //   return;
  // }

  void HAndEImageAnalysisFilter::_segmentHematoxylinImage()
  {
    HematoxylinImageSegmentationFilter imageAnalyzer;
    imageAnalyzer.setHematoxylinImage( m_hematoxylinImage );
    imageAnalyzer.setNucleusSizeThreshold(m_nucleusSizeThreshold);
    imageAnalyzer.update();

    m_nucleusBinaryMask = imageAnalyzer.getNucleiBinaryMaskImage();

    return;
  }



  // void HAndEImageAnalysisFilter::_segmentHematoxylinImage()
  // {
  //   typedef itk::CastImageFilter< itkUCharImage2DType, itkFloatImage2DType> CastFilterType_UC_to_F;
  //   CastFilterType_UC_to_F::Pointer castFilter1 = CastFilterType_UC_to_F::New();
  //   castFilter1->SetInput(m_hematoxylinImage);
  //   castFilter1->Update();

  //   itkFloatImage2DType::Pointer originalImage = castFilter1->GetOutput();;

  //   itkShortImage2DType::Pointer mask = thresholdHematoxylinImage(originalImage);
  //   itkShortImage2DType::PixelType* maskBufferPointer = mask->GetBufferPointer();

  //   itkFloatImage2DType::Pointer dist = gth818n::bw2dist(mask);
  //   dist = GaussianSmoothing(dist, 60);

  //   //gth818n::writeImage<itkFloatImage2DType>(dist, outputImageName, true);

  //   typedef itk::GradientMagnitudeImageFilter<itkFloatImage2DType, itkFloatImage2DType >  GradientMagnitudeImageFilterType;
  //   GradientMagnitudeImageFilterType::Pointer gradientMagnitudeImageFilter = GradientMagnitudeImageFilterType::New();
  //   gradientMagnitudeImageFilter->SetInput(originalImage);
  //   gradientMagnitudeImageFilter->Update();
  //   itkFloatImage2DType::Pointer gradImage = gradientMagnitudeImageFilter->GetOutput();

  //   itkUCharImage2DType::Pointer localMinMask = gth818n::localMinimaMaskFilter(dist);
  //   typedef itk::ConnectedComponentImageFilter <itkUCharImage2DType, itkUIntImage2DType > ConnectedComponentImageFilterType;
  //   ConnectedComponentImageFilterType::Pointer connected1 = ConnectedComponentImageFilterType::New();
  //   connected1->SetInput(localMinMask);
  //   connected1->Update();
  //   itkUIntImage2DType::Pointer seedImage = connected1->GetOutput();

  //   long numPixels = gradImage->GetLargestPossibleRegion().GetNumberOfPixels();

  //   itkUIntImage2DType::Pointer scResults = itkUIntImage2DType::New();

  //   std::vector<long> imSize(2);
  //   imSize[0] = gradImage->GetLargestPossibleRegion().GetSize()[0];
  //   imSize[1] = gradImage->GetLargestPossibleRegion().GetSize()[1];

  //   std::vector<itkFloatImage2DType::PixelType> m_imSrcVec(numPixels);
  //   std::vector<itkUIntImage2DType::PixelType> m_imSeedVec(numPixels);

  //   itkFloatImage2DType::PixelType* gradImageBufferPointer = gradImage->GetBufferPointer();
  //   itkUIntImage2DType::PixelType* seedImageBufferPointer = seedImage->GetBufferPointer();

  //   for (long it = 0; it < numPixels; ++it)
  //     {
  //       m_imSrcVec[it] = gradImageBufferPointer[it];
  //       m_imSeedVec[it] = seedImageBufferPointer[it];
  //     }

  //   FGC::FastGrowCut<itkFloatImage2DType::PixelType, itkUIntImage2DType::PixelType> *m_fastGC = new FGC::FastGrowCut<itkFloatImage2DType::PixelType, itkUIntImage2DType::PixelType>();
  //   m_fastGC->SetSourceImage(m_imSrcVec);
  //   m_fastGC->SetSeedlImage(m_imSeedVec);
  //   m_fastGC->SetImageSize(imSize);
  //   m_fastGC->SetWorkMode(false);

  //   // Do Segmentation
  //   m_fastGC->DoSegmentation();

  //   std::vector<itkUIntImage2DType::PixelType> m_imLabVec(numPixels); ///< to store the output
  //   m_fastGC->GetLabeImage(m_imLabVec);


  //   //itkUIntImage2DType::Pointer labelImage = itkUIntImage2DType::New();
  //   scResults->SetRegions(gradImage->GetLargestPossibleRegion() );
  //   scResults->Allocate();
  //   scResults->FillBuffer(0);

  //   itkUIntImage2DType::PixelType* scResultsBufferPointer = scResults->GetBufferPointer();

  //   for (long it = 0; it < numPixels; ++it)
  //     {
  //       scResultsBufferPointer[it] = m_imLabVec[it];
  //     }

  //   gth818n::writeImage<itkUIntImage2DType>(scResults, "scResults.nrrd", true);


  //   delete m_fastGC;

  //   std::cout<<"before edgesOfDifferentLabelRegion\n"<<std::flush;
  //   itkUCharImage2DType::Pointer edgeOfScResults = gth818n::edgesOfDifferentLabelRegion(scResults);
  //   itkUCharImage2DType::PixelType* edgeOfScResultsBufferPointer = edgeOfScResults->GetBufferPointer();
  //   std::cout<<"after edgesOfDifferentLabelRegion\n"<<std::flush;

  //   gth818n::writeImage<itkUCharImage2DType>(edgeOfScResults, "edgeOfSC.nrrd", true);

  //   bool maxOfBufferIs0 = true;
  //   for (long it = 0; it < numPixels; ++it)
  //     {
  //       maskBufferPointer[it] *= (1 - edgeOfScResultsBufferPointer[it]);
  //       if (maskBufferPointer[it] > 0)
  //         {
  //           maxOfBufferIs0 = false;
  //         }
  //     }

  //   std::cout<<"before casting\n"<<std::flush;

  //   if (!maxOfBufferIs0)
  //     {
  //       typedef itk::CastImageFilter< itkShortImage2DType, itkUCharImage2DType> CastFilterType;
  //       CastFilterType::Pointer castFilter = CastFilterType::New();
  //       castFilter->SetInput(mask);
  //       castFilter->Update();

  //       std::cout<<"after casting, before Chan Vese\n"<<std::flush;

  //       //  itkUIntImage2DType::Pointer sizeLabel = gth818n::labelImageToSizeLabeledImage(castFilter->GetOutput(), 5);

  //       int numiter = 10;
  //       CSFLSLocalChanVeseSegmentor2D< itkFloatImage2DType::PixelType > cv;
  //       cv.setImage(originalImage);
  //       cv.setMask( castFilter->GetOutput() );
  //       cv.setNumIter(numiter);
  //       float lambda = 0.3;
  //       cv.setCurvatureWeight(lambda);
  //       cv.doSegmenation();

  //       std::cout<<"after Chan Vese\n"<<std::flush;


  //       CSFLSLocalChanVeseSegmentor2D< itkFloatImage2DType::PixelType >::LSImageType::Pointer phi = cv.mp_phi;

  //       itkUCharImage2DType::Pointer lsMask = itkUCharImage2DType::New();
  //       lsMask->SetRegions(phi->GetLargestPossibleRegion() );
  //       lsMask->Allocate();
  //       lsMask->FillBuffer(0);

  //       itkUCharImage2DType::PixelType* lsMaskBufferPtr = lsMask->GetBufferPointer();
  //       CSFLSLocalChanVeseSegmentor2D< itkFloatImage2DType::PixelType >::LSImageType::PixelType* phiBufferPointer = phi->GetBufferPointer();

  //       for (long it = 0; it < numPixels; ++it)
  //         {
  //           lsMaskBufferPtr[it] = phiBufferPointer[it]<=0?1:0;
  //         }

  //       m_nucleusBinaryMask = lsMask;
  //     }
  //   else
  //     {
  //       m_nucleusBinaryMask = itkUCharImage2DType::New();
  //       m_nucleusBinaryMask->SetRegions(mask->GetLargestPossibleRegion() );
  //       m_nucleusBinaryMask->Allocate();
  //       m_nucleusBinaryMask->FillBuffer(0);
  //     }

  //   std::cout<<"after Chan Vese post labeling\n"<<std::flush;

  //   return;
  //   // typedef itk::ConnectedComponentImageFilter <itkUCharImage2DType, itkUIntImage2DType > ConnectedComponentImageFilterType;
  //   // ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();
  //   // //connected->SetInput(thresholder->GetOutput());
  //   // connected2->SetInput(lsMask);
  //   // connected2->Update();

  //   // itkUIntImage2DType::Pointer sizeLabelImg = gth818n::labelImageToSizeLabeledImage(connected2->GetOutput(), 10);

  //   //gth818n::writeImage<itkUIntImage2DType>(sizeLabelImg, outputImageName, true);
  // }

  // void HAndEImageAnalysisFilter::_segment()
  // {
  //   cv::Mat imageGray(m_HAndEImage.size(), CV_8UC1, cv::Scalar(0));
  //   cv::Mat labelImage(m_HAndEImage.size(), CV_8UC1, cv::Scalar(0));

  //   int64_t numOfPixelPerTile = m_HAndEImage.total();

  //   //#pragma omp parallel for
  //   for (int64_t it = 0; it < numOfPixelPerTile; ++it)
  //     {
  //       // Threshold
  //       float b = static_cast<float>(m_HAndEImage.at<cv::Vec3b>(it)[0]);
  //       float g = static_cast<float>(m_HAndEImage.at<cv::Vec3b>(it)[1]);
  //       float r = static_cast<float>(m_HAndEImage.at<cv::Vec3b>(it)[2]);

  //       float v = (b + g + r)/3.0;

  //       labelImage.at<unsigned char>(it) = v<m_threshold?1:0;

  //       imageGray.at<unsigned char>(it) = v;
  //     }

  //   ////////////////////////////////////////////////////////////////////////////////
  //   /// Convert cv::Mat image to itk image
  //   itkUCharImage2DType::Pointer itkImageTmp = itk::OpenCVImageBridge::CVMatToITKImage< itkUCharImage2DType >( imageGray );

  //   typedef itk::CastImageFilter< itkUCharImage2DType, itkFloatImage2DType > CastFilterType;
  //   CastFilterType::Pointer castFilter = CastFilterType::New();
  //   castFilter->SetInput(itkImageTmp);
  //   castFilter->Update();
  //   itkFloatImage2DType::Pointer itkImageGrayCurrentTile = castFilter->GetOutput();

  //   m_nucleusBinaryMask = itk::OpenCVImageBridge::CVMatToITKImage< itkUCharImage2DType >( labelImage );
  //   /// Convert cv::Mat image to itk image
  //   ////////////////////////////////////////////////////////////////////////////////

  //   typedef itk::DanielssonDistanceMapImageFilter< itkUCharImage2DType, itkFloatImage2DType >  DanielssonDistanceMapImageFilterType;
  //   DanielssonDistanceMapImageFilterType::Pointer sdfFilter = DanielssonDistanceMapImageFilterType::New();
  //   sdfFilter->SetInput( m_nucleusBinaryMask );
  //   sdfFilter->InputIsBinaryOn();
  //   sdfFilter->UseImageSpacingOn();
  //   sdfFilter->Update();

  //   typedef itk::ThresholdSegmentationLevelSetImageFilter< itkFloatImage2DType, itkFloatImage2DType > ThresholdSegmentationLevelSetImageFilterType;
  //   ThresholdSegmentationLevelSetImageFilterType::Pointer thresholdSegmentation = ThresholdSegmentationLevelSetImageFilterType::New();
  //   thresholdSegmentation->SetCurvatureScaling( 5.0 );
  //   thresholdSegmentation->SetPropagationScaling( 0.1 ); ///< this param and the SetCurvatureScaling one, as long as having the same ratio, (5, 1) and (50, 10) seem to be the same.
  //   thresholdSegmentation->SetMaximumRMSError( 1e-6 );
  //   thresholdSegmentation->SetNumberOfIterations( 100 ); ///< 1000 may be good
  //   thresholdSegmentation->SetUpperThreshold( m_threshold ); ///< change to 75 will significantly shrink to darker area
  //   thresholdSegmentation->SetLowerThreshold( 0 ); ///< change this to -10 does not affect result
  //   thresholdSegmentation->SetIsoSurfaceValue(1.0); ///< the output of sdfFilter has 1.0-isocontour the same as sdf's input binary
  //   thresholdSegmentation->SetInput( sdfFilter->GetOutput() );
  //   thresholdSegmentation->SetFeatureImage( itkImageGrayCurrentTile );
  //   thresholdSegmentation->Update();

  //   typedef itk::BinaryThresholdImageFilter<itkFloatImage2DType, itkUCharImage2DType> ThresholdingFilterType;
  //   ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  //   thresholder->SetLowerThreshold( -1000.0 );
  //   thresholder->SetUpperThreshold(     0.0 );

  //   thresholder->SetOutsideValue(  0  );
  //   thresholder->SetInsideValue(  25 );
  //   thresholder->SetInput( thresholdSegmentation->GetOutput() );
  //   thresholder->Update();

  //   m_nucleusBinaryMask = thresholder->GetOutput();

  //   return;
  // }


  void HAndEImageAnalysisFilter::_computeNucleiLabelImage()
  {
    ////////////////////////////////////////////////////////////////////////////////
    /// Connected component to separate each nucleus
    ///
    /// The max of <unsigned int>, as tested in
    /// ../test/mainNumericTraits.cxx, is 4G. Should be large enough for
    /// a tile. Enough for a WSI?
    /// typedef unsigned int LabelType;
    /// typedef itk::Image< LabelType, Dimension > LabelImageType;
    typedef itk::ConnectedComponentImageFilter <itkUCharImage2DType, itkUIntImage2DType > ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
    connected->SetInput(m_nucleusBinaryMask);
    connected->Update();
    /// Connected component to separate each nucleus
    ////////////////////////////////////////////////////////////////////////////////

    m_totalNumberOfConnectedComponents = connected->GetObjectCount();

    m_nucleiLabelImage = connected->GetOutput();

    return;
  }

  void HAndEImageAnalysisFilter::_computeAllNucleiSizes()
  {
    //std::cout<<"m_totalNumberOfConnectedComponents = "<<m_totalNumberOfConnectedComponents<<std::endl<<std::flush;

    m_sizesOfAllNuclei.resize(m_totalNumberOfConnectedComponents + 1);

    //std::cout<<"----------- 1 ----------------\n"<<std::endl;

    const itkUIntImage2DType::PixelType* labelImageOfCurrentTilePtr = static_cast<itkUIntImage2DType::PixelType*>(m_nucleiLabelImage->GetBufferPointer());

    //std::cout<<"----------- 2 ----------------\n"<<std::endl;

    unsigned long numberOfPixels = m_nucleiLabelImage->GetLargestPossibleRegion().GetNumberOfPixels();

    //std::cout<<"numberOfPixels = "<<numberOfPixels<<std::endl<<std::flush;

    float magnificationAdjustRatio = static_cast<float>(m_magnification)/20.0; ///< So that the size will the equivalent to under the magnification of 20x

    //std::cout<<"magnificationAdjustRatio = "<<magnificationAdjustRatio<<std::endl<<std::flush;

    for (unsigned long it = 0; it < numberOfPixels; ++it)
      {
        ++m_sizesOfAllNuclei[labelImageOfCurrentTilePtr[it]];
      }

    for (unsigned long it = 0; it < m_sizesOfAllNuclei.size(); ++it)
      {
        m_sizesOfAllNuclei[it] /= magnificationAdjustRatio;
      }

    return;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// The roundness of an object is computed as 1 - FA where FA is the
  /// fractional anisotropy computed from the two eigen values of the
  /// covariance matrix of all the points in this object.
  ///
  /// Should be from 0 (non-round) to 1 (pure round)
  void HAndEImageAnalysisFilter::_computeAllNucleiRoundness()
  {
    /// go thru each pixel, put its coordinate to the corresponding list
    std::vector< std::vector<itk2DIndexType> > listOfCoordinatesOfEachObject(m_totalNumberOfConnectedComponents + 1);

    const itkUIntImage2DType::PixelType* labelImageOfCurrentTilePtr = static_cast<itkUIntImage2DType::PixelType*>(m_nucleiLabelImage->GetBufferPointer());

    unsigned long nx = m_nucleiLabelImage->GetLargestPossibleRegion().GetSize()[0];
    unsigned long ny = m_nucleiLabelImage->GetLargestPossibleRegion().GetSize()[1];

    {
      unsigned long it = 0;
      itk2DIndexType idx;

      for (unsigned long iy = 0; iy < ny; ++iy)
        {
          idx[1] = iy;
          for (unsigned long ix = 0; ix < nx; ++ix)
            {
              idx[0] = ix;
              itkUIntImage2DType::PixelType objectId = labelImageOfCurrentTilePtr[it++];
              listOfCoordinatesOfEachObject[objectId].push_back(idx);
            }
        }
    }

    /// For each object, if it contains more than pixels in [10, 200],
    /// compute the fractional anisotrpoy. If we have 3 pixels, we can
    /// compute the co-variance matrix and then the eigen
    /// decomposition. But only 10 to 200 pixel area may be
    /// meaningful.
    m_roundnessOfAllNuclei.resize(m_totalNumberOfConnectedComponents + 1);

    for (itkUIntImage2DType::PixelType itObject = 0; itObject < m_totalNumberOfConnectedComponents; ++itObject)
      {
        if (10 > m_sizesOfAllNuclei[itObject] || 200 < m_sizesOfAllNuclei[itObject])
          {
            m_roundnessOfAllNuclei[itObject] = -1;
            continue;
          }

        const std::vector<itk2DIndexType>& thePointsInThisObject = listOfCoordinatesOfEachObject[itObject];
        int numOfPointsInThisObject = thePointsInThisObject.size();

        float meanX = 0;
        float meanY = 0;

        for (int itPoint = 0; itPoint < numOfPointsInThisObject; ++itPoint)
          {
            meanX += thePointsInThisObject[itPoint][0];
            meanY += thePointsInThisObject[itPoint][1];
          }

        meanX /= static_cast<float>(numOfPointsInThisObject);
        meanY /= static_cast<float>(numOfPointsInThisObject);

        vnlFloatMatrixType M(2, numOfPointsInThisObject);
        for (int itPoint = 0; itPoint < numOfPointsInThisObject; ++itPoint)
          {
            M(0, itPoint) = thePointsInThisObject[itPoint][0] - meanX;
            M(1, itPoint) = thePointsInThisObject[itPoint][1] - meanY;
          }

        vnlFloatMatrixType covMatrix = M*( M.transpose() );

        vnl_symmetric_eigensystem<float> eig(covMatrix);
        float ev0 = eig.get_eigenvalue(0);
        float ev1 = eig.get_eigenvalue(1);

        float fa = fabs(ev0 - ev1)/sqrt(ev0*ev0 + ev1*ev1);

        m_roundnessOfAllNuclei[itObject] = 1.0 - fa;
      }

    return;
  }




  HAndEImageAnalysisFilter::HAndEImageAnalysisFilter()
  {
    m_nucleusSizeThreshold = 0;

    m_allDone = false;

    m_threshold = 50; ///< 75 is small, miss a lot of nuclei. But of course the false positive is low. To make the false positive even lower, let me use 50
    //m_threshold = 150; ///< 150 till now over seg all luad tiles. But i think it safer to put to 175 if we want small false negative

    return;
  }



  // cv::Mat HAndEImageAnalysisFilter::getNormalizedImage()
  // {
  //   if (!m_allDone)
  //     {
  //       std::cerr<<"Error: computation not performed.\n";
  //       abort();
  //     }

  //   return m_HAndEImageNormalized;
  // }


  HAndEImageAnalysisFilter::itkIntImage2DType::Pointer
  HAndEImageAnalysisFilter::getNucleiLabelImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_nucleiLabelImage;
  }

  HAndEImageAnalysisFilter::itkBinaryMaskImage2DType::Pointer
  HAndEImageAnalysisFilter::getNucleiBinaryMaskImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_nucleusBinaryMask;
  }

  std::vector<long>
  HAndEImageAnalysisFilter::getSizesOfAllNuclei()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_sizesOfAllNuclei;
  }


  std::vector<float>
  HAndEImageAnalysisFilter::getRoundnessOfAllNuclei()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_roundnessOfAllNuclei;
  }

  int64_t
  HAndEImageAnalysisFilter::getTotalNumberOfConnectedComponents()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not performed.\n";
        abort();
      }

    return m_totalNumberOfConnectedComponents;
  }

  void
  HAndEImageAnalysisFilter::_extractHematoxylinImage()
  {
    m_hematoxylinImage = gth818n::ExtractHematoxylinChannel(m_HAndEImage);

    return;
  }




}// namespace gth818n
