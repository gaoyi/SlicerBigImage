#include <map>

// itk
#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"


// local
#include "utilitiesImage.h"
#include "BinaryMaskAnalysisFilter.h"
//#include "MeanshiftClusteringFilter.h"
#include "HierarchicalMeanshiftClusteringFilter.h"


namespace gth818n
{
  BinaryMaskAnalysisFilter::BinaryMaskAnalysisFilter()
  {
    m_featureColoredImage = 0;

    m_connectedComponentLabelImage = 0;

    m_objectAreas.clear();
    m_objectPerimeters.clear();
    m_objectEquivalentSphericalRadius.clear();

    m_micronPerPixel = -1.0;

    m_binaryMask = 0;
    m_inputImage = 0;

    m_numberOfObjects = 0;

    m_objectSizeThreshold = 8; ///< smallest cell, sperm, nutraphils, platlets (may be considered as cell fragments) are around 3 micron in one dim. So set 8 um^2 as lower limit
    m_objectSizeUpperThreshold = 200;

    m_allDone = false;

    return;
  }


  void BinaryMaskAnalysisFilter::update()
  {
    if (m_micronPerPixel < 0)
      {
        std::cerr<<"ERROR: mpp not set yet.\n";
        abort();
      }

    ///< spacing AFFECT the ShapeLabelObject filter in things like
    ///EquivalentSphericalRadius, but not those features that are
    ///pixel based.
    m_binaryMask->SetSpacing(m_micronPerPixel/1000.0);  //SetSpacing is mm based, so /1000.

    _computeConnectedComponentsLabelImage();

    _computeLabelMap();

    _computeObjectFeatures();

    _findObjectsToBreak();

    _breakRegion();

    m_allDone = true;

    return;
  }

  void BinaryMaskAnalysisFilter::_findObjectsToBreak()
  {
    m_objectToBreak.assign( m_numberOfObjects, 0 ); ///< same as m_objectAreas.size()

    /*--------------------------------------------------------------------------------*/
    // /// Criteria 1. If one area is too large, break it
    // for (std::size_t objectId = 0; objectId < m_objectAreas.size(); ++objectId)
    //   {
    //     if (m_objectAreas[objectId] > m_objectSizeUpperThreshold)
    //       {
    //         m_objectToBreak[objectId] = 1;
    //       }
    //     else
    //       {
    //         m_objectToBreak[objectId] = 0;
    //       }
    //   }


    // /*--------------------------------------------------------------------------------*/
    // /// Criteria 2. If one perimeter/ratio is too high, break
    // for (std::size_t objectId = 0; objectId < m_objectAreas.size(); ++objectId)
    //   {
    //     double measure = m_objectPerimeters[objectId]*m_objectPerimeters[objectId]/m_objectAreas[objectId];

    //     if (measure > 50) ///< For COAD A6-2671, 10 is too low and break some single nucleus, 100 is too large and can't break some that should break
    //       {
    //         m_objectToBreak[objectId] = 1;
    //       }
    //     else
    //       {
    //         m_objectToBreak[objectId] = 0;
    //       }
    //   }


    /*--------------------------------------------------------------------------------*/
    /// Criteria 3. If one perimeter/ratio is too high, break
    for (std::size_t objectId = 0; objectId < m_objectAreas.size(); ++objectId)
      {
        double measure = m_objectPerimeters[objectId]*m_objectPerimeters[objectId]/m_objectAreas[objectId];

        if (m_objectAreas[objectId] > m_objectSizeUpperThreshold)
          {
            m_objectToBreak[objectId] = 1;
            continue;
          }

        if (measure > 30) ///< For COAD A6-2671, 10 is too low and break some single nucleus, 100 is too large and can't break some that should break
          {
            m_objectToBreak[objectId] = 1;
            continue;
          }
      }

    return;
  }

  void BinaryMaskAnalysisFilter::_breakRegion()
  {
    //typedef gth818n::MeanshiftClusteringFilter<float, 2> MeanshiftClusteringFilterType;
    typedef gth818n::HierarchicalMeanshiftClusteringFilter<float, 2> MeanshiftClusteringFilterType;

    typedef MeanshiftClusteringFilterType::VectorType VectorType;
    typedef itk::Statistics::ListSample< VectorType > SampleType;

    /// Step 10. Go through all labels, find ones that are larger than 400um^2
    /// Step 20. Get all index of the region with this label, save as listSamples
    long nx = m_connectedComponentLabelImage->GetLargestPossibleRegion().GetSize()[0];
    long ny = m_connectedComponentLabelImage->GetLargestPossibleRegion().GetSize()[1];

    std::vector<itkIntImage2DType::PixelType> labelsWithLargeAreas;
    std::vector<SampleType::Pointer> indexOfThatLabel;

    itk2DIndexType idx;
    VectorType idxVec;

    itkIntImage2DType::PixelType currentLargestLabel = 0;

    std::map<itkIntImage2DType::PixelType, SampleType::Pointer> labelToIndexListMap;

    for (long iy = 0; iy < ny; ++iy)
      {
        idx[1] = iy;
        idxVec[1] = static_cast<VectorType::ValueType>(iy);
        for (long ix = 0; ix < nx; ++ix)
          {
            idx[0] = ix;
            idxVec[0] = static_cast<VectorType::ValueType>(ix);

            itkIntImage2DType::PixelType thisLabel = m_connectedComponentLabelImage->GetPixel(idx);
            currentLargestLabel = currentLargestLabel>thisLabel?currentLargestLabel:thisLabel;

            if (thisLabel)
              {
                long objectId = thisLabel - 1; ///< This "- 1" is the way the itkLabelImageToShapeLabelMapFilter works
                if (1 == m_objectToBreak[objectId]  && labelToIndexListMap.find(thisLabel) == labelToIndexListMap.end() )
                  {
                    labelToIndexListMap[thisLabel] = SampleType::New();
                  }
                else if (1 == m_objectToBreak[objectId] )
                  {
                    labelToIndexListMap[thisLabel]->PushBack(idxVec);
                  }
              }// thisLabel != 0
          }// ix
      }// iy


    /// Step 30. Use mean shift to cluster those index
    /// Step 40. For each label, assign to the largest-current-label + 1, then largest-current-label add by 1
    for (std::map<itkIntImage2DType::PixelType, SampleType::Pointer>::iterator it = labelToIndexListMap.begin(); it != labelToIndexListMap.end(); ++it)
      {
        //dbg
        std::cout<<"MS cluster "<<it->first<<" label. It has "<<it->second->Size()<<" points......... "<<std::flush;
        //dbg, end
        MeanshiftClusteringFilterType ms;
        ms.setRadius(20.0);
        ms.setInputPointSet(it->second);
        ms.update();
        std::vector<long> label = ms.getLabelOfPoints();
        long maxLabel = label[0];

        itk2DIndexType idx;
        for (std::size_t ip = 0; ip < label.size(); ++ip)
          {
            idx[0] = static_cast<itk2DIndexType::IndexValueType>(it->second->GetMeasurementVector( ip )[0]);
            idx[1] = static_cast<itk2DIndexType::IndexValueType>(it->second->GetMeasurementVector( ip )[1]);

            maxLabel = maxLabel>label[ip]?maxLabel:label[ip];

            m_connectedComponentLabelImage->SetPixel(idx, 1 + label[ip] + currentLargestLabel);
          }

        //dbg
        std::cout<<" It has "<<maxLabel + 1<<" centers"<<std::endl<<std::flush;
        //dbg, end

        currentLargestLabel = currentLargestLabel + maxLabel + 1;
      }

    typedef itk::RelabelComponentImageFilter<itkIntImage2DType, itkIntImage2DType> FilterType;
    FilterType::Pointer relabelFilter = FilterType::New();
    relabelFilter->SetInput(m_connectedComponentLabelImage);
    //relabelFilter->SetMinimumObjectSize(static_cast<FilterType::ObjectSizeType>(m_objectSizeThreshold/m_micronPerPixel/m_micronPerPixel)); // This command takes number of pixels as input
    relabelFilter->Update();
    m_connectedComponentLabelImage = relabelFilter->GetOutput();

    return;
  }



  void BinaryMaskAnalysisFilter::_colorObjectByFeature(unsigned char featureType)
  {
    m_featureColoredImage = itkFloatImage2DType::New();
    m_featureColoredImage->SetRegions(m_binaryMask->GetLargestPossibleRegion() );
    m_featureColoredImage->Allocate();
    m_featureColoredImage->CopyInformation(m_binaryMask);
    m_featureColoredImage->FillBuffer(0.0);

    m_featureColoredImage->SetSpacing(m_micronPerPixel/1000.); //SetSpacing is mm based, so /1000.

    long numPixels = m_featureColoredImage->GetLargestPossibleRegion().GetNumberOfPixels();

    //if (0 == m_connectedComponentLabelImage)
    if (!m_connectedComponentLabelImage)
      {
        _computeConnectedComponentsLabelImage();
      }

    const itkIntImage2DType::PixelType* m_connectedComponentLabelImageBufferPtr = m_connectedComponentLabelImage->GetBufferPointer();
    itkFloatImage2DType::PixelType* m_featureColoredImageBufferPtr = m_featureColoredImage->GetBufferPointer();


    if (featureType == 1)
      {
        if (m_objectAreas.empty())
          {
            _computeObjectFeatures();
          }

        /// color by size
        for (long it = 0; it < numPixels; ++it)
          {
            itkIntImage2DType::PixelType thisLabel = m_connectedComponentLabelImageBufferPtr[it];
            if (thisLabel)
              {
                long objectId = thisLabel - 1; ///< This "- 1" is the way the itkLabelImageToShapeLabelMapFilter works
                m_featureColoredImageBufferPtr[it] = m_objectAreas[objectId];
              }
          }
      }
    else if (featureType == 2)
      {
        if (m_objectPerimeters.empty())
          {
            _computeObjectFeatures();
          }

        /// color by perimeter^2/size
        for (long it = 0; it < numPixels; ++it)
          {
            itkIntImage2DType::PixelType thisLabel = m_connectedComponentLabelImageBufferPtr[it];
            if (thisLabel)
              {
                long objectId = thisLabel - 1; ///< This "- 1" is the way the itkLabelImageToShapeLabelMapFilter works
                m_featureColoredImageBufferPtr[it] = m_objectPerimeters[objectId]*m_objectPerimeters[objectId]/m_objectAreas[objectId];
              }
          }
      }

    return;
  }


  BinaryMaskAnalysisFilter::itkIntImage2DType::Pointer BinaryMaskAnalysisFilter::getConnectedComponentLabelImage()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not done.\n";
        abort();
      }

    return m_connectedComponentLabelImage;
  }

  itkFloatImage2DType::Pointer BinaryMaskAnalysisFilter::getFeatureColoredImage(unsigned char featureType)
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not done.\n";
        abort();
      }

    if (!m_featureColoredImage)
      {
        _colorObjectByFeature(featureType);
      }

    return m_featureColoredImage;
  }


  void BinaryMaskAnalysisFilter::_computeObjectFeatures()
  {
    // std::cout << "File " << "\"" << fileName << "\""
    //           << " has " << labelMap->GetNumberOfLabelObjects() << " labels." << std::endl;

    // Retrieve all attributes
    m_objectAreas.resize(m_numberOfObjects);
    m_objectPerimeters.resize(m_numberOfObjects);
    m_objectEquivalentSphericalRadius.resize(m_numberOfObjects);
    //m_objectNecessityOfBreaking.resize(m_numberOfObjects);

    for (unsigned int n = 0; n < m_numberOfObjects; ++n)
      {
        ShapeLabelObjectType *labelObject = m_labelMap->GetNthLabelObject(n);

        // std::cout<<n<<"-th object's" << " Label: "
        //           << itk::NumericTraits<LabelMapType::LabelType>::PrintType(labelObject->GetLabel()) << std::endl;
        // std::cout << "    BoundingBox: "
        //           << labelObject->GetBoundingBox() << std::endl;

        //m_objectAreas[n] = m_micronPerPixel*m_micronPerPixel*static_cast<float>(labelObject->GetNumberOfPixels());
        m_objectAreas[n] = labelObject->GetPhysicalSize();
        //m_objectAreas[n] = static_cast<float>(labelObject->GetNumberOfPixels());
        m_objectPerimeters[n]  = labelObject->GetPerimeter();
        m_objectEquivalentSphericalRadius[n] = labelObject->GetEquivalentSphericalRadius();
        //m_objectNecessityOfBreaking[n] = m_objectPerimeters[n]*m_objectPerimeters[n]/m_objectAreas[n];

        // m_fileForOutputNucluesFeatures << labelObject->GetNumberOfPixels() << ",";
        // m_fileForOutputNucluesFeatures << m_micronPerPixel*m_micronPerPixel*static_cast<float>(labelObject->GetNumberOfPixels()) << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetNumberOfPixelsOnBorder() << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetFeretDiameter() << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetPrincipalMoments()[0] << "," << labelObject->GetPrincipalMoments()[1] << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetElongation() << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetRoundness() << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetEquivalentSphericalRadius() << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetEquivalentSphericalPerimeter() << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetEquivalentEllipsoidDiameter()[0] << "," << labelObject->GetEquivalentEllipsoidDiameter()[1] << ",";
        // m_fileForOutputNucluesFeatures << labelObject->GetFlatness() << std::endl << std::flush;
      }

    ////////////////////////////////////////////////////////////////////////////////
    /// Now save the resuts:
    /// Only save on average 1/100 of tiles to save space and perform random check if needed
    // cv::RNG rng;
    // rng.state = cv::getTickCount();
    // double saveOrNot = rng.uniform(0.0, 1.0);

    // //if (saveOrNot < 0.01)
    //   {
    //     // char outputName[1000];
    //     // sprintf(outputName, "%s_appMag_%d_%ld_%ld-tile.png", m_outputPrefix.c_str(), m_magnification, m_currentTileIndexW, m_currentTileIndexH);
    //     // //cv::imwrite(outputName, m_currentTile); ///< Save the raw-tile before normalization
    //     // std::cout<<outputName<<std::endl<<std::flush;
    //     // writeImage<itkRGBImage2DType>(m_currentTile, outputName);

    //     char outputLabelName[1000];
    //     sprintf(outputLabelName, "%s_appMag_%d_%ld_%ld-label.nrrd", m_outputPrefix.c_str(), m_magnification, m_currentTileIndexW, m_currentTileIndexH);
    //     //cv::imwrite(outputLabelName, labelImage);
    //     std::cout<<outputLabelName<<std::endl<<std::flush;
    //     //writeImage<itkBinaryMaskImage2DType>(m_binaryMaskOfCurrentTile, outputLabelName);
    //     writeImage<itkIntImage2DType>(m_labelImageOfCurrentTile, outputLabelName);
    //     std::cout<<"done"<<std::endl<<std::flush;
    //   }

    //   std::cout<<"Done\n"<<std::endl<<std::flush;


    // /// Save the nuclei sizes
    // FILE* pFile;

    // char outputNucleiSizeName[1000];
    // sprintf(outputNucleiSizeName, "%s_appMag_%d_%ld_%ld.features", m_outputPrefix.c_str(), m_magnification, m_currentTileIndexW, m_currentTileIndexH);

    // pFile = fopen(outputNucleiSizeName, "w");

    // // start with size 1 coz do not output the size of background

    // for (itkUIntImage2DType::PixelType it = 1; it < m_totalNumberOfConnectedComponentsInThisTile; ++it)
    //   {
    //     fprintf(pFile, "%u, %ld, %f\n", it, m_sizesOfNucleiCurrentTile[it], m_roundnessOfNucleiCurrentTile[it]);
    //   }

    // fclose (pFile);
    // /// Now save the resuts
    // ////////////////////////////////////////////////////////////////////////////////

    return;
  }

  void BinaryMaskAnalysisFilter::setMaskImage(const itkBinaryMaskImage2DType* img)
  {
    m_inputImage = img;

    typedef itk::BinaryThresholdImageFilter<itkBinaryMaskImage2DType, itkBinaryMaskImage2DType> ThresholdingFilterType;
    ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
    thresholder->SetLowerThreshold( 1 ); ///< the filter's rule is <=, i.e. lower <= x <= upper
    thresholder->SetUpperThreshold( itk::NumericTraits< itkBinaryMaskImage2DType::PixelType >::max() );

    thresholder->SetOutsideValue(  0  );
    thresholder->SetInsideValue(  1 );
    thresholder->SetInput( m_inputImage );
    thresholder->Update();

    m_binaryMask = thresholder->GetOutput();

    return;
  }


  void BinaryMaskAnalysisFilter::setMPP(float mpp)
  {
    if (mpp > 0)
      {
        m_micronPerPixel = mpp;
      }
    else
      {
        std::cerr<<"Error: mpp should be > 0. But got "<<mpp<<std::endl;
      }

    return;
  }

  void BinaryMaskAnalysisFilter::_computeConnectedComponentsLabelImage()
  {
    typedef itk::ConnectedComponentImageFilter <itkBinaryMaskImage2DType, itkIntImage2DType > ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
    connected->SetInput(m_binaryMask);
    connected->Update();

    typedef itk::RelabelComponentImageFilter<itkIntImage2DType, itkIntImage2DType> FilterType;
    FilterType::Pointer relabelFilter = FilterType::New();
    relabelFilter->SetInput(connected->GetOutput());
    //std::cout<<static_cast<FilterType::ObjectSizeType>(m_objectSizeThreshold/m_micronPerPixel/m_micronPerPixel)<<std::endl;
    relabelFilter->SetMinimumObjectSize(static_cast<FilterType::ObjectSizeType>(m_objectSizeThreshold/m_micronPerPixel/m_micronPerPixel)); // This command takes number of pixels as input
    relabelFilter->Update();

    m_connectedComponentLabelImage = relabelFilter->GetOutput();

    return;
  }

  void BinaryMaskAnalysisFilter::_computeLabelMap()
  {
    I2LType::Pointer i2l = I2LType::New();
    i2l->SetInput( m_connectedComponentLabelImage );
    i2l->SetComputePerimeter(true);
    i2l->Update();

    m_labelMap = i2l->GetOutput();

    m_numberOfObjects = m_labelMap->GetNumberOfLabelObjects();

    return;
  }

  // void BinaryMaskAnalysisFilter::_computeObjectNecessityOfBreakingValues()
  // {
  //   m_objectAreas.resize(m_numberOfObjects);
  //   m_objectPerimeters.resize(m_numberOfObjects);
  //   m_objectEquivalentSphericalRadius.resize(m_numberOfObjects);

  //   for (unsigned int n = 0; n < m_numberOfObjects; ++n)
  //     {

  //     }

  //   //m_objectNecessityOfBreaking.
  // }

  const std::vector<double>& BinaryMaskAnalysisFilter::getObjectAreas()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not done.\n";
        abort();
      }

    return m_objectAreas;
  }

  const std::vector<double>& BinaryMaskAnalysisFilter::getObjectPerimeters()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not done.\n";
        abort();
      }

    return m_objectPerimeters;
  }

  const std::vector<double>& BinaryMaskAnalysisFilter::getObjectEquivalentSphericalRadius()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: computation not done.\n";
        abort();
      }

    return m_objectEquivalentSphericalRadius;
  }





}// namespace gth818n
