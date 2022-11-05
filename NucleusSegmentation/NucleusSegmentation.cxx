#include <cstdio>
#include <iostream>
#include <string>
#include <algorithm>

// local
#include "HAndEImageAnalysisFilter.h"

#include "itkTypedefs.h"
#include "Image/gth818nImage.h"
#include "utilitiesImage.h"
#include "ColorDecomposition.h"

#include "NucleusSegmentationCLP.h"

#include "BinaryMaskAnalysisFilter.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  //const int ImageDimension = 2;

  typedef gth818n::itkRGBImage2DType itkRGBImage2DType;
  itkRGBImage2DType::Pointer img = gth818n::readImage<itkRGBImage2DType>(inputfileName.c_str());

  typedef itk::Image<unsigned char, 2> itkUCharImage2DType;
  itkUCharImage2DType::Pointer eosinChannelImage = gth818n::ExtractEosinChannel(img);
  gth818n::writeImage<itkUCharImage2DType>(eosinChannelImage, outputEImage.c_str(), true);

  itkUCharImage2DType::Pointer hematoxylinChannelImage = gth818n::ExtractHematoxylinChannel(img);
  gth818n::writeImage<itkUCharImage2DType>(hematoxylinChannelImage, outputHImage.c_str(), true);


  int magnification = 20;

  gth818n::HAndEImageAnalysisFilter tileAnalyzer;
  tileAnalyzer.setInputHAndEImage(img);
  tileAnalyzer.setMagnification(magnification);
  tileAnalyzer.setThreshold(threshold);
  tileAnalyzer.update();

  //gth818n::writeImage<gth818n::HAndEImageAnalysisFilter::itkIntImage2DType>(tileAnalyzer.getNucleiLabelImage(), outputImageName.c_str());

  //std::cout<<tileAnalyzer.getTotalNumberOfConnectedComponents()<<std::endl;

  gth818n::BinaryMaskAnalysisFilter binaryMaskAnalyzer;
  binaryMaskAnalyzer.setMaskImage( tileAnalyzer.getNucleiBinaryMaskImage() );
  binaryMaskAnalyzer.setMPP(mpp);
  binaryMaskAnalyzer.update();

  unsigned char featureType = 1; ///< area of the objects;
  gth818n::writeImage<itk::Image<short, 2> >(gth818n::castItkImage<gth818n::BinaryMaskAnalysisFilter::itkFloatImage2DType, itk::Image<short, 2> >(binaryMaskAnalyzer.getFeatureColoredImage(featureType)), outputImageName.c_str());

  //binaryMaskAnalyzer.update();


  return EXIT_SUCCESS;
}
