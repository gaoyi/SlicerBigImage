#include <cstdio>
#include <iostream>
#include <string>
#include <algorithm>

// local
#include "HAndEImageAnalysisFilter.h"

#include "itkTypedefs.h"
//#include "Image/gth818nImage.h"
#include "include/gth818nImageIO.h"
#include "utilitiesImage.h"
#include "ColorDecomposition.h"

#include "NucleusSegmentationCLP.h"

//#include "BinaryMaskAnalysisFilter.h"

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

  // itkUCharImage2DType::Pointer eosinChannelImage = gth818n::ExtractEosinChannel(img);
  // gth818n::writeImage<itkUCharImage2DType>(eosinChannelImage, outputEImage.c_str(), true);

  // itkUCharImage2DType::Pointer hematoxylinChannelImage = gth818n::ExtractHematoxylinChannel(img);
  // gth818n::writeImage<itkUCharImage2DType>(hematoxylinChannelImage, outputHImage.c_str(), true);

  //int magnification = 20; //hack should NOT be fixed

  gth818n::HAndEImageAnalysisFilter tileAnalyzer;
  tileAnalyzer.setInputHAndEImage(img);
  //tileAnalyzer.setMagnification(magnification);
  tileAnalyzer.setMPP(mpp);
  tileAnalyzer.setThreshold(threshold);
  tileAnalyzer.update();


  itkUCharImage2DType::Pointer seg = tileAnalyzer.getNucleiBinaryMaskImage();
  seg->SetSpacing(mpp/1000.);

  //gth818n::writeImage<itkUCharImage2DType>(seg, outputImageName.c_str());

  gth818n::writeImage<gth818n::HAndEImageAnalysisFilter::itkIntImage2DType>(tileAnalyzer.getNucleiLabelImage(), outputImageName.c_str());

  //std::cout<<tileAnalyzer.getTotalNumberOfConnectedComponents()<<std::endl;


  return EXIT_SUCCESS;
}
