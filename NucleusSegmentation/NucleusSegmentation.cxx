#include <cstdio>
#include <iostream>
#include <string>
#include <algorithm>

// local
#include "include/itkTypedefs.h"
#include "include/ColorDecomposition.h"
#include "include/utilitiesIO.h"
#include "include/utilitiesImage.h"
#include "include/HAndEImageAnalysisFilter.h"

// CLI
#include "NucleusSegmentationCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const int ImageDimension = 2;

  typedef gth818n::itkRGBImage2DType itkRGBImage2DType;
  itkRGBImage2DType::Pointer img = gth818n::readImage<itkRGBImage2DType>(inputfileName.c_str());

  typedef itk::Image<unsigned char, ImageDimension> itkUCharImage2DType;

  gth818n::HAndEImageAnalysisFilter tileAnalyzer;
  tileAnalyzer.setInputHAndEImage(img);
  tileAnalyzer.setMPP(mpp);
  tileAnalyzer.setThreshold(threshold);
  tileAnalyzer.update();


  itkUCharImage2DType::Pointer seg = tileAnalyzer.getNucleiBinaryMaskImage();
  seg->SetSpacing(mpp/1000.);

  //gth818n::writeImage<itkUCharImage2DType>(seg, outputImageName.c_str()); // output binary mask for all nuclei
  gth818n::writeImage<gth818n::HAndEImageAnalysisFilter::itkIntImage2DType>(tileAnalyzer.getNucleiLabelImage(), outputImageName.c_str());

  return EXIT_SUCCESS;
}
