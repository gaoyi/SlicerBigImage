#ifndef ColorDecomposition_h_
#define ColorDecomposition_h_


#include <string>

// itk
#include "itkImage.h"

// local
#include "itkTypedefs.h"

namespace gth818n
{
  ////////////////////////////////////////////////////////////////////////////////
  /// This is an implementation of the paper:
  ///
  /// Anal Quant Cytol Histol. 2001 Aug;23(4):291-9. Quantification of
  /// histochemical staining by color deconvolution. Ruifrok AC,
  /// Johnston DA.
  ///
  /// The implementation is adopted from the java ImageJ plugin at
  /// http://www.mecourse.com/landinig/software/cdeconv/cdeconv.html
  ///
  /// The ImageJ plugin takes an input image, then a pop-up box is
  /// given to ask user to determine the type of stain. It has many
  /// different choices. It output several channels, such as
  /// Hematoxylin, Eosin, etc., depending on the input stains.
  ///
  /// In this implementation, however, we only take the H&E stain as
  /// input type, and only output the Hematoxylin channel as output.
  itkUCharImage2DType::Pointer ExtractHematoxylinChannel(itkRGBImage2DType::Pointer HAndEImage);
  itkUCharImage2DType::Pointer ExtractEosinChannel(itkRGBImage2DType::Pointer HAndEImage);



}// namespace gth818n


#endif
