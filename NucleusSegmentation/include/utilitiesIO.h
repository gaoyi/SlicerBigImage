#ifndef gth818nImage_h_
#define gth818nImage_h_


#include <string>

// itk
#include "itkImage.h"

// // vnl
// #include "vnl/vnl_matrix.h"


namespace gth818n
{
  /**********************************************************************************
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName);

  /************************************************************************************
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName, bool compress = true);


}// namespace gth818n


#include "utilitiesIO.hxx"

#endif
