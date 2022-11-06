#ifndef SFLSSegmentor2D_h_
#define SFLSSegmentor2D_h_

#include "SFLS.h"

#include <list>

//itk
#include "itkImage.h"

template< typename TPixel >
class CSFLSSegmentor2D : public CSFLS
{
public:
  typedef CSFLSSegmentor2D< TPixel > Self;

  typedef CSFLS SuperClassType;
  typedef SuperClassType::NodeType NodeType;
  typedef SuperClassType::CSFLSLayer CSFLSLayer;

  typedef itk::Image<TPixel, 2> ImageType;
  typedef itk::Image<double, 2> LSImageType;
  typedef itk::Image<char, 2> LabelImageType;
  typedef itk::Image<unsigned char, 2> MaskImageType;

  CSFLSSegmentor2D();
  virtual ~CSFLSSegmentor2D() {}

  /* ============================================================
   * functions         */
  void basicInit();

  void setNumIter(unsigned long n);

  void setImage(typename ImageType::Pointer img);
  void setMask(typename MaskImageType::Pointer mask);

  virtual void computeForce() = 0;

  void normalizeForce();

  //     double maxPhi(long ix, long iy, long iz, double level);
  //     double minPhi(long ix, long iy, long iz, double level);
  bool getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(long ix, long iy, double& thePhi);

  void oneStepLevelSetEvolution();

  void getSFLSFromPhi();

  void initializeSFLS() { initializeSFLSFromMask(); }
  void initializeSFLSFromMask();

  void initializeLabel();
  void initializePhi();

  virtual void doSegmenation() = 0;

  typename MaskImageType::Pointer getSegmentationMask(typename LSImageType::PixelType levelsetThreshold = 0);

  // geometry
  double computeKappa(long ix, long iy);

  void setCurvatureWeight(double a);


  /* ============================================================
   * data     */
  //    CSFLS::Pointer mp_sfls;

  typename ImageType::Pointer mp_img;
  typename LabelImageType::Pointer mp_label;
  typename MaskImageType::Pointer mp_mask;
  typename LSImageType::Pointer mp_phi;

protected:
  double m_curvatureWeight;

  long m_nx;
  long m_ny;


  std::vector< double > m_force;

  double m_timeStep;

  unsigned long m_numIter;


  inline bool doubleEqual(double a, double b, double eps = 1e-10)
  {
    return (a-b < eps && b-a < eps);
  }


  /*----------------------------------------------------------------------
    These two record the pts which change status 

    Because they are created and visited sequentially, and when not
    needed, are clear-ed as a whole. No random insertion or removal is
    needed. So use vector is faster than list.  */ 
  CSFLSLayer m_lIn2out;
  CSFLSLayer m_lOut2in;


  //     //debug//
  //     void labelsCoherentCheck();
  //     void labelsCoherentCheck1();
  //     ////debug////
};



#include "SFLSSegmentor2D.hxx"


#endif
