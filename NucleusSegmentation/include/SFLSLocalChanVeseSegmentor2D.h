#ifndef SFLSLocalChanVeseSegmentor2D_h_
#define SFLSLocalChanVeseSegmentor2D_h_

#include "SFLSSegmentor2D.h"

#include <list>


template< typename TPixel >
class CSFLSLocalChanVeseSegmentor2D : public CSFLSSegmentor2D< TPixel >
{
public:
  typedef CSFLSSegmentor2D< TPixel > SuperClassType;

  //    typedef boost::shared_ptr< CSFLSLocalChanVeseSegmentor2D< TPixel > > Pointer;

  typedef typename SuperClassType::NodeType NodeType;
  typedef typename SuperClassType::CSFLSLayer CSFLSLayer; 


 /*================================================================================
    ctor */
  CSFLSLocalChanVeseSegmentor2D() : CSFLSSegmentor2D< TPixel >()
  {
    basicInit();
  }

  void basicInit();

  void setNBHDSize(long nbx, long nby)
  {
    m_nbx = nbx;
    m_nby = nby;
  }


  //     /* ============================================================
  //        New    */
  //     static Pointer New() 
  //     {
  //       return Pointer(new CSFLSLocalChanVeseSegmentor2D< TPixel >);
  //     }

  // data
  double m_areaIn;
  double m_areaOut;

  double m_meanIn;
  double m_meanOut;

  long m_nbx, m_nby;

  /* ============================================================
   * functions
   * ============================================================*/
  void computeMeans();
  void computeMeansAt(long ix, long iy);
  //  void updateMeans();

  //void doChanVeseSegmenation();
  void doSegmenation();


  /* ============================================================
     computeForce    */
  void computeForce();

  void setInflation(float f) {m_globalInflation = f;}


private:
  float m_globalInflation;

};


#include "SFLSLocalChanVeseSegmentor2D.hxx"

#endif
