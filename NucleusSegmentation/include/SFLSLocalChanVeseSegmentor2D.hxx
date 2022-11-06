#ifndef SFLSLocalChanVeseSegmentor2D_hpp_
#define SFLSLocalChanVeseSegmentor2D_hpp_

#include "SFLSLocalChanVeseSegmentor2D.h"

#include <algorithm>


/* ============================================================
   basicInit    */
template< typename TPixel >
void
CSFLSLocalChanVeseSegmentor2D< TPixel >
::basicInit()
{
  SuperClassType::basicInit();

  m_nbx = 5;
  m_nby = 5;

  m_globalInflation = 0.0; // pos: inflation; neg: contraction

}


/* ============================================================
   computeForce    */
template< typename TPixel >
void
CSFLSLocalChanVeseSegmentor2D< TPixel >
::computeForce()
{
    this->m_force.clear();

    double fmax = -1e10;

    long n = this->m_lz.size();
    double* kappaOnZeroLS = new double[ n ];
    double* cvForce = new double[ n ];

  {
    long i = 0;
    for (typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz, ++i)
      {
        long ix = (*itz)[0];
        long iy = (*itz)[1];

        typename itk::Image<TPixel, 2>::IndexType idx = {{ix, iy}};

        kappaOnZeroLS[i] = this->computeKappa(ix, iy);

        computeMeansAt(ix, iy);

        double I = this->mp_img->GetPixel(idx);
        double a = (I - m_meanIn)*(I - m_meanIn) - (I - m_meanOut)*(I - m_meanOut) - m_globalInflation;

        fmax = fabs(a)>fmax?fabs(a):fmax;

        cvForce[i] = a;
      }
  }


  for (long i = 0; i < n; ++i)
    {
      this->m_force.push_back( cvForce[i]/(fmax + 1e-10) +  (this->m_curvatureWeight)*kappaOnZeroLS[i]);
    }


  delete[] kappaOnZeroLS;
  delete[] cvForce;
}


/* ============================================================
   doSegmenation    */
template< typename TPixel >
void
CSFLSLocalChanVeseSegmentor2D< TPixel >
::doSegmenation()
{
  /*============================================================
   * From the initial mask, generate: 1. SFLS, 2. mp_label and
   * 3. mp_phi.
   */
  this->initializeSFLS();

  //computeMeans();

  //gth818n::saveAsImage2< double >(mp_phi, "initPhi.nrrd");
  for (unsigned int it = 0; it < this->m_numIter; ++it)
    {
      computeForce();

      this->normalizeForce();

      this->oneStepLevelSetEvolution();
    }
}


/* ============================================================
   computeMeansAt    */
template< typename TPixel >
void
CSFLSLocalChanVeseSegmentor2D< TPixel >
::computeMeansAt(long ix, long iy)
{
  /*----------------------------------------------------------------------
    Compute the local meanIn/Out areaIn/Out at this pixel. */

  m_areaIn = 0;
  m_areaOut = 0;

  m_meanIn = 0;
  m_meanOut = 0;

  for (long iix = ix-m_nbx; iix <= ix+m_nbx; ++iix)
    {
      for (long iiy = iy-m_nby; iiy <= iy+m_nby; ++iiy)
        {
          if (iix >= 0 && iix < this->m_nx && iiy >= 0 && iiy < this->m_ny)
            {
              typename itk::Index<2> idx = {{iix, iiy}};

              TPixel imgVal = (this->mp_img)->GetPixel(idx);
              double phi = (this->mp_phi)->GetPixel(idx);

              if (phi <= 0)
                {
                  // in
                  ++m_areaIn;
                  m_meanIn += imgVal;

                }
              else 
                {
                  ++m_areaOut;
                  m_meanOut += imgVal;;
                }
            }
        }
    }

  m_meanIn /= (m_areaIn + vnl_math::eps);
  m_meanOut /= (m_areaOut + vnl_math::eps);

  return;
}


/* ============================================================
   computeMeans    */
template< typename TPixel >
void
CSFLSLocalChanVeseSegmentor2D< TPixel >
::computeMeans()
{
  m_areaIn = 0;
  m_areaOut = 0;

  m_meanIn = 0;
  m_meanOut = 0;

  for (long ix = 0; ix < this->m_nx; ++ix)
    {
      for (long iy = 0; iy < this->m_ny; ++iy)
        {
          typename itk::Image<TPixel, 2>::IndexType idx = {{ix, iy}};

          double phi = this->mp_phi->GetPixel(idx);
          double v = this->mp_img->GetPixel(idx);

          if (phi <= 0)
            {
              ++m_areaIn;
              m_meanIn += v;
            }
          else
            {
              ++m_areaOut;
              m_meanOut += v;
            }

        }
    }

  m_meanIn /= (m_areaIn + vnl_math::eps);
  m_meanOut /= (m_areaOut + vnl_math::eps);
}


// /* ============================================================
//    updateMeans    */
// template< typename TPixel >
// void
// CSFLSLocalChanVeseSegmentor2D< TPixel >
// ::updateMeans()
// {
//   double sumIn = m_meanIn*m_areaIn;
//   double sumOut = m_meanOut*m_areaOut;

//   for (typename CSFLSLayer::const_iterator it = this->m_lIn2out.begin(); it != this->m_lIn2out.end(); ++it)
//     {
//       long ix = (*it)[0];
//       long iy = (*it)[1];

//       typename itk::Image<TPixel, 2>::IndexType idx = {{ix, iy}};

//       sumIn  -= this->mp_img->GetPixel(idx);
//       --m_areaIn;

//       sumOut += this->mp_img->GetPixel(idx);
//       ++m_areaOut;
//     }

//   for (typename CSFLSLayer::const_iterator it = this->m_lOut2in.begin(); it != this->m_lOut2in.end(); ++it)
//     {
//       long ix = (*it)[0];
//       long iy = (*it)[1];

//       typename itk::Image<TPixel, 2>::IndexType idx = {{ix, iy}};

//       sumIn  += this->mp_img->GetPixel(idx);
//       ++m_areaIn;

//       sumOut -= this->mp_img->GetPixel(idx);
//       --m_areaOut;
//     }


//   m_meanIn = (sumIn/m_areaIn + vnl_math::eps);
//   m_meanOut = (sumOut/m_areaOut + vnl_math::eps);
// }




#endif
