#ifndef SFLSSegmentor2D_hpp_
#define SFLSSegmentor2D_hpp_

#include "SFLSSegmentor2D.h"

#include <algorithm>
#include <cmath>


#include <csignal>


template< typename TPixel >
CSFLSSegmentor2D< TPixel >
::CSFLSSegmentor2D() : CSFLS()
{
  basicInit(); 
}


/* ============================================================
   basicInit    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::basicInit()
{
  m_numIter = 100;
  m_timeStep = 1.0;

  m_curvatureWeight = 0.0;


  m_nx = 0;
  m_ny = 0;
}

/* ============================================================
   setImage    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::setNumIter(unsigned long n)
{
  m_numIter = n;
}

/* ============================================================
   setCurvatureWeight    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::setCurvatureWeight(double a) 
{
  if (a < 0)
    {
      std::cerr<<"Error: curvature weight < 0\n";
      raise(SIGABRT);
    }

  m_curvatureWeight = a;

  return;
}

/* ============================================================
   setImage    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::setImage(typename ImageType::Pointer img)
{
  mp_img = img;

  typename ImageType::IndexType start = mp_img->GetLargestPossibleRegion().GetIndex();
  typename ImageType::IndexType origin = {{0, 0}};
  if (start != origin)
    {
      std::cout<<"Warrning: Force image start to be (0, 0)\n";

      typename ImageType::RegionType region = mp_img->GetLargestPossibleRegion();
      region.SetIndex(origin);

      mp_img->SetRegions(region);
    }


  typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();

  if (m_nx + m_ny == 0)
    {
      m_nx = size[0];
      m_ny = size[1];
    }
  else if ( m_nx != (long)size[0] || m_ny != (long)size[1] )
    {
      std::cerr<<"image sizes do not match. abort\n";
      raise(SIGABRT);
    }
}

/* ============================================================
   setMask    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::setMask(typename MaskImageType::Pointer mask)
{
  mp_mask = mask;

  typename ImageType::IndexType start = mp_mask->GetLargestPossibleRegion().GetIndex();
  typename ImageType::IndexType origin = {{0, 0}};
  if (start != origin)
    {
      std::cout<<"Warrning: Force mask start to be (0, 0)\n";

      typename ImageType::RegionType region = mp_mask->GetLargestPossibleRegion();
      region.SetIndex(origin);

      mp_mask->SetRegions(region);
    }


  typename ImageType::SizeType size = mask->GetLargestPossibleRegion().GetSize();

  if (m_nx + m_ny == 0)
    {
      m_nx = size[0];
      m_ny = size[1];
    }
  else if ( m_nx != (long)size[0] || m_ny != (long)size[1] )
    {
      std::cerr<<"image sizes do not match. abort\n";
      raise(SIGABRT);
    }
}


template< typename TPixel >
//typename CSFLSSegmentor2D< TPixel >::MaskImageType::Pointer
typename itk::Image<unsigned char, 2>::Pointer
CSFLSSegmentor2D< TPixel >::getSegmentationMask(typename LSImageType::PixelType levelsetThreshold)
{
  typename MaskImageType::Pointer seg = MaskImageType::New();
  seg->SetRegions(mp_phi->GetLargestPossibleRegion() );
  seg->Allocate();
  seg->CopyInformation(mp_phi);
  seg->FillBuffer(0);

  std::size_t numPixels = mp_phi->GetLargestPossibleRegion().GetNumberOfPixels();
  typename MaskImageType::PixelType* segBufferPointer = seg->GetBufferPointer();
  const typename LSImageType::PixelType* phiBufferPointer = mp_phi->GetBufferPointer();

  for (std::size_t it = 0; it < numPixels; ++it)
    {
      segBufferPointer[it] = phiBufferPointer[it] <= levelsetThreshold?1:0;
    }

  return seg;
}

template< typename TPixel >
bool
CSFLSSegmentor2D< TPixel >
::getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(long ix, long iy, double& thePhi)
{
  /*--------------------------------------------------
   *
   * Look in all the neighbors, to find the phi value of the nbhd:
   * this nbhd should satisfy: 1. its layer is strictly closer to
   * the zero layer hence its value is thought to be updated. 2. If
   * there are several nbhd's belonging to the same layer, choose
   * the one whose phi value has the smallest abs value.  If (ix,
   * iy) is outside, go through all nbhd who is in the layer of
   * label = mylevel-1 pick the SMALLEST phi. If (ix, iy) is inside,
   * go through all nbhd who is in the layer of label = mylevel+1
   * pick the LARGEST phi.
   */

  typename ImageType::IndexType idx = {{ix, iy}};
  char mylevel = mp_label->GetPixel(idx);

  bool foundNbhd = false;

  if (mylevel > 0)
    {
      // find the SMALLEST phi
      thePhi = 10000;

      typename ImageType::IndexType idx3 = {{ix, iy+1}};
      if ( (iy+1 < m_ny) && (mp_label->GetPixel(idx3) == mylevel - 1) )
        {
          double itsPhi = mp_phi->GetPixel(idx3);
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};
      if ( ((iy-1) >= 0 ) && (mp_label->GetPixel(idx4) == mylevel - 1) ) 
        {
          double itsPhi = mp_phi->GetPixel(idx4);
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      typename ImageType::IndexType idx1 = {{ix+1, iy}};
      if ( (ix+1 < m_nx) && (mp_label->GetPixel(idx1) == mylevel - 1) )
        {
          double itsPhi = mp_phi->GetPixel(idx1);
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      typename ImageType::IndexType idx2 = {{ix-1, iy}};
      if ( (ix-1 >= 0 ) && (mp_label->GetPixel(idx2) == mylevel - 1) )
        {
          double itsPhi = mp_phi->GetPixel(idx2);
          thePhi = thePhi<itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }
    }
  else
    {
      // find the LARGEST phi
      thePhi = -10000;

      typename ImageType::IndexType idx3 = {{ix, iy+1}};
      if ( (iy+1 < m_ny) && (mp_label->GetPixel(idx3) == mylevel + 1) )
        {
          double itsPhi = mp_phi->GetPixel(idx3);
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};
      if ( ((iy-1) >= 0 ) && (mp_label->GetPixel(idx4) == mylevel + 1) ) 
        {
          double itsPhi = mp_phi->GetPixel(idx4);
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }


      typename ImageType::IndexType idx1 = {{ix+1, iy}};
      if ( (ix+1 < m_nx) && (mp_label->GetPixel(idx1) == mylevel + 1) )
        {
          double itsPhi = mp_phi->GetPixel(idx1);
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }


      typename ImageType::IndexType idx2 = {{ix-1, iy}};
      if ( (ix-1 >= 0  ) && (mp_label->GetPixel(idx2) == mylevel + 1) )
        {
          double itsPhi = mp_phi->GetPixel(idx2);
          thePhi = thePhi>itsPhi?thePhi:itsPhi;

          foundNbhd = true;
        }
    }



  return foundNbhd;
}


/* ============================================================
   normalizeForce   
   Normalize m_force s.t. max(abs(m_force)) < 0.5 */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::normalizeForce()
{
  unsigned long nLz = m_lz.size();

  if (m_force.size() != nLz )
    {
      std::cerr<<"m_force.size() = "<<m_force.size()<<std::endl;
      std::cerr<<"nLz = "<<nLz<<std::endl;

      std::cerr<<"m_force.size() != nLz, abort.\n";
      raise(SIGABRT);
    }

  double fMax = fabs( m_force.front() );

  {
    long nf = m_force.size();
    for (long itf = 0; itf < nf; ++itf)
      {
        double v = fabs(m_force[itf]);
        fMax = fMax>v?fMax:v;
      }
  }

  fMax /= 0.49;

  {
    long nf = m_force.size();

    for (long itf = 0; itf < nf; ++itf)
      {
        m_force[itf] /= (fMax + 1e-10);
      }
  }
}

/* ============================================================
   oneStepLevelSetEvolution    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::oneStepLevelSetEvolution()
{
  // create 'changing status' lists
  CSFLSLayer Sz;
  CSFLSLayer Sn1;
  CSFLSLayer Sp1;
  CSFLSLayer Sn2;
  CSFLSLayer Sp2;

  m_lIn2out.clear();
  m_lOut2in.clear();

  /*--------------------------------------------------
    1. add F to phi(Lz), create Sn1 & Sp1 
    scan Lz values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========                */
  {
    std::vector<double>::const_iterator itf = m_force.begin();
    for (CSFLSLayer::iterator itz = m_lz.begin(); \
         itz != m_lz.end(); \
         ++itf)
      {
        long ix = (*itz)[0];
        long iy = (*itz)[1];

        typename ImageType::IndexType idx = {{ix, iy}};

        double phi_old = mp_phi->GetPixel(idx);
        double phi_new = phi_old + (*itf);

        /*----------------------------------------------------------------------
          Update the lists of pt who change the state, for faster
          energy fnal computation. */
        if ( phi_old <= 0 && phi_new > 0 )
          {
            m_lIn2out.push_back(NodeType(ix, iy, 0));
          }

        if( phi_old>0  && phi_new <= 0)
          {
            m_lOut2in.push_back(NodeType(ix, iy, 0));
          }


        //           // DEBUG
        //           if (phi_new > 3.1 || phi_new < -3.1)
        //             {
        //               std::cout<<"phi_old = "<<phi_old<<std::endl;
        //               std::cout<<"its lbl = "<<(int)mp_label->get(ix, iy)<<std::endl;

        //               std::cerr<<"phi_new > 3.1 || phi_new < -3.1\n";
        //               raise(SIGABRT);
        //             }

        mp_phi->SetPixel(idx, phi_new);

        if(phi_new > 0.5)
          {
            Sp1.push_back(*itz);
            itz = m_lz.erase(itz);
          }
        else if (phi_new < -0.5)
          {
            Sn1.push_back(*itz);
            itz = m_lz.erase(itz);
          }
        else
          {
            ++itz;              
          }

        /*--------------------------------------------------
          NOTE, mp_label are (should) NOT update here. They should
          be updated with Sz, Sn/p's
          --------------------------------------------------*/
      }
  }


  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    2. update Ln1,Lp1,Lp2,Lp2, ****in that order****

    2.1 scan Ln1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5] 
    ==========                     */
  for (CSFLSLayer::iterator itn1 = m_ln1.begin(); itn1 != m_ln1.end(); )
    {
      long ix = (*itn1)[0];
      long iy = (*itn1)[1];

      typename ImageType::IndexType idx = {{ix, iy}};

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, thePhi);

      if (found)
        {
          double phi_new = thePhi-1;
          mp_phi->SetPixel(idx, phi_new);

          if (phi_new >= -0.5)
            {
              Sz.push_back(*itn1);
              itn1 = m_ln1.erase(itn1);
            }
          else if (phi_new < -1.5)
            {
              Sn2.push_back(*itn1);
              itn1 = m_ln1.erase(itn1);
            }
          else
            {
              ++itn1;
            }
        }
      else 
        {
          /*--------------------------------------------------
            No nbhd in inner (closer to zero contour) layer, so
            should go to Sn2. And the phi shold be further -1
          */
          Sn2.push_back(*itn1);
          itn1 = m_ln1.erase(itn1);

          mp_phi->SetPixel(idx, mp_phi->GetPixel(idx) - 1);
        }
    }



  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    2.2 scan Lp1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5] 
    ========          */
  for (CSFLSLayer::iterator itp1 = m_lp1.begin(); itp1 != m_lp1.end();)
    {
      long ix = (*itp1)[0];
      long iy = (*itp1)[1];

      typename ImageType::IndexType idx = {{ix, iy}};

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, thePhi);

      if (found)
        {
          double phi_new = thePhi+1;
          mp_phi->SetPixel(idx, phi_new);

          if (phi_new <= 0.5)
            {
              Sz.push_back(*itp1);
              itp1 = m_lp1.erase(itp1);
            }
          else if (phi_new > 1.5)
            {
              Sp2.push_back(*itp1);
              itp1 = m_lp1.erase(itp1);
            }
          else
            {
              ++itp1;
            }
        }
      else 
        {
          /*--------------------------------------------------
            No nbhd in inner (closer to zero contour) layer, so
            should go to Sp2. And the phi shold be further +1
          */

          Sp2.push_back(*itp1);
          itp1 = m_lp1.erase(itp1);

          mp_phi->SetPixel(idx, mp_phi->GetPixel(idx) + 1);  
        }
    }


  //     // debug
  //     labelsCoherentCheck1();



  /*--------------------------------------------------
    2.3 scan Ln2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5] 
    ==========                                      */
  for (CSFLSLayer::iterator itn2 = m_ln2.begin(); itn2 != m_ln2.end(); )
    {
      long ix = (*itn2)[0];
      long iy = (*itn2)[1];

      typename ImageType::IndexType idx = {{ix, iy}};

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, thePhi);

      if (found)
        {
          double phi_new = thePhi-1;
          mp_phi->SetPixel(idx, phi_new);

          if (phi_new >= -1.5)
            {
              Sn1.push_back(*itn2);
              itn2 = m_ln2.erase(itn2);
            }
          else if (phi_new < -2.5)
            {
              itn2 = m_ln2.erase(itn2);
              mp_phi->SetPixel(idx, -3);
              mp_label->SetPixel(idx, -3);
            }
          else
            {
              ++itn2;
            }
        }
      else 
        {
          itn2 = m_ln2.erase(itn2);
          mp_phi->SetPixel(idx, -3);
          mp_label->SetPixel(idx, -3);
        }
    }


  //     // debug
  //     labelsCoherentCheck1();



  /*--------------------------------------------------
    2.4 scan Lp2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========= */
  for (CSFLSLayer::iterator itp2 = m_lp2.begin(); itp2 != m_lp2.end(); )
    {
      long ix = (*itp2)[0];
      long iy = (*itp2)[1];

      typename ImageType::IndexType idx = {{ix, iy}};

      double thePhi;
      bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, thePhi);

      if (found)
        {
          double phi_new = thePhi+1;
          mp_phi->SetPixel(idx, phi_new);

          if (phi_new <= 1.5)
            {
              Sp1.push_back(*itp2);
              itp2 = m_lp2.erase(itp2);
            }
          else if (phi_new > 2.5)
            {
              itp2 = m_lp2.erase(itp2);
              mp_phi->SetPixel(idx, 3);
              mp_label->SetPixel(idx, 3);
            }
          else
            {
              ++itp2;
            }
        }
      else 
        {
          itp2 = m_lp2.erase(itp2);
          mp_phi->SetPixel(idx, 3);
          mp_label->SetPixel(idx, 3);
        }
    }


  //     // debug
  //     labelsCoherentCheck1();



  /*--------------------------------------------------
    3. Deal with S-lists Sz,Sn1,Sp1,Sn2,Sp2 
    3.1 Scan Sz */
  for (CSFLSLayer::iterator itSz = Sz.begin(); itSz != Sz.end(); ++itSz)
    {
      long ix = (*itSz)[0];
      long iy = (*itSz)[1];
        
      typename ImageType::IndexType idx = {{ix, iy}};
      
      m_lz.push_back(*itSz);
      mp_label->SetPixel(idx, 0);
    }


  //     // debug
  //     labelsCoherentCheck1();



  /*--------------------------------------------------
    3.2 Scan Sn1     */
  for (CSFLSLayer::iterator itSn1 = Sn1.begin(); itSn1 != Sn1.end(); ++itSn1)
    {
      long ix = (*itSn1)[0];
      long iy = (*itSn1)[1];
        
      typename ImageType::IndexType idx = {{ix, iy}};

      m_ln1.push_back(*itSn1);
      // itSn1 = Sn1.erase(itSn1);

      mp_label->SetPixel(idx, -1);

      typename ImageType::IndexType idx1 = {{ix+1, iy}};
      if ( (ix+1 < m_nx) && doubleEqual(mp_phi->GetPixel(idx1), -3.0) )
        {
          Sn2.push_back(NodeType(ix+1, iy, 0));
          mp_phi->SetPixel(idx1, mp_phi->GetPixel(idx) - 1 ); 
        }

      typename ImageType::IndexType idx2 = {{ix-1, iy}};
      if ( (ix-1 >= 0) && doubleEqual(mp_phi->GetPixel(idx2), -3.0) )
        {
          Sn2.push_back(NodeType(ix-1, iy, 0));
          mp_phi->SetPixel(idx2, mp_phi->GetPixel(idx) - 1); 
        }

      typename ImageType::IndexType idx3 = {{ix, iy+1}};
      if ( (iy+1 < m_ny) && doubleEqual(mp_phi->GetPixel(idx3), -3.0) )
        {
          Sn2.push_back(NodeType(ix, iy+1, 0));
          mp_phi->SetPixel(idx3, mp_phi->GetPixel(idx) - 1 ); 
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};
      if ( (iy-1>=0) && doubleEqual(mp_phi->GetPixel(idx4), -3.0) )
        {
          Sn2.push_back(NodeType(ix, iy-1, 0) );
          mp_phi->SetPixel(idx4, mp_phi->GetPixel(idx) - 1 );
        }
    }


  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    3.3 Scan Sp1     */
  for (CSFLSLayer::iterator itSp1 = Sp1.begin(); itSp1 != Sp1.end(); ++itSp1)
    {
      long ix = (*itSp1)[0];
      long iy = (*itSp1)[1];

      typename ImageType::IndexType idx = {{ix, iy}};
        
      m_lp1.push_back(*itSp1);
      mp_label->SetPixel(idx, 1);

      typename ImageType::IndexType idx3 = {{ix, iy+1}};
      if ( (iy+1 < m_ny) && doubleEqual(mp_phi->GetPixel(idx3), 3.0) )
        {
          Sp2.push_back(NodeType(ix, iy+1, 0));
          mp_phi->SetPixel(idx3, mp_phi->GetPixel(idx) + 1 ); 
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};
      if ( (iy-1>=0) && doubleEqual(mp_phi->GetPixel(idx4), 3.0) )
        {
          Sp2.push_back(NodeType(ix, iy-1, 0) );
          mp_phi->SetPixel(idx4, mp_phi->GetPixel(idx) + 1 );
        }

      typename ImageType::IndexType idx1 = {{ix+1, iy}};
      if ( (ix+1 < m_nx) && doubleEqual(mp_phi->GetPixel(idx1), 3.0) )
        {
          Sp2.push_back(NodeType(ix+1, iy, 0));
          mp_phi->SetPixel(idx1, mp_phi->GetPixel(idx) + 1 ); 
        }

      typename ImageType::IndexType idx2 = {{ix-1, iy}};
      if ( (ix-1 >= 0) && doubleEqual(mp_phi->GetPixel(idx2), 3.0) )
        {
          Sp2.push_back(NodeType(ix-1, iy, 0));
          mp_phi->SetPixel(idx2, mp_phi->GetPixel(idx) + 1); 
        }
    } 


  //     // debug
  //     labelsCoherentCheck1();


  /*--------------------------------------------------
    3.4 Scan Sn2     */
  {
    //debug 
    int aaa = 0;
    for (CSFLSLayer::iterator itSn2 = Sn2.begin(); itSn2 != Sn2.end(); ++itSn2, ++aaa)
      {
        long ix = (*itSn2)[0];
        long iy = (*itSn2)[1];
        
        typename ImageType::IndexType idx = {{ix, iy}};

        m_ln2.push_back(*itSn2);

        mp_label->SetPixel(idx, -2);

        //           // debug
        //           labelsCoherentCheck1();
      }
  }



  /*--------------------------------------------------
    3.5 Scan Sp2     */
  for (CSFLSLayer::iterator itSp2 = Sp2.begin(); itSp2 != Sp2.end(); ++itSp2)
    {
      long ix = (*itSp2)[0];
      long iy = (*itSp2)[1];
        
      typename ImageType::IndexType idx = {{ix, iy}};

      m_lp2.push_back(*itSp2);

      mp_label->SetPixel(idx, 2);
    }


  //     // debug
  //     labelsCoherentCheck1();

}

/*================================================================================
  initializeLabel*/
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::initializeLabel()
{
  if (m_nx + m_ny == 0)
    {
      std::cerr<<"set mp_img first.\n";
      raise(SIGABRT);
    }

  //find interface and mark as 0, create Lz
  char defaultLabel = 0;

  mp_label = LabelImageType::New();

  LabelImageType::RegionType region = mp_img->GetLargestPossibleRegion();
  
  mp_label->SetRegions( region );
  mp_label->Allocate();
  mp_label->CopyInformation(mp_img);


  mp_label->FillBuffer(defaultLabel);

  return;
}


/*================================================================================
  initializePhi*/
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::initializePhi()
{
  if (m_nx + m_ny == 0)
    {
      std::cerr<<"set mp_img first.\n";
      raise(SIGABRT);
    }

  double arbitraryInitPhi = 1000;

  mp_phi = LSImageType::New();
  LSImageType::IndexType start = {{0, 0}};

  LSImageType::SizeType size = {{static_cast<LSImageType::SizeValueType>(m_nx), static_cast<LSImageType::SizeValueType>(m_ny)}};

  LSImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  
  mp_phi->SetRegions( region );
  mp_phi->Allocate();
  mp_phi->CopyInformation(mp_img);


  mp_phi->FillBuffer(arbitraryInitPhi);


  return;
}


/* ============================================================
   initializeSFLSFromMask    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::initializeSFLSFromMask()
{
  if (!mp_mask)
    {
      std::cerr<<"set mp_mask first.\n";
      raise(SIGABRT);
    }


  initializePhi();
  initializeLabel();

  for (long ix = 0; ix < m_nx; ++ix) 
    {
      for (long iy = 0; iy < m_ny; ++iy) 
        {
          typename ImageType::IndexType idx = {{ix, iy}};
          typename ImageType::IndexType idx1 = {{ix-1, iy}};
          typename ImageType::IndexType idx2 = {{ix+1, iy}};
          typename ImageType::IndexType idx3 = {{ix, iy-1}};
          typename ImageType::IndexType idx4 = {{ix, iy+1}};

          //mark the inside and outside of label and phi
          if( mp_mask->GetPixel(idx) == 0 )
            { 
              mp_label->SetPixel(idx, 3); 
              mp_phi->SetPixel(idx, 3); 
            }
          else
            { 
              mp_label->SetPixel(idx, -3); 
              mp_phi->SetPixel(idx, -3); 

              if ( (iy+1 < m_ny && mp_mask->GetPixel(idx4) == 0)	\
                   || (iy-1 >= 0 && mp_mask->GetPixel(idx3) == 0)	\
                   || (ix+1 < m_nx && mp_mask->GetPixel(idx2) == 0)   \
                   || (ix-1 >= 0 && mp_mask->GetPixel(idx1) == 0)	)
                {
                  m_lz.push_back(NodeType(ix, iy, 0)); // z-idx = 0 for 2D image

                  mp_label->SetPixel(idx, 0);
                  mp_phi->SetPixel(idx, 0.0);
                }
            }
        }
    }


  //scan Lz to create Ln1 and Lp1
  for (CSFLSLayer::const_iterator it = m_lz.begin(); it != m_lz.end(); ++it)
    {
      long ix = (*it)[0];
      long iy = (*it)[1];
    
      if(iy+1 < m_ny)
        {// up
          typename ImageType::IndexType idx = {{ix, iy+1}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back(NodeType(ix, iy+1, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix, iy+1, 0) );
            }
        }

      if(iy-1 >= 0)
        {// down
          typename ImageType::IndexType idx = {{ix, iy-1}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back( NodeType(ix, iy-1, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix, iy-1, 0) );
            }
        }


      if(ix+1 < m_nx)
        {// right
          typename ImageType::IndexType idx = {{ix+1, iy}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back( NodeType(ix+1, iy, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix+1, iy, 0) );
            }
        }

      if(ix-1 >= 0)
        {// left
          typename ImageType::IndexType idx = {{ix-1, iy}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back( NodeType(ix-1, iy, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix-1, iy, 0) );
            }
        }
    }



  //scan Ln1 to create Ln2
  for (CSFLSLayer::const_iterator it = m_ln1.begin(); it != m_ln1.end(); ++it)
    {
      long ix = (*it)[0];
      long iy = (*it)[1];

        
      typename ImageType::IndexType idx3 = {{ix, iy+1}};
      if(iy+1 < m_ny && mp_label->GetPixel(idx3) == -3 )
        {
          mp_label->SetPixel(idx3, -2); 
          mp_phi->SetPixel(idx3, -2); 
            
          m_ln2.push_back( NodeType(ix, iy+1, 0) );
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};
      if(iy-1 >= 0 && mp_label->GetPixel(idx4) == -3 )
        {
          mp_label->SetPixel(idx4, -2); 
          mp_phi->SetPixel(idx4, -2); 
            
          m_ln2.push_back( NodeType(ix, iy-1, 0) );
        }

      typename ImageType::IndexType idx1 = {{ix+1, iy}};
      if(ix+1 < m_nx && mp_label->GetPixel(idx1) == -3 )
        {
          mp_label->SetPixel(idx1, -2); 
          mp_phi->SetPixel(idx1, -2); 
            
          m_ln2.push_back( NodeType(ix+1, iy, 0) );
        }

      typename ImageType::IndexType idx2 = {{ix-1, iy}};
      if(ix-1 >= 0 && mp_label->GetPixel(idx2) == -3 )
        {
          mp_label->SetPixel(idx2, -2); 
          mp_phi->SetPixel(idx2, -2); 
            
          m_ln2.push_back( NodeType(ix-1, iy, 0) );
        }
    }

  //scan Lp1 to create Lp2
  for (CSFLSLayer::const_iterator it = m_lp1.begin(); it != m_lp1.end(); ++it)
    {
      long ix = (*it)[0];
      long iy = (*it)[1];

        
      typename ImageType::IndexType idx3 = {{ix, iy+1}};
      if(iy+1 < m_ny && mp_label->GetPixel(idx3) == 3 )
        {
          mp_label->SetPixel(idx3, 2); 
          mp_phi->SetPixel(idx3, 2); 
            
          m_lp2.push_back( NodeType(ix, iy+1, 0) );
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};
      if(iy-1 >= 0 && mp_label->GetPixel(idx4) == 3 )
        {
          mp_label->SetPixel(idx4, 2); 
          mp_phi->SetPixel(idx4, 2); 
            
          m_lp2.push_back( NodeType(ix, iy-1, 0) );
        }

      typename ImageType::IndexType idx1 = {{ix+1, iy}};
      if(ix+1 < m_nx && mp_label->GetPixel(idx1) == 3 )
        {
          mp_label->SetPixel(idx1, 2); 
          mp_phi->SetPixel(idx1, 2); 
            
          m_lp2.push_back( NodeType(ix+1, iy, 0) );
        }

      typename ImageType::IndexType idx2 = {{ix-1, iy}};
      if(ix-1 >= 0 && mp_label->GetPixel(idx2) == 3 )
        {
          mp_label->SetPixel(idx2, 2); 
          mp_phi->SetPixel(idx2, 2); 
            
          m_lp2.push_back( NodeType(ix-1, iy, 0) );
        }
    }
}


/* ============================================================
   getSFLSFromPhi    */
template< typename TPixel >
void
CSFLSSegmentor2D< TPixel >
::getSFLSFromPhi()
{
  initializePhi();
  initializeLabel();

  for (long ix = 0; ix < m_nx; ++ix) 
    {
      for (long iy = 0; iy < m_ny; ++iy) 
        {
          typename ImageType::IndexType idx = {{ix, iy}};
          typename ImageType::IndexType idx1 = {{ix-1, iy}};
          typename ImageType::IndexType idx2 = {{ix+1, iy}};
          typename ImageType::IndexType idx3 = {{ix, iy-1}};
          typename ImageType::IndexType idx4 = {{ix, iy+1}};

          //mark the inside and outside of label and phi
          if( mp_mask->GetPixel(idx) == 0 )
            { 
              mp_label->SetPixel(idx, 3); 
              mp_phi->SetPixel(idx, 3); 
            }
          else
            { 
              mp_label->SetPixel(idx, -3); 
              mp_phi->SetPixel(idx, -3); 

              if ( (iy+1 < m_ny && mp_mask->GetPixel(idx4) == 0)	\
                   || (iy-1 >= 0 && mp_mask->GetPixel(idx3) == 0)	\
                   || (ix+1 < m_nx && mp_mask->GetPixel(idx2) == 0) \
                   || (ix-1 >= 0 && mp_mask->GetPixel(idx1) == 0)	)
                {
                  m_lz.push_back(NodeType(ix, iy, 0)); // z-idx = 0 for 2D image

                  mp_label->SetPixel(idx, 0);
                  mp_phi->SetPixel(idx, 0.0);
                }
            }
        }
    }


  //scan Lz to create Ln1 and Lp1
  for (CSFLSLayer::const_iterator it = m_lz.begin(); it != m_lz.end(); ++it)
    {
      long ix = (*it)[0];
      long iy = (*it)[1];
    
      if(iy+1 < m_ny)
        {// up
          typename ImageType::IndexType idx = {{ix, iy+1}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back(NodeType(ix, iy+1, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix, iy+1, 0) );
            }
        }

      if(iy-1 >= 0)
        {// down
          typename ImageType::IndexType idx = {{ix, iy-1}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back( NodeType(ix, iy-1, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix, iy-1, 0) );
            }
        }


      if(ix+1 < m_nx)
        {// right
          typename ImageType::IndexType idx = {{ix+1, iy}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back( NodeType(ix+1, iy, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix+1, iy, 0) );
            }
        }

      if(ix-1 >= 0)
        {
          typename ImageType::IndexType idx = {{ix-1, iy}};

          if ( mp_label->GetPixel(idx) == 3 )
            {
              mp_label->SetPixel(idx, 1); 
              mp_phi->SetPixel(idx, 1); 
            
              m_lp1.push_back( NodeType(ix-1, iy, 0) );
            }
          else if ( mp_label->GetPixel(idx) == -3 )
            {
              mp_label->SetPixel(idx, -1); 
              mp_phi->SetPixel(idx, -1); 
            
              m_ln1.push_back( NodeType(ix-1, iy, 0) );
            }
        }
    }



  //scan Ln1 to create Ln2
  for (CSFLSLayer::const_iterator it = m_ln1.begin(); it != m_ln1.end(); ++it)
    {
      long ix = (*it)[0];
      long iy = (*it)[1];

        
      typename ImageType::IndexType idx3 = {{ix, iy+1}};
      if(iy+1 < m_ny && mp_label->GetPixel(idx3) == -3 )
        {
          mp_label->SetPixel(idx3, -2); 
          mp_phi->SetPixel(idx3, -2); 
            
          m_ln2.push_back( NodeType(ix, iy+1, 0) );
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};
      if(iy-1 >= 0 && mp_label->GetPixel(idx4) == -3 )
        {
          mp_label->SetPixel(idx4, -2); 
          mp_phi->SetPixel(idx4, -2); 
            
          m_ln2.push_back( NodeType(ix, iy-1, 0) );
        }

      typename ImageType::IndexType idx1 = {{ix+1, iy}};
      if(ix+1 < m_nx && mp_label->GetPixel(idx1) == -3 )
        {
          mp_label->SetPixel(idx1, -2); 
          mp_phi->SetPixel(idx1, -2); 
            
          m_ln2.push_back( NodeType(ix+1, iy, 0) );
        }

      typename ImageType::IndexType idx2 = {{ix-1, iy}};
      if(ix-1 >= 0 && mp_label->GetPixel(idx2) == -3 )
        {
          mp_label->SetPixel(idx2, -2); 
          mp_phi->SetPixel(idx2, -2); 
            
          m_ln2.push_back( NodeType(ix-1, iy, 0) );
        }
    }

  //scan Lp1 to create Lp2
  for (CSFLSLayer::const_iterator it = m_lp1.begin(); it != m_lp1.end(); ++it)
    {
      long ix = (*it)[0];
      long iy = (*it)[1];

        
      typename ImageType::IndexType idx3 = {{ix, iy+1}};        
      if(iy+1 < m_ny && mp_label->GetPixel(idx3) == 3 )
        {
          mp_label->SetPixel(idx3, 2); 
          mp_phi->SetPixel(idx3, 2); 
            
          m_lp2.push_back( NodeType(ix, iy+1, 0) );
        }

      typename ImageType::IndexType idx4 = {{ix, iy-1}};        
      if(iy-1 >= 0 && mp_label->GetPixel(idx4) == 3 )
        {
          mp_label->SetPixel(idx4, 2); 
          mp_phi->SetPixel(idx4, 2); 

          m_lp2.push_back( NodeType(ix, iy-1, 0) );
        }

      typename ImageType::IndexType idx1 = {{ix+1, iy}};        
      if(ix+1 < m_nx && mp_label->GetPixel(idx1) == 3 )
        {
          mp_label->SetPixel(idx1, 2); 
          mp_phi->SetPixel(idx1, 2); 
            
          m_lp2.push_back( NodeType(ix+1, iy, 0) );
        }

      typename ImageType::IndexType idx2 = {{ix-1, iy}};        
      if(ix-1 >= 0 && mp_label->GetPixel(idx2) == 3 )
        {
          mp_label->SetPixel(idx2, 2); 
          mp_phi->SetPixel(idx2, 2); 
            
          m_lp2.push_back( NodeType(ix-1, iy, 0) );
        }
    }
}


// /* ============================================================
//    doSegmenation    */
// template< typename TPixel >
// void
// CSFLSSegmentor2D< TPixel >
// ::doSegmenation()
// {
//   // gth818n::saveAsImage2< double >(mp_phi, "init0.nrrd");


//   /*============================================================
//    * From the initial mask, generate: 1. SFLS, 2. mp_label and
//    * 3. mp_phi.      
//    */
//   initializeSFLS();

//   //gth818n::saveAsImage2< double >(mp_phi, "initPhi.nrrd");

//   for (unsigned int it = 0; it < m_numIter; ++it)
//     {
//       /*--------------------------------------------------
//         Compute the force on the zero level set, NOT on the whole domain.
//         This is NOT implemented in this base class.    

//         This function will compute the m_force. m_force has the same
//         size as the m_ln, indicating the change at each pixel on the
//         zero level set.
//       */
//       computeForce(); 


//       normalizeForce();

//       //         // debug
//       //         for (std::list< double >::const_iterator itf = m_force.begin(); itf != m_force.end(); ++itf)
//       //           {
//       //             std::cout<<(*itf)<<", ";
//       //           }
//       //         std::cout<<std::endl<<it<<std::endl<<std::endl;



//       //         //debug//
//       //         labelsCoherentCheck1();

//       oneStepLevelSetEvolution();


//       //         //debug//
//       //         std::cout<<"-----------------------"<<it<<"---------------------------"<<std::endl;
//       //         std::cout<<"lz \t ln1 \t ln2 \t lp1 \t lp2 \n";
//       //         std::cout<<m_lz.size()<<"\t"<<m_ln1.size()<<"\t"<<m_ln2.size()<<"\t"<<m_lp1.size()<<"\t"<<m_lp2.size()<<std::endl;
//       //         std::cout<<"--------------------------------------------------"<<std::endl;


//       //         // debug
//       //         labelsCoherentCheck1();


//       //        gth818n::saveAsImage2< double >(mp_phi, "temp.nrrd");
//     }
// }

//   /* ============================================================
//      labelsCoherentCheck    */
//   template< typename TPixel >
//   void
//   CSFLSSegmentor2D< TPixel >
//   ::labelsCoherentCheck()
//   {
//     // check all in m_lz has the label 0
//     for (CSFLSLayer::const_iterator it = m_lz.begin(); it != m_lz.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_label->get(ix, iy) != 0)
//           {
//             std::cerr<<"mp_label->get(ix, iy) != 0\n";
//             raise(SIGABRT);
//           }

//         if (mp_phi->get(ix, iy) > 0.5 || mp_phi->get(ix, iy) < -0.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > 0.5 || mp_phi->get(ix, iy) < -0.5\n";
//             raise(SIGABRT);
//           }
//       }


//     // check all in m_lz has the label -1
//     for (CSFLSLayer::const_iterator it = m_ln1.begin(); it != m_ln1.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_label->get(ix, iy) != -1)
//           {
//             std::cerr<<"mp_label->get(ix, iy) != -1\n";
//             raise(SIGABRT);
//           }

//         if (mp_phi->get(ix, iy) > -0.5 || mp_phi->get(ix, iy) < -1.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > -0.5 || mp_phi->get(ix, iy) < -1.5\n";
//             raise(SIGABRT);
//           }
//       }

//     // check all in m_lz has the label -2
//     for (CSFLSLayer::const_iterator it = m_ln2.begin(); it != m_ln2.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_label->get(ix, iy) != -2)
//           {
//             std::cerr<<"mp_label->get(ix, iy) != -2\n";
//             raise(SIGABRT);
//           }

//         if (mp_phi->get(ix, iy) > -1.5 || mp_phi->get(ix, iy) < -2.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > -1.5 || mp_phi->get(ix, iy) < -2.5\n";
//             raise(SIGABRT);
//           }
//       }


//     // check all in m_lz has the label 2
//     for (CSFLSLayer::const_iterator it = m_lp2.begin(); it != m_lp2.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_label->get(ix, iy) != 2)
//           {
//             std::cerr<<"mp_label->get(ix, iy) != 2\n";
//             raise(SIGABRT);
//           }

//         if (mp_phi->get(ix, iy) > 2.5 || mp_phi->get(ix, iy) < 1.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > 2.5 || mp_phi->get(ix, iy) < 1.5\n";
//             raise(SIGABRT);
//           }
//       }


//     // check all in m_lz has the label 1
//     for (CSFLSLayer::const_iterator it = m_lp1.begin(); it != m_lp1.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_label->get(ix, iy) != 1)
//           {
//             std::cerr<<"mp_label->get(ix, iy) != 1\n";
//             raise(SIGABRT);
//           }

//         if (mp_phi->get(ix, iy) > 1.5 || mp_phi->get(ix, iy) < 0.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > 1.5 || mp_phi->get(ix, iy) < 0.5\n";
//             raise(SIGABRT);
//           }
//       }
//   }



/*============================================================

  computeKappa 

  Compute kappa at a point in the zero level set  */
template< typename TPixel >
double
CSFLSSegmentor2D< TPixel >
::computeKappa(long ix, long iy)
{
  //    double kappa;

  double dx = 0;
  double dy = 0;
  double dxx = 0;
  double dyy = 0;
  double dx2 = 0;
  double dy2 = 0;
  double dxy = 0;

  char xok = 0;
  char yok = 0;

  typename ImageType::IndexType idx = {{ix, iy}};
  typename ImageType::IndexType idx1 = {{ix-1, iy}};
  typename ImageType::IndexType idx2 = {{ix+1, iy}};
  typename ImageType::IndexType idx3 = {{ix, iy-1}};
  typename ImageType::IndexType idx4 = {{ix, iy+1}};



  if( ix+1 < m_nx && ix-1 >=0 )
    { 
      xok = 1;
    }

  if( iy+1 < m_ny && iy-1 >=0 ) 
    {
      yok = 1;
    }

  if (xok)
    {
      dx  = (mp_phi->GetPixel(idx2) - mp_phi->GetPixel(idx1) )/2.0;
      dxx = mp_phi->GetPixel(idx2) - 2.0*(mp_phi->GetPixel(idx)) + mp_phi->GetPixel(idx1);
      dx2 = dx*dx;
    }

  if (yok)
    {
      dy  = (mp_phi->GetPixel(idx4) - mp_phi->GetPixel(idx3) )/2.0;
      dyy = mp_phi->GetPixel(idx4) - 2*(mp_phi->GetPixel(idx)) + mp_phi->GetPixel(idx3);
      dy2 = dy*dy;
    }

  if(xok && yok)
    {// (ul+dr-ur-dl)/4
      typename ImageType::IndexType idx_1 = {{ix+1, iy+1}};
      typename ImageType::IndexType idx_2 = {{ix-1, iy-1}};
      typename ImageType::IndexType idx_3 = {{ix+1, iy-1}};
      typename ImageType::IndexType idx_4 = {{ix-1, iy+1}};

      dxy = 0.25*(mp_phi->GetPixel(idx_1) + mp_phi->GetPixel(idx_2) - mp_phi->GetPixel(idx_3) - mp_phi->GetPixel(idx_4));
    }
  
  return (dxx*dy2 + dyy*dx2 - 2*dx*dy*dxy)/(dx2 + dy2 + vnl_math::eps);
}


//   /* ============================================================
//      labelsCoherentCheck    */
//   template< typename TPixel >
//   void
//   CSFLSSegmentor2D< TPixel >
//   ::labelsCoherentCheck1()
//   {
//     // check all in m_lz has the label 0
//     for (CSFLSLayer::const_iterator it = m_lz.begin(); it != m_lz.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_phi->get(ix, iy) > 0.5 || mp_phi->get(ix, iy) < -0.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > 0.5 || mp_phi->get(ix, iy) < -0.5\n";
//             raise(SIGABRT);
//           }
//       }


//     // check all in m_lz has the label -1
//     for (CSFLSLayer::const_iterator it = m_ln1.begin(); it != m_ln1.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_phi->get(ix, iy) > -0.5 || mp_phi->get(ix, iy) < -1.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > -0.5 || mp_phi->get(ix, iy) < -1.5\n";
//             raise(SIGABRT);
//           }
//       }

//     // check all in m_lz has the label -2
//     for (CSFLSLayer::const_iterator it = m_ln2.begin(); it != m_ln2.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_phi->get(ix, iy) > -1.5 || mp_phi->get(ix, iy) < -2.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > -1.5 || mp_phi->get(ix, iy) < -2.5\n";
//             raise(SIGABRT);
//           }
//       }


//     // check all in m_lz has the label 2
//     for (CSFLSLayer::const_iterator it = m_lp2.begin(); it != m_lp2.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_phi->get(ix, iy) > 2.5 || mp_phi->get(ix, iy) < 1.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > 2.5 || mp_phi->get(ix, iy) < 1.5\n";
//             raise(SIGABRT);
//           }
//       }


//     // check all in m_lz has the label 1
//     for (CSFLSLayer::const_iterator it = m_lp1.begin(); it != m_lp1.end(); ++it)
//       {
//         long ix = (*it)[0];
//         long iy = (*it)[1];

//         if (mp_phi->get(ix, iy) > 1.5 || mp_phi->get(ix, iy) < 0.5)
//           {
//             std::cerr<<"mp_phi->get(ix, iy) = "<<mp_phi->get(ix, iy)<<std::endl;
//             std::cerr<<"mp_phi->get(ix, iy) > 1.5 || mp_phi->get(ix, iy) < 0.5\n";
//             raise(SIGABRT);
//           }
//       }
//   }


#endif
