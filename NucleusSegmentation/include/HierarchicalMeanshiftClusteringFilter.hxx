#ifndef HierarchicalMeanshiftClusteringFilter_hxx_
#define HierarchicalMeanshiftClusteringFilter_hxx_

#include "itkNumericTraits.h"
#include "vnl/vnl_random.h"

#include "HierarchicalMeanshiftClusteringFilter.h"
#include "MeanshiftClusteringFilter.h"

namespace gth818n
{

  template< typename TCoordRep, unsigned int NPointDimension >
  HierarchicalMeanshiftClusteringFilter<TCoordRep, NPointDimension>::HierarchicalMeanshiftClusteringFilter()
  {
    m_inputPointSet = 0;

    m_subsampleRatio = 0.1;

    m_radius = 3.0;
    m_allDone = false;
  }


  template< typename TCoordRep, unsigned int NPointDimension >
  void HierarchicalMeanshiftClusteringFilter<TCoordRep, NPointDimension>::update()
  {
    /// Sub-sample input point set
    typename VectorSampleType::Pointer subsample = VectorSampleType::New();

    vnl_random rg;
    for (unsigned int i = 0 ; i < m_inputPointSet->Size() ; ++i )
      {
        if (rg.drand64() <= m_subsampleRatio)
          {
            subsample->PushBack(m_inputPointSet->GetMeasurementVector(i) );
          }
      }


    /// Run mean shift on sub-sample
    typedef gth818n::MeanshiftClusteringFilter<TCoordRep, NPointDimension> MeanshiftClusteringFilterType;
    MeanshiftClusteringFilterType ms;
    ms.setInputPointSet(subsample);
    ms.setRadius(m_radius);
    ms.update();

    /// Get label back to input point set by closest point criteria
    std::vector<long> sublabel = ms.getLabelOfPoints();
    m_centers = ms.getCenters();

    typedef typename itk::Statistics::KdTreeGenerator< VectorSampleType > TreeGeneratorType;
    typedef typename TreeGeneratorType::KdTreeType TreeType;

    typename TreeGeneratorType::Pointer subTreeGenerator = TreeGeneratorType::New();
    subTreeGenerator->SetSample( subsample );
    subTreeGenerator->SetBucketSize( 16 );
    subTreeGenerator->Update();
    typename TreeType::Pointer subtree = subTreeGenerator->GetOutput();

    unsigned int numberOfNeighbors = 1;
    typename TreeType::InstanceIdentifierVectorType neighbors;

    m_labelOfPoints.resize(m_inputPointSet->Size());

    for (unsigned int i = 0 ; i < m_inputPointSet->Size() ; ++i )
      {
        const VectorType& queryPoint = m_inputPointSet->GetMeasurementVector(i);
        subtree->Search( queryPoint, numberOfNeighbors, neighbors ) ;
        m_labelOfPoints[i] = static_cast<long>(sublabel[neighbors[0]]);
      }






    m_allDone = true;

    return;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  void HierarchicalMeanshiftClusteringFilter<TCoordRep, NPointDimension>::setRadius(RealType rad)
  {
    if (rad <= 0.0)
      {
        std::cerr<<"Error: rad should > 0, but got "<<rad<<std::endl;
        abort();
      }

    m_radius = rad;

    return;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  void HierarchicalMeanshiftClusteringFilter<TCoordRep, NPointDimension>::setSubsampleRatio(RealType ratio)
  {
    if (ratio <= 0.0 || ratio > 1.0)
      {
        std::cerr<<"Error: ratio should in (0, 1], but got "<<ratio<<std::endl;
        abort();
      }

    m_subsampleRatio = ratio;

    return;
  }


  template< typename TCoordRep, unsigned int NPointDimension >
  typename HierarchicalMeanshiftClusteringFilter<TCoordRep, NPointDimension>::VectorSampleType::Pointer
  HierarchicalMeanshiftClusteringFilter<TCoordRep, NPointDimension>::getCenters()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_centers;
  }


  template< typename TCoordRep, unsigned int NPointDimension >
  std::vector<long>
  HierarchicalMeanshiftClusteringFilter<TCoordRep, NPointDimension>::getLabelOfPoints()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_labelOfPoints;
  }


}// namespace gth818n


#endif
