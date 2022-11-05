#ifndef MeanshiftClusteringFilter_hxx_
#define MeanshiftClusteringFilter_hxx_

#include "itkNumericTraits.h"


#include "MeanshiftClusteringFilter.h"


namespace gth818n
{

  template< typename TCoordRep, unsigned int NPointDimension >
  MeanshiftClusteringFilter<TCoordRep, NPointDimension>::MeanshiftClusteringFilter()
  {
    m_epoch = 2;
    m_inputPointSet = 0;
    m_seedPoints = 0;

    m_numberOfMSIteration = 100;

    m_radius = 3.0;
    m_allDone = false;
  }


  template< typename TCoordRep, unsigned int NPointDimension >
  void MeanshiftClusteringFilter<TCoordRep, NPointDimension>::update()
  {
    // //dbg
    // std::cout<<"_constructKdTree..."<<std::flush;
    // //dbg, end
    _constructKdTree();
    // //dbg
    // std::cout<<"done"<<std::endl<<std::flush;
    // //dbg, end

    if (!m_seedPoints)
      {
        _constructSeedPoints();
      }

    _meanshiftIteration();

    _findUniqueCenters();

    _findLabelOfPoints();

    m_allDone = true;

    return;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  void MeanshiftClusteringFilter<TCoordRep, NPointDimension>::setRadius(RealType rad)
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
  void MeanshiftClusteringFilter<TCoordRep, NPointDimension>::setEpoch(int epoch)
  {
    if (epoch < 0)
      {
        std::cerr<<"Error: epoch should >= 0, but got "<<epoch<<std::endl;
        abort();
      }

    m_epoch = epoch;

    return;
  }




  template< typename TCoordRep, unsigned int NPointDimension >
  void MeanshiftClusteringFilter<TCoordRep, NPointDimension>::_constructKdTree()
  {
    m_treeGenerator = TreeGeneratorType::New();
    m_treeGenerator->SetSample( m_inputPointSet );
    m_treeGenerator->SetBucketSize( 16 );
    m_treeGenerator->Update();

    m_tree = m_treeGenerator->GetOutput();

    // VectorType queryPoint;
    // queryPoint[0] = 10.0;
    // queryPoint[1] = 7.0;

    // // K-Neighbor search
    // std::cout << "K-Neighbor search:" << std::endl;
    // unsigned int numberOfNeighbors = 3;
    // typename TreeType::InstanceIdentifierVectorType neighbors;
    // m_tree->Search( queryPoint, numberOfNeighbors, neighbors ) ;

    // for ( unsigned int i = 0 ; i < neighbors.size() ; ++i )
    //   {
    //     std::cout << m_tree->GetMeasurementVector( neighbors[i] ) << std::endl;
    //   }

    // // Radius search
    // std::cout << "Radius search:" << std::endl;
    // double radius = 4.0;
    // m_tree->Search( queryPoint, radius, neighbors ) ;
    // std::cout << "There are " << neighbors.size() << " neighbors." << std::endl;
    // for ( unsigned int i = 0 ; i < neighbors.size() ; ++i )
    //   {
    //     std::cout << m_tree->GetMeasurementVector( neighbors[i] ) << std::endl;
    //   }

    return;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  void MeanshiftClusteringFilter<TCoordRep, NPointDimension>::_computeInputPointRange()
  {
    VectorType inputPointSetMin;
    inputPointSetMin.Fill(itk::NumericTraits< RealType >::max());

    VectorType inputPointSetMax;
    inputPointSetMax.Fill(itk::NumericTraits< RealType >::min());

    for (long itp = 0; itp < m_inputPointSet->Size(); ++itp)
      {
        VectorType thisPoint = m_inputPointSet->GetMeasurementVector(itp);

        for (unsigned int idim = 0; idim < NPointDimension; ++idim)
          {
            inputPointSetMin[idim] = inputPointSetMin[idim]<thisPoint[idim]?inputPointSetMin[idim]:thisPoint[idim];
            inputPointSetMax[idim] = inputPointSetMax[idim]<thisPoint[idim]?inputPointSetMax[idim]:thisPoint[idim];
          }
      }

    m_inputPointSetRange = inputPointSetMax - inputPointSetMin;

    return;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  void MeanshiftClusteringFilter<TCoordRep, NPointDimension>::_meanshiftIteration()
  {
    VectorType queryPoint;
    typename TreeType::InstanceIdentifierVectorType neighbors;

    for (long it = 0; it < m_numberOfMSIteration; ++it)
      {
        for (long itp = 0; itp < m_seedPoints->Size(); ++itp)
          {
            queryPoint = m_seedPoints->GetMeasurementVector(itp);
            m_tree->Search( queryPoint, m_radius, neighbors ) ;

            VectorType newPosition;
            newPosition.Fill(0);

            for ( unsigned int i = 0 ; i < neighbors.size() ; ++i )
              {
                newPosition += m_tree->GetMeasurementVector( neighbors[i] );
                //std::cout << m_tree->GetMeasurementVector( neighbors[i] ) << std::endl;
              }

            newPosition /= static_cast<RealType>(neighbors.size());

            m_seedPoints->SetMeasurementVector(itp, newPosition);

            /// If relative increamental is small enough, break
            VectorType del = queryPoint - newPosition;
            for (unsigned int idim = 0; idim < NPointDimension; ++idim)
              {
                del[idim] /= m_inputPointSetRange[idim];
              }

            if (del.GetNorm() < 1e-1)
              {
                break;
              }
          }
      }

    if (m_epoch)
      {
        MeanshiftClusteringFilter<TCoordRep, NPointDimension> ms;
        ms.setInputPointSet(m_seedPoints);
        ms.setEpoch(--m_epoch);
        ms.update();
        m_seedPoints = ms._getSeedPoints();
      }

    return;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  typename MeanshiftClusteringFilter<TCoordRep, NPointDimension>::VectorSampleType::Pointer
  MeanshiftClusteringFilter<TCoordRep, NPointDimension>::_getSeedPoints()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_seedPoints;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  typename MeanshiftClusteringFilter<TCoordRep, NPointDimension>::VectorSampleType::Pointer
  MeanshiftClusteringFilter<TCoordRep, NPointDimension>::getCenters()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_centers;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  void
  MeanshiftClusteringFilter<TCoordRep, NPointDimension>::_findUniqueCenters()
  {
    m_centers = VectorSampleType::New();
    m_centers->PushBack( m_seedPoints->GetMeasurementVector(0) );

    for (unsigned int i = 1 ; i < m_seedPoints->Size() ; ++i )
      {
        const VectorType& newCenterCandidate = m_seedPoints->GetMeasurementVector(i);
        bool IAmNotCloseToAnyExistingCenter = true;

        for (unsigned int ii = 0 ; ii < m_centers->Size() ; ++ii )
          {
            VectorType distVector = newCenterCandidate - m_centers->GetMeasurementVector(ii);
            if (distVector.GetNorm() < m_radius )
              {
                IAmNotCloseToAnyExistingCenter = false;
                break;
              }
          }

        if (IAmNotCloseToAnyExistingCenter)
          {
            m_centers->PushBack( newCenterCandidate );
          }
      }

    m_numberOfModes = m_centers->Size();

    return;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  std::vector<long>
  MeanshiftClusteringFilter<TCoordRep, NPointDimension>::getLabelOfPoints()
  {
    if (!m_allDone)
      {
        std::cerr<<"Error: not done.\n";
      }

    return m_labelOfPoints;
  }

  template< typename TCoordRep, unsigned int NPointDimension >
  void
  MeanshiftClusteringFilter<TCoordRep, NPointDimension>::_findLabelOfPoints()
  {
    typename TreeGeneratorType::Pointer newTreeGen = TreeGeneratorType::New();
    newTreeGen->SetSample( m_centers );
    newTreeGen->SetBucketSize( 16 );
    newTreeGen->Update();
    typename TreeType::Pointer newTree = newTreeGen->GetOutput();

    m_labelOfPoints.resize(m_inputPointSet->Size());

    // K-Neighbor search
    //std::cout << "K-Neighbor search:" << std::endl;
    unsigned int numberOfNeighbors = 1;
    typename TreeType::InstanceIdentifierVectorType neighbors;

    for (unsigned int i = 0 ; i < m_seedPoints->Size() ; ++i )
      {
        const VectorType& queryPoint = m_seedPoints->GetMeasurementVector(i);
        newTree->Search( queryPoint, numberOfNeighbors, neighbors ) ;
        m_labelOfPoints[i] = static_cast<long>(neighbors[0]);
      }

    return;
  }


  template< typename TCoordRep, unsigned int NPointDimension >
  void MeanshiftClusteringFilter<TCoordRep, NPointDimension>::_constructSeedPoints()
  {
    /// Duplicate input points as seed points
    m_seedPoints = VectorSampleType::New();
    m_seedPoints->Resize(m_inputPointSet->Size() );

    for (unsigned int i = 0 ; i < m_inputPointSet->Size() ; ++i )
      {
        m_seedPoints->SetMeasurementVector( i, m_inputPointSet->GetMeasurementVector(i) );
      }

    return;
  }


}// namespace gth818n


#endif
