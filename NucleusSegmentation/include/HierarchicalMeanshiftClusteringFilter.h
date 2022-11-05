#ifndef HierarchicalMeanshiftClusteringFilter_h_
#define HierarchicalMeanshiftClusteringFilter_h_


// itk
#include "itkVector.h"
#include "itkListSample.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistanceMetric.h"


namespace gth818n
{
  template< typename TCoordRep, unsigned int NPointDimension >
  class HierarchicalMeanshiftClusteringFilter
  {
  public:
    //------------------------------------------------------------------------------
    /// typedef

    /** Standard class typedefs. */
    typedef HierarchicalMeanshiftClusteringFilter Self;

    /** ValueType can be used to declare a variable that is the same type
     * as a data element held in an Point.   */
    typedef TCoordRep ValueType;
    typedef TCoordRep CoordRepType;

    typedef typename itk::NumericTraits< ValueType >::RealType RealType;

    typedef itk::Vector< TCoordRep, NPointDimension > VectorType;

    typedef typename itk::Statistics::ListSample< VectorType > VectorSampleType;
    /// typedef, end
    //------------------------------------------------------------------------------


    //------------------------------------------------------------------------------
    /// ctor
    HierarchicalMeanshiftClusteringFilter();
    ~HierarchicalMeanshiftClusteringFilter() {}
    /// ctor, end
    //------------------------------------------------------------------------------


    //------------------------------------------------------------------------------
    /// public fn
    /** Get the dimension (size) of the point. */
    static unsigned int GetPointDimension() { return m_PointDimension; }

    void setInputPointSet(typename VectorSampleType::Pointer inputPointSet) {m_inputPointSet = inputPointSet;}
    void setRadius(RealType rad);
    void setSubsampleRatio(RealType ratio);
    void update();

    typename VectorSampleType::Pointer getCenters();
    std::vector<long> getLabelOfPoints();
    /// public fn, end
    //------------------------------------------------------------------------------


  private:
    static const unsigned int m_PointDimension = NPointDimension;

    typename VectorSampleType::Pointer m_inputPointSet;
    typename VectorSampleType::Pointer m_centers;

    RealType m_radius;

    RealType m_subsampleRatio;

    std::vector<long> m_labelOfPoints;

    bool m_allDone;
    /// private data, end
    //------------------------------------------------------------------------------


    //------------------------------------------------------------------------------
    /// private fn
    void _subsample();
    void _runMeanshiftOnSubsample();
    void _getLabelFromResultsOfSubsample();
    /// private fn
    //------------------------------------------------------------------------------
  };

}// namespace gth818n

#include "HierarchicalMeanshiftClusteringFilter.hxx"

#endif // HierarchicalMeanshiftClusteringFilter_h_
