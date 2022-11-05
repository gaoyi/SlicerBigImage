#ifndef MeanshiftClusteringFilter_h_
#define MeanshiftClusteringFilter_h_


// itk
#include "itkVector.h"
#include "itkListSample.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistanceMetric.h"


namespace gth818n
{
  template< typename TCoordRep, unsigned int NPointDimension >
  class MeanshiftClusteringFilter
  {
  public:
    //------------------------------------------------------------------------------
    /// typedef

    /** Standard class typedefs. */
    typedef MeanshiftClusteringFilter Self;

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
    MeanshiftClusteringFilter();
    ~MeanshiftClusteringFilter() {}
    /// ctor, end
    //------------------------------------------------------------------------------


    //------------------------------------------------------------------------------
    /// public fn
    /** Get the dimension (size) of the point. */
    static unsigned int GetPointDimension() { return m_PointDimension; }

    void setInputPointSet(typename VectorSampleType::Pointer inputPointSet) {m_inputPointSet = inputPointSet;}
    void setRadius(RealType rad);
    void setEpoch(int epoch);
    void update();

    typename VectorSampleType::Pointer getCenters();
    std::vector<long> getLabelOfPoints();
    /// public fn, end
    //------------------------------------------------------------------------------


  private:
    typedef typename itk::Statistics::KdTreeGenerator< VectorSampleType > TreeGeneratorType;
    typedef typename TreeGeneratorType::KdTreeType TreeType;
    typedef typename TreeType::NearestNeighbors NeighborsType;
    typedef typename TreeType::KdTreeNodeType NodeType;

    typename TreeGeneratorType::Pointer m_treeGenerator;
    typename TreeType::Pointer m_tree;

    static const unsigned int m_PointDimension = NPointDimension;

    typename VectorSampleType::Pointer m_inputPointSet;
    typename VectorSampleType::Pointer m_seedPoints;
    typename VectorSampleType::Pointer m_centers;

    VectorType m_inputPointSetRange;


    RealType m_radius;
    long m_numberOfMSIteration;
    long m_numberOfModes;

    int m_epoch; ///< after mean shift on input data, will recursively
                 ///run mean shift on the obtained data. This give the
                 ///number of recusion. If 0, just run MS for once:
                 ///there may be some close but not too close mode
                 ///centers.

    std::vector<long> m_labelOfPoints;

    bool m_allDone;
    /// private data, end
    //------------------------------------------------------------------------------


    //------------------------------------------------------------------------------
    /// private fn
    void _computeInputPointRange();
    void _constructKdTree();
    void _meanshiftIteration();
    void _constructSeedPoints();
    void _findUniqueCenters();
    void _findLabelOfPoints();
    typename VectorSampleType::Pointer _getSeedPoints();
    /// private fn
    //------------------------------------------------------------------------------
  };

}// namespace gth818n

#include "MeanshiftClusteringFilter.hxx"

#endif // MeanshiftClusteringFilter_h_
