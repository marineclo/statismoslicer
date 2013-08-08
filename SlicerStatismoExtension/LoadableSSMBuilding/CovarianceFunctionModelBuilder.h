#ifndef __COVARIANCE_MODEL_BUILDER_H
#define __COVARIANCE_MODEL_BUILDER_H
#include "itkVectorImageRepresenter.h"
#include "itkMeshRepresenter.h"
#include "statismo/Config.h"
#include "statismo/ModelInfo.h"
#include "statismo/ModelBuilder.h"
#include "statismo/DataManager.h"
#include "statismo/StatisticalModel.h"
#include "statismo/CommonTypes.h"
#include "Kernel.h"
#include <vector>
#include <cmath>
#include <memory>
#include <future>
#include <thread>
#include "RandSVD.h"
#include <Eigen/QR>
#include <Eigen/SVD>
namespace statismo {


  struct ResEigenfunctionPointComputations { 
  
  ResEigenfunctionPointComputations(unsigned _lowerInd, unsigned _upperInd, const MatrixType& _resMat) : 
    lowerInd(_lowerInd), upperInd(_upperInd), resMatrix(_resMat) {}

    MatrixType resMatrix;
    unsigned lowerInd;
    unsigned upperInd;
  };


  template<typename Representer>
    class CovarianceFunctionModelBuilder: public ModelBuilder<Representer> {

  public:

    typedef ModelBuilder<Representer> Superclass;
    typedef typename Superclass::StatisticalModelType StatisticalModelType;

    typedef typename Representer::PointType PointType;

    typedef statismo::Domain<typename Representer::PointType> DomainType;
    typedef typename DomainType::DomainPointsListType DomainPointsListType;

    /**
     * Factory method to create a new PCAModelBuilder
     */
    static CovarianceFunctionModelBuilder* Create(
						  const Representer* representer) {
      return new CovarianceFunctionModelBuilder(representer);
    }

    /**
     * Destroy the object.
     * The same effect can be achieved by deleting the object in the usual
     * way using the c++ delete keyword.
     */
    void Delete() {
      delete this;
    }

    /**
     * The desctructor
     */
    virtual ~CovarianceFunctionModelBuilder() {
    }



    ResEigenfunctionPointComputations computeEigenfunctionRow(const Kernel* kernel,  unsigned numComponents, unsigned numDomainPoints, const  std::vector<VectorType>& xs, const  std::vector<VectorType> & domainPts, const MatrixType& M, unsigned lowerInd, unsigned upperInd) const { 

	unsigned n = numDomainPoints;
	unsigned m = xs.size();

	if (upperInd >= domainPts.size()) { 
	  std::cout << "somethings terribly wrong, exiting" << std::endl;
	  exit(-1);
	}
      MatrixType resMat = MatrixType::Zero(upperInd - lowerInd, numComponents);


      for (unsigned i = lowerInd; i < upperInd; i++) { 
	
	VectorType Ei = VectorType::Zero(m);

	for (unsigned j = 0; j < m; j++) {
	  Ei(j) = (*kernel)(domainPts[i], xs[j]);
	}
	resMat.row(i - lowerInd) =  M  * Ei;

	if ((i - lowerInd) % 100000 == 0) { 
	  std::cout << "processed 100'000 items out of " << upperInd - lowerInd << " "  << "(id = " << std::this_thread::get_id() << ")" << std::endl;
	}
      }
	return ResEigenfunctionPointComputations(lowerInd, upperInd, resMat);
    }


    /**
     * Build a new model using a zero mean GAussian Process with covariance function given by the kernel.
     * \param numPointsForNystrom  The number of points used in the Nystrom Approximation
     * \param numComponents The number of principal components to be computed
     *
     * \return a new statistical model representing the given gaussian process
     */
    StatisticalModelType* BuildNewModel(const Kernel* kernel,
					unsigned numPointsForNystrom, unsigned numComponents) const {

      DomainType domain = m_representer->GetDomain();
      unsigned n = domain.GetNumberOfPoints();
      unsigned p = n;
      unsigned outputDimension = Representer::GetDimensions();
	
	//std::cout<< "n = "<<n<< " p = "<<p<<" outputDimension = "<<outputDimension<<std::endl;
 
      // convert all the points in the domain to a vector representation
      std::vector < VectorType > domainPoints;
      for (typename DomainPointsListType::const_iterator it =
	     domain.GetDomainPoints().begin();
	   it != domain.GetDomainPoints().end(); ++it) {
	// TODO this should be generic, and not assume that the point type is an ITK point
	domainPoints.push_back(convertITKPointToVectorType(*it));
      }
      assert(domainPoints.size() == n); //check if all points have been converted (if not return an error)

      // select a subset of the points for computing the nystrom approximation
      std::vector<VectorType> xs;
      int step = domainPoints.size() / numPointsForNystrom;

	//std::cout<<"domainPoints.size()= "<<domainPoints.size()<<" step= "<<step<<std::endl;
    
	for (typename std::vector<VectorType>::const_iterator it =
	     domainPoints.begin(); it < domainPoints.end(); it += step) {
	xs.push_back(*it);

      }

      int m = xs.size();

	//std::cout<<"m= "<<m<<std::endl;

      // compute a eigenvalue decomposition of the kernel matrix, evaluated at the points used for the
      // nystrom approximation
      MatrixType U = MatrixType::Zero(m, numComponents); // eigenvectors (principal components)
      VectorType D = RowVectorType::Zero(numComponents); // eigenvalues (variance)
      computeKernelMatrixDecomposition(kernel, xs, numComponents, U, D);

      // Compute an m times n extension matrix E and multiply it by u, to extend 
      // the estimation of the eigenvector to the full function (i.e. all image domain points). 
      // To avoid having to compute the huge matrix E, we implement the matrix multiplication 
      // by n vector x matrix multiplications.

      std::vector < std::future<ResEigenfunctionPointComputations> > resvec;


      MatrixType M =  sqrt(m / float(n)) * (U.leftCols(numComponents)* D.topRows(numComponents).asDiagonal().inverse()).transpose();

      MatrixType pcaBasis = MatrixType::Zero(n, numComponents);
      unsigned numChunks = 24;
      for (unsigned i = 0; i < numChunks ; i++) { 
	int chunkSize = n / numChunks;
	int lowerInd = i * chunkSize;
	int upperInd = std::min(n - 1, (i + 1) * chunkSize);

	resvec.push_back(
			 std::async(std::launch::async,
				    &CovarianceFunctionModelBuilder<Representer>::computeEigenfunctionRow,
				    this, kernel, numComponents, n, xs, domainPoints, M, lowerInd, upperInd));

      }

      for (unsigned i = 0; i < resvec.size(); i++) { 
	ResEigenfunctionPointComputations res = resvec[i].get();
	pcaBasis.block(res.lowerInd, 0, res.upperInd - res.lowerInd, pcaBasis.cols()) = res.resMatrix;
	std::cout << "\rcomputed chunk " << i << " of " << numChunks << "chunks" << std::flush;
      }
      std::cout << std::endl;
    
      std::cout << " now computing the pca Basis Matrix for the ND case" << std::endl;
      
      MatrixType pcaVariance = n / float(m) * D.topRows(numComponents);

      std::cout << "pcaVariance " << pcaVariance << std::endl;

      // all the computations were now for the 1d case. now we have several output directions.
      // We just assumed that they are all independent.
      MatrixType pcaBasisND = MatrixType::Zero(n * outputDimension,
					       numComponents * outputDimension);
      VectorType pcaVarianceND = VectorType::Zero(
						  numComponents * outputDimension);

      for (unsigned k = 0; k < numComponents; k++) {
	for (unsigned i = 0; i < n; i++) {
	  for (unsigned d = 0; d < outputDimension; d++) {
	    pcaBasisND(i * outputDimension + d, k * outputDimension + d) =
	      pcaBasis(i, k);
	  }
	}

	for (unsigned d = 0; d < outputDimension; d++) {
	  pcaVarianceND(k * outputDimension + d) = pcaVariance(k);
	}

      }


      
      
      // depending on the model we build (shape or deformation model) the model mean might not be zero.
      // the following call accounts for it.
      RowVectorType mu = getMean();
      StatisticalModelType* model = StatisticalModelType::Create(
								 m_representer, mu, pcaBasisND, pcaVarianceND, 0);

      // there are no scores in our case. We just store an empty list to satisfy the compiler
      MatrixType scores;

      // finally add meta data to the model info
      typename ModelInfo::BuilderInfoList bi;
      bi.push_back(
		   ModelInfo::KeyValuePair("BuilderName ",
					   "CovarianceFunctionModelBuilder"));
      bi.push_back(
		   ModelInfo::KeyValuePair("NoiseVariance ",
					   Utils::toString(0)));

      typename ModelInfo::DataInfoList dataInfo;

      ModelInfo info(scores, dataInfo, bi);
      model->SetModelInfo(info);

      return model;
    }

  private:
    RowVectorType getMean() const;


    /**
     * Compute the kernel matrix for all points given in xs and
     * return a matrix U with the first numComponents eigenvectors and a vector D with
     * the corresponding eigenvalues of this kernel matrix
     */
    void computeKernelMatrixDecomposition(const Kernel* kernel,
					  const std::vector<VectorType>& xs, unsigned numComponents,
					  MatrixType& U, VectorType& D) const {
      double n = xs.size();
      MatrixTypeDoublePrecision K = MatrixTypeDoublePrecision::Zero(n, n);
      for (unsigned i = 0; i < n; ++i) {
	for (unsigned j = i; j < n; ++j) {
	  K(i, j) = (*kernel)(xs[i], xs[j]);
	  K(j, i) = K(i, j);
	}
      }


      typedef RandSVD<double> SVDType;
      SVDType svd(K, numComponents);
      U = svd.matrixU().cast<ScalarType>();
      D = svd.singularValues().cast<ScalarType>();
      std::cout << "DL: " << D.topRows(numComponents) << std::endl;

      std::cout << "nystroem computed " << std::endl;
    }

    /**
     * little helper routing to convert itk point types to a statismo vector type.
     */
    statismo::VectorType convertITKPointToVectorType(
						     const typename Representer::PointType& itkPt) const {
      unsigned dim = Representer::PointType::Dimension;
      statismo::VectorType pt = statismo::VectorType::Zero(dim, 1);

      for (unsigned i = 0; i < dim; ++i) {
	pt(i) = itkPt[i];
      }
      return pt;
    }


    /**
     * constructor - only used internally
     */
  CovarianceFunctionModelBuilder(const Representer* representer) :
    m_representer(representer) {
    }


    // purposely not implemented
    CovarianceFunctionModelBuilder(const CovarianceFunctionModelBuilder& orig);
    CovarianceFunctionModelBuilder& operator=(
					      const CovarianceFunctionModelBuilder& rhs);

    const Representer* m_representer;

  };
  template <class RepresenterType>
    RowVectorType CovarianceFunctionModelBuilder<RepresenterType>::getMean() const {throw std::runtime_error("specialize get mean for different representers"); }


  // deformation fields are modeled as zero mean guassian processes
  template <>
    RowVectorType CovarianceFunctionModelBuilder<itk::VectorImageRepresenter<float, 3, 3> >::getMean() const {
    return RowVectorType::Zero(this->m_representer->GetDomain().GetNumberOfPoints() * 3);
  }

  template <>
    RowVectorType CovarianceFunctionModelBuilder<itk::VectorImageRepresenter<float, 2, 2> >::getMean() const {
    return RowVectorType::Zero(this->m_representer->GetDomain().GetNumberOfPoints() * 2);
  }

  // for mesh models, we choose the reference as the mean.
  template <>
    RowVectorType CovarianceFunctionModelBuilder<itk::MeshRepresenter<float, 2> >::getMean() const { return this->m_representer->SampleToSampleVector(this->m_representer->GetReference()); }

  template <>
    RowVectorType CovarianceFunctionModelBuilder<itk::MeshRepresenter<float, 3> >::getMean() const { return this->m_representer->SampleToSampleVector(this->m_representer->GetReference()); }


} // namespace statismo

#endif // __COVARIANCE_MODEL_BUILDER_H
