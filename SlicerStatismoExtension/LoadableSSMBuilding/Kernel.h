/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianKernel.h,v $
  Language:  C++
  Date:      $Date: 2009-11-13 13:22:00 $
  Version:   $Revision: 16784 $
  Author:	 $Author: Marcel L�thi
  
=========================================================================*/
#ifndef KERNEL_H
#define KERNEL_H

#include <vector>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_diag_matrix.h"
#include "vnl/vnl_math.h"
#include <map>
#include <cmath>
#include "statismo/CommonTypes.h"
#include <exception>

#include <unistd.h>


namespace statismo {

/**
 * \brief Creates StatisticalModel using arbitrary covariance matrix
 *
 */

class Kernel {
public:



	virtual ~Kernel() {
	}
	virtual double operator()(const VectorType& x,	const VectorType& y) const = 0;
};



class GaussianKernel: public Kernel {

public:

	GaussianKernel(double sigma2) :
			m_sigma2(sigma2) {
	}
	virtual ~GaussianKernel() {
	}

	inline double operator()(const VectorType& x, const VectorType& y) const {
		VectorType r = x - y;
		return exp(-r.dot(r) / m_sigma2);
	}

private:
	// private, to prevent use
	GaussianKernel(const GaussianKernel& orig);
	GaussianKernel& operator=(const GaussianKernel& rhs);

	double m_sigma2;
};


class MultiscaleGaussianKernel: public Kernel {

public:

	MultiscaleGaussianKernel(double sigma2Lower, double sigma2Upper, unsigned numLevels)
	 {
		double step = (sigma2Upper - sigma2Lower) / numLevels;
		for (unsigned i = 0;i < numLevels; i++) {
			m_sigma2s.push_back(sigma2Lower + i * step);
		}
		std::cout<<"MultiscaleGaussianKernel"<<std::endl;
	}
	virtual ~MultiscaleGaussianKernel() {
	}

	inline double operator()(const VectorType& x, const VectorType& y) const {
		VectorType r = x - y;
		double sum = 0;
		for(std::list<double>::const_iterator it = m_sigma2s.begin(); it != m_sigma2s.end(); ++it) {
			sum += exp(-r.dot(r) / *it);
		}
		//std::cout<<" test "<<std::endl;
		//sleep(10);
		return sum;
	}

private:
	// private, to prevent use
	MultiscaleGaussianKernel(const MultiscaleGaussianKernel& orig);
	MultiscaleGaussianKernel& operator=(const MultiscaleGaussianKernel& rhs);

	std::list<double> m_sigma2s;
};


class BSplineKernel : public Kernel {
public:

	BSplineKernel(int scaleLevel)  {
		m_support = (3+1)*0.5; // (degree + 1)*0.5

		m_factor = scaleLevel;

	}

	virtual ~BSplineKernel() {}

	double operator()(const VectorType& x, const VectorType& t) const{
		unsigned n = 3;
		unsigned K = n + 1;
		unsigned dims = x.rows();

		switch(dims) {
		case 1:
			//double k1 = ceil(x - (n + 1) / 2.);
			//k += BSpline3(scaled_x[0] - z0)*BSpline3(scaled_t[0] - z0);
//			for (unsigned k = k1; k < K - 1; k++ ) {
//
//			}
//			break;
		case 2:
			double scale = m_factor;
			double k1 = ceil(x(0)/scale - (n + 1) / 2.);
			double l1 = ceil(x(1)/scale - (n + 1) / 2.);

			double val = 0;
			for (unsigned k = k1 - 1; k <= k1 +  K - 1 + 1; k++) {
				for (unsigned l = l1 - 1 ; l <= l1 +  K - 1 + 1; l++) {
					val += BSpline3(x(0) / scale - k ) * BSpline3(x(1)/scale - l)  * BSpline3(t(0)/scale - k ) * BSpline3(t(1)/scale - l);;
				}
			}

			return val;

		}


	}
	double BSpline32D(const VectorType& x) {
		return BSpline3(x(0)) * (x(1));
	}
		double BSpline3(double x_)const{
			double x = std::fabs(x_);

			if(x < 1){
				return 2.0/3.0 - x*x + (x*x*x)/2.0;
			}
			else if(x < 2){
				return ((2-x)*(2-x)*(2-x))/6.0;
			}
			else{
				return 0;
			}
		}


private:


	double m_factor;
	double m_support;
};


//class BSplineKernel : public Kernel {
//public:
//
//	BSplineKernel(int scaleLevel)  {
//		m_support = (3+1)*0.5; // (degree + 1)*0.5
//
//		if(scaleLevel >= 0){
//			m_factor = pow(2,scaleLevel);
//		}
//		else{
//			m_factor = 1.0/pow(2,-scaleLevel);
//		}
//
//	}
//
//	virtual ~BSplineKernel() {}
//
//	double operator()(const VectorType& x, const VectorType& t) const{
//		double k = 0;
//		unsigned dims = x.rows();
//		VectorType scaled_x = VectorType::Zero(dims);
//		VectorType scaled_t = VectorType::Zero(dims);
//
//		for(unsigned i=0; i<dims; i++){
//			scaled_x[i] = m_factor*x(i);
//			scaled_t[i] = m_factor*t(i);
//		}
//
//		VectorType diff = scaled_t - scaled_x;
//
//		bool isInSupport = true;
//		for(unsigned i=0; i<dims; i++){
//			if(diff[i] - m_support  > 2 || diff[i] + m_support < -2){
//				isInSupport = false;
//			}
//		}
//
//		if(isInSupport){
//
//			// k = floor(2^j*x) - (degree+1)*0.5 : ceil(2^j*x) + (degree+1)*0.5
//			std::vector<int> u_limit_x(dims);
//			std::vector<int> l_limit_x(dims);
//			for(unsigned i=0; i<dims; i++){
//				u_limit_x[i] = ceil(scaled_x[i] + m_support);
//				l_limit_x[i] = floor(scaled_x[i] - m_support);
//			}
//
//			// sum k_1,k_2 elem Z^d phi(2^j*x(1)-k_1)*phi(2^j*t(1)-k_1)*phi(2^j*x(2)-k_2)*phi(2^j*t(2)-k_2)
//			switch(dims){
//			case 1:
//				for(signed z0 = l_limit_x[0]; z0 <= u_limit_x[0]; z0++){
//					k += BSpline3(scaled_x[0] - z0)*BSpline3(scaled_t[0] - z0);
//				}
//				k *= m_factor; // * spacing*spacing
//				break;
//			case 2:
//				for(signed z0 = l_limit_x[0]; z0 <= u_limit_x[0]; z0++){
//					for(signed z1 = l_limit_x[1]; z1 <= u_limit_x[1]; z1++){
//						k += BSpline3(scaled_x[0] - z0)*BSpline3(scaled_x[1] - z1) *
//								BSpline3(scaled_t[0] - z0)*BSpline3(scaled_t[1] - z1);
//					}
//				}
//				k *= m_factor*m_factor; // * spacing*spacing
//				break;
//			case 3:
//				for(signed z0 = l_limit_x[0]; z0 <= u_limit_x[0]; z0++){
//					for(signed z1 = l_limit_x[1]; z1 <= u_limit_x[1]; z1++){
//						for(signed z2 = l_limit_x[2]; z2 <= u_limit_x[2]; z2++){
//							k += BSpline3(scaled_x[0] - z0)*BSpline3(scaled_x[1] - z1)*BSpline3(scaled_x[2] - z2) *
//									BSpline3(scaled_t[0] - z0)*BSpline3(scaled_t[1] - z1)*BSpline3(scaled_t[2] - z2);
//						}
//					}
//				}
//				k *= m_factor*m_factor*m_factor; // * spacing*spacing
//				break;
//			}
//		}
//		return k;
//	}
//
//	float BSpline3(double x)const{
//		x = std::fabs(x);
//
//		if(x < 1){
//			return 2.0/3.0 - x*x + (x*x*x)/2.0;
//		}
//		else if(x < 2){
//			return ((2-x)*(2-x)*(2-x))/6.0;
//		}
//		else{
//			return 0;
//		}
//	}
//
//
//
//private:
//
//
//	double m_factor;
//	double m_support;
//};


class MultiScaleKernel : public Kernel {
public:

	MultiScaleKernel(int lowerScaleLevel, int upperScaleLevel, std::vector<double> lambdas) : m_l(lowerScaleLevel), m_u(upperScaleLevel) , m_lambda(lambdas) {


		if(m_u < m_l) throw std::runtime_error("upper scale level must be greater or equal the lower scale level.");
		if(m_lambda.empty()) throw std::runtime_error("no lambda defined.");
		if((int)m_lambda.size() != m_u-m_l+1){
			std::cout << "#lambdas: " << m_lambda.size() << " should be: " << m_u-m_l+1 << std::endl;
			throw std::runtime_error("too few lambdas defined. ");
		}

		// set up the kernels for each scale elvel
	    for(int j=m_l; j<=m_u; j++){
	    	kernels[j] =  new BSplineKernel(j);
	    }


	}

	virtual ~MultiScaleKernel() {
	    for(int j=m_l; j<=m_u; j++){
	    	delete kernels.at(j);
	    }
	}


	double operator()(const VectorType& x, const VectorType& t) const{
		float k = 0.0;
	    unsigned i=0;
	    for(int j=m_l; j<=m_u; j++){
	    	k += (*kernels.at(j))(x,t)*m_lambda[i++];
	    }
	    return k;
	}

	bool IsSeparable() const{ return true; }

	void Initialize(){
	}



private:
	MultiScaleKernel(const MultiScaleKernel&); //purposely not implemented
	MultiScaleKernel& operator=(const MultiScaleKernel&); //purposely not implemented

	int m_u; // upper scale level
	int m_l; // lower scale level
	std::vector<double> m_lambda;
	std::map<int, BSplineKernel*> kernels;
};




}

#endif /* inclde guard */
