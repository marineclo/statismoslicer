/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
 *
 * Copyright (c) 2011 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#ifndef ITKCOVARIANCEFUNCIONMODELBUILDER_H_
#define ITKCOVARIANCEFUNCIONMODELBUILDER_H_

#include "itkObject.h"
#include "itkObjectFactory.h"

#include "statismoITKConfig.h"
#include "itkDataManager.h"
#include "itkStatisticalModel.h"
#include "CovarianceFunctionModelBuilder.h"

namespace itk
{

/**
 * \brief ITK Wrapper for the statismo::PCAModelBuilder class.
 * \see statismo::PCAModelBuilder for detailed documentation.
 */
template <class Representer>
class CovarianceFunctionModelBuilder : public Object {
public:

	typedef CovarianceFunctionModelBuilder            Self;
	typedef Object	Superclass;
	typedef SmartPointer<Self>                Pointer;
	typedef SmartPointer<const Self>          ConstPointer;

	itkNewMacro( Self );
	itkTypeMacro( CovarianceFunctionModelBuilder, Object );


	typedef statismo::CovarianceFunctionModelBuilder<Representer> ImplType;

	CovarianceFunctionModelBuilder() : m_impl(0){}

	virtual ~CovarianceFunctionModelBuilder() {
		if (m_impl) {
			delete m_impl;
			m_impl = 0;
		}
	}


	ImplType* GetstatismoImplObj() const { return m_impl; }

  	void SetstatismoImplObj(ImplType* impl) {
  		if (m_impl) {
  			delete m_impl;
  		}
  		m_impl = impl;
  	}

	void SetRepresenter(const Representer* representer) {
		SetstatismoImplObj(ImplType::Create(representer));
	}



	template <class F>
	typename std::tr1::result_of<F()>::type callstatismoImpl(F f) const {
		try {
			  return f();
		}
		 catch (statismo::StatisticalModelException& s) {
			itkExceptionMacro(<< s.what());
		}
	}



	typename StatisticalModel<Representer>::Pointer BuildNewModel(const statismo::Kernel* kernel, unsigned numPointsForNystrom, unsigned numComponents) const {
		statismo::StatisticalModel<Representer>* model_statismo = callstatismoImpl(std::tr1::bind(&ImplType::BuildNewModel, this->m_impl, kernel, numPointsForNystrom, numComponents));
		typename StatisticalModel<Representer>::Pointer model_itk = StatisticalModel<Representer>::New();
		model_itk->SetstatismoImplObj(model_statismo);
		return model_itk;
	}


private:
	CovarianceFunctionModelBuilder(const CovarianceFunctionModelBuilder& orig);
	CovarianceFunctionModelBuilder& operator=(const CovarianceFunctionModelBuilder& rhs);

	ImplType* m_impl;
};


}

#endif /* ITKCOVARIANCEFUNCIONMODELBUILDER_H_ */
