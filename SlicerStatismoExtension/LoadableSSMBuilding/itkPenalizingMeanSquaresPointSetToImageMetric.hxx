/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkPenalizingMeanSquaresPointSetToImageMetric_hxx
#define __itkPenalizingMeanSquaresPointSetToImageMetric_hxx

#include "itkPenalizingMeanSquaresPointSetToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
namespace itk
{


/**
 * Constructor
 */
template <class TFixedPointSet, class TMovingImage>
PenalizingMeanSquaresPointSetToImageMetric<TFixedPointSet, TMovingImage>
::PenalizingMeanSquaresPointSetToImageMetric() : m_regularizationParameter(0)
{
}

/**
 * Get the match Measure
 */
template <class TFixedPointSet, class TMovingImage>
typename PenalizingMeanSquaresPointSetToImageMetric<TFixedPointSet, TMovingImage>::MeasureType
PenalizingMeanSquaresPointSetToImageMetric<TFixedPointSet, TMovingImage>
::GetValue(const TransformParametersType & parameters) const
{
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if( !fixedPointSet )
    {
    itkExceptionMacro(<< "Fixed point set has not been assigned");
    }

  PointIterator pointItr = fixedPointSet->GetPoints()->Begin();
  PointIterator pointEnd = fixedPointSet->GetPoints()->End();

  PointDataIterator pointDataItr = fixedPointSet->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedPointSet->GetPointData()->End();

  MeasureType measure = NumericTraits<MeasureType>::Zero;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters(parameters);

  typedef  typename NumericTraits<MeasureType>::AccumulateType AccumulateType;

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    OutputPointType transformedPoint =
      this->m_Transform->TransformPoint(inputPoint);

    if( this->m_Interpolator->IsInsideBuffer(transformedPoint) )
      {
      const RealType movingValue  = this->m_Interpolator->Evaluate(transformedPoint);
      const RealType fixedValue   = pointDataItr.Value();
      const RealType diff = movingValue - fixedValue;
      measure += diff * diff;
      this->m_NumberOfPixelsCounted++;
      }

    ++pointItr;
    ++pointDataItr;
    }

  double regValue = 0;
  for( unsigned int par = 0; par < parameters.Size(); par++ )
    {
      regValue += parameters[par] * parameters[par];
    }



  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
  else
    {
    measure /= this->m_NumberOfPixelsCounted;
    }




  return measure + m_regularizationParameter * regValue;
}

/**
 * Get the Derivative Measure
 */
template <class TFixedPointSet, class TMovingImage>
void
PenalizingMeanSquaresPointSetToImageMetric<TFixedPointSet, TMovingImage>
::GetDerivative(const TransformParametersType & parameters,
                DerivativeType & derivative) const
{
  if( !this->GetGradientImage() )
    {
    itkExceptionMacro(<< "The gradient image is null, maybe you forgot to call Initialize()");
    }

  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if( !fixedPointSet )
    {
    itkExceptionMacro(<< "Fixed image has not been assigned");
    }

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters(parameters);

  typedef  typename NumericTraits<MeasureType>::AccumulateType AccumulateType;

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType(ParametersDimension);
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::Zero);

  PointIterator pointItr = fixedPointSet->GetPoints()->Begin();
  PointIterator pointEnd = fixedPointSet->GetPoints()->End();

  PointDataIterator pointDataItr = fixedPointSet->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedPointSet->GetPointData()->End();

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    OutputPointType transformedPoint =
      this->m_Transform->TransformPoint(inputPoint);

    if( this->m_Interpolator->IsInsideBuffer(transformedPoint) )
      {
      const RealType movingValue  = this->m_Interpolator->Evaluate(transformedPoint);
      const RealType fixedValue   = pointDataItr.Value();

      this->m_NumberOfPixelsCounted++;
      const RealType diff = movingValue - fixedValue;

      // Now compute the derivatives
      TransformJacobianType jacobian;
      this->m_Transform->ComputeJacobianWithRespectToParameters(inputPoint, jacobian);

      // Get the gradient by NearestNeighboorInterpolation:
      // which is equivalent to round up the point components.
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef ContinuousIndex<CoordRepType, MovingImageType::ImageDimension>
      MovingImageContinuousIndexType;

      MovingImageContinuousIndexType tempIndex;
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex(transformedPoint, tempIndex);

      typename MovingImageType::IndexType mappedIndex;
      mappedIndex.CopyWithRound(tempIndex);

      const GradientPixelType gradient =
        this->GetGradientImage()->GetPixel(mappedIndex);
      for( unsigned int par = 0; par < ParametersDimension; par++ )
        {
        RealType sum = NumericTraits<RealType>::Zero;
        for( unsigned int dim = 0; dim < Self::FixedPointSetDimension; dim++ )
          {
          sum += 2.0 *diff *jacobian(dim, par) * gradient[dim];
          }
        derivative[par] += sum;
        }
      }

    ++pointItr;
    ++pointDataItr;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
  else
    {
    for( unsigned int i = 0; i < ParametersDimension; i++ )
      {
      derivative[i] /= this->m_NumberOfPixelsCounted;

      derivative[i] += parameters[i] * 2 * m_regularizationParameter;
      }
    }


}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <class TFixedPointSet, class TMovingImage>
void
PenalizingMeanSquaresPointSetToImageMetric<TFixedPointSet, TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
  if( !this->GetGradientImage() )
    {
    itkExceptionMacro(<< "The gradient image is null, maybe you forgot to call Initialize()");
    }

  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if( !fixedPointSet )
    {
    itkExceptionMacro(<< "Fixed image has not been assigned");
    }

  this->m_NumberOfPixelsCounted = 0;
  MeasureType measure = NumericTraits<MeasureType>::Zero;

  this->SetTransformParameters(parameters);

  typedef  typename NumericTraits<MeasureType>::AccumulateType AccumulateType;

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType(ParametersDimension);
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::Zero);

  PointIterator pointItr = fixedPointSet->GetPoints()->Begin();
  PointIterator pointEnd = fixedPointSet->GetPoints()->End();

  PointDataIterator pointDataItr = fixedPointSet->GetPointData()->Begin();
  PointDataIterator pointDataEnd = fixedPointSet->GetPointData()->End();

  while( pointItr != pointEnd && pointDataItr != pointDataEnd )
    {
    InputPointType inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    OutputPointType transformedPoint =
      this->m_Transform->TransformPoint(inputPoint);

    if( this->m_Interpolator->IsInsideBuffer(transformedPoint) )
      {
      const RealType movingValue  = this->m_Interpolator->Evaluate(transformedPoint);
      const RealType fixedValue   = pointDataItr.Value();

      this->m_NumberOfPixelsCounted++;

      // Now compute the derivatives
      TransformJacobianType jacobian;
      this->m_Transform->ComputeJacobianWithRespectToParameters(inputPoint, jacobian);

      const RealType diff = movingValue - fixedValue;

      measure += diff * diff;

      // Get the gradient by NearestNeighboorInterpolation:
      // which is equivalent to round up the point components.
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef ContinuousIndex<CoordRepType, MovingImageType::ImageDimension>
      MovingImageContinuousIndexType;

      MovingImageContinuousIndexType tempIndex;
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex(transformedPoint, tempIndex);

      typename MovingImageType::IndexType mappedIndex;
      mappedIndex.CopyWithRound(tempIndex);

      const GradientPixelType gradient =
        this->GetGradientImage()->GetPixel(mappedIndex);
      for( unsigned int par = 0; par < ParametersDimension; par++ )
        {
        RealType sum = NumericTraits<RealType>::Zero;
        for( unsigned int dim = 0; dim < Self::FixedPointSetDimension; dim++ )
          {
          sum += 2.0 *diff *jacobian(dim, par) * gradient[dim];
          }
        derivative[par] += sum;
        }
      }

    ++pointItr;
    ++pointDataItr;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
  else
    {
    for( unsigned int i = 0; i < ParametersDimension; i++ )
      {
      derivative[i] /= this->m_NumberOfPixelsCounted;
      derivative[i] += 2 * m_regularizationParameter * parameters[i];
      }
    measure /= this->m_NumberOfPixelsCounted;
    }



  double regValue = 0;
  for( unsigned int par = 0; par < parameters.Size(); par++ )
    {
      regValue += parameters[par] * parameters[par];
    }


  value = measure + m_regularizationParameter * regValue;
}

} // end namespace itk

#endif
