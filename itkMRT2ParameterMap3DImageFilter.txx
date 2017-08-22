/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMRT2ParameterMap3DImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/01/03 22:21:12 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMRT2ParameterMap3DImageFilter_txx
#define __itkMRT2ParameterMap3DImageFilter_txx

#include "itkMRT2ParameterMap3DImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkArray.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_least_squares_function.h"
#include "vnl/algo/vnl_levenberg_marquardt.h"
#include <math.h>

namespace itk
{

class vnl_exponential_function : public vnl_least_squares_function
{
public:
  vnl_exponential_function(bool with_grad, unsigned int n,
    vnl_vector<double> t, vnl_vector<double> s)
    : vnl_least_squares_function(2,n,with_grad ? use_gradient : no_gradient)
    {this->m_NSignals=n; this->m_Time=t; this->m_Signal=s;}

  static double compute(double t, double a, double b) {
    return (a * exp( -t * b ));
    }
  static double compute_a(double t, double a, double b) {
    return (exp( -t * b ));
    }
  static double compute_b(double t, double a, double b) {
    return (- t * a * exp( -t * b ));
    }

  void f(vnl_vector<double> const& x, vnl_vector<double>& y) {
    for (unsigned int i=0; i<this->m_NSignals; ++i)       {
      y[i] = compute(this->m_Time[i], x(0), x(1) ) - this->m_Signal[i];
      }
    }

  void gradf(vnl_vector<double> const& x, vnl_matrix<double> &J) {
    for (unsigned int i=0; i<this->m_NSignals; ++i) {
      J(i,0) = compute_a(this->m_Time[i], x(0), x(1) );
      }
    for (unsigned int i=0; i<this->m_NSignals; ++i) {
      J(i,1) = compute_b(this->m_Time[i], x(0), x(1) );
      }
    }

private:
  unsigned int        m_NSignals;
  vnl_vector<double>  m_Time;
  vnl_vector<double>  m_Signal;
};

class vnl_exponential_with_constant_function : public vnl_least_squares_function
{
public:
  vnl_exponential_with_constant_function(bool with_grad, unsigned int n,
    vnl_vector<double> t, vnl_vector<double> s)
  : vnl_least_squares_function(3, n,with_grad ? use_gradient : no_gradient)
  {this->m_NSignals=n; this->m_Time=t; this->m_Signal=s;}

  static double compute(double t, double a, double b, double c) {
    return (a * exp( -t * b ) + c);
    }
  static double compute_a(double t, double a, double b, double c) {
    return (exp( -t * b ));
    }
  static double compute_b(double t, double a, double b, double c) {
    return (- t * a * exp( -t * b ));
    }
  static double compute_c(double t, double a, double b, double c) {
    return ( 1.0 );
    }

  void f(vnl_vector<double> const& x, vnl_vector<double>& y) {
    for (unsigned int i=0; i<this->m_NSignals; ++i)       {
      y[i] = compute(this->m_Time[i], x(0), x(1), x(2) ) - this->m_Signal[i];
      }
    }

  void gradf(vnl_vector<double> const& x, vnl_matrix<double> &J) {
    for (unsigned int i=0; i<this->m_NSignals; ++i) {
      J(i,0) = compute_a(this->m_Time[i], x(0), x(1), x(2) );
      }
    for (unsigned int i=0; i<this->m_NSignals; ++i) {
      J(i,1) = compute_b(this->m_Time[i], x(0), x(1), x(2) );
      }
    for (unsigned int i=0; i<this->m_NSignals; ++i) {
      J(i,2) = compute_c(this->m_Time[i], x(0), x(1), x(2) );
      }
    }

private:
  unsigned int        m_NSignals;
  vnl_vector<double>  m_Time;
  vnl_vector<double>  m_Signal;
};

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>::MRT2ParameterMap3DImageFilter()
{
  // At least 1 input is necessary for a vector image.
  // For images added one at a time we need at least 2
  this->SetNumberOfRequiredInputs( 1 );
  this->m_NumberOfEchoImages = 0;
  this->m_MaxT2Time = 10.0f;
  this->m_PerformR2Mapping = false;
  this->m_MREchoImageTypeEnumeration = Else;
  this->m_EchoTimeContainer = NULL;
  this->m_Algorithm = LINEAR;
}


template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "PerformR2Mapping: " << (this->m_PerformR2Mapping ? "On"
    : "Off") << std::endl;
  if ( this->m_EchoTimeContainer )
    {
    os << indent << "EchoTimeContainer: " << this->m_EchoTimeContainer
      << std::endl;
    }
  else
    {
    os << indent << "EchoTimeContainer: (Echo times not set)" << std::endl;
    }
  os << indent << "NumberOfEchoImages: " << this->m_NumberOfEchoImages
    << std::endl;
  os << indent << "Maximum T2 time: " << this->m_MaxT2Time << std::endl;
  if ( this->m_MREchoImageTypeEnumeration == MREchoIsInASingleImage )
    {
    os << indent << "MR echo images haven been supplied " << std::endl;
    }
  else if ( this->m_MREchoImageTypeEnumeration == MREchoIsInManyImages )
    {
    os << indent << "A multicomponent MR echo image has been supplied"
      << std::endl;
    }
  if ( this->m_Algorithm == LINEAR )
    {
    os << indent << "The LINEAR algorithm is being used for the T2 fitting"
      << std::endl;
    }
  else if ( this->m_Algorithm == NON_LINEAR )
    {
    os << indent << "The NON_LINEAR algorithm is being used for the T2 fitting"
      << std::endl;
    }
  else
    {
    os << indent << "The NON_LINEAR_WITH_CONSTANT algorithm is being used for "
      "the T2 fitting" << std::endl;
    }
}

//----------------------------------------------------------------------------
template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::GenerateOutputInformation()
{
  // Override the method to set vector length
  Superclass::GenerateOutputInformation();

  typename Superclass::OutputImagePointer output = this->GetOutput();
  if( !output )
    {
    return;
    }
  // Vector length is always 4.
  output->SetVectorLength(4);
}

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::BeforeThreadedGenerateData()
{
  const unsigned int numberOfInputs = this->GetNumberOfInputs();

  // There need to be at least 2 echo images to be able to compute the
  // T2 map.
  if( this->m_NumberOfEchoImages < 2 )
    {
    itkExceptionMacro( << "At least 2 echo images are required" );
    }

  // If there is only 1 echo image, it must be an itk::VectorImage. Otherwise we
  // must have a container of (numberOfInputs-1) itk::Image. Check to make sure
  if ( (numberOfInputs == 1)
      && (this->m_MREchoImageTypeEnumeration != MREchoIsInASingleImage) )
    {
    std::string echoImageClassName(
      this->ProcessObject::GetInput(0)->GetNameOfClass());
    if ( strcmp(echoImageClassName.c_str(),"VectorImage") != 0 )
      {
      itkExceptionMacro( <<
          "There is only one echo image. It should be a VectorImage. "
          << "But its of type: " << echoImageClassName );
      }
    }
}

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::FitLinearExponential(ExponentialFitType X, ExponentialFitType Y,
  unsigned int num, MRParameterMapPixelType &output)
{
  EchoTimeType Sumxy=0, Sumx=0, Sumy=0, Sumx2=0, Sumy2=0, b=0, denom=0;

  for(unsigned int i=0; i<num; i++)
    {
    Sumxy += X[i]*log(Y[i]);
    Sumx += X[i];
    Sumy += log(Y[i]);
    Sumy2 += log(Y[i])*log(Y[i]);
    Sumx2 += X[i]*X[i];
    }
  denom = Sumx2-(Sumx*Sumx/static_cast<EchoTimeType>(num));
  if( denom == 0 )
    {
    b = NumericTraits< EchoTimeType >::max() *
      ((Sumxy-(Sumx*Sumy/static_cast<EchoTimeType>(num))) < 0)?-1.0f:1.0f;
    }
  else
    {
    b = (Sumxy-(Sumx*Sumy/static_cast<EchoTimeType>(num)))/denom;
    }
  if( b == 0 )
    {
    b = NumericTraits< EchoTimeType >::max();
    }
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(-b); // T2
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>
    (exp((Sumy-b*Sumx)/static_cast<EchoTimeType>(num))); // Constant
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>
    ((Sumxy*Sumxy)/(Sumy2*Sumx2)); // R-squared
}

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::FitNonlinearExponential(ExponentialFitType X, ExponentialFitType Y,
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp;
  EchoTimeType averageY = 0;
  EchoTimeType SSE = 0;
  EchoTimeType SST = 0;

  // Find linear fit for A & B as initial estimation.
  FitLinearExponential(X,Y,num,output);

  vnl_exponential_function f(true,num,X,Y);

  vnl_levenberg_marquardt lm(f);

  vnl_vector<EchoTimeType> x1(2);
  x1[0] = output[1];
  x1[1] = output[0];
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);

  // Only use the output if there were no failures.
  switch(lm.get_failure_code())
    {
    case vnl_levenberg_marquardt::CONVERGED_FTOL:
    case vnl_levenberg_marquardt::CONVERGED_XTOL:
    case vnl_levenberg_marquardt::CONVERGED_XFTOL:
    case vnl_levenberg_marquardt::CONVERGED_GTOL:
      output[1] = static_cast<typename MRParameterMapPixelType::ValueType>
        (x1[0]);
      output[0] = static_cast<typename MRParameterMapPixelType::ValueType>
        (x1[1]);
      break;
    default:
      break;
    }

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_exponential_function::compute(X[i],output[1],output[0])
      -Y[i];
    SSE += (err*err);
    averageY += Y[i];
    }
  averageY /= (double)num;
  temp = Y - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::FitNonlinearExponentialWithConstant(ExponentialFitType X, ExponentialFitType Y,
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp;
  EchoTimeType averageY = 0;
  EchoTimeType SSE = 0;
  EchoTimeType SST = 0;

  // Find linear fit for A & B as initial estimation.
  FitLinearExponential(X,Y,num,output);

  vnl_exponential_with_constant_function f(true,num,X,Y);

  vnl_levenberg_marquardt lm(f);

  vnl_vector<EchoTimeType> x1(3);
  x1[0] = output[1];
  x1[1] = output[0];
  x1[2] = 0; // Set initial constant to zero.
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);

  // Only use the output if there were no failures.
  switch(lm.get_failure_code())
    {
    case vnl_levenberg_marquardt::CONVERGED_FTOL:
    case vnl_levenberg_marquardt::CONVERGED_XTOL:
    case vnl_levenberg_marquardt::CONVERGED_XFTOL:
    case vnl_levenberg_marquardt::CONVERGED_GTOL:
      output[1] = static_cast<typename MRParameterMapPixelType::ValueType>
        (x1[0]);
      output[0] = static_cast<typename MRParameterMapPixelType::ValueType>
        (x1[1]);
      output[2] = static_cast<typename MRParameterMapPixelType::ValueType>
        (x1[2]);
      break;
    default:
      break;
    }

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_exponential_with_constant_function::compute(X[i],output[1],
      output[0],output[2])-Y[i];
    SSE += (err*err);
    averageY += Y[i];
    }
  averageY /= (double)num;
  temp = Y - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
  ThreadIdType threadId)
{
  typename OutputImageType::Pointer outputImage =
    static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));

  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

  ProgressReporter progress(this, threadId,
    outputRegionForThread.GetNumberOfPixels(), 100);

  // Two cases here .
  // 1. If the echoes have been specified in multiple images, we will create
  // 'n' iterators for each of the echo images and fit the T2 curve for each
  // voxel.
  // 2. If the echo images have been specified in a single multi-component image,
  // one iterator will suffice to do the same.

  if( this->m_MREchoImageTypeEnumeration == MREchoIsInManyImages )
    {

    typedef ImageRegionConstIterator< MREchoImageType > MREchoIteratorType;
    std::vector< MREchoIteratorType * > echoItContainer;
    int nonzeroCount = 0;

    for( unsigned int i = 0; i< this->m_NumberOfEchoImages; i++ )
      {
      typename MREchoImageType::Pointer echoImagePointer = NULL;

      echoImagePointer = static_cast< MREchoImageType * >(
        this->ProcessObject::GetInput(i) );

      MREchoIteratorType *eit = new MREchoIteratorType( echoImagePointer,
        outputRegionForThread );
      eit->GoToBegin();
      echoItContainer.push_back(eit);
      }

    // Iterate over the echo images and fit the T2 curve.
    while( !oit.IsAtEnd() )
      {

      MRParameterMapPixelType map;
      map.SetSize(4);
      map.Fill(0);
      // Create array for T2 calculation.
      nonzeroCount = -1;
      vnl_vector<EchoTimeType> echoImageValues(this->m_NumberOfEchoImages);
      vnl_vector<EchoTimeType> echoTimes(this->m_NumberOfEchoImages);
      for( unsigned int i = 0; i< this->m_NumberOfEchoImages; i++ )
        {
        echoImageValues[i] = static_cast<EchoTimeType>(
          echoItContainer[i]->Get());
        echoTimes[i] = this->m_EchoTimeContainer->ElementAt(i);
        // The echo image values should never be less than or equal
        // to zero.  This little loop will count the number of contiguous
        // non-zero values, stopping at the first zero reached.  Only the
        // first contiguous non-zero values will be used for the calculation.
        if( (echoImageValues[i] <= 0) && (nonzeroCount < 0) )
          {
          nonzeroCount = i;
          }
        ++(*echoItContainer[i]);
        }
      if( nonzeroCount < 0 )
        {
        nonzeroCount = this->m_NumberOfEchoImages;
        }

      // Only do the calculation if we have at least 2 contiguous non-zero
      // values.
      if( nonzeroCount >= 2 )
        {
        switch( this->m_Algorithm )
          {
          case LINEAR:
            FitLinearExponential(echoTimes,echoImageValues,nonzeroCount,map);
            break;
          case NON_LINEAR:
            FitNonlinearExponential(echoTimes,echoImageValues,nonzeroCount,map);
            break;
          case NON_LINEAR_WITH_CONSTANT:
            FitNonlinearExponentialWithConstant(echoTimes,echoImageValues,
              nonzeroCount,map);
            break;
          default:
            itkExceptionMacro( << "Unknown fit type: " << this->m_Algorithm );
          }
        // Do inverse of R2 if not
        // performing R2 mapping.
        if( !this->m_PerformR2Mapping )
          {
          if( map[0] != 0 )
            {
            map[0] = 1.0f/map[0];
            }
          if( static_cast<EchoTimeType>(map[0]) > this->m_MaxT2Time )
            {
            map[0] = static_cast<MRParameterPixelType>(this->m_MaxT2Time);
            }
          }
        else
          {
          if( this->m_MaxT2Time == 0 )
            {
            this->m_MaxT2Time = 1.0f/NumericTraits< EchoTimeType >::max();
            }
          if( static_cast<EchoTimeType>(map[0]) < (1.0f/this->m_MaxT2Time) )
            {
            map[0] = static_cast<MRParameterPixelType>(1.0f/this->m_MaxT2Time);
            }
          }

        // Should never have a value less than zero,
        // so automatically filter that out.
        if( map[0] < 0 )
          {
          map.Fill(0);
          }
        }

        oit.Set( map );
        ++oit;
        progress.CompletedPixel();
      }

    for( unsigned int i = 0; i< echoItContainer.size(); i++ )
      {
      delete echoItContainer[i];
      }
    }
  // The echoes are specified in a single multi-component image
  else if( this->m_MREchoImageTypeEnumeration == MREchoIsInASingleImage )
    {
    typedef ImageRegionConstIterator< MREchoImagesType > MREchoIteratorType;
    typedef typename MREchoImagesType::PixelType         MREchoVectorType;
    typename MREchoImagesType::Pointer echoImagePointer = NULL;
    int nonzeroCount = 0;

    echoImagePointer = static_cast< MREchoImagesType * >(
      this->ProcessObject::GetInput(0) );

    MREchoIteratorType eit(echoImagePointer, outputRegionForThread );
    eit.GoToBegin();

    while( !eit.IsAtEnd() )
      {
      MREchoVectorType echoValues = eit.Get();

      MRParameterMapPixelType map;
      map.SetSize(4);
      map.Fill(0);
      // Create array for T2 calculation.
      nonzeroCount = -1;
      vnl_vector<EchoTimeType> echoImageValues(this->m_NumberOfEchoImages);
      vnl_vector<EchoTimeType> echoTimes(this->m_NumberOfEchoImages);
      for( unsigned int i = 0; i< this->m_NumberOfEchoImages; i++ )
        {
        echoImageValues[i] = static_cast<EchoTimeType>(echoValues[i]);
        echoTimes[i] = this->m_EchoTimeContainer->ElementAt(i);
        // The echo image values should never be less than or equal
        // to zero.  This little loop will count the number of contiguous
        // non-zero values, stopping at the first zero reached.  Only the
        // first contiguous non-zero values will be used for the calculation.
        if( (echoImageValues[i] <= 0) && (nonzeroCount < 0) )
          {
          nonzeroCount = i;
          }
        }
      if( nonzeroCount < 0 )
        {
        nonzeroCount = this->m_NumberOfEchoImages;
        }

      // Only do the calculation if we have at least 2 contiguous non-zero
      // values.
      if( nonzeroCount >= 2 )
        {
        switch( this->m_Algorithm )
          {
          case LINEAR:
            FitLinearExponential(echoTimes,echoImageValues,nonzeroCount,map);
            break;
          case NON_LINEAR:
            FitNonlinearExponential(echoTimes,echoImageValues,nonzeroCount,map);
            break;
          case NON_LINEAR_WITH_CONSTANT:
            FitNonlinearExponentialWithConstant(echoTimes,echoImageValues,
              nonzeroCount,map);
            break;
          default:
            itkExceptionMacro( << "Unknown fit type: " << this->m_Algorithm );
          }

        // Do inverse of R2 if not
        // performing R2 mapping.
        if( !this->m_PerformR2Mapping )
          {
          if( map[0] != 0 )
            {
            map[0] = 1.0f/map[0];
            }
          if( static_cast<EchoTimeType>(map[0]) > this->m_MaxT2Time )
            {
            map[0] = static_cast<MRParameterPixelType>(this->m_MaxT2Time);
            }
          }
        else
          {
          if( this->m_MaxT2Time == 0 )
            {
            this->m_MaxT2Time = 1.0f/NumericTraits< EchoTimeType >::max();
            }
          if( static_cast<EchoTimeType>(map[0]) < (1.0f/this->m_MaxT2Time) )
            {
            map[0] = static_cast<MRParameterPixelType>(1.0f/this->m_MaxT2Time);
            }
          }

        // Should never have a value less than zero,
        // so automatically filter that out.
        if( map[0] < 0 )
          {
          map.Fill(0);
          }
        }

      oit.Set( map );
      ++oit; // Output (fitted image parameters) iterator
      ++eit; // Echo image iterator
      progress.CompletedPixel();
      }
  }
}

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::AddMREchoImage( EchoTimeType echoTime, const MREchoImageType *image)
{
  // Make sure crazy users did not call both AddMREchoImage and
  // SetMREchoImage
  if( this->m_MREchoImageTypeEnumeration == MREchoIsInASingleImage)
    {
    itkExceptionMacro( << "Cannot call both methods:"
    << "AddMREchoImage and SetMREchoImage. Please call only one of them.");
    }

  // If the container to hold the echo times hasn't been allocated
  // yet, allocate it.
  if( !this->m_EchoTimeContainer )
    {
    this->m_EchoTimeContainer = EchoTimeContainerType::New();
    }

  m_EchoTimeContainer->InsertElement( this->m_NumberOfEchoImages, echoTime );
  this->ProcessObject::SetNthInput( this->m_NumberOfEchoImages,
      const_cast< MREchoImageType* >(image) );
  ++this->m_NumberOfEchoImages;
  this->m_MREchoImageTypeEnumeration = MREchoIsInManyImages;
}

template< class TMREchoImagePixelType, class TMRParameterMapImagePixelType >
void
MRT2ParameterMap3DImageFilter<TMREchoImagePixelType,
  TMRParameterMapImagePixelType>
::SetMREchoImage( EchoTimeContainerType *echoContainer,
  const MREchoImagesType *image )
{
  // Make sure crazy users did not call both AddMREchoImage and
  // SetMREchoImage
  if( this->m_MREchoImageTypeEnumeration == MREchoIsInManyImages )
    {
    itkExceptionMacro( << "Cannot call both methods:"
    << "AddMREchoImage and SetMREchoImage. Please call only one of them.");
    }

  this->m_EchoTimeContainer = echoContainer;

  this->m_NumberOfEchoImages = echoContainer->Size();

  // ensure that the echo image we received has as many components as
  // the number of echo times
  if( image->GetVectorLength() != this->m_NumberOfEchoImages )
    {
    itkExceptionMacro( << this->m_NumberOfEchoImages <<
      " echo times specified but image has " << image->GetVectorLength()
      << " components.");
    }

  this->ProcessObject::SetNthInput( 0,
      const_cast< MREchoImagesType* >(image) );
  this->m_MREchoImageTypeEnumeration = MREchoIsInASingleImage;
}


} // end namespace itk

#endif
