/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMRT2ParameterMap3DImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007/10/03 18:00:00 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMRT2ParameterMap3DImageFilter_h
#define __itkMRT2ParameterMap3DImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include "vnl/vnl_vector.h"

namespace itk
{

/** \class MRT2ParameterMap3DImageFilter
 * \brief Generate a 3D MR T2 parameter map using multiple echo images.
 *
 * This filter is templated over the input pixel type and the output pixel
 * type. The input T2-weighted MR images must all have the same size and
 * dimensions.  The 3D \c VectorImage output will have four components.
 * \c VectorImageToImageAdaptor should be used to extract the various
 * components for viewing or saving to file.
 *
 * Given a multi-echo image data set of T2-weighted MR images and a vector
 * of TE (echo) times, generate the T2 or R2 parameter map.  The exponential
 * fitting function may be changed using the \c SetAlgorithm function.
 *
 * \par Inputs and Usage
 * There are two ways to use this class. When you have multiple echo images
 * you would use the class as
 * \code
 *       filter->AddMREchoImage( echoTime1, image1 );
 *       filter->AddMREchoImage( echoTime2, image2 );
 *   ...
 * \endcode
 *
 * \par
 * When you have the 'n' echo images in a single multi-component image
 * (VectorImage), you can specify the images simply as
 * \code
 *       filter->SetMREchoImage( echoTimeContainer, vectorImage );
 * \endcode
 *
 * \par Outputs
 * The output image is a vector image containing 4 values:
 * The first component of the output will be the T2 time in seconds (R2 in Hz),
 * the second component will be the constant A as shown in the function below,
 * and the fourth component will be the r-squared value from the curve fitting.
 * The third component will vary depending on the type of T2 fitting selected.
 * For LINEAR and NON_LINEAR the third component will be zero.  For
 * NON_LINEAR_WITH_CONSTANT the third component will be the value C as shown
 * below.
 *
 * LINEAR (Linear least squares):
 * S = Ae^-(TE/T2)
 *
 * NON_LINEAR (Non-linear least squares using Levenberg-Marquardt):
 * S = Ae^-(TE/T2)
 *
 * NON_LINEAR_WITH_CONSTANT (Non-linear least squares using Levenberg-Marquardt):
 * S = Ae^-(TE/T2) + C
 *
 * \code
 *       VectorImage< TMRParameterMapImagePixelType, 3 >
 * \endcode
 *
 * \par Parameters
 * \li Algorithm -  Set/Get the T2 fitting algorithm used
 * (LINEAR, NON_LINEAR, and NON_LINEAR_WITH_CONSTANT).
 * \li MaxT2Time - Set the maximum T2 time (T2 times greater than or equal to
 * to this value will be set to this value).
 * \li PerformR2Mapping - If On R2 (in Hz) will be calculated instead of T2.
 *
 * \ingroup Multithreaded
 *
 * \author Don Bigler, Center for Nuclear Magnetic Resonance Research,
 * Penn State Milton S. Hershey Medical Center
 *
 */
template< class TMREchoImagePixelType,
  class TMRParameterMapImagePixelType=double >
class ITK_EXPORT MRT2ParameterMap3DImageFilter:
    public ImageToImageFilter<Image< TMREchoImagePixelType, 3 >,
                              VectorImage< TMRParameterMapImagePixelType, 3 > >
{
public:
  /** Standard class typedefs. */
  typedef MRT2ParameterMap3DImageFilter                 Self;
  typedef ImageToImageFilter<Image< TMREchoImagePixelType, 3 >,
    VectorImage< TMRParameterMapImagePixelType, 3 > >   Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MRT2ParameterMap3DImageFilter, ImageToImageFilter);


  typedef TMREchoImagePixelType                         MREchoPixelType;

  typedef typename VectorImage< TMRParameterMapImagePixelType, 3 >::PixelType
                                                        MRParameterMapPixelType;
  typedef TMRParameterMapImagePixelType                 MRParameterPixelType;

  /** Multi-echo T2-weighted MR image data. */
  typedef typename Superclass::InputImageType           MREchoImageType;

  /** An alternative typedef defining one (of the many) echo images.
   * It will be assumed that the vectorImage has a vector length
   * parameter of \c n (number of echo times) */
  typedef VectorImage< MREchoPixelType, 3 >             MREchoImagesType;

  typedef typename Superclass::OutputImageType          MRParameterMapImageType;
  typedef MRParameterMapImageType                       OutputImageType;
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;

  /** Data type of echo times */
  typedef double                                        EchoTimeType;
  typedef vnl_vector<EchoTimeType>                      ExponentialFitType;

  /** Container to hold echo times of the 'n' MR echo image measurements */
  typedef VectorContainer< unsigned int,EchoTimeType >  EchoTimeContainerType;

  /** Set method to add an echo time (in seconds) and its corresponding image. */
  void AddMREchoImage( EchoTimeType, const MREchoImageType *image);

  /** Another set method to add an echo time (in seconds) and its corresponding
   * image. The image here is a VectorImage. The user is expected to pass the
   * echo times in a container. The ith element of the container corresponds
   * to the echo time of the ith component image of the VectorImage. */
  void SetMREchoImage( EchoTimeContainerType *, const MREchoImagesType *image);

  /** define values used to determine which algorithm to use */
  static const int LINEAR = 0;
  static const int NON_LINEAR = 1;
  static const int NON_LINEAR_WITH_CONSTANT = 2;

  /** Set/Get the T2 fitting algorithm (default is LINEAR). */
  itkSetMacro(Algorithm, int);
  itkGetMacro(Algorithm, int);

  /** Set/Get the maximum T2 time (default is 10.0). */
  itkSetMacro(MaxT2Time, EchoTimeType);
  itkGetMacro(MaxT2Time, EchoTimeType);

  /** Set/Get whether to perform R2 mapping instead of T2 mapping
  * (default is false). */
  itkSetMacro( PerformR2Mapping, bool );
  itkGetConstReferenceMacro( PerformR2Mapping, bool );
  itkBooleanMacro( PerformR2Mapping );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(MREchoPixelConvertibleToDoubleCheck,
    (Concept::Convertible<MREchoPixelType, double>));
  itkConceptMacro(MRParameterMapPixelConvertibleToDoubleCheck,
    (Concept::Convertible<MRParameterPixelType, double>));
 /** End concept checking */
#endif

protected:
  MRT2ParameterMap3DImageFilter();
  ~MRT2ParameterMap3DImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void FitLinearExponential(ExponentialFitType X, ExponentialFitType Y,
    unsigned int num, MRParameterMapPixelType &output);
  void FitNonlinearExponential(ExponentialFitType X, ExponentialFitType Y,
    unsigned int num, MRParameterMapPixelType &output);
  void FitNonlinearExponentialWithConstant(ExponentialFitType X,
    ExponentialFitType Y, unsigned int num, MRParameterMapPixelType &output);

  /** Setup vector image vector length. */
  virtual void GenerateOutputInformation();

  void BeforeThreadedGenerateData();
  void ThreadedGenerateData( const
      OutputImageRegionType &outputRegionForThread, ThreadIdType);

  /** enum to indicate if the MR echo image is specified as a single multi-
   * component image or as several separate images */
  typedef enum
    {
    MREchoIsInASingleImage = 1,
    MREchoIsInManyImages,
    Else
    } MREchoImageTypeEnumeration;

private:
  MRT2ParameterMap3DImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  int                             m_Algorithm;

  EchoTimeType                    m_MaxT2Time;

  bool                            m_PerformR2Mapping;

  /** container to hold echo times */
  EchoTimeContainerType::Pointer  m_EchoTimeContainer;

  /** Number of echo images */
  unsigned int                    m_NumberOfEchoImages;

  /** MR echo image was specified in a single image or in multiple images */
  MREchoImageTypeEnumeration      m_MREchoImageTypeEnumeration;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRT2ParameterMap3DImageFilter.txx"
#endif

#endif
