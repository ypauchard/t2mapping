#include "itkWin32Header.h"
#include <iostream>
#include <fstream>
//#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMRT2ParameterMap3DImageFilter.h"
//#include "itkVectorContainer.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "vnl/vnl_vector_fixed.h"
#include <algorithm>

int main(int argc, char **argv)
{
  if( argc < 8 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " outBaseName algorithm threshold image1 TE1 image2 TE2 [image3 TE3] [image4 TE4] etc.";
    std::cerr << " " << std::endl;
    return 1;
  }

//build/t2mapping data/mojo_0033 2 0 data/mojo_0033_18ms.mha 18.0 data/mojo_0033_45ms.mha 45.0 data/mojo_0033_99ms.mha 99.0 data/mojo_0033_126ms.mha 126.0

  std::string outBaseName(argv[1]);
  std::string outputT2Filename = outBaseName + "_t2.mha";
  std::string outputExpConstFilename = outBaseName + "_A.mha";
  std::string outputConstFilename = outBaseName + "_C.mha";
  std::string outputRSquaredFilename = outBaseName + "_R2.mha";

  int algorithm = atoi( argv[2] );

  float threshold = atof(argv[3]);

  typedef itk::Image< float, 4 > ImageType4D;
  typedef itk::Image< float, 3 > ImageType;

  typedef itk::ThresholdImageFilter< ImageType >
    ThresholdImageFilterType;

  typedef itk::MRT2ParameterMap3DImageFilter< ImageType::PixelType >
    MRT2ParameterMap3DImageFilterType;

  typedef itk::VectorIndexSelectionCastImageFilter<
    MRT2ParameterMap3DImageFilterType::OutputImageType, ImageType >
      VectorIndexSelectionCastImageFilterType;

  typedef itk::ImageFileReader< ImageType >
        ReaderType;

  typedef itk::ImageFileWriter< ImageType >
      WriterType;

  bool r2Mapping = false;
  //int algorithm = MRT2ParameterMap3DImageFilterType::NON_LINEAR_WITH_CONSTANT;
  double maxT2Time = 1000.0f; //in milli-seconds




  // Create T2 mapping class.
  MRT2ParameterMap3DImageFilterType::Pointer t2Map =
  MRT2ParameterMap3DImageFilterType::New();
  // Select the fit type.
  switch(algorithm)
  {
    case MRT2ParameterMap3DImageFilterType::LINEAR:
    t2Map->SetAlgorithm(MRT2ParameterMap3DImageFilterType::LINEAR);
    break;
    case MRT2ParameterMap3DImageFilterType::NON_LINEAR:
    t2Map->SetAlgorithm(MRT2ParameterMap3DImageFilterType::NON_LINEAR);
    break;
    case MRT2ParameterMap3DImageFilterType::NON_LINEAR_WITH_CONSTANT:
    t2Map->SetAlgorithm(
      MRT2ParameterMap3DImageFilterType::NON_LINEAR_WITH_CONSTANT);
      break;
      default:
      std::cerr << "In valid algorithm = " << algorithm << std::endl;
      return 1;
    }
    t2Map->SetMaxT2Time(maxT2Time);
    if( r2Mapping )
    {
      t2Map->PerformR2MappingOn();
    }




    //foreach image
    for(int i = 4; i<argc-1; i += 2 )
    {
      std::cout << "Reading " << argv[i] << std::endl;
      ReaderType::Pointer imageIO = ReaderType::New();
      imageIO->SetFileName(argv[i]);
      try
        {
        imageIO->Update();
        }
      catch( itk::ExceptionObject &err )
        {
        std::cerr << "ExceptionObject caught";
        std::cerr << " : "  << err.GetDescription();
        return 1;
        }

        if (threshold > 0)
        {
          // Threshold the image.
          std::cout << "Applying threshold " << threshold << std::endl;
          ThresholdImageFilterType::Pointer imageThreshold = ThresholdImageFilterType::New();
          imageThreshold->SetOutsideValue(0.0f);
          imageThreshold->SetInput(imageIO->GetOutput());
          imageThreshold->ThresholdBelow(threshold);
          imageThreshold->Update();

          t2Map->AddMREchoImage(atoi(argv[i+1]), //time in milli-seconds
          imageThreshold->GetOutput());
        }else
        {
          t2Map->AddMREchoImage(atoi(argv[i+1]), //time in milli-seconds
          imageIO->GetOutput());
        }

    }

    std::cout<< t2Map <<std::endl;
    t2Map->Update();
    // Extract each output component and write to disk.
    VectorIndexSelectionCastImageFilterType::Pointer extractComp =
    VectorIndexSelectionCastImageFilterType::New();
    extractComp->SetInput(t2Map->GetOutput());
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(extractComp->GetOutput());

    // T2/R2 map.
    writer->SetFileName( outputT2Filename.c_str() );
    extractComp->SetIndex(0);
    try
    {
      writer->Update();
    }
    catch(...)
    {
      std::cerr << "Error during write of " << outputT2Filename << std::endl;
      return 1;
    }
    // Exponent Constant map.
    extractComp->SetIndex(1);
    writer->SetFileName( outputExpConstFilename.c_str() );
    try
    {
      writer->Update();
    }
    catch(...)
    {
      std::cerr << "Error during write of " << outputExpConstFilename << std::endl;
      return 1;
    }
    // Constant map.
    extractComp->SetIndex(2);
    writer->SetFileName( outputConstFilename.c_str() );
    try
    {
      writer->Update();
    }
    catch(...)
    {
      std::cerr << "Error during write of " << outputConstFilename << std::endl;
      return 1;
    }
    // Rsquared map.
    extractComp->SetIndex(3);
    writer->SetFileName( outputRSquaredFilename.c_str() );
    try
    {
      writer->Update();
    }
    catch(...)
    {
      std::cerr << "Error during write of " << outputRSquaredFilename << std::endl;
      return 1;
    }

    return 0;
  }
