#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "ColorDecompositionCLP.h"


//--------------------------------------------------------------------------------
// I added theses. Yi Gao

// itk
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"



// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

  const unsigned int ImageDimension2D = 2;

  typedef itk::RGBPixel<unsigned char> itkRGBPixelType;
  typedef itk::Image<itkRGBPixelType, ImageDimension2D> itkRGBImage2DType;

  typedef unsigned char OutputPixelType;
  typedef itk::Image<OutputPixelType, ImageDimension2D> itkUCharImage2DType;

  itkUCharImage2DType::Pointer ExtractChannel(itkRGBImage2DType::Pointer rgbImage, int outputChannelId, std::string stainType);

  int DoIt( int argc, char * argv[] )
  {
    PARSE_ARGS;

    typedef itk::ImageFileReader<itkRGBImage2DType>  ReaderType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputVolume.c_str() );
    reader->Update();
    itkRGBImage2DType::Pointer rgbImage = reader->GetOutput();

    itkUCharImage2DType::Pointer outputImage = ExtractChannel(rgbImage, outputChannel - 1, stainChoice);
    // "outputChannel - 1" because outputChannel is 1-based for easier UI

    typedef itk::ImageFileWriter<itkUCharImage2DType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputVolume.c_str() );
    writer->SetInput( outputImage );
    writer->SetUseCompression(1);
    writer->Update();

    return EXIT_SUCCESS;
  }



  itkUCharImage2DType::Pointer ExtractChannel(itkRGBImage2DType::Pointer rgbImage, int outputChannelId, std::string stainType)
  {
    if (outputChannelId < 0 || outputChannelId > 2)
      {
        std::cerr<<"outputChannelId should be in {0, 1, 2}. But got "<<outputChannelId<<std::endl;
        abort();
      }

    int width = rgbImage->GetLargestPossibleRegion().GetSize()[0];
    int height = rgbImage->GetLargestPossibleRegion().GetSize()[1];

    double leng, A, V, C, log255=log(255.0);
    int i,j;
    //int outputChannelId = 1; ///< 0: Hematoxylin; 1: Eosin

    double MODx[3];
    double MODy[3];
    double MODz[3];
    double cosx[3];
    double cosy[3];
    double cosz[3];
    double len[3];
    double q[9];
    // unsigned char rLUT[256];
    // unsigned char gLUT[256];
    // unsigned char bLUT[256];


    if (0 == stainType.compare("H-E"))
      {
        // GL Haem matrix
        MODx[0]= 0.644211; //0.650;
        MODy[0]= 0.716556; //0.704;
        MODz[0]= 0.266844; //0.286;
        // GL Eos matrix
        MODx[1]= 0.092789; //0.072;
        MODy[1]= 0.954111; //0.990;
        MODz[1]= 0.283111; //0.105;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("H-E 2"))
      {
        // GL Haem matrix
        MODx[0]= 0.49015734;
        MODy[0]= 0.76897085;
        MODz[0]= 0.41040173;
        // GL Eos matrix
        MODx[1]= 0.04615336;
        MODy[1]= 0.8420684;
        MODz[1]= 0.5373925;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("H-DAB"))
      {
        // 3,3-diamino-benzidine tetrahydrochloride
        // Haem matrix
        MODx[0]= 0.650;
        MODy[0]= 0.704;
        MODz[0]= 0.286;
        // DAB matrix
        MODx[1]= 0.268;
        MODy[1]= 0.570;
        MODz[1]= 0.776;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("Feulgen Light Green"))
      {
        //GL Feulgen & light green
        //Feulgen
        MODx[0]= 0.46420921;
        MODy[0]= 0.83008335;
        MODz[0]= 0.30827187;
        // light green
        MODx[1]= 0.94705542;
        MODy[1]= 0.25373821;
        MODz[1]= 0.19650764;
        // Zero matrix
        MODx[2]= 0.0; // 0.0010000
        MODy[2]= 0.0; // 0.47027777
        MODz[2]= 0.0; //0.88235928
      }
    else if (0 == stainType.compare("Giemsa"))
      {
        // GL  Methylene Blue and Eosin
        MODx[0]= 0.834750233;
        MODy[0]= 0.513556283;
        MODz[0]= 0.196330403;
        // GL Eos matrix
        MODx[1]= 0.092789;
        MODy[1]= 0.954111;
        MODz[1]= 0.283111;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("FastRed FastBlue DAB"))
      {
        //fast red
        MODx[0]= 0.21393921;
        MODy[0]= 0.85112669;
        MODz[0]= 0.47794022;
        // fast blue
        MODx[1]= 0.74890292;
        MODy[1]= 0.60624161;
        MODz[1]= 0.26731082;
        // dab
        MODx[2]= 0.268;
        MODy[2]= 0.570;
        MODz[2]= 0.776;
      }
    else if (0 == stainType.compare("Methyl Green DAB"))
      {
        // MG matrix (GL)
        MODx[0]= 0.98003;
        MODy[0]= 0.144316;
        MODz[0]= 0.133146;
        // DAB matrix
        MODx[1]= 0.268;
        MODy[1]= 0.570;
        MODz[1]= 0.776;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("H-E DAB"))
      {
        // Haem matrix
        MODx[0]= 0.650;
        MODy[0]= 0.704;
        MODz[0]= 0.286;
        // Eos matrix
        MODx[1]= 0.072;
        MODy[1]= 0.990;
        MODz[1]= 0.105;
        // DAB matrix
        MODx[2]= 0.268;
        MODy[2]= 0.570;
        MODz[2]= 0.776;
      }
    else if (0 == stainType.compare("H-AEC"))
      {
        // 3-amino-9-ethylcarbazole
        // Haem matrix
        MODx[0]= 0.650;
        MODy[0]= 0.704;
        MODz[0]= 0.286;
        // AEC matrix
        MODx[1]= 0.2743;
        MODy[1]= 0.6796;
        MODz[1]= 0.6803;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("Azan-Mallory"))
      {
        //Azocarmine and Aniline Blue (AZAN)
        // GL Blue matrix Anilline Blue
        MODx[0]= .853033;
        MODy[0]= .508733;
        MODz[0]= .112656;
        // GL Red matrix Azocarmine
        MODx[1]=0.09289875;
        MODy[1]=0.8662008;
        MODz[1]=0.49098468;
        //GL  Orange matrix Orange-G
        MODx[2]=0.10732849;
        MODy[2]=0.36765403;
        MODz[2]=0.9237484;
      }
    else if (0 == stainType.compare("Masson Trichrome"))
      {
        // GL Methyl blue
        MODx[0]=0.7995107;
        MODy[0]=0.5913521;
        MODz[0]=0.10528667;
        // GL Ponceau Fuchsin has 2 hues, really this is only approximate
        MODx[1]=0.09997159;
        MODy[1]=0.73738605;
        MODz[1]=0.6680326;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
        // GL Iron Haematoxylin, but this does not seem to work well because it gets confused with the other 2 components
        //			MODx[2]=0.6588232;
        //			MODy[2]=0.66414213;
        //			MODz[2]=0.3533655;
      }
    else if (0 == stainType.compare("Alcian blue-H"))
      {
        // GL Alcian Blue matrix
        MODx[0]= 0.874622;
        MODy[0]= 0.457711;
        MODz[0]= 0.158256;
        // GL Haematox after PAS matrix
        MODx[1]= 0.552556;
        MODy[1]= 0.7544;
        MODz[1]= 0.353744;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("H-PAS"))
      {
        // GL Haem matrix
        MODx[0]= 0.644211; //0.650;
        MODy[0]= 0.716556; //0.704;
        MODz[0]= 0.266844; //0.286;
        // GL PAS matrix
        MODx[1]= 0.175411;
        MODy[1]= 0.972178;
        MODz[1]= 0.154589;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("RGB"))
      {
        //R
        MODx[0]= 0.0;
        MODy[0]= 1.0;
        MODz[0]= 1.0;
        //G
        MODx[1]= 1.0;
        MODy[1]= 0.0;
        MODz[1]= 1.0;
        //B
        MODx[2]= 1.0;
        MODy[2]= 1.0;
        MODz[2]= 0.0;
      }
    else if (0 == stainType.compare("CMY"))
      {
        //C
        MODx[0]= 1.0;
        MODy[0]= 0.0;
        MODz[0]= 0.0;
        //M
        MODx[1]= 0.0;
        MODy[1]= 1.0;
        MODz[1]= 0.0;
        //Y
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 1.0;
      }
    else
      {
        std::cout << "Unknown stain type, use H&E by default.\n"<<std::flush;
        // GL Haem matrix
        MODx[0]= 0.644211; //0.650;
        MODy[0]= 0.716556; //0.704;
        MODz[0]= 0.266844; //0.286;
        // GL Eos matrix
        MODx[1]= 0.092789; //0.072;
        MODy[1]= 0.954111; //0.990;
        MODz[1]= 0.283111; //0.105;
        // Zero matrix
        MODx[2]= 0.0;
        MODy[2]= 0.0;
        MODz[2]= 0.0;
      }



    /*-------------------- Now let's start --------------------*/
    for (i=0; i<3; i++)
      {
        //normalise vector length
        cosx[i]=cosy[i]=cosz[i]=0.0;
        len[i]=sqrt(MODx[i]*MODx[i] + MODy[i]*MODy[i] + MODz[i]*MODz[i]);
        if (len[i] != 0.0)
          {
            cosx[i]= MODx[i]/len[i];
            cosy[i]= MODy[i]/len[i];
            cosz[i]= MODz[i]/len[i];
          }
      }


    // translation matrix
    if (cosx[1]==0.0)
      { //2nd colour is unspecified
        if (cosy[1]==0.0)
          {
            if (cosz[1]==0.0)
              {
                cosx[1]=cosz[0];
                cosy[1]=cosx[0];
                cosz[1]=cosy[0];
              }
          }
      }

    if (cosx[2]==0.0)
      { // 3rd colour is unspecified
        if (cosy[2]==0.0)
          {
            if (cosz[2]==0.0)
              {
                if ((cosx[0]*cosx[0] + cosx[1]*cosx[1])> 1)
                  {
                    cosx[2]=0.0;
                  }
                else
                  {
                    cosx[2]=sqrt(1.0-(cosx[0]*cosx[0])-(cosx[1]*cosx[1]));
                  }

                if ((cosy[0]*cosy[0] + cosy[1]*cosy[1])> 1)
                  {
                    cosy[2]=0.0;
                  }
                else
                  {
                    cosy[2]=sqrt(1.0-(cosy[0]*cosy[0])-(cosy[1]*cosy[1]));
                  }

                if ((cosz[0]*cosz[0] + cosz[1]*cosz[1])> 1)
                  {
                    cosz[2]=0.0;
                  }
                else
                  {
                    cosz[2]=sqrt(1.0-(cosz[0]*cosz[0])-(cosz[1]*cosz[1]));
                  }
              }
          }
      }

    leng=sqrt(cosx[2]*cosx[2] + cosy[2]*cosy[2] + cosz[2]*cosz[2]);

    cosx[2]= cosx[2]/leng;
    cosy[2]= cosy[2]/leng;
    cosz[2]= cosz[2]/leng;

    for (i=0; i<3; i++)
      {
        if (cosx[i] == 0.0) cosx[i] = 0.001;
        if (cosy[i] == 0.0) cosy[i] = 0.001;
        if (cosz[i] == 0.0) cosz[i] = 0.001;
      }

    //matrix inversion
    A = cosy[1] - cosx[1] * cosy[0] / cosx[0];
    V = cosz[1] - cosx[1] * cosz[0] / cosx[0];
    C = cosz[2] - cosy[2] * V/A + cosx[2] * (V/A * cosy[0] / cosx[0] - cosz[0] / cosx[0]);
    q[2] = (-cosx[2] / cosx[0] - cosx[2] / A * cosx[1] / cosx[0] * cosy[0] / cosx[0] + cosy[2] / A * cosx[1] / cosx[0]) / C;
    q[1] = -q[2] * V / A - cosx[1] / (cosx[0] * A);
    q[0] = 1.0 / cosx[0] - q[1] * cosy[0] / cosx[0] - q[2] * cosz[0] / cosx[0];
    q[5] = (-cosy[2] / A + cosx[2] / A * cosy[0] / cosx[0]) / C;
    q[4] = -q[5] * V / A + 1.0 / A;
    q[3] = -q[4] * cosy[0] / cosx[0] - q[5] * cosz[0] / cosx[0];
    q[8] = 1.0 / C;
    q[7] = -q[8] * V / A;
    q[6] = -q[7] * cosy[0] / cosx[0] - q[8] * cosz[0] / cosx[0];

    // initialize 3 output colour stacks
    itkUCharImage2DType::Pointer hematoxylinChannel = itkUCharImage2DType::New();
    hematoxylinChannel->SetRegions(rgbImage->GetLargestPossibleRegion());
    hematoxylinChannel->Allocate();
    hematoxylinChannel->CopyInformation(rgbImage);

    itkUCharImage2DType::PixelType* newpixels = hematoxylinChannel->GetBufferPointer();

    // for (j=0; j<256; j++)
    //   {
    //     rLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosx[outputChannelId]);
    //     gLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosy[outputChannelId]);
    //     bLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosz[outputChannelId]);
    //   }

    // translate ------------------
    int imagesize = width * height;

    itkRGBImage2DType::PixelType* rgbImagePointer = rgbImage->GetBufferPointer();

    for (j=0;j<imagesize;j++)
      {
        // log transform the RGB data
        unsigned char R = rgbImagePointer[j].GetRed();
        unsigned char G = rgbImagePointer[j].GetGreen();
        unsigned char B = rgbImagePointer[j].GetBlue();
        double Rlog = -((255.0*log((static_cast<double>(R)+1)/255.0))/log255);
        double Glog = -((255.0*log((static_cast<double>(G)+1)/255.0))/log255);
        double Blog = -((255.0*log((static_cast<double>(B)+1)/255.0))/log255);

        // rescale to match original paper values
        double Rscaled = Rlog * q[outputChannelId*3];
        double Gscaled = Glog * q[outputChannelId*3+1];
        double Bscaled = Blog * q[outputChannelId*3+2];
        double output = exp(-((Rscaled + Gscaled + Bscaled) - 255.0) * log255 / 255.0);

        if(output>255)
          {
            output=255;
          }

        unsigned char idx = static_cast<unsigned char>(0xff & static_cast<int>(floor(output + .5)));

        ////////////////////////////////////////////////////////////////////////////////
        /// Index to rgb to gray
        ///
        /// The original ImageJ plugin output the colorful H
        /// channel as a indexed image using a LUT. Here I first
        /// use the LUT to convert the index to rgb and then
        /// convert to gray. The formula below is what we have
        /// in matlab rgb2gray
        //newpixels[j] = static_cast<itkUCharImage2DType::PixelType>(0.2989*rLUT[idx] + 0.5870*gLUT[idx] + 0.1140*bLUT[idx]);
        newpixels[j] = idx;
        /// Index to rgb to gray, end
        ////////////////////////////////////////////////////////////////////////////////
      }


    return hematoxylinChannel;
  }

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // itk::ImageIOBase::IOPixelType     pixelType;
  // itk::ImageIOBase::IOComponentType componentType;

  try
    {
      return DoIt( argc, argv );
      // itk::GetImageType(inputVolume, pixelType, componentType);

      // // This filter handles all types on input, but only produces
      // // signed types

      // switch( componentType )
      //   {
      //   case itk::ImageIOBase::RGB:

      //     break;
      //   case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      //   default:
      //     std::cerr << "Unknown input image pixel component type: ";
      //     std::cerr << itk::ImageIOBase::GetComponentTypeAsString( componentType );
      //     std::cerr << "Input image should be 2D RGB image.";
      //     std::cerr << std::endl;
      //     return EXIT_FAILURE;
      //     break;
      //   }
    }

  catch( itk::ExceptionObject & excep )
    {
      std::cerr << argv[0] << ": exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
