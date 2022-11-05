#include <string>

// itk
#include "itkImage.h"

// local
#include "itkTypedefs.h"

#include "ColorDecomposition.h"

namespace gth818n
{
  ////////////////////////////////////////////////////////////////////////////////
  /// This is an implementation of the paper:
  ///
  /// Anal Quant Cytol Histol. 2001 Aug;23(4):291-9. Quantification of
  /// histochemical staining by color deconvolution. Ruifrok AC,
  /// Johnston DA.
  ///
  /// The implementation is adopted from the java ImageJ plugin at
  /// http://www.mecourse.com/landinig/software/cdeconv/cdeconv.html
  ///
  /// The ImageJ plugin takes an input image, then a pop-up box is
  /// given to ask user to determine the type of stain. It has many
  /// different choices. It output several channels, such as
  /// Hematoxylin, Eosin, etc., depending on the input stains.
  ///
  /// In this implementation, however, we only take the H&E stain as
  /// input type, and only output the Hematoxylin channel as output.
  itkUCharImage2DType::Pointer ExtractHematoxylinChannel(itkRGBImage2DType::Pointer HAndEImage)
  {
    int width = HAndEImage->GetLargestPossibleRegion().GetSize()[0];
    int height = HAndEImage->GetLargestPossibleRegion().GetSize()[1];

    std::string stainType = "H&E";

    double leng, A, V, C, log255=log(255.0);
    int i,j;
    int channelIndex = 0; ///< 0: Hematoxylin; 1: Eosin
    double MODx[3];
    double MODy[3];
    double MODz[3];
    double cosx[3];
    double cosy[3];
    double cosz[3];
    double len[3];
    double q[9];
    unsigned char rLUT[256];
    unsigned char gLUT[256];
    unsigned char bLUT[256];

    //if (!stainType.compare("H&E"))
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


    // start
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
    hematoxylinChannel->SetRegions(HAndEImage->GetLargestPossibleRegion());
    hematoxylinChannel->Allocate();
    hematoxylinChannel->CopyInformation(HAndEImage);

    itkUCharImage2DType::PixelType* newpixels = hematoxylinChannel->GetBufferPointer();

    for (j=0; j<256; j++)
      {
        rLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosx[channelIndex]);
        gLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosy[channelIndex]);
        bLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosz[channelIndex]);
      }

    // translate ------------------
    int imagesize = width * height;

    itkRGBImage2DType::PixelType* rgbImagePointer = HAndEImage->GetBufferPointer();

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
        double Rscaled = Rlog * q[channelIndex*3];
        double Gscaled = Glog * q[channelIndex*3+1];
        double Bscaled = Blog * q[channelIndex*3+2];
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


  itkUCharImage2DType::Pointer ExtractEosinChannel(itkRGBImage2DType::Pointer HAndEImage)
  {
    int width = HAndEImage->GetLargestPossibleRegion().GetSize()[0];
    int height = HAndEImage->GetLargestPossibleRegion().GetSize()[1];

    std::string stainType = "H&E";

    double leng, A, V, C, log255=log(255.0);
    int i,j;
    int channelIndex = 1; ///< 0: Hematoxylin; 1: Eosin
    double MODx[3];
    double MODy[3];
    double MODz[3];
    double cosx[3];
    double cosy[3];
    double cosz[3];
    double len[3];
    double q[9];
    unsigned char rLUT[256];
    unsigned char gLUT[256];
    unsigned char bLUT[256];

    //if (!stainType.compare("H&E"))
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


    // start
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
    hematoxylinChannel->SetRegions(HAndEImage->GetLargestPossibleRegion());
    hematoxylinChannel->Allocate();
    hematoxylinChannel->CopyInformation(HAndEImage);

    itkUCharImage2DType::PixelType* newpixels = hematoxylinChannel->GetBufferPointer();

    for (j=0; j<256; j++)
      {
        rLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosx[channelIndex]);
        gLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosy[channelIndex]);
        bLUT[255-j]=static_cast<unsigned char>(255.0 - (double)j * cosz[channelIndex]);
      }

    // translate ------------------
    int imagesize = width * height;

    itkRGBImage2DType::PixelType* rgbImagePointer = HAndEImage->GetBufferPointer();

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
        double Rscaled = Rlog * q[channelIndex*3];
        double Gscaled = Glog * q[channelIndex*3+1];
        double Bscaled = Blog * q[channelIndex*3+2];
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

}// namespace gth818n
