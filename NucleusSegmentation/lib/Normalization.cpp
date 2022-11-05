/*
 * Normalization.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: gteodor
 */

#include "Normalization.h"


namespace nscale
{
  /*% Segment foreground from background using discriminant functions
    %inputs:
    %I - M x N x 3 rgb image for color normalization.
    %M - a 4x2 matrix containing the Fisher's linear discriminant function coefficients.
    % DiscriminantF_BG =  M(1,1)*R + M(1,2)*G + M(1,3)*B + M(1,4);
    % DiscriminantF_FG =  M(2,1)*R + M(2,2)*G + M(2,3)*B + M(2,4);
    %ouput:
    %I_seg_fg - segmented binary image*/
  cv::Mat Normalization::segFG(cv::Mat I, cv::Mat M)
  {
    std::vector<cv::Mat> bgr;
    //	cv::imwrite("I.pm", I);
    split(I, bgr);

    // Convert RG image to float
    Size s = bgr[0].size();
    Mat bd(s, CV_32FC1);
    Mat gd(s, CV_32FC1);
    Mat rd(s, CV_32FC1);
    bgr[0].convertTo(bd, bd.type(), 1.0, 0.0);
    bgr[1].convertTo(gd, gd.type(), 1.0, 0.0);
    bgr[2].convertTo(rd, rd.type(), 1.0, 0.0);

    //	discriminant1 = M(1,1)*double(I(:,:,1)) + M(1,2)*double(I(:,:,2)) + M(1,3)*double(I(:,:,3)) + M(1,4);
    cv::Mat discriminant1 = M.at<float>(0,0) * rd + M.at<float>(0,1) * gd + M.at<float>(0,2) * bd + M.at<float>(0,3);

    //	discriminant2 = M(2,1)*double(I(:,:,1)) + M(2,2)*double(I(:,:,2)) + M(2,3)*double(I(:,:,3)) + M(2,4);
    cv::Mat discriminant2 = M.at<float>(1,0) * rd + M.at<float>(1,1) * gd + M.at<float>(1,2) * bd + M.at<float>(1,3);

    cv::Mat I_seg_fg = discriminant1 < discriminant2;
    // true is 255, whereas it is 1 in matlab
    //I_seg_fg /= 255;

    bgr[0].release();
    bgr[1].release();
    bgr[2].release();
    bd.release();
    gd.release();
    rd.release();

    discriminant1.release();
    discriminant2.release();

    return I_seg_fg;
  }

  /*function LAB = rgb2lab(RGB)
    % Convert  RGB to lab(l,alpha,beta) to RGB
    % Reinhard et al. Color Transfer between Images, IEEE Computer Graphics and Application,2001
    %input:
    %RGB - RGB sginals
    %output:
    %LAB - lab(l,alpha,beta)signals*/
  cv::Mat Normalization::bgr2Lab(cv::Mat I)
  {
    cv::Mat LMS(I.size(), CV_32FC3);
    cv::Mat LAB(I.size(), CV_32FC3);
    /*
      % transformation between RGB and LMS cone space
      Matrix1 = [0.3811 0.5783 0.0402;0.1967 0.7244 0.0782;0.0241 0.1288 0.8444];
      LMS = [Matrix1 * RGB']';*/
    for (int i=0; i<I.rows; i++){
      for (int j=0; j<I.cols; j++){
        float b = I.ptr<float>(i)[j*3];
        float g = I.ptr<float>(i)[j*3+1];
        float r = I.ptr<float>(i)[j*3+2];
        LMS.ptr<float>(i)[j*3] = 0.3811*r+0.5783*g+0.0402*b;
        LMS.ptr<float>(i)[j*3+1] = 0.1967*r+0.7244*g+0.0782*b;
        LMS.ptr<float>(i)[j*3+2] = 0.0241*r+0.1288*g+0.8444*b;
      }
    }

    /*% converting the data to logarithmic space
      log_LMS = log10(LMS);*/
    cv::Mat log_LMS(I.size(), CV_32FC3);

    for (int i=0; i<I.rows; i++){
      // get pointer to beginning of each line
      float *LMS_ptr = LMS.ptr<float>(i);
      float *log_LMS_ptr = log_LMS.ptr<float>(i);
      for (int j=0; j<I.cols; j++){
        log_LMS_ptr[j*3] = log10f(LMS_ptr[j*3]);
        log_LMS_ptr[j*3+1] = log10f(LMS_ptr[j*3+1]);
        log_LMS_ptr[j*3+2] = log10f(LMS_ptr[j*3+2]);
        /*			if(i==0  && j < 10){
				std::cout << "log_LMS(0,0): "<< log_LMS_ptr[j*3] <<" (0,1):"<<  log_LMS_ptr[j*3+1] <<" (0,2):"<< log_LMS_ptr[j*3+2] << " LMS(0,2):"<<LMS_ptr[j*3+2] << std::endl;
                                }*/
      }
    }
    /*	% transformation between LMS and lab
	Matrix2 = diag([1/sqrt(3) 1/sqrt(6) 1/sqrt(2)]) * [1 1 1; 1 1 -2; 1 -1 0];*/
    cv::Mat Matrix2 = cv::Mat(3,3,CV_32FC1, 0.0);
    Matrix2.at<float>(0,0)= 1/sqrtf(3.0);
    Matrix2.at<float>(1,1)= 1/sqrtf(6.0);
    Matrix2.at<float>(2,2)= 1/sqrtf(2.0);

    float auxData[9] = {1, 1, 1, 1, 1, -2, 1, -1, 0};
    cv::Mat aux = cv::Mat(3, 3, CV_32FC1, &auxData);
    Matrix2 *=aux;

    // Print Matrix2
    /*	for (int i=0; i<Matrix2.rows; i++){
    // get pointer to beginning of each line
    std::cout << Matrix2.ptr<float>(i)[0] <<", " << Matrix2.ptr<float>(i)[1] <<", "<< Matrix2.ptr<float>(i)[2] << std::endl;

    }*/

    /*	LAB = [Matrix2 * log_LMS']';*/
    // Attribute matrix values to variables to avoid accessing matrix several times (expensive).
    float m11 = Matrix2.at<float>(0,0), m12 = Matrix2.at<float>(0,1), m13 = Matrix2.at<float>(0,2),
      m21 = Matrix2.at<float>(1,0), m22 = Matrix2.at<float>(1,1), m23 = Matrix2.at<float>(1,2),
      m31 = Matrix2.at<float>(2,0), m32 = Matrix2.at<float>(2,1), m33 = Matrix2.at<float>(2,2);
    for (int i=0; i<I.rows; i++){
      // get pointer to beginning of each line
      float *LAB_ptr = LAB.ptr<float>(i);
      float *log_LMS_ptr = log_LMS.ptr<float>(i);
      for (int j=0; j<I.cols; j++){
        float l = log_LMS_ptr[j*3];
        float m = log_LMS_ptr[j*3+1];
        float s = log_LMS_ptr[j*3+2];

        LAB_ptr[j*3] = l * m11 + m * m12 + s * m13;
        LAB_ptr[j*3+1] = l * m21 + m * m22 + s * m23;
        LAB_ptr[j*3+2] = l * m31 + m * m32 + s * m33;
        //			if(i==0  && j == 0){
        //				std::cout << "LAB(0,0): "<< LAB_ptr[j*3] <<" (0,1):"<<  LAB_ptr[j*3+1] <<" (0,2):"<< LAB_ptr[j*3+2] << std::endl;
        //			}
      }
    }

    LMS.release();
    log_LMS.release();
    Matrix2.release();
    aux.release();
    return LAB;
  }

  // Assuming n > 0
  int Normalization::rndint(float n)//round float to the nearest integer
  {
    int ret = (int)floor(n);
    float t;
    t=n-floor(n);
    if (t>=0.5)
      {
        ret = (int)floor(n) + 1;
      }
    return ret;
  }




  /*input: LAB - lab(l,alpha,beta)signals
    output:RGB - RGB signals */
  cv::Mat Normalization::lab2BGR(cv::Mat LAB)
  {
    /*from lab to LMS */
    //Matrix1 =  [1 1 1; 1 1 -1; 1 -2 0] * diag([sqrt(3)/3 sqrt(6)/6 sqrt(2)/2]);
    float auxData[9] = {1, 1, 1, 1, 1, -1, 1, -2, 0};
    cv::Mat Matrix1 = cv::Mat(3, 3, CV_32FC1, &auxData);
    cv::Mat diag = cv::Mat(3,3,CV_32FC1, 0.0);
    diag.at<float>(0,0)= 1/sqrtf(3.0);
    diag.at<float>(1,1)= 1/sqrtf(6.0);
    diag.at<float>(2,2)= 1/sqrtf(2.0);
    Matrix1 *=diag;

    //log_LMS = [Matrix1 * LAB']';
    float m11 = Matrix1.at<float>(0,0), m12 = Matrix1.at<float>(0,1), m13 = Matrix1.at<float>(0,2),
      m21 = Matrix1.at<float>(1,0), m22 = Matrix1.at<float>(1,1), m23 = Matrix1.at<float>(1,2),
      m31 = Matrix1.at<float>(2,0), m32 = Matrix1.at<float>(2,1), m33 = Matrix1.at<float>(2,2);

    //	// Print Matrix1
    //	for (int i=0; i<Matrix1.rows; i++){
    //		// get pointer to beginning of each line
    //		std::cout << Matrix1.ptr<float>(i)[0] <<", " << Matrix1.ptr<float>(i)[1] <<", "<< Matrix1.ptr<float>(i)[2] << std::endl;
    //
    //	}
    cv::Mat log_LMS(LAB.size(), CV_32FC3);

    for (int i=0; i<LAB.rows; i++){
      // get pointer to beginning of each line
      float *LAB_ptr = LAB.ptr<float>(i);
      float *log_LMS_ptr = log_LMS.ptr<float>(i);

      for (int j=0; j<LAB.cols; j++){
        float l = LAB_ptr[j*3];
        float a = LAB_ptr[j*3+1];
        float b = LAB_ptr[j*3+2];

        log_LMS_ptr[j*3]   = l * m11 + a * m12 + b * m13;
        log_LMS_ptr[j*3+1] = l * m21 + a * m22 + b * m23;
        log_LMS_ptr[j*3+2] = l * m31 + a * m32 + b * m33;
        //			if(i==0  && j == 0){
        //				std::cout << "l: "<< l <<" a:"<<a << " b:"<<b<<std::endl;
        //					std::cout << "lab2BGR: log_LMS(0,0): "<< log_LMS_ptr[j*3] <<" (0,1):"<<  log_LMS_ptr[j*3+1] <<" (0,2):"<< log_LMS_ptr[j*3+2] << std::endl;
        //			}
      }
    }
    /*	% conver back from log space to linear space
        LMS = 10.^log_LMS; */
    cv::Mat LMS(LAB.size(), CV_32FC3);
    for (int i=0; i<LMS.rows; i++){
      // get pointer to beginning of each line
      float *LMS_ptr = LMS.ptr<float>(i);
      float *log_LMS_ptr = log_LMS.ptr<float>(i);

      for (int j=0; j<LMS.cols; j++){
        LMS_ptr[j*3] = pow(10.0, log_LMS_ptr[j*3]);
        LMS_ptr[j*3+1] = pow(10.0, log_LMS_ptr[j*3+1]);
        LMS_ptr[j*3+2] = pow(10.0, log_LMS_ptr[j*3+2]);
        //			if(i==0  && j <2 ){
        //				std::cout << "pow: " << pow(10.0, log_LMS_ptr[j*3+2]) << std::endl;
        //				std::cout << "lab2BGR: log_LMS(0,0): "<< log_LMS_ptr[j*3] <<" (0,1):"<<  log_LMS_ptr[j*3+1] <<" (0,2):"<< log_LMS_ptr[j*3+2] << std::endl;
        //				std::cout << "lab2BGR: LMS(0,0): "<< LMS_ptr[j*3] <<" (0,1):"<<  LMS_ptr[j*3+1] <<" (0,2):"<< LMS_ptr[j*3+2] << std::endl;
        //			}
      }
    }

    /*	LMS(find(LMS==-Inf)) = 0;
        LMS(isnan(LMS)) = 0;*/
    for (int i=0; i<LMS.rows; i++){
      // get pointer to beginning of each line
      float *LMS_ptr = LMS.ptr<float>(i);
      for (int j=0; j<LMS.cols; j++){
        for(int c = 0; c < LMS.channels(); c++){
          float x = LMS_ptr[j*3+c];
          if(!(x<= DBL_MAX && x >= -DBL_MAX)){
            LMS_ptr[j*3+c] = 0.0;
          }
        }
      }
    }

    /*% from linear LMS to RGB
      Matrix2 = [4.4687 -3.5887 0.1196;-1.2197 2.3831 -0.1626;0.0585 -0.2611 1.2057];
      RGB = [Matrix2 * LMS']';*/
    // OpenCV stores RGB as BGR, so we have the channels inverted as compared to Matlab
    cv::Mat BGRF(LMS.size(), LMS.type());
    float Matrix2[9] = {4.4687, -3.5887, 0.1196, -1.2197, 2.3831, -0.1626, 0.0585, -0.2611, 1.2057};
    m11 = Matrix2[0], m12 = Matrix2[1], m13 = Matrix2[2],
      m21 = Matrix2[3], m22 = Matrix2[4], m23 = Matrix2[5],
      m31 = Matrix2[6], m32 = Matrix2[7], m33 = Matrix2[8];

    for (int i=0; i<LMS.rows; i++){
      // get pointer to beginning of each line
      float *LMS_ptr = LMS.ptr<float>(i);
      float *BGRF_ptr = BGRF.ptr<float>(i);

      for (int j=0; j<LMS.cols; j++){
        float l = LMS_ptr[3*j];
        float m = LMS_ptr[3*j+1];
        float s = LMS_ptr[3*j+2];

        float b = l*m31 + m*m32 + s *m33;
        float g = l*m21 + m*m22 + s *m23;
        float r = l*m11 + m*m12 + s *m13;
        BGRF_ptr[3*j] = b;
        BGRF_ptr[3*j+1] = g;
        BGRF_ptr[3*j+2] = r;
        //			if(i==0  && j < 2){
        //				std::cout << "lab2BGR: BGR(0,0): "<< BGRF_ptr[j*3] <<" (0,1):"<<  BGRF_ptr[j*3+1] <<" (0,2):"<< BGRF_ptr[j*3+2] << std::endl;
        //			}
      }
    }

    // ColorNormI = reshape(uint8(new_im*255), r,c,d);
    cv::Mat BGR(LMS.size(), CV_8UC3);
    for(int i = 0; i < BGR.rows; i++){
      float* BGRF_ptr = BGRF.ptr<float>(i);
      unsigned char* BGR_ptr = BGR.ptr<unsigned char>(i);
      for(int j = 0; j < BGR.cols; j++){
        for(int c=0; c< BGR.channels();c++){
          // *255
          float bgrfScaled = BGRF_ptr[j*3+c] *255.0;
          // uint8
          if(bgrfScaled < 0.0) bgrfScaled = 0.0;
          if(bgrfScaled > 255.0) bgrfScaled = 255.0;
          BGR_ptr[j*3+c] = (unsigned char)Normalization::rndint(bgrfScaled);
        }
      }
    }

    return BGR;

  }

  /*%fg - RGB sginals(foreground)
    %bg - RGB signals(backgorund)
    %output:
    %idx_fg - indices of foreground pixels
    %idx_bg - indices of background pixels
    %fg_lab - converted lab sginals(foreground)
    %bg_lab - converted lab signals(backgorund)*/
  void nscale::Normalization::PixelClass(cv::Mat I, cv::Mat o_fg, cv::Mat o_bg, cv::Mat& o_fg_lab, cv::Mat& o_bg_lab)
  {

    //	 % foreground : tissue
    //	    idx_fg=find(fg == 1);
    //	    [r,c,d] = size(Img);
    //	    Img1 = reshape(double(Img), r*c, d);
    //	cv::Mat Img1 =
    //	    Img_fg=Img1(idx_fg,:);
    //	    fg_rgb = Img_fg/255;
    //	    fg_lab = rgb2lab(fg_rgb);

    // copy BRG foreground image to fg_bgr, and convert it to float
    Size s = I.size();
    Mat float_bgr(s, CV_32FC1);
    I.convertTo(float_bgr, float_bgr.type(), 1.0, 0.0);
    float_bgr /= 255;
    cv::Mat LAB = nscale::Normalization::bgr2Lab(float_bgr);

    //		float *LAB_ptr = LAB.ptr<float>(0);
    //		int j=0;
    //	std::cout << "return: LAB(0,0): "<< LAB_ptr[j*3] <<" (0,1):"<<  LAB_ptr[j*3+1] <<" (0,2):"<< LAB_ptr[j*3+2] << std::endl;
    LAB.copyTo(o_fg_lab, o_fg);

    /*		cv::Mat labImg;
		cv::cvtColor(I, labImg,CV_BGR2Lab);
		labImg.copyTo(o_fg_lab, o_fg);*/

    //	    % background :  non-tissue
    //	    idx_bg=find(bg == 1);
    //	    Img_bg=Img1(idx_bg,:);
    //	    bg_rgb = Img_bg/255;
    //	    bg_lab = rgb2lab(bg_rgb);

    //		labImg.copyTo(o_bg_lab, o_bg);
    LAB.copyTo(o_bg_lab, o_bg);
  }

  cv::Mat nscale::Normalization::TransferI(cv::Mat fg_lab, cv::Mat fg_mask, float meanT[3], float stdT[3])
  {
    float meanLAB[3], stdLAB[3];
    cv::Mat transferred(fg_lab.size(), CV_32FC3);;
    std::vector<cv::Mat> lab;

    split(fg_lab, lab);

    /*	mask1 = ~isnan(im) & (im~=-Inf);
	mask1 = mask1(:,1).*mask1(:,2).*mask1(:,3);
	valid_im = im(find(mask1>0),:);*/
    cv::Mat multres = lab[0].mul(lab[1]);
    multres = multres.mul(lab[2]);
    cv::Mat mask1 = (multres > 0);
    cv::bitwise_and(mask1, fg_mask, mask1);

    /*	meanI = mean(valid_im);
	stdI= std(valid_im);
    */

    // calculates mean and std for LAB channels
    for(int n = 0; n < 3; n++){
      double sum = 0.0;
      double sq_sum = 0.0;
      int count = 0;
      for(int i = 0; i < lab[n].rows; i++){
        float *labi_ptr = lab[n].ptr<float>(i);
        unsigned char *mask1_ptr = mask1.ptr<unsigned char>(i);
        for(int j =0; j < lab[n].cols; j++){
          if(mask1_ptr[j]){
            sum += labi_ptr[j];
            sq_sum += labi_ptr[j] * labi_ptr[j];
            count++;
          }
        }
      }
      double mean = sum / count;
      double variance = sq_sum / count - mean * mean;
      double std = sqrt(variance);

      //		std::cout <<  "Mean: "<< mean << " std: "<< std<< " sq_sum: "<< sq_sum<<std::endl;
      meanLAB[n] = mean;
      stdLAB[n] = std;
    }
    //	std::cout << "After mean" << std::endl;

    /*
      % scale the data
      [N,d] = size(valid_im);
      new_im = valid_im - repmat(meanI, N, 1);
      new_im = repmat(stdT./stdI, N, 1).*new_im;
      new_im = new_im + repmat(meanT, N, 1); */
    //	std::cout << "transfer: LAB(0,0): "<< lab[0].ptr<float>(0)[0] <<" (0,1):"<<  lab[1].ptr<float>(0)[0] <<" (0,2):"<< lab[2].ptr<float>(0)[0] << std::endl;
    for(int i = 0; i < 3; i++){
      lab[i] -= meanLAB[i];
      lab[i] *= stdT[i]/stdLAB[i];
      lab[i] += meanT[i];

    }

    //	std::cout << "transfer: LAB(0,0): "<< lab[0].ptr<float>(0)[0] <<" (0,1):"<<  lab[1].ptr<float>(0)[0] <<" (0,2):"<< lab[2].ptr<float>(0)[0] << std::endl;


    // merge lab channels into a single multi-channel image
    cv::merge(lab, transferred);
    //std::cout << "transfer after merge: LAB(0,0): "<< transferred.ptr<float>(0)[0] <<" (0,1):"<<  transferred.ptr<float>(0)[1] <<" (0,2):"<< transferred.ptr<float>(0)[2] << std::endl;

    // release temporary data. It would be released automatically, buy why not doing it explicitly?
    lab[0].release();
    lab[1].release();
    lab[2].release();
    multres.release();
    mask1.release();

    return transferred;

  }

  /*Performs color normalization on a single tile
    % Inputs:
    % OriginalI - M x N x 3 rgb original image to be normalized.
    % target_mean - mapping parameter, output of the TargetParameters function
    % target_std - mapping parameter, output of the TargetParameters function
    % M - a 4x2 matrix containing the Fisher's linear discriminant function coefficients.
    % M =[-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];
    % DiscriminantF_BG =  M(1,1)*R + M(1,2)*G + M(1,3)*B + M(1,4);
    % DiscriminantF_FG =  M(2,1)*R + M(2,2)*G + M(2,3)*B + M(2,4);
    % Output:
    % ColorNormI - normalized image */

  cv::Mat Normalization::normalization(const cv::Mat& originalI, float targetMean[3], float targetStd[3])
  {
    // Output normalized image
    cv::Mat ColorNormI;
    //Segmenting foreground from background using discriminant functions
    //[r,c,d] = size(OriginalI);
    int r = originalI.rows;
    int c = originalI.cols;
    int d = originalI.channels();
    assert(d == 3);

    //% M =[-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];
    float mData[8] = {-0.154, 0.035, 0.549, -45.718, -0.057, -0.817, 1.170, -49.887};
    cv::Mat M = cv::Mat(2, 4, CV_32FC1, &mData);

    //	o_fg = SegFG(OriginalI,M);
    cv::Mat o_fg = segFG(originalI, M);
    //cv::imwrite("segOut.ppm", o_fg);

    //	o_bg = imcomplement(o_fg);
    cv::Mat o_bg = nscale::PixelOperations::invert<unsigned char>(o_fg);

    //cv::imwrite("o_bg.ppm", o_bg);

    //Applying Reinhard's method
    //Converting RGB signals to lab(l,alpha,beta)
    //[o_idx_fg,o_idx_bg,o_fg_lab,o_bg_lab] = PixelClass(OriginalI,o_fg,o_bg);
    cv::Mat o_fg_lab = cv::Mat::zeros(c, r, CV_32FC3);
    cv::Mat o_bg_lab = cv::Mat::zeros(c, r, CV_32FC3);

    nscale::Normalization::PixelClass(originalI, o_fg, o_bg, o_fg_lab, o_bg_lab);

    //	cv::imwrite("o_fg_lab.ppm", o_fg_lab);
    //	cv::imwrite("o_bg_lab.ppm", o_bg_lab);


    //	   % Mapping the color distribution of an image to that of the target image
    //	     new_fg = transferI(o_fg_lab, target_mean, target_std);
    //	     new_im=zeros(length(o_idx_fg)+length(o_idx_bg),3);
    //	     new_im(o_idx_fg,:)=new_fg;
    cv::Mat fg_transferred = nscale::Normalization::TransferI(o_fg_lab, o_fg, targetMean, targetStd);

    //	     new_im(o_idx_bg,:)=o_bg_lab;
    o_bg_lab.copyTo(fg_transferred, o_bg);

    //
    //	   % Converting lab to RGB
    //	     new_im = lab2rgb(new_im);
    //	     ColorNormI = reshape(uint8(new_im*255), r,c,d);
    //	     % figure, imshow(ColorNormI)
    //	%       imwrite(ColorNormI,'ColorNormI.tif','Compression','lzw');
    ColorNormI = nscale::Normalization::lab2BGR(fg_transferred);

    return ColorNormI;
  }

  /*function [Mean Std] = TargetParameters(TargetI, M)
    % Calculates mapping parameters for use in color normalization.
    %inputs:
    %TargetI - M x N x 3 rgb target image for color normalization.
    %M - a 4x2 matrix containing the Fisher's linear discriminant function coefficients.
    % DiscriminantF_BG =  M(1,1)*R + M(1,2)*G + M(1,3)*B + M(1,4);
    % DiscriminantF_FG =  M(2,1)*R + M(2,2)*G + M(2,3)*B + M(2,4);
    %ouputs:
    %Mean - scalar mean parameter for mapping.
    %Std - scalar variance parameter for mapping.*/
  void Normalization::targetParameters(const cv::Mat& originalI, float (&targetMean)[3], float (&targetStd)[3])
  {
    for(int i = 0; i < 3; i++)
      {
        targetMean[i] = i;
        targetStd[i] = i;
      }
    int r = originalI.rows;
    int c = originalI.cols;
    int d = originalI.channels();
    assert(d == 3);

    /*% Target image processing
      % Segmenting foreground from background using discriminant functions
      if nargin == 1
      M =[-0.154 0.035 0.549 -45.718; -0.057 -0.817 1.170 -49.887];*/

    float mData[8] = {-0.154, 0.035, 0.549, -45.718, -0.057, -0.817, 1.170, -49.887};
    cv::Mat M = cv::Mat(2, 4, CV_32FC1, &mData);

    /*
      t_fg = SegFG(TargetI,M);
      t_fg = im2bw(t_fg);*/
    cv::Mat o_fg = segFG(originalI, M);
    cv::Mat o_bg = nscale::PixelOperations::invert<unsigned char>(o_fg);

    cv::Mat o_fg_lab = cv::Mat::zeros(c, r, CV_32FC3);
    cv::Mat o_bg_lab = cv::Mat::zeros(c, r, CV_32FC3);

    nscale::Normalization::PixelClass(originalI, o_fg, o_bg, o_fg_lab, o_bg_lab);

    std::vector<cv::Mat> lab;
    split(o_fg_lab, lab);

    /*	mask1 = ~isnan(im) & (im~=-Inf);
	mask1 = mask1(:,1).*mask1(:,2).*mask1(:,3);
	valid_im = im(find(mask1>0),:);*/
    cv::Mat multres = lab[0].mul(lab[1]);
    multres = multres.mul(lab[2]);
    cv::Mat mask1 = (multres > 0);
    cv::bitwise_and(mask1, o_fg, mask1);

    // calculates mean and std for LAB channels
    for(int n = 0; n < 3; n++){
      double sum = 0.0;
      double sq_sum = 0.0;
      int count = 0;
      for(int i = 0; i < lab[n].rows; i++){
        float *labi_ptr = lab[n].ptr<float>(i);
        unsigned char *mask1_ptr = mask1.ptr<unsigned char>(i);
        for(int j =0; j < lab[n].cols; j++){
          if(mask1_ptr[j]){
            sum += labi_ptr[j];
            sq_sum += labi_ptr[j] * labi_ptr[j];
            count++;
          }
        }
      }
      double mean = sum / count;
      double variance = sq_sum / count - mean * mean;
      double std = sqrt(variance);

      targetMean[n] = mean;
      targetStd[n] = std;
    }
    o_fg.release();
    o_bg.release();
    o_fg_lab.release();
    o_bg_lab.release();
    lab[0].release(); lab[1].release();lab[2].release();
    multres.release();
    mask1.release();
  }

}// nscale namespace


