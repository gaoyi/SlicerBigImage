/*
 * ScanlineOperations.cpp
 *
 *  Created on: Aug 2, 2011
 *      Author: tcpan
 */

#include "PixelOperations.h"
#include <limits>
#include "Logger.h"
#include "TypeUtils.h"

namespace nscale {

using namespace cv;

template <typename T>
Mat PixelOperations::invert(const Mat& img) {
	// write the raw image
	CV_Assert(img.channels() == 1);

	if (std::numeric_limits<T>::is_integer) {

		if (std::numeric_limits<T>::is_signed) {
			Mat output;
			bitwise_not(img, output);
			return output + 1;
		} else {
			// unsigned int
			return std::numeric_limits<T>::max() - img;
		}

	} else {
		// floating point type
		return -img;
	}


}

template <typename T>
Mat PixelOperations::mod(Mat& img, T mod) {
	// write the raw image
	CV_Assert(img.channels() == 1);
	CV_Assert(std::numeric_limits<T>::is_integer);

	Mat result(img.size(), img.type());
	T *ptr, *res;
	for(int y=0; y< img.rows; y++){
		ptr = img.ptr<T>(y);
		res = result.ptr<T>(y);
		for (int x = 0; x < img.cols; ++x) {
			res[x] = ptr[x] % mod;
		}
	}

	return result;
}

Mat PixelOperations::bgr2gray(const ::cv::Mat& img){
	int imageChannels = img.channels();
	assert(imageChannels == 3);
	assert(img.type() == CV_8UC3);

	Mat gray = Mat(img.size(), CV_8UC1);

	// Same constants as used by Matlab
	double r_const = 0.298936021293776;
	double g_const = 0.587043074451121;
	double b_const = 0.114020904255103;

	int nr = img.rows, nc = img.cols;

	for(int i=0; i<nr; i++){
		const unsigned char* data_in = img.ptr<unsigned char>(i);
		unsigned char* data_out = gray.ptr<unsigned char>(i);
		for(int j=0; j<nc; j++){
			unsigned char b = data_in[j * imageChannels];
			unsigned char g = data_in[j * imageChannels + 1];
			unsigned char r = data_in[j * imageChannels + 2];
			double grayPixelValue = r_const * (double)r + g_const * (double)g + b_const * (double)b;
			data_out[j] = cci::common::type::double2uchar(grayPixelValue);
		}
	}
	return gray;
}

Mat PixelOperations::ComputeInverseStainMatrix(const Mat& M, const Mat& b) {
	assert(M.type() == CV_64FC1);
	assert(b.rows == 1);  // should be a row vector

	long t1 = cci::common::event::timestampInUS();
	//initialize normalized stain deconvolution matrix
	Mat normal_M(M);

	//stain normalization
	double col_Norm;
	for(int i=0; i< M.cols; i++){
		col_Norm = norm(M.col(i));
		if(  col_Norm > (double) 0 ){
			normal_M.col(i) = M.col(i) / col_Norm;
		}
	}
	//showMat(normal_M, "normal_M--stain normalization: ");

	//find last column of the normalized stain deconvolution matrix
	int last_col_index = M.cols-1;
//	Mat last_col = normal_M.col(last_col_index); //or//Mat last_col = normal_M(Range::all(), Range(2,3));
	//showMat( last_col, "last column " );

	//normalize the stain deconvolution matrix again
	if(norm(normal_M.col(last_col_index)) == (double) 0){
		double nrm;

		for(int i=0; i< normal_M.rows; i++){
			nrm = norm(normal_M.row(i));
			if( nrm > 1 ){
				normal_M.at<double>(i,last_col_index) = 0;
			}
			else{
				normal_M.at<double>(i,last_col_index) =  sqrt( 1 - nrm * nrm );
			}
		}
		normal_M.col(last_col_index) = normal_M.col(last_col_index) / norm(normal_M.col(last_col_index));
	}
	//showMat(normal_M, "normal_M");

	//take the inverse of the normalized stain deconvolution matrix
	Mat Q = normal_M.inv();
	normal_M.release();

	//select rows in Q with a true marker in b
	Mat T(1,Q.cols,Q.type());
	for (int i=0; i<b.cols; i++){
		if( b.at<char>(0,i) == 1 )
			T.push_back(Q.row(i));
	}
	Q.release();
	Mat T2;
	T.rowRange(Range(1,T.rows)).copyTo(T2);
	T.release();
	Mat output;
	T2.convertTo(output, CV_32FC1);
	T2.release();

	long t2 = cci::common::event::timestampInUS();

	//cout << "  Before normalized = "<< t2-t1 <<endl;

	return output;
}

vector<float> PixelOperations::ComputeLookupTable() {
	vector<float> precomp_res;
	float temp;
	const float k = -255.0 / log(255.0);
	for(int i=0; i < 256; i++){
		//temp = -(255.0* log(((double)i +1.0)/255.0))/log(255.0);
		temp = k * log((float)i + 1.0) + 255.0;
		precomp_res.push_back(temp);
	}
	return precomp_res;
}

void PixelOperations::ColorDeconv( const Mat& image, const Mat& Q, const vector<float> &lut, Mat& H, Mat& E, bool BGR2RGB){
	assert(image.channels() == 3); // Image must have 3 channels for color deconvolution;
	assert(Q.isContinuous());  // Q should be continuous.  use copyTo to ensure
	assert(H.rows == image.rows);
	assert(E.rows == image.rows);
	assert(H.cols == image.cols);
	assert(E.cols == image.cols);

	//normalized image
	long t2 = cci::common::event::timestampInUS();

	int nr = image.rows, nc = image.cols, step1 = image.step1(); // step1: number of channel elements in a row.
	//Mat dn = Mat::zeros(nr, nc, CV_64FC3);
	Mat dn = Mat::zeros(nr, nc, CV_32FC3);
	LUT(image, lut, dn);
//	for(int i=0; i<nr; i++){
//		const unsigned char* data_in = image.ptr<unsigned char>(i);
//		//double* data_out = dn.ptr<double>(i);
//		float* data_out = dn.ptr<float>(i);
//		for(int j=0; j<step1; j++){
//			data_out[j] = (float)lut[data_in[j]];
//		}
//	}

	long t1loop = cci::common::event::timestampInUS();
	//cout << "  After first loop = "<< t1loop - t2 <<endl;

	//channel deconvolution
	Mat dn2 = dn.reshape(1, nr*nc);  // col = 3, rows = numpixels;
	Mat Q4;
	if (BGR2RGB) {
		flip(Q.t(), Q4, 0); // flip vertically of the transpose
	} else {
		Q4 = Q.t();
	}
	Mat cn = dn2 * Q4;  // matrix multiply.  result = col=2, rows=numpixels, channels = 1;
	Q4.release();
	dn.release();

////	Mat cn = Mat::zeros(nr, nc, CV_64FC2);
//	Mat cn = Mat::zeros(nr, nc, CV_32FC2);
//	int dn_channels = dn.channels();
//	int cn_channels = cn.channels();

//	const double *Q_ptr = Q.ptr<double>(0);
//
//	for(int i=0; i<nr; i++){
////		const double *dn_ptr = dn.ptr<double>(i);
////		double *cn_ptr = cn.ptr<double>(i);
//		const float *dn_ptr = dn.ptr<float>(i);
//		float *cn_ptr = cn.ptr<float>(i);
//
//		for(int j=0; j<nc; j++){
//			for(int k=0; k<dn_channels; k++){
//				for(int Q_i=0; Q_i<Q.rows; Q_i++)
//					if( BGR2RGB ){
//						cn_ptr[j * cn_channels + Q_i] += (float)(Q_ptr[Q_i * Q.cols + k]  * (double)dn_ptr[ j * dn_channels + dn_channels-1-k]);
//					}else{
//						cn_ptr[j * cn_channels + Q_i] += (float)(Q_ptr[Q_i * Q.cols + k]  * (double)dn_ptr[ j * dn_channels + k]);
//					}
//			}
//		}
//	}
//	dn.release();

	long t2loop = cci::common::event::timestampInUS();
	//cout << "  After 2 loop = "<< t2loop - t1loop <<endl;

	//denormalized H and E channels
	float log255div255 = -log(255.0)/255.0;
	Mat deconved(cn.size(), cn.type());
	exp((cn - (float)255.0) * log255div255, deconved);
	cn.release();
	Mat HE(nr, nc, CV_8UC2);
	deconved.reshape(2, nr).convertTo(HE, CV_8UC2);  // convert back to a 2D multi channel image
	deconved.release();
	vector<Mat> output;
	split(HE, output);  // split back to H and E.
	HE.release();

	output[0].copyTo(H);
	output[1].copyTo(E);
	output[0].release();
	output[1].release();
//	double temp;
//	for(int i=0; i<nr; i++){
//		unsigned char *E_ptr = E.ptr<unsigned char>(i);
//		unsigned char *H_ptr = H.ptr<unsigned char>(i);
//		//const double *cn_ptr = cn.ptr<double>(i);
//		const float *cn_ptr = cn.ptr<float>(i);
//
//		for(int j=0; j<nc; j++){
//			temp = exp(-((double)cn_ptr[j * cn_channels]-255.0)*log255div255);
//			H_ptr[j] = cci::common::type::double2uchar(temp);
//
//			temp = exp(-((double)cn_ptr[j * cn_channels + 1]-255.0)*log255div255);
//
//			E_ptr[j] = cci::common::type::double2uchar(temp);
//		}
//	}

	long t3 = cci::common::event::timestampInUS();
	//cout << "  Rest = "<< t3-t2loop<<endl;
}

template <typename T>
Mat PixelOperations::replace(const Mat& img, T oldval, T newval) {
	// write the raw image
	CV_Assert(img.channels() == 1);

	Mat result(img.size(), img.type());
	const T *ptr;
	T *res;
	for(int y=0; y< img.rows; y++){
		ptr = img.ptr<T>(y);
		res = result.ptr<T>(y);
		for (int x = 0; x < img.cols; ++x) {
			res[x] = (ptr[x] == oldval ? newval : ptr[x]);
		}
	}

	return result;
}




template Mat PixelOperations::invert<unsigned char>(const Mat&);
template Mat PixelOperations::invert<float>(const Mat&);
template Mat PixelOperations::invert<int>(const Mat&);  // for imfillholes


template Mat PixelOperations::mod<unsigned char>(Mat&, unsigned char mod);
template Mat PixelOperations::mod<int>(Mat&, int mod);

template Mat PixelOperations::replace<unsigned char>(const Mat&, unsigned char oldval, unsigned char newval);
template Mat PixelOperations::replace<int>(const Mat&, int oldval, int newval);

}


