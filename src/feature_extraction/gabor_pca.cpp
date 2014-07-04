/*
 * gabor_pca.cpp
 *
 *  Created on: Jan 26, 2013
 *      Author: jmendoza
 */

#include "gabor_pca.h"

GaborPCA::GaborPCA( const GaborWavelet& wavelets, ImageSet &images, int image_rows, int image_cols, int dim ) {
	wavelets_ = wavelets;
	image_rows_ = image_rows;
	image_cols_ = image_cols;
	images_ptr_ = &images;
	proy_space_dim_ = dim;
	performPCA();
}
cv::Mat GaborPCA::roughGaborFeatureExtraction ( const cv::Mat &image ) {
	cv::Mat temp;
	cvtColor ( image, temp, CV_BGR2HSV );
	std::vector<cv::Mat> images = wavelets_.getMagnitudes(temp);
	cv::Size size ( image_rows_, image_cols_ );
	cv::Mat fv ( image_rows_ * image_cols_ * images.size(), 1, CV_64FC1 );
	for ( int i = 0; i < images.size(); ++ i ) {
		images[i] = cv::abs ( images[i] );
		cv::resize( images[i], images[i], size );
	}
	//grid(images, 3, 15);
	for ( int i = 0; i < images.size(); ++ i ) {
		for ( int j = 0; j < image_rows_; ++ j )
			for ( int k = 0; k < image_cols_; ++ k )
				fv.ptr<double>(i*image_rows_*image_cols_ + j*image_cols_ + k)[0] = images[i].ptr<double>(j)[k];
	}
	return fv;
}
void GaborPCA::performPCA() {
	std::cout << "perform pca" << std::endl;
	std::vector<cv::Mat> features;
	int len = images_ptr_->size();
	features.resize(len);
	for ( int i = 0; i < len; ++ i ) {
		features[i] = roughGaborFeatureExtraction( images_ptr_->getImage(i) );
	}
	cv::Mat data ( len, features[0].rows , CV_64FC1 );
	for ( int i = 0; i < len; ++ i )
		for ( int j = 0; j < features[0].rows; ++ j )
			data.ptr<double>(i)[j] = features[i].ptr<double>(j)[0];
	cv::transpose(data, data);
	pca_( data, cv::Mat(), CV_PCA_DATA_AS_COL, proy_space_dim_ );
	std::cout << "end pca" << std::endl;
}

cv::Mat GaborPCA::extractFeatures ( cv::Mat &image ) {
	cv::Mat features = roughGaborFeatureExtraction(image);
	cv::Mat ans = pca_.project(features);
	//std::cout << "Size " << ans.rows << " " << ans.cols << std::endl;
	return ans;
}
