/*
 * gabor_pca.h
 *
 *  Created on: Jan 26, 2013
 *      Author: jmendoza
 */

#ifndef GABOR_PCA_H_
#define GABOR_PCA_H_

#include "feature_extractor.h"
#include "gabor_wavelet.h"
#include "../new/data.h"
class GaborPCA : public FeatureExtractorInterface {
public:
	GaborPCA( const GaborWavelet& wavelets, ImageSet &images, int rows, int cols, int dim );
	virtual cv::Mat extractFeatures ( cv::Mat &image );
private:
	void performPCA();
	cv::Mat roughGaborFeatureExtraction( const cv::Mat &image );
	GaborWavelet wavelets_;
	cv::PCA pca_;
	ImageSet* images_ptr_;
	int image_rows_;
	int image_cols_;
	int proy_space_dim_;
};

#endif /* GABOR_PCA_H_ */
