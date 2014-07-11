/*
 * gabor_wavelet.h
 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */

#ifndef GABOR_WAVELET_H_
#define GABOR_WAVELET_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <boost/shared_ptr.hpp>

#include "../ipcore/util.h"
#include "../new/data.h"
#include "feature_extractor.h"
using namespace std;
using namespace cv;
//namespace ipcore {

struct FilterBankParams : public AbstractParams {
	int num_orientations;
	int num_scales;
	double low_frequency_bound;
	double high_frequency_bound;
	int kernel_side;
	bool remove_dc;
};
struct GaborWaveletParams : public FilterBankParams {
	// Gabor filter bank params
	int metrics;
	int response;
	// Gabor feature extraction params
	enum MaskType { NO_MASK, BACKGROUND_REMOVAL_MASK, MID_REGION_MASK, LARGEST_REGION_MASK, CUSTOM_MASK };
	enum Metrics { STD = 1, MEAN = 2, ENTROPY = 4 };
	enum Responses { REAL = 1, IMAGINARY = 2, MAGNITUDE = 4 };
	int mask_type;
	GaborWaveletParams() {
		num_orientations = 5;
		num_scales = 3;
		low_frequency_bound = 0.1;
		high_frequency_bound = 0.4;
		kernel_side = 30;
		remove_dc = false;
		mask_type = NO_MASK;
		metrics = STD | MEAN;
		response = MAGNITUDE;
	}
};

class FilterBankBuilderInterface {
public:
	virtual std::vector<cv::Mat> getFilterBank ( FilterBankParams params_ ) = 0;
};

typedef boost::shared_ptr<FilterBankBuilderInterface> FilterBankBuilderPtr;

class ManjunathFilterBank : public FilterBankBuilderInterface {
	FilterBankParams params_;
public:
	ManjunathFilterBank ( FilterBankParams params_ );
	void gabor_filter ( cv::Mat &real, cv::Mat &image, int scale, int orientation );
	virtual std::vector<cv::Mat> getFilterBank ( FilterBankParams );
};

class LogGaborFilterBank : public FilterBankBuilderInterface {
	virtual std::vector<cv::Mat> getFilterBank ( FilterBankParams );
};

class GaborWavelet : public FeatureExtractorInterface {
	private:
	//Gabor filter bank params
	std::vector<cv::Mat> filter_bank;
	//Gabor wavelet feature extractor params
	cv::Mat mask_;
	int mask_type;
	GaborWaveletParams params_;
	boost::shared_ptr<FilterBankBuilderInterface> filter_bank_builder_ptr_;

	//private methods
	void update();
	void gabor_filter ( cv::Mat&, cv::Mat&, int, int);
	cv::Mat extract_features_one_channel ( const cv::Mat & im );

	public:
	GaborWavelet ( int sc, int ors, double u_l, double u_h, int side, bool rdc );
	GaborWavelet ( GaborWaveletParams params );
	GaborWavelet ();
	bool persist ( std::string filename );
	bool load ( std::string filename );
	// IFeatureExtractor Methods
	virtual cv::Mat extractFeatures ( cv::Mat &image );
	std::vector<cv::Mat> convolveOneChannel ( const cv::Mat &image ) const;
	std::vector<cv::Mat> convolve ( const cv::Mat &image ) const;
	std::vector<cv::Mat> getMagnitudes( const std::vector<cv::Mat> &real_imag ) const;
	std::vector<cv::Mat> getMagnitudes( const cv::Mat &image ) const;
	std::vector<cv::Mat> getResponses ( const cv::Mat &image ) const;
	cv::Mat measure ( std::vector<cv::Mat> inputs );
	static void generate_mask ( const cv::Mat& input, cv::Mat& mask, int mask_type );
	void set_mask( const cv::Mat &mask );
	const cv::Mat& getFrequencyProfile() const;
	const int getNumScales() const;
	const int getNumOrientations() const;
	const double getHighFrequencyBound() const;
	const double getLowFrequencyBound() const;
	const int getNumMetrics() const;
	double setHighFrequencyBound(double value);
	double setLowFrequencyBound(double value);
	void setFilterBank(std::vector<cv::Mat>);
	std::vector<cv::Mat> getFilterBank();
	void setMaskType( int type );
	void setMask( const cv::Mat &mask );
	void setFilterBankBuilderPtr( FilterBankBuilderPtr );
	FilterBankBuilderPtr getFilterBankBuilderPtr();
};

class PGabor : public GaborWavelet {
private:
	int deep;
public:
	PGabor ( GaborWaveletParams params );
	// override
	virtual cv::Mat extractFeatures ( cv::Mat &image );
};

#endif /* GABOR_WAVELET_H_ */
