/*
 * cdh.h
 *
 *  Created on: Jan 20, 2013
 *      Author: jmendoza
 */

#ifndef CDH_H_
#define CDH_H_

#include <vector>
#include <opencv/cv.h>
#include "../ipcore/util.h"
#include "feature_extractor.h"
using namespace ipcore;
using namespace std;
using namespace cv;


class CDHFeatureExtractor : public FeatureExtractorInterface {
	private:
	static double lab_coeff[3][3];
	static double illuminant_coeff[3];
	static double channel_ranges[3][2];
	double fl ( double x );
	double f ( double x );
	double fab ( double x );
	Mat rgb2lab ( const Mat &image );
	vector < vector < double > > quantize_lab_colors ( const Mat& lab_im, int n_l, int n_a, int n_b );
	vector < vector < double > > get_edge_orientations ( const Mat& edge_im );
	vector<double> build_cdh ( const Mat& lab_im, const vector < vector < double > > &ori );
	inline int quantize_orientation ( double theta );
	inline int quantize_color ( double l, double a, double b );

	public:
	int q_channels[3]; // quantizations factors for channels L*a*b
	int q_phi; // quantization factor for edge orientation angle
	int distance; // neighborhood distance
	CDHFeatureExtractor ( );
	virtual cv::Mat extractFeatures ( cv::Mat &image );
};

#endif /* CDH_H_ */
