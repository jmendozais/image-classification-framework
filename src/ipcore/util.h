/*
 * util.h
 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <iostream>
#include <map>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/ml.h>

namespace ipcore {
	// MACROS
#define PI (2*acos(0))
#define CUANTIZE(var,var_min,var_max,q_factor) ((int)(((var-var_min)*1.0/(var_max-var_min))*q_factor+0.5))
#define NORMALIZE(x,x_min,x_max,y_min,y_max) (((x-x_min)*1.0/(x_max-x_min))*(y_max-y_min) + y_min)
#define EPS (1e-10)
#define two(x) (x)*(x)

	// UTILITIES
double normalize ( double x, double x_min, double x_max, double y_min, double y_max );
cv::Mat normalizeC1( cv::Mat &in );
cv::Mat normalize( cv::Mat &in );
cv::Mat normalize( cv::Mat &in, double min, double max );

/* STATISTICAL FUNCTIONS */
double mean ( const cv::Mat &in );
double mean ( const cv::Mat &in, const cv::Mat &mask );
double standard_deviation ( const cv::Mat &in, double _mean );
double standard_deviation ( const cv::Mat &in, const cv::Mat &mask, double _mean );
double entropy ( const cv::Mat &in );
double entropy ( const cv::Mat &in, const cv::Mat &mask );
double variance ( const cv::Mat &in, double _mean);

double f_test ( const cv::Mat &x, const cv::Mat& h );
cv::Mat covariance ( const cv::Mat &data );
cv::Mat correlation ( const cv::Mat &data );

// Parzen windows method to esticv::Mate P(x) in continous variables
double kernel_density_function ( const cv::Mat& X, double x, int N, double h );

// CRAPPY THINGS -> IMPORTANT: This functions should be cleaned in the future
std::vector<double> Mat1D_64F_to_vector ( const cv::Mat & some );

}
#endif /* UTIL_H_ */
