/*
 * util.cpp

 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */
#include "util.h"
namespace ipcore {
	double normalize ( double x, double x_min, double x_max, double y_min, double y_max ) {
		if ( x_min == x_max )
			return y_min;
		else return NORMALIZE( x, x_min, x_max, y_min, y_max );
	}
	cv::Mat normalizeC1( cv::Mat &in ) {
		double rmax = -1e10;
		double rmin = 1e10;
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j ) {
				if ( in.at<double>(i,j) > rmax )
					rmax = in.at<double>(i,j);
				if ( in.at<double>(i,j) < rmin )
					rmin = in.at<double>(i,j);
			}
		cv::Mat ans ( in.rows, in.cols, CV_64FC1 );
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j )
				ans.at<double>(i,j) = NORMALIZE(in.at<double>(i,j), rmin, rmax, 0, 1);
		return ans;
	}
	cv::Mat normalize( cv::Mat &in ) {
		std::vector<cv::Mat> vin;
		cv::split ( in, vin );
		std::vector<cv::Mat> resp ( vin.size() );
		for ( int i = 0; i < vin.size(); ++ i )
			resp[i] = normalizeC1 ( vin[i] );
		cv::Mat ans;
		cv::merge ( resp, ans );
		return ans;
	}
	cv::Mat normalizeC1( cv::Mat &in, double min_, double max_) {
		double rmax = -1e10;
		double rmin = 1e10;
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j ) {
				if ( in.at<double>(i,j) > rmax )
					rmax = in.at<double>(i,j);
				if ( in.at<double>(i,j) < rmin )
					rmin = in.at<double>(i,j);
			}
		cv::Mat ans ( in.rows, in.cols, CV_64FC1 );
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j )
				ans.at<double>(i,j) = NORMALIZE(in.at<double>(i,j), rmin, rmax, min_, max_);
		return ans;
	}
	cv::Mat normalize( cv::Mat &in, double min, double max ) {
		std::vector<cv::Mat> vin;
		cv::split ( in, vin );
		std::vector<cv::Mat> resp ( vin.size() );
		for ( int i = 0; i < vin.size(); ++ i )
			resp[i] = normalizeC1 ( vin[i], min, max );
		cv::Mat ans;
		cv::merge ( resp, ans );
		return ans;
	}
	double mean ( const cv::Mat &in ) {
		double sum = 0;
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j )
				sum += in.ptr<double>(i)[j];
		return sum/(in.rows * in.cols );
	}
	double mean ( const cv::Mat &in, const cv::Mat &mask ) {
		assert ( mask.type() == cv::Mat(1,1,CV_8UC1).type() );
		assert ( in.type() == cv::Mat(1,1,CV_64FC1).type() );
		assert ( in.rows == mask.rows && in.cols == mask.cols );

		double sum = 0;
		int counter = 0;
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j )
				if ( mask.ptr<uchar>(i)[j] != 0 ) {
					++ counter;
					sum += in.ptr<double>(i)[j];
				}
		return sum/counter;
	}
	double standard_deviation ( const cv::Mat &in, double _mean = 0 ) {
		if ( _mean == 0 )
			_mean = mean ( in );
		double sq_diffs = 0, temp;
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j ) {
				temp = in.ptr<double>(i)[j] - _mean;
				sq_diffs += temp * temp;
			}
		return sqrt ( sq_diffs / ( in.rows * in.cols - 1 ) );
	}
	double standard_deviation ( const cv::Mat &in, const cv::Mat &mask, double _mean = 0 ) {

		assert ( mask.type() == cv::Mat(1,1,CV_8UC1).type() );
		assert ( in.type() == cv::Mat(1,1,CV_64FC1).type() );
		assert ( in.rows == mask.rows && in.cols == mask.cols );

		if ( _mean == 0 )
			_mean = mean ( in );
		double sq_diffs = 0, temp;
		int counter = 0;
		for ( int i = 0; i < in.rows; ++ i )
			for ( int j = 0; j < in.cols; ++ j ) {
				if ( mask.ptr<uchar>(i)[j] != 0 ) {
					++ counter;
					temp = in.ptr<double>(i)[j] - _mean;
					sq_diffs += temp * temp;
				}
			}
		return sqrt ( sq_diffs / ( counter - 1 ) );
	}
	double entropy ( const cv::Mat &in ) {
			double max_, min_;
			max_ = -1e100;
			min_ = -max_;

			for ( int i = 0; i < in.rows; ++ i )
				for ( int j = 0; j < in.cols; ++ j ) {
					double val = in.ptr<double>(i)[j];
					max_ = std::max ( max_, val);
					min_ = std::min ( min_, val );
				}
			double p[256];
			for ( int i = 0; i < 256; ++ i ) p[i] = 0;

			for ( int i = 0; i < in.rows; ++ i )
				for ( int j = 0; j < in.cols; ++ j ) {
					double x = in.ptr<double>(i)[j];
					int x256 = ipcore::normalize( x, min_, max_, 0, 255.0 );
					p[x256] += 1;
				}
			double sum = 0;
			for ( int i = 0; i < 256; ++ i ) {
				p[i] /= (in.rows*in.cols);
				sum += p[i];
			}
			double E = 0;
			for ( int i = 0; i < 256; ++ i ) {
				if ( p[i] == 0 || p[i] == 1 ) continue;
					E += p[i] * (log(p[i])/log(2) );
			}
			return -E;
	}

	double entropy ( const cv::Mat &in, const cv::Mat &mask ) {
			double max_, min_;
			max_ = -1e100;
			min_ = -max_;

			for ( int i = 0; i < in.rows; ++ i )
				for ( int j = 0; j < in.cols; ++ j ) {
					double val = in.ptr<double>(i)[j];
					max_ = std::max ( max_, val);
					min_ = std::min ( min_, val );
				}

			double p[256];
			for ( int i = 0; i < 256; ++ i )
				p[i] = 0;
			int counter = 0;
			for ( int i = 0; i < in.rows; ++ i )
				for ( int j = 0; j < in.cols; ++ j ) {
					if ( mask.ptr<uchar>(i)[j] == 0 ) continue;
					double x = in.ptr<double>(i)[j];
					int x256 = (int) NORMALIZE(x,min_,max_,0,255);
					p[x256] += 1;
					++counter;
				}
			for ( int i = 0; i < 256; ++ i )
				p[i] /= counter;

			double E = 0;
			for ( int i = 0; i < 256; ++ i ) {
				if ( p[i] == 0 || p[i] == 1 ) continue;
				E += p[i] * (log(p[i])/log(2) );
			}
			return -E;
	}

	double variance ( const cv::Mat &in, double _mean = 0 ) {
		double s = standard_deviation ( in );
		return s * s;
	}
	double f_test ( const cv::Mat &x, const cv::Mat& h ) {
		double xm = mean ( x );
		std::map < double, int > cc;
		for ( int i = 0; i < h.rows; ++ i )
			cc[h.ptr<double>(i)[0]] ++;
		int K = cc.size();
		int N = x.rows;
		double num = 0;
		double den = 0;
		for ( std::map < double, int >::iterator it = cc.begin(); it != cc.end(); ++ it ) {
			std::pair < double, int > c = (*it);
			double xk = 0;
			int nk = c.second;
			double var_k = 0;
			for ( int i = 0; i < h.rows; ++ i )
				if ( h.ptr<double>(i)[0] == c.first ) {
					xk += x.ptr<double>(i)[0];
				}
			// class mean
			xk /= nk;
			// class variance
			for ( int i = 0; i < h.rows; ++ i )
				if ( h.ptr<double>(i)[0] == c.first ) {
					double diff =  ( x.ptr<double>(i)[0] - xk );
					var_k += diff * diff;
				}
			// variance of class k
			var_k /= ( nk - 1 );

			//num += nk * ( xk - xm );
			num += nk * abs( xk - xm );

			den += nk * var_k;
		}
		num /= (K-1);
		den /= (N-K);
		return num/den;
	}
	cv::Mat covariance ( const cv::Mat &data ) {
		int n = data.rows;
		int m = data.cols;
		// calc means
		cv::Mat means ( m, 1, CV_64F );
		for ( int i = 0; i < m; ++ i ) {
			means.ptr<double>(i)[0] = mean ( data ( cv::Range::all(), cv::Range(i,i+1) ) );
		}
		cv::Mat cov ( m, m, CV_64F );
		for ( int i = 0; i < m; ++ i )
			for ( int j = 0; j < m; ++ j ) {
				cov.ptr<double>(i)[j] = 0;
				for ( int k = 0; k < n; ++ k ) {
					if ( i != j )
					cov.ptr<double>(i)[j] += ( ( data.ptr<double>(k)[i] - means.ptr<double>(i)[0] ) *
																		 ( data.ptr<double>(k)[j] - means.ptr<double>(j)[0] ) );
				}
				cov.ptr<double>(i)[j] /= (n-1);
			}
		return cov;
	}
	cv::Mat correlation ( const cv::Mat &data ) {
		int n = data.rows;
		int m = data.cols;
		// calc means
		cv::Mat means ( m, 1, CV_64F );
		for ( int i = 0; i < m; ++ i ) {
			means.ptr<double>(i)[0] = mean ( data ( cv::Range::all(), cv::Range(i,i+1) ) );
		}
		cv::Mat corr ( m, m, CV_64F );
		for ( int i = 0; i < m; ++ i )
			for ( int j = 0; j < m; ++ j ) {
				corr.ptr<double>(i)[j] = 0;
				for ( int k = 0; k < n; ++ k ) {
					if ( i != j )
					corr.ptr<double>(i)[j] += ( ( data.ptr<double>(k)[i] - means.ptr<double>(i)[0] ) *
																		 ( data.ptr<double>(k)[j] - means.ptr<double>(j)[0] ) );
				}
				corr.ptr<double>(i)[j] /= (n-1);
				corr.ptr<double>(i)[j] /= ( means.ptr<double>(i)[0] * means.ptr<double>(j)[0] );
			}
		return corr;
	}


	// Parzen windows method to estimate P(x) in continous variables
	double kernel_density_function ( const cv::Mat& X, double x, int N, double h ) {
		return -1;
	}

	// CRAPPY THINGS -> IMPORTANT: This functions should be cleaned in the future
	std::vector<double> Mat1D_64F_to_vector ( const cv::Mat & some ) {
		std::vector< double > ans ( some.rows );
		for ( int i = 0; i < ans.size(); ++ i )
			ans[i] = *some.ptr<double>(i);
		return ans;
	}
}



