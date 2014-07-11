/*
 * cdh.cpp
 *
 *  Created on: Jan 20, 2013
 *      Author: jmendoza
 */

#include "cdh.h"
#include <fstream>

double sobel2d_kernel[2][3][3] = { { { -1 -2 -1 }, { 0, 0, 0 }, { 1, 2, 1 } }, { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } } };
	Mat sobel2d ( Mat image ) {
		Mat result = image.clone();
		for ( int i = 0; i + 2 < image.rows; ++ i )
			for ( int j = 0; j + 2 < image.cols; ++ j ) {
				Vec3d spixel[2];
				for ( int x = 0; x < 3; ++ x )
					for ( int y = 0; y < 3; ++ y ) {
						Vec3b pixel = image.at<Vec3b>(i+x, j+y);
						for ( int k = 0; k < 3; ++ k )
							for ( int m = 0; m < 2; ++ m )
								spixel[m][k] += pixel[k] * sobel2d_kernel[m][x][y];
					}
				for ( int k = 0; k < 3; ++ k )
					result.at<Vec3b>(i+1, j+1)[k] = sqrt ( two ( spixel[0][k] ) + two ( spixel[1][k] ) );
			}
		return result;
	}

	pair<Mat,Mat> sobel2d_gradients ( const Mat &image ) {
		Mat resultx ( image.rows, image.cols, CV_64FC3 );
		Mat resulty ( image.rows, image.cols, CV_64FC3 );
		for ( int i = 0; i < image.rows; ++ i )
			for ( int j = 0; j < image.cols; ++ j )
				for ( int k = 0; k < 3; ++ k )
					resultx.at<Vec3d>(i,j)[k] = resulty.at<Vec3d>(i,j)[k] = image.at<Vec3d>(i,j)[k];

		for ( int i = 0; i + 2 < image.rows; ++ i ) {
			for ( int j = 0; j + 2 < image.cols; ++ j ) {
				Vec3d spixel[2];
				for ( int x = 0; x < 3; ++ x )
					for ( int y = 0; y < 3; ++ y ) {
						Vec3d pixel = image.at<Vec3d>(i+x, j+y);
						for ( int k = 0; k < 3; ++ k )
							for ( int m = 0; m < 2; ++ m )
								spixel[m][k] += pixel[k] * sobel2d_kernel[m][x][y];
				for ( int k = 0; k < 3; ++ k ) {
					resultx.at<Vec3d>(i+1, j+1)[k] = (spixel[0][k]);
					resulty.at<Vec3d>(i+1, j+1)[k] = (spixel[1][k]);
				}
			}
		}
		return make_pair(resultx, resulty);
	}
}


	CDHFeatureExtractor::CDHFeatureExtractor() {
		q_channels[0] = 10;
		q_channels[1] = q_channels[2] = 3;
		q_phi = 18;
		distance = 1;
	}
	double CDHFeatureExtractor::lab_coeff[3][3] = {{0.412453, 0.357580, 0.180423}, {0.212671, 0.715160, 0.072169}, {0.019334, 0.119193, 0.950227}};
	double CDHFeatureExtractor::illuminant_coeff[3] = {0.950450, 1, 1.088754};
	double CDHFeatureExtractor::channel_ranges[3][2] = {{0,100},{-60,60},{-60,60}};
	double CDHFeatureExtractor::fl ( double u ) {
		if ( u > 0.008856 )
			return 116 * pow ( u, 1.0/3 ) - 16;
		else
			return 903.3 * pow ( u, 1.0/3 );
	}
	double CDHFeatureExtractor::fab ( double u ) {
		if ( u > 0.008856 )
			return pow ( u, 1.0/3 );
		else
			return 7.787 * u + 16.0/116;
	}
	double CDHFeatureExtractor::f ( double u ) {
		if ( u > pow ( 6.0 / 29 , 3 ) )
			return pow ( u, 1.0 / 3 );
		else
			return ( 1.0/3 ) * pow ( 29.0/6, 2 ) * u + 4.0/29;
	}

	Mat CDHFeatureExtractor::rgb2lab ( const Mat &image ) {
		Mat result ( image.rows, image.cols, CV_64FC3 );

		Mat result_showable ( image.rows, image.cols, CV_64FC3 );

		double range[3][2];
		for ( int i = 0 ; i < 3; ++ i ) {
			range[i][0] = 1e30;
			range[i][1] = -1e30;
		}
		for ( int i = 0; i < result.rows; ++ i ) {
			for ( int j = 0; j < result.cols; ++ j ) {
				Vec3b rgb = image.at<Vec3b>(i,j);
				Vec3d x;
				for ( int c = 0; c < 3; ++ c ) {
					x[c] = 0;
					for ( int k = 0; k < 3; ++ k )
						x[c] += ((int)rgb[k])/255.0 * lab_coeff[c][k];
				}
				result.at<Vec3d>(i,j)[0] = 116 * f ( x[1]/illuminant_coeff[1] ) - 16;
				result.at<Vec3d>(i,j)[1] = 500.0 * ( f ( x[0]/illuminant_coeff[0] ) - f ( x[1]/illuminant_coeff[1] ) );
				result.at<Vec3d>(i,j)[2] = 200.0 * ( f ( x[0]/illuminant_coeff[0] ) - f ( x[2]/illuminant_coeff[2] ) );

				result_showable.at<Vec3d>(i,j)[0] = NORMALIZE(result.at<Vec3d>(i,j)[0],0,100,0,1);
				result_showable.at<Vec3d>(i,j)[1] = NORMALIZE(result.at<Vec3d>(i,j)[1],-60,60,0,1);
				result_showable.at<Vec3d>(i,j)[2] = NORMALIZE(result.at<Vec3d>(i,j)[2],-60,60,0,1);

				for ( int k = 0; k < 3; ++ k ) {
					range[k][0] = min ( range[k][0], result.at<Vec3d>(i,j)[k] );
					range[k][1] = max ( range[k][1], result.at<Vec3d>(i,j)[k] );
				}
				/*
				cout << result.at<Vec3d>(i,j)[0] << " ";
				cout << result.at<Vec3d>(i,j)[1] << " ";// <<  500.0 * ( fab ( x[0]/illuminant_coeff[0] ) - fab ( x[1]/illuminant_coeff[1] ) ) << endl;
				cout << result.at<Vec3d>(i,j)[2] << " ";// <<  200.0 * ( fab ( x[0]/illuminant_coeff[0] ) - fab ( x[2]/illuminant_coeff[2] ) ) << endl;
				cout << endl;
				*/
			}
			//cout << endl;
		}
		/*
		cout << "RANGES" << endl;
		for ( int i = 0; i < 3; ++ i )
			cout << range[i][0] << " " << range[i][1] << endl;
		/*
		namedWindow( "Display Image 2", CV_WINDOW_FREERATIO);
		imshow( "Display Image 2", result_showable );
		*/
		return result;
	}
	vector < vector < double > > CDHFeatureExtractor::get_edge_orientations ( const Mat& image ) {
		int N = image.rows;
		int M = image.cols;
		pair<Mat,Mat> sobel_gradients = sobel2d_gradients ( image );
		Mat mg_1 = image.clone();
		Mat mg_2 = image.clone();

		double g_xx, g_yy, g_xy, tempx, tempy, phi, g_1, g_2;
		double range[2] = { 1e30, -1e30 };
		double phirange[2] = { 1e30, -1e30 };

		vector < vector < pair < double, double > > > mat_gradients ( N );
		vector < vector < double > > orientations ( N ); // Max gradients and orientations
		for ( int i = 0; i < sobel_gradients.first.rows; ++ i ) {
			mat_gradients[i] = vector < pair < double, double > > ( M );
			orientations[i] = vector < double > ( M );
		}
		for ( int i = 0; i < sobel_gradients.first.rows; ++ i )
			for ( int j = 0; j < sobel_gradients.first.cols; ++ j ) {
				g_xx = g_yy = g_xy = 0;
				for ( int k = 0; k < 3; ++ k ) {
					tempx = sobel_gradients.first.at<Vec3d>(i,j)[k];
					g_xx += two ( tempx );
					tempy = sobel_gradients.second.at<Vec3d>(i,j)[k];
					g_yy += two ( tempy );
					g_xy += tempx * tempy;
				}
				phi = 0.5 * ( atan ( 2.0 * g_xy / ( g_xx - g_yy + EPS ) ) );
				//cout << "END READ" << endl;
				g_1 = sqrt ( 0.5 * ( ( g_xx + g_yy ) + ( g_xx - g_yy ) * cos ( 2 * phi ) + 2 * g_xy * sin ( 2 * phi ) ) );
				g_2 = sqrt ( 0.5 * ( ( g_xx + g_yy ) + ( g_xx - g_yy ) * cos ( 2 * ( phi + PI/2 ) ) + 2 * g_xy * sin ( 2 * ( phi + PI/2 ) ) ) );

				// Temporal
				if ( g_1 > g_2 ) orientations[i][j] = phi;
				else orientations[i][j] = phi + PI/2;
				orientations[i][j] += PI/4;
				//cout << ( g_xy / ( g_xx - g_yy + EPS ) ) << " " <<  orientations[i][j] << endl;

				if ( phirange[0] > orientations[i][j] )
					phirange[0] = orientations[i][j];
				if ( phirange[1] < orientations[i][j] )
					phirange[1] = orientations[i][j];

				mat_gradients[i][j] = make_pair ( min ( g_1, g_2 ), max ( g_1, g_2 ) );
				range[0] = min ( range[0], mat_gradients[i][j].first );
				range[1] = max ( range[1], mat_gradients[i][j].second );

			}
		/*
		cout << "RANGE " << range[0] << " " << range[1] << endl;
		cout << "RANGE FOR PHI " << phirange[0] << " " << phirange[1] << endl;
		*/
		for ( int i = 0; i < sobel_gradients.first.rows; ++ i )
			for ( int j = 0; j < sobel_gradients.first.cols; ++ j ) {
				//cout << g_1 << " " << g_2 << endl;
				for ( int k = 0; k < 3; ++ k ) {
					mg_1.at<Vec3d>(i,j)[k] = NORMALIZE(mat_gradients[i][j].first,(int)(range[0]),(int)(range[1]),0,1);
					mg_2.at<Vec3d>(i,j)[k] = NORMALIZE(mat_gradients[i][j].second,(int)(range[0]),(int)(range[1]),0,1);
				}
			}
		/*
		namedWindow( "Minimum Gradient Image", CV_WINDOW_FREERATIO);
		namedWindow( "Maximum Gradient Image", CV_WINDOW_FREERATIO);
	  imshow( "Minimum Gradient Image", mg_1 );
	  imshow( "Maximum Gradient Image", mg_2 );
		*/
		return orientations;
	}

	inline int CDHFeatureExtractor::quantize_orientation ( double theta ) {
		return CUANTIZE(theta, 0, PI, CDHFeatureExtractor::q_phi );
	}
	inline int CDHFeatureExtractor::quantize_color ( double l, double a, double b ) {
		return CUANTIZE(l,channel_ranges[0][0],channel_ranges[0][1],q_channels[0]) * q_channels[1] * q_channels[2] +
										CUANTIZE(a,channel_ranges[1][0],channel_ranges[1][1],q_channels[1]) * q_channels[2] +
										CUANTIZE(b,channel_ranges[2][0],channel_ranges[2][1],q_channels[2]);
	}
	vector<double> CDHFeatureExtractor::build_cdh ( const Mat &image, const vector < vector < double > > &orientations ) {
		vector<double> cdh ( q_channels[0] * q_channels[1] * q_channels[2] + q_phi + 1 );
		int N = image.rows;
		int M = image.cols;
		for ( int i = 0; i < N; ++ i )
			for ( int j = 0; j < M; ++ j ) {
				Vec3d pixel = image.at<Vec3d>(i,j);
				for ( int x = 0; x <= distance; ++ x )
					for ( int y = 0; y <= distance; ++ y ) {
						if ( x == 0 && y == 0 ) continue;
						int nx = i + x;
						int ny = j + y;
						if ( nx >= N || ny >= M ) continue;
						Vec3d npixel = image.at<Vec3d>(nx,ny);

						// H(color)
						int q_theta_1 = CDHFeatureExtractor::quantize_orientation(orientations[i][j]);
						int q_theta_2 = CDHFeatureExtractor::quantize_orientation(orientations[nx][ny]);
						int q_color_1 = CDHFeatureExtractor::quantize_color(pixel[0], pixel[1], pixel[2]);
						int q_color_2 = CDHFeatureExtractor::quantize_color(npixel[0], npixel[1], npixel[2]);

						double difference = 0;
						for ( int k = 0; k < 3; ++ k )
							difference += two(pixel[k] - npixel[k]);
						difference = sqrt ( difference );

						if ( q_theta_1 == q_theta_2 )
							cdh[ q_color_1 ] += difference;
						if ( q_color_1 == q_color_2 ) {
							//cout << 90 + q_theta_1 << " " << difference << endl;
							cdh[ 90 + q_theta_1 ] += difference;
						}
					}
			}
		return cdh;
	}
	cv::Mat CDHFeatureExtractor::extractFeatures ( Mat& image ) {
		int N = image.rows;
		int M = image.cols;
		Mat lab = rgb2lab ( image );
		vector < vector < double > > edge_orientations = CDHFeatureExtractor::get_edge_orientations(lab);
		cout << "FEATURES BEGIN" << endl;
		vector < double > features = CDHFeatureExtractor::build_cdh(lab,edge_orientations);
		Mat ans ( features.size(), 1, CV_64F );
		for ( int i = 0; i < features.size(); ++ i )
			ans.ptr<double>(i)[0] = features[i];
		//waitKey();
		cout << "FEATURES END" << endl;
		// print features
		ofstream out("histogram.m");
		out << "Y = [ " << endl;
		for ( int i = 0; i < features.size(); ++ i )
			out << features[i] << endl;
		out << " ] " << endl;
		// end features
		return ans;
	}
