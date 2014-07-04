/*
 * image_transform.h
 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */

#ifndef IMAGE_TRANSFORM_H_
#define IMAGE_TRANSFORM_H_

#include <vector>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/ml.h>
class ImageTransform {
public:
	virtual bool transform ( const cv::Mat& input, cv::Mat& single_output ) = 0;
};

class OrientationNormalization : public ImageTransform {
public:
	virtual bool transform ( const cv::Mat& in, cv::Mat& output ) {
		// otsu thresholding
		cv::Mat input = in.clone();
		cv::Mat gray_input, u_input, temp, mask;
		cvtColor( input, gray_input, CV_BGR2GRAY );
		gray_input.convertTo( u_input, CV_8UC1 );
		cv::threshold(u_input, temp, 100, 255, cv::THRESH_OTSU );

		// invert scalars
		for ( int i = 0; i < temp.rows; ++ i )
			for ( int j = 0; j < temp.cols; ++ j )
				if ( temp.ptr<uchar>(i)[j] == 0 )
					temp.ptr<uchar>(i)[j] = 255;
				else temp.ptr<uchar>(i)[j] = 0;

		// closing holes ( max rad of holes 5px )
		int side = 5;
		int dif = abs ( temp.rows - temp.cols );
		if ( temp.rows < temp.cols ) {
			cv::copyMakeBorder(temp, temp, dif/2 + 1, dif/2 + 1, 0, 0, cv::BORDER_CONSTANT, cv::Scalar(0) );
			cv::copyMakeBorder(input, input, dif/2 + 1, dif/2 + 1, 0, 0, cv::BORDER_CONSTANT, cv::Scalar(255,255,255) );
		} else {
			cv::copyMakeBorder(temp, temp, 0, 0, dif/2 + 1, dif/2 + 1, cv::BORDER_CONSTANT, cv::Scalar(0) );
			cv::copyMakeBorder(input, input, 0, 0, dif/2 + 1, dif/2 + 1, cv::BORDER_CONSTANT, cv::Scalar(255,255,255) );
		}

		cv::copyMakeBorder(temp, temp, side + 1, side + 1, side + 1, side + 1, cv::BORDER_CONSTANT );
		cv::Mat element = cv::getStructuringElement( 2, cv::Size( 2*side + 1, 2*side+1 ), cv::Point( side, side ) );
		cv::morphologyEx ( temp, temp, cv::MORPH_CLOSE, element );
		cv::Mat cropped ( temp, cv::Rect ( side + 1, side + 1, input.cols, input.rows ) );
		mask = cropped;

		cv::Mat threshold_output;
		std::vector<std::vector<cv::Point> > contours;
		std::vector<cv::Vec4i> hierarchy;
		cv::findContours( mask, contours, hierarchy, CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE );
		if ( contours[0].size() > 5 ) {
			cv::RotatedRect ellipse = cv::fitEllipse( cv::Mat(contours[0]) );
			cv::Mat rot_mat = cv::getRotationMatrix2D( ellipse.center, ellipse.angle, 1.0 );
			cv::warpAffine ( input, output, rot_mat, input.size(), 0, cv::BORDER_CONSTANT, cv::Scalar(255,255,255) );
		} else {
			std::cerr << "Error al hallar el contorno en la imagen" << std::endl;
			output = input;
			return false;
		}
		return true;
	}
};
class GaborTextureSegmentator : public ImageTransform {
public:
	virtual bool transform ( const cv::Mat& in, cv::Mat& output ) {

	}
};



#endif /* IMAGE_TRANSFORM_H_ */
