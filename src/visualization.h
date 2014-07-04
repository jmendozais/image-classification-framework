/*
 * visualization.h
 *
 *  Created on: Jan 31, 2013
 *      Author: jmendoza
 */

#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_
#include "ipcore/util.h"
class GaborWaveletVisualizer {
public:
	GaborWaveletVisualizer( const GaborWavelet &wavelets );

	void show() const;

	const cv::Mat& getImage() const;
	void setImage(const cv::Mat& image);

	const GaborWavelet& getGaborWavelets() const {
		return wavelets_;
	}
	void setGaborWavelets( GaborWavelet& wavelets ) {
		wavelets_ = wavelets;
	}

	// statics functions
	static void on_trackbar( int, void* );
	//static void on_trackbar_uf( int, void* );
	//static void update_view( int, void* );

	// slider
	int lf_slider;
	int lf_slider_max;
	int uf_slider;
	int uf_slider_max;

private:
	cv::Mat image_;
	GaborWavelet wavelets_;
};

inline GaborWaveletVisualizer::GaborWaveletVisualizer(
		const GaborWavelet& wavelets) {
	wavelets_ = wavelets;
	double lf = wavelets_.getLowFrequencyBound();
	double uf = wavelets_.getHighFrequencyBound();
	uf_slider_max = 100;
	lf_slider_max = 100;
	uf_slider = uf * 100;
	lf_slider = lf * 100;
	cv::namedWindow("Gabor");
	cv::createTrackbar( "Lower freq bound ", "Gabor", &lf_slider, lf_slider_max, GaborWaveletVisualizer::on_trackbar, this );
	cv::createTrackbar( "High freq bound ", "Gabor", &uf_slider, uf_slider_max, GaborWaveletVisualizer::on_trackbar, this );
}

inline void GaborWaveletVisualizer::show() const {
	std::vector<cv::Mat> mags = wavelets_.getMagnitudes(image_);
	GaborWaveletVisualizer::on_trackbar(lf_slider, (void*)this);
}

inline const cv::Mat& GaborWaveletVisualizer::getImage() const {
	return image_;
}

inline void GaborWaveletVisualizer::setImage(const cv::Mat& single_channel_image) {

	cv::Mat gray;
	cvtColor ( single_channel_image, gray, CV_BGR2GRAY );
	cv::Mat_<double> gray_double ( gray );
	image_ = ipcore::normalize(gray_double);
}

inline void GaborWaveletVisualizer::on_trackbar(int int1, void* object) {

	GaborWaveletVisualizer* visualizer_ptr_ = (GaborWaveletVisualizer*) object;
	int lf_slider = visualizer_ptr_->lf_slider;
	int lf_slider_max = visualizer_ptr_->lf_slider_max;
	int uf_slider = visualizer_ptr_->uf_slider;
	int uf_slider_max = visualizer_ptr_->uf_slider_max;

	cv::Mat query_image_ = visualizer_ptr_->getImage();
	GaborWavelet wavelets_ = visualizer_ptr_->getGaborWavelets();
	int rows = wavelets_.getNumScales();
	int cols = wavelets_.getNumOrientations();

	double lf = lf_slider*1.0/lf_slider_max;
	double uf = uf_slider*1.0/uf_slider_max;
	wavelets_.setLowFrequencyBound(lf);
	wavelets_.setHighFrequencyBound(uf);
	cv::vector<cv::Mat> mod_fb, fb = wavelets_.getFilterBank();
	for ( int i = 0; i < rows; ++ i ) {
		int filter_id = (i*cols+2);
		mod_fb.push_back ( fb[2*(filter_id)] );
		mod_fb.push_back ( fb[2*(filter_id)+1] );
		/*
		cv::namedWindow("filter");
		cv::imshow("filter", ipcore::normalize(fb[2*(filter_id)]));
		cv::waitKey();
		*/
	}
	wavelets_.setFilterBank(mod_fb);

	std::vector<cv::Mat> images = wavelets_.getMagnitudes(query_image_);
	visualizer_ptr_->setGaborWavelets(wavelets_);

 	for ( int i = 0; i < images.size(); ++ i ) {
		cv::resize(images[i], images[i], cv::Size(200,200));
		images[i] = ipcore::normalize(images[i]);
	}
	std::cout << "frequency bounds " << wavelets_.getLowFrequencyBound() << " " << wavelets_.getHighFrequencyBound() << images.size() << std::endl;
	int max_rows, max_cols;
	max_rows = max_cols = -1;
	for ( int i = 0; i < images.size(); ++ i ) {
		max_rows = std::max ( max_rows, images[i].rows );
		max_cols = std::max ( max_cols, images[i].cols );
	}
	rows = 3; cols = 3;
	cv::Mat resp ( rows * max_rows, cols * max_cols, CV_64FC1 );
	for ( int i = 0; i < rows; ++ i )
		for ( int j = 0; j < cols; ++ j ) {
			if ( i*cols + j >= images.size() ) continue;
			cv::Rect roi ( j*max_cols, i*max_rows, max_cols, max_rows );
			cv::Mat image_roi ( resp, roi );
			ipcore::normalize(images[i*cols+j]).copyTo(image_roi);
		}
	std::stringstream ss;
	ss << images.size();
	std::string tok;
	ss >> tok;
	cv::imshow("Gabor", resp);
}



class GaborFilterVisualizer {
public:
	double lambda;
	double psi;
	double ganma;
	double theta;
	double sigma;
	cv::Mat createFilter() {
		double sigma_x = sigma;
		double sigma_y = sigma/ganma;
		double step = 3;
	}
};
class GaborClustererVisualizer {
private:
	GaborWavelet wavelets_;
	cv::Mat image_;

public:
	GaborClustererVisualizer( GaborWavelet wavelets ) {
		wavelets_ = wavelets;
	}

	const cv::Mat& getImage() const {
		return image_;
	}

	void setImage(const cv::Mat& image) {
		image_ = image;
	}

	const GaborWavelet& getWavelets() const {
		return wavelets_;
	}

	void setWavelets(const GaborWavelet& wavelets) {
		wavelets_ = wavelets;
	}
	cv::Mat clusterize( cv::Mat input, int k ) {
		// cast
		cv::Mat gray;
		cvtColor ( input, gray, CV_BGR2GRAY );
		cv::Mat_<double> gray_double ( gray );
		input = gray_double;
		// create data
		std::vector<cv::Mat> magnitudes = wavelets_.getMagnitudes(input);

		std::cout << "frequency band " << wavelets_.getLowFrequencyBound() << " "
				<< wavelets_.getHighFrequencyBound() <<  " magnitudes size " << magnitudes.size() << std::endl;
		cv::Mat data ( input.rows * input.cols, magnitudes.size(), CV_32FC1 );
		for ( int i = 0; i < input.rows; ++ i ) {
			for ( int j = 0; j < input.cols; ++ j )
				for ( int k = 0; k < magnitudes.size(); ++k )
					data.ptr<float>(i*input.cols + j)[k] = (float)magnitudes[k].ptr<double>(i)[j];
		}
		// clustering
		cv::Mat labels, centers;
		cv::kmeans(data, k, labels, cv::TermCriteria( CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 0.0001, 10000), 5, KMEANS_PP_CENTERS, centers );
		image_ = input.clone();
		std::set<int> slabels;
		for ( int i = 0; i < input.rows; ++ i ) {
			for ( int j = 0; j < input.cols; ++ j ) {
				slabels.insert(labels.ptr<int>(i*input.cols + j)[0]);
				double gray_level = labels.ptr<int>(i*input.cols + j)[0]*1.0/labels.rows;
				image_.ptr<double>(i)[j] = gray_level;
			}
		}
		cv::namedWindow("TEST");
		cv::imshow("TEST", ipcore::normalize(image_));
		std::cout << "end clustering" << std::endl;
		return image_;
	}
};

#endif /* VISUALIZATION_H_ */
