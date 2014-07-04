/*
 * knn.cpp
 *
 *  Created on: Jan 19, 2013
 *      Author: jmendoza
 */

#include "knn.h"
#include <opencv/cv.h>
#include <opencv/ml.h>
#include <cassert>

KNNClassifier::KNNClassifier( ){
	k_ = 10;
}
KNNClassifier::KNNClassifier( int k ){
	k_ = k;
}
void KNNClassifier::train ( const std::vector< std::vector<double> > features, const std::vector<int> classes ) {
	cv::Mat training_data ( features.size() , features[0].size() , CV_32FC1 );
	cv::Mat	training_responses ( features.size(), 1, CV_32FC1 );
	int size = features.size(); int dim = features[0].size();
	for ( int i = 0; i < size; ++ i ) {
		training_responses.ptr<float>(i)[0] = (float) classes[i];
		for ( int j = 0; j < dim; ++ j )
			training_data.ptr<float>(i)[j] = features[i][j];
	}
	dim_ = features[0].size();
	knn_ = KNearestPtr( new cv::KNearest() );
	//std::cout << "BEGIN TRAIN " << features.size() << " " << features[0].size() << std::endl;
	knn_->train(training_data, training_responses);
	//std::cout << "END TRAIN " << features.size() << " " << features[0].size() << std::endl;
}
int KNNClassifier::predict ( const std::vector<double>& feature ) {
	//std::cout << feature.size() << "  " << dim_ << std::endl;
	//std::assert ( feature.size() == dim_ );
	cv::Mat sample ( 1, dim_, CV_32FC1 );
	for ( int i = 0; i < dim_; ++ i )
		sample.ptr<float>(0)[i] = feature[i];
	cv::Mat nearests ( 1, k_, CV_32FC1 ), a, b;
	float resp = (int)knn_->find_nearest ( sample, k_, a, nearests, b );
	/*
	for ( int i = 0; i < k_; ++ i )
		std::cout << nearests.ptr<float>(0)[i] << " ";
	std::cout << std::endl;
	*/
	return (int)resp;
}
