/*
 * ClassificationModel.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: jmendoza
 */

#include "classification_model.h"

ClassificationModel::ClassificationModel() {
	// TODO Auto-generated constructor stub

}

ClassificationModel::~ClassificationModel() {
	// TODO Auto-generated destructor stub
}

void ClassificationModel::set_image_transform ( ImageTransform* transform ) {
	transform_ = transform;
}
void ClassificationModel::set_feature_extractor ( FeatureExtractorInterface *feature_extractor ) {
	feature_extractor_ = feature_extractor;
}
void ClassificationModel::set_feature_selector ( FilterSelectorInterface *feature_selector ) {
	feature_selector_ = feature_selector;
}
void ClassificationModel::set_classifier ( ClassifierInterface *classifier ) {
	classifier_ = classifier;
}

FeatureSet ClassificationModel::extractFeatures ( ClassificationImageSet image_set ) {
	cv::Mat temp = image_set.get_image(0);
	int dim = feature_extractor_->extractFeatures ( temp ).rows;
	FeatureSet ans ( image_set.size(), dim );
	for ( int i = 0; i < image_set.size(); ++ i ) {
		cv::Mat image = image_set.get_image(i);
		cv::Mat fv = feature_extractor_->extractFeatures(image) ;
		std::cout << i << ".- FV -> " << fv << std::endl;
		for ( int j = 0; j < dim; ++ j )
			ans.set_value(i, j, fv.ptr<double>(j)[0] );
	}
	return ans;
}
void ClassificationModel::transformImages ( const ClassificationImageSet &input, ClassificationImageSet &output ) {
	for ( int i = 0; i < input.size(); ++ i ) {
		cv::Mat out;
		if ( !transform_->transform( input.get_image(i), out ) )
			std::cerr << "Image " << i << " no procesada" << std::endl;
		output.set_image( out, i );
	}
}




