/*
 * batch.cpp
 *
 *  Created on: Jan 24, 2013
 *      Author: jmendoza
 */
#include "batch.h"
#include <iostream>
ImageTransformBatch::ImageTransformBatch ( ImageSet &source, ImageSet &target, ImageTransform &transform ) {
	source_ = &source;
	target_ = &target;
}
void ImageTransformBatch::perform() {
	std::cerr << "Not implemented yet" << std::endl;
}

FeatureExtractorBatch::FeatureExtractorBatch ( ImageSet &source, FeatureSet &target, FeatureExtractorInterface &extractor ) {
	source_ = &source;
	target_ = &target;
	extractor_ = &extractor;
}
void FeatureExtractorBatch::perform() {
	cv::Mat temp = source_->getImage(0);
	int dim = extractor_->extractFeatures ( temp ).rows;

	target_->resize ( source_->size(), dim );
	for ( int i = 0; i < source_->size(); ++ i ) {
		cv::Mat image = source_->getImage(i);
		cv::Mat fv = extractor_->extractFeatures(image) ;
		std::cout << "FV " << i << " -> " << fv << std::endl;
		for ( int j = 0; j < dim; ++ j )
			target_->setValue(i, j, fv.ptr<double>(j)[0] );
	}
}




