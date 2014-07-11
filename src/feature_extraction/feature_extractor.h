/*
 * FeatureExtractorInterface.h
 *
 *  Created on: Jan 14, 2013
 *      Author: jmendoza
 */

#ifndef FEATURE_EXTRACTOR_H_
#define FEATURE_EXTRACTOR_H_

#include <vector>
#include <opencv/cv.h>
//#include "../data_model/data_model.h";

class FeatureExtractorInterface {
public:
	virtual cv::Mat extractFeatures ( cv::Mat &image ) = 0;
};

#endif /* FEATUREEXTRACTORINTERFACE_H_ */
