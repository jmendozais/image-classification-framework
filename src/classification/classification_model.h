/*
 * ClassificationModel.h
 *
 *  Created on: Jan 14, 2013
 *      Author: jmendoza
 */

#ifndef CLASSIFICATION_MODEL_H_
#define CLASSIFICATION_MODEL_H_

#include "../image_transformation/image_transform.h"
#include "../feature_extraction/feature_extractor.h"
#include "../feature_selection/feature_selector.h"
#include "classifier.h"
//#include "data_model/data_model.h";


class ClassificationModel {
public:
	ClassificationModel();
	virtual ~ClassificationModel();

	enum TransformType { SINGLE_OUTPUT_TRANSFORM, MULTIPLE_OUTPUT_TRANSFORM };
	// Accessors & Mutators
	void set_image_set ( ClassificationImageSet image_set );
	void set_image_transform ( ImageTransform* transform );
	void set_feature_extractor ( FeatureExtractorInterface *feature_extractor );
	void set_feature_selector ( FilterSelectorInterface *feature_selector );
	void set_classifier ( ClassifierInterface* classifier);

	void transformImages ( const ClassificationImageSet &input, ClassificationImageSet &output );
	FeatureSet extractFeatures ( ClassificationImageSet image_set );
	FeatureSet selectFeatures ( FeatureSet feature_set );

	void kFoldCrossValidation (  ClassificationDataSet, FeatureSet, int );
	void trainClassifier ( FeatureSet feature_set );
	// SomeOutput test_model ();

private:
	// DataTransform *transform_;
	ImageTransform *transform_;
	TransformType transform_type_;
	FeatureExtractorInterface *feature_extractor_;
	FilterSelectorInterface *feature_selector_;
	ClassifierInterface *classifier_;
};

#endif /* CLASSIFICATIONMODEL_H_ */
