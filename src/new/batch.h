/*
 * batch.h
 *
 *  Created on: Jan 24, 2013
 *      Author: jmendoza
 */

#ifndef BATCH_H_
#define BATCH_H_

#include "data.h"
#include "../image_transformation/image_transform.h"
#include "../feature_extraction/gabor_wavelet.h"
class Functor {
public:
	virtual void perform() = 0;
};
class Batch : public Functor {
private:
	Functor *f;
public:
	void setFunctor ( Functor &func ) {
		f = &func;
	}
	virtual void perform() = 0;
};
class ImageTransformBatch : public Batch {
private:
	ImageSet * source_;
	ImageSet * target_;
	// bad idea
	ImageTransform * transform_;
public:
	ImageTransformBatch ( ImageSet &source, ImageSet &target, ImageTransform &transform );
	virtual void perform();
};
class FeatureExtractorBatch : public Batch {
private:
	ImageSet * source_;
	FeatureSet * target_;
	// bad idea
	FeatureExtractorInterface * extractor_;
public:
	FeatureExtractorBatch ( ImageSet &source, FeatureSet &target, FeatureExtractorInterface &extractor );
	virtual void perform();
};


#endif /* BATCH_H_ */
