/*
 * mrMr.h
 *
 *  Created on: Jan 16, 2013
 *      Author: jmendoza
 */

#ifndef MRMR_H_
#define MRMR_H_

#include "feature_selector.h"
#include "../new/data.h"

#include <opencv/cv.h>

struct MRMR_params : AbstractParams {
	int target_dim;
	ClassificationDataSet data_set;
};
class GreedyMinRedundanceMaxRelevanceSelector : public FilterSelectorInterface {
private:
	MRMR_params params_;
public:
	GreedyMinRedundanceMaxRelevanceSelector ( MRMR_params params );
	virtual FeatureSet selectFeatures ( const FeatureSet& features );
};

typedef GreedyMinRedundanceMaxRelevanceSelector GreedyMRMRSelector;

#endif /* MRMR_H_ */
