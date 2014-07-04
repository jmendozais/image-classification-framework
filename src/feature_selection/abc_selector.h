/*
 * abc_selector.h
 *
 *  Created on: Jan 19, 2013
 *      Author: jmendoza
 */

#ifndef ABC_SELECTOR_H_
#define ABC_SELECTOR_H_
#include "feature_selector.h";
class ABCWrapperSelector : public AbstractWrapperFeatureSelector {
public:
		ABCWrapperSelector ( ClassifierInterface * classifier, AbstractClassifierEvaluator * evaluator ) : AbstractWrapperFeatureSelector ( classifier, evaluator ) {}
		virtual std::vector< std::vector< double > > selectFeatures ( const std::vector<std::vector<double> > &data, const std::vector<int> &classes );
};

#endif /* ABC_SELECTOR_H_ */
