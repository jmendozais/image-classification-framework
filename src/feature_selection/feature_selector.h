/*
 * feature_selector.h
 *
 *  Created on: Jan 14, 2013
 *      Author: jmendoza
 */

#ifndef FEATURE_SELECTOR_H_
#define FEATURE_SELECTOR_H_

//#include "../data_model/data_model.h";
#include "../evolutive_computation.h"
#include "../classification/classifier.h"
#include "../validation.h"
#include <cmath>
class FilterSelectorInterface {
public:
	virtual FeatureSet selectFeatures ( const FeatureSet& ) = 0;
};

class ClassifierFitnessFunction : public ObjectiveFunction {
public:
	ClassifierFitnessFunction( ClassifierInterface * classifier, AbstractClassifierEvaluator* evaluator, const std::vector< std::vector<double> > &data, const std::vector<int> &responses  ) {
		classifier_ = classifier;
		data_ = data;
		responses_ = responses;
		evaluator_ = evaluator;
	}
	virtual double perform( double *v ) {
		//std::cout << "FUNCTIONOID INHERITANCE ... OK " << data_[0].size() << std::endl;
		int new_dim = 0;
		for ( int i = 0; i < data_[0].size(); ++ i )
			if ( v[i] > 0.5 ) ++ new_dim;
		//std::cout << "new dim " << new_dim << std::endl;
		std::vector< std::vector<double> > selected_data;
		selected_data.resize(data_.size());
		for ( int i = 0; i < data_.size(); ++ i ) {
			int idx = 0;
			selected_data[i].resize(new_dim);
			for ( int j = 0; j < data_[i].size(); ++ j )
				if ( v[j] > 0.5 )
					selected_data[i][idx++] = data_[i][j];
		}
		evaluator_->changeData(selected_data, responses_);
		ClassifierEvaluatorResult result = evaluator_->evaluate();
		//std::cout << "output " << result.global_accuracy << std::endl;
		return pow(2, 100*(1-result.global_accuracy));
	}
private:
	ClassifierInterface * classifier_;
	AbstractClassifierEvaluator * evaluator_;
	std::vector< std::vector<double> > data_;
	std::vector< int > responses_;
};
class AbstractWrapperFeatureSelector {
protected:
	ClassifierInterface *classifier_;
	AbstractClassifierEvaluator *evaluator_;
public:
	AbstractWrapperFeatureSelector ( ClassifierInterface *classifier, AbstractClassifierEvaluator * evaluator ) {
		classifier_ = classifier;
		evaluator_ = evaluator;
	}
	virtual std::vector< std::vector< double > > selectFeatures ( const std::vector<std::vector<double> > &data, const std::vector<int> &classes ) = 0;
};

#endif /* FEATURE_SELECTOR_H_ */
