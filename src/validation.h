/*
 * validation.h
 *
 *  Created on: Jan 20, 2013
 *      Author: jmendoza
 */

#ifndef VALIDATION_H_
#define VALIDATION_H_

#include "classification/classifier.h";
#include "new/data.h"
#include <algorithm>

class AbstractClassifierEvaluator {
public:
	AbstractClassifierEvaluator ( ClassifierInterface * classifier, ImageSet* images, FeatureSet* features ) {
		classifier_ = classifier;
		image_set_ = images;
		feature_set_ = features;
	}
	// virtual methods
	virtual ClassifierEvaluatorResult evaluate() = 0;
	// implemented methods
	void changeData ( const std::vector<std::vector<double> > &data, const std::vector<int> &responses ) {
		feature_set_->setData(data);
		feature_set_->setDataClassIds(responses);
	}
	std::vector< std::vector< std::vector<int> > > getBelongMatrix() {
		return belong_matrix_;
	}

	void setBelongMatrix(
			const std::vector<std::vector<std::vector<int> > >& belongMatrix) {
		belong_matrix_ = belongMatrix;
	}

	ClassifierInterface* getClassifier() const {
		return classifier_;
	}

	void setClassifier(ClassifierInterface* classifier) {
		classifier_ = classifier;
	}

	const std::vector<std::vector<int> >& getConfussionMatrix() const {
		return confussion_matrix_;
	}

	void setConfussionMatrix(
			const std::vector<std::vector<int> >& confussionMatrix) {
		confussion_matrix_ = confussionMatrix;
	}

	FeatureSet* getFeatureSet() const {
		return feature_set_;
	}

	void setFeatureSet(FeatureSet* featureSet) {
		feature_set_ = featureSet;
	}

	ImageSet* getImageSet() const {
		return image_set_;
	}

	void setImageSet(ImageSet* imageSet) {
		image_set_ = imageSet;
	}

protected:
	ClassifierInterface* classifier_;
	ImageSet* image_set_;
	FeatureSet* feature_set_;
	std::vector< std::vector<int> > confussion_matrix_;
	std::vector< std::vector< std::vector<int> > > belong_matrix_;
};

class LOOCV :public AbstractClassifierEvaluator {
public:
	LOOCV ( ClassifierInterface * classifier, ImageSet* images, FeatureSet* features) : AbstractClassifierEvaluator(classifier, images, features) {}
	virtual ClassifierEvaluatorResult evaluate() {
		int num_classes_ = feature_set_->num_classes();
		int data_size_ = feature_set_->size();
		std::vector<std::vector<double> > data_ = feature_set_->std_mat();
		std::vector<int> responses_ = feature_set_->getDataClassIds();

		confussion_matrix_.resize(num_classes_);
		belong_matrix_.resize(num_classes_);
		for ( int i = 0; i < num_classes_; ++ i ) {
			confussion_matrix_[i].resize(num_classes_);
			belong_matrix_[i].resize(num_classes_);
		}
		for ( int i = 0; i < num_classes_; ++ i )
			for ( int j = 0; j < num_classes_; ++ j ) {
				confussion_matrix_[i][j] = 0;
				belong_matrix_[i][j] = std::vector<int>();
			}

		for ( int i = 0; i < data_size_; ++ i ) {
			std::cout << "EVALUATOR " << i << " " << data_size_ << std::endl;
			std::vector<double> current = data_[i];
			int response = responses_[i];
			data_[i] = data_[(i+1)%data_.size()];
			responses_[i] = responses_[(i+1)%data_.size()];
			classifier_->train(data_,responses_);
			int outcome = classifier_->predict(current);
			std::cout << "response, outcome " << response << " " << outcome << std::endl;
			confussion_matrix_[response][outcome]++;
			data_[i] = current;
			responses_[i] = response;
			std::cout << "NEXT" << std::endl;
		}
		ClassifierEvaluatorResult result;
		result.confussion_matrix = confussion_matrix_;
		result.accuracy.resize( num_classes_ );
		result.precision.resize( num_classes_ );
		result.sensibility.resize( num_classes_ );
		result.specificity.resize( num_classes_ );
		result.belong_matrix = belong_matrix_;
		for ( int i = 0; i < num_classes_; ++ i ) {
			int tp, tn, fp, fn;
			tp = tn = fp = fn = 0;
			for ( int j = 0; j < num_classes_; ++ j )
				for ( int k = 0; k < num_classes_; ++ k ) {
					if ( i != j && i != k )
						tn += confussion_matrix_[j][k];
					else if ( i == j && i != k )
						fn += confussion_matrix_[j][k];
					else if ( i != j && i == k )
						fp += confussion_matrix_[j][k];
					else if ( i == j && i == k )
						tp += confussion_matrix_[j][k];
				}
			result.precision[i] = tp * 1.0 / ( tp + fp );
			result.sensibility[i] = tp * 1.0 / ( tp + fn );
			result.specificity[i] = tn * 1.0 / ( tn + fp );
			result.accuracy[i] = ( tp + tn ) * 1.0 / ( tp + tn + fp + fn );
		}
		result.global_accuracy = 1e-9;
		for ( int i = 0; i < num_classes_; ++ i)
			result.global_accuracy += confussion_matrix_[i][i];
		result.global_accuracy /= data_.size();
		return result;
	}
};

class KFOLDCV : public AbstractClassifierEvaluator {
public:
	KFOLDCV ( ClassifierInterface * classifier, ImageSet* image_set, FeatureSet* feature_set, int k = 10 ) : AbstractClassifierEvaluator(classifier, image_set, feature_set) {
		k_ = k;
	}
	virtual ClassifierEvaluatorResult evaluate() {
		int num_classes_ = feature_set_->num_classes();
		int data_size_ = feature_set_->size();
		std::cout << "num classes " << num_classes_ << std::endl;
		std::vector<std::vector<double> > data_ = feature_set_->std_mat();
		std::vector<int> responses_ = feature_set_->getDataClassIds();

		confussion_matrix_.resize(num_classes_);
		belong_matrix_.resize(num_classes_);
		for ( int i = 0; i < num_classes_; ++ i ) {
			confussion_matrix_[i].resize(num_classes_);
			belong_matrix_[i].resize(num_classes_);
		}
		for ( int i = 0; i < num_classes_; ++ i )
			for ( int j = 0; j < num_classes_; ++ j ) {
				confussion_matrix_[i][j] = 0;
				belong_matrix_[i][j] = std::vector<int>();
			}

		std::vector<int> idx ( data_.size() );
		for ( int i = 0; i < data_.size(); ++ i )
			idx[i] = i;
		std::random_shuffle( idx.begin(), idx.end() );
		int len = data_.size();
		int fold_len = ceil ( len * 1.0 / k_ );
		// data
		std::vector< std::vector<double> > train_data, test_data;
		std::vector< int > train_responses, test_responses;
		int klen;
		for ( int k = 0; k < k_; ++ k ) {
			std::cout << "K " << k << std::endl;
			klen = std::min( (k+1)*fold_len, len ) - k*fold_len;
			train_data.resize(len - klen);
			train_responses.resize(len - klen);
			test_data.resize(klen);
			test_responses.resize(klen);
			int	c1, c2; c1 = c2 = 0;
			std::vector<int> test_ids;
			for ( int i = 0; i < len; ++ i ) {
				if ( i >= k*fold_len && i < std::min( (k+1)*fold_len, len ) ) {
					test_data[c1] = data_[idx[i]];
					test_responses[c1++] = responses_[idx[i]];
					test_ids.push_back(idx[i]);
				} else {
					train_data[c2] = data_[idx[i]];
					train_responses[c2++] = responses_[idx[i]];
				}
			}
			classifier_->train(train_data, train_responses);
			for ( int i = 0; i < test_data.size(); ++ i ) {
				int outcome = classifier_->predict(test_data[i]);
				//std::cout << "response, outcome " << test_responses[i] << " " << outcome << std::endl;
				confussion_matrix_[ test_responses[i] ][ outcome ]++;
				belong_matrix_[ test_responses[i] ][ outcome ].push_back( test_ids[i] );
			}
		}
		ClassifierEvaluatorResult result;
		result.confussion_matrix = confussion_matrix_;
		result.belong_matrix = belong_matrix_;
		result.accuracy.resize( num_classes_ );
		result.precision.resize( num_classes_ );
		result.sensibility.resize( num_classes_ );
		result.specificity.resize( num_classes_ );

		for ( int i = 0; i < num_classes_; ++ i ) {
			int tp, tn, fp, fn;
			tp = tn = fp = fn = 0;
			for ( int j = 0; j < num_classes_; ++ j )
				for ( int k = 0; k < num_classes_; ++ k ) {
					if ( i != j && i != k )
						tn += confussion_matrix_[j][k];
					else if ( i == j && i != k )
						fn += confussion_matrix_[j][k];
					else if ( i != j && i == k )
						fp += confussion_matrix_[j][k];
					else if ( i == j && i == k )
						tp += confussion_matrix_[j][k];
				}
			result.precision[i] = tp * 1.0 / ( tp + fp );
			result.sensibility[i] = tp * 1.0 / ( tp + fn );
			result.specificity[i] = tn * 1.0 / ( tn + fp );
			result.accuracy[i] = ( tp + tn ) * 1.0 / ( tp + tn + fp + fn );
		}
		result.global_accuracy = 1e-9;
		for ( int i = 0; i < num_classes_; ++ i)
			result.global_accuracy += confussion_matrix_[i][i];
		result.global_accuracy /= data_.size();
		return result;
	}
private:
	int k_;
};

#endif /* VALIDATION_H_ */
