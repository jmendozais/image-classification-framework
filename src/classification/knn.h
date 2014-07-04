/*
 * knn.h
 *
 *  Created on: Jan 19, 2013
 *      Author: jmendoza
 */

#ifndef KNN_H_
#define KNN_H_

#include "classifier.h"
#include <opencv/ml.h>
#include <boost/shared_ptr.hpp>

typedef boost::shared_ptr<cv::KNearest> KNearestPtr;
class KNNClassifier: public  ClassifierInterface {
public:
	KNNClassifier();
	KNNClassifier( int K );
	virtual void train ( const std::vector< std::vector<double> >, const std::vector<int> );
	virtual int predict ( const std::vector<double>& );
private:
	int k_;
	int dim_;
	KNearestPtr knn_;
};

#endif /* KNN_H_ */
