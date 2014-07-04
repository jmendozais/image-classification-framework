/*
 * classifier.h
 *
 *  Created on: Jan 14, 2013
 *      Author: jmendoza
 */

#ifndef CLASSIFIER_H_
#define CLASSIFIER_H_

#include <vector>
class ClassifierInterface {
public:
	virtual void train ( const std::vector< std::vector<double> >, const std::vector<int> ) = 0;
	virtual int predict ( const std::vector<double>& ) = 0;
};

#endif /* CLASSIFIER_H_ */
