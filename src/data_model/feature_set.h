/*
 * feature_set.h
 *
 *  Created on: Jan 14, 2013
 *      Author: jmendoza
 */

#ifndef FEATURE_SET_H_
#define FEATURE_SET_H_

#include "../new/data.h"


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <opencv/cv.h>
class FeatureSet {
private:
	std::vector<int> data_ids_;
	std::vector<std::vector<double> > mat_;
	int size_;
	int dim_;
	std::string path_;
public:
	FeatureSet ();
	FeatureSet ( int size, int dim );
	void write( std::string path );
	void read( std::string path );
	void change_feature_vectors ( cv::Mat fv );
	void normalize ();

	void set_value( int i, int j, double val );
	cv::Mat cv_mat() const;
	std::vector< std::vector<double> > std_mat() const;
	double get_value( int i, int j) const;


	int size() const { return size_; }
	int dim() const { return dim_; }
};

#endif /* FEATURE_SET_H_ */
