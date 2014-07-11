/*
 * feature_set.cpp
 *
 *  Created on: Jan 14, 2013
 *      Author: jmendoza
 */

#include "feature_set.h"
#include "../ipcore/util.h"
FeatureSet::FeatureSet ( int size , int dims ) {
	size_ = size;
	dim_ = dims;
	data_ids_.resize(size_);
	mat_.resize(size_);
	for ( int i = 0; i < size_; ++ i )
		mat_[i].resize(dim_);
}
FeatureSet::FeatureSet () {
	size_ = 0;
	dim_ = 0;
	data_ids_.resize(size_);
	mat_.resize(size_);
	for ( int i = 0; i < size_; ++ i )
		mat_[i].resize(dim_);
}
void FeatureSet::write( std::string path ) {
	path_ = path;
	std::ofstream fout( path_.c_str() );
	fout << size_ << "\n";
	fout << dim_ << "\n";
	for ( int i = 0; i < size_; ++ i ) {
		fout << data_ids_[i];
		for ( int j = 0; j < dim_; ++ j )
			fout << " " << mat_[i][j];
		fout << "\n";
	}
	fout.close();
}
void FeatureSet::read(std::string path ) {
	path_ = path;
	std::ifstream fin ( path_.c_str() );
	fin >> size_;
	fin >> dim_;
	data_ids_.resize(size_);
	mat_.resize(size_);
	for ( int i = 0; i < size_; ++ i ) {
		fin >> data_ids_[i];
		mat_[i].resize(dim_);
		for ( int j = 0; j < dim_; ++ j )
			fin >> mat_[i][j];
	}
}
cv::Mat FeatureSet::cv_mat() const {
	cv::Mat ans ( size_, dim_, CV_64F);
	for ( int i = 0; i < size_; ++ i )
		for ( int j = 0; j < dim_; ++ j )
			ans.ptr<double>(i)[j] = mat_[i][j];
	return ans;
}
std::vector< std::vector<double> > FeatureSet::std_mat() const {
	return mat_;
}
#include <iostream>
void FeatureSet::change_feature_vectors ( cv::Mat fv ) {
	if ( fv.rows != mat_.size() ) {
		std::cerr << "changing feature vectors: Not compatible sizes" << std::endl;
		return;
	}
	dim_ = fv.cols;
	for ( int i = 0; i < fv.rows; ++ i ) {
		mat_[i].resize(dim_);
		for ( int j = 0; j < fv.cols; ++ j )
			mat_[i][j] = fv.ptr<double>(i)[j];
	}
}
void FeatureSet::normalize() {
	for ( int i = 0; i < dim_; ++ i ) {
		double min_ = 1e30;
		double max_ = -min_;
		double val;
		for ( int j = 0; j < size_; ++ j ) {
			val = get_value (j, i);
			if ( min_ > val ) min_ = val;
			if ( max_ < val ) max_ = val;
		}
		std::cout << min_ << ", " << max_ << std::endl;
		for ( int j = 0; j < size_; ++ j ) {
			val = get_value (j, i);
			val = ipcore::normalize ( val, min_, max_, -1, 1 );
			set_value (j, i, val );
		}
	}
}
double FeatureSet::get_value( int i, int j ) const { return mat_[i][j]; }
void FeatureSet::set_value( int i, int j, double val ) { mat_[i][j] = val; }

