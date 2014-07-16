/*
 * data.cpp
 *
 *  Created on: Jan 24, 2013
 *      Author: jmendoza
 */

#include "data.h"

#include "../ipcore/util.h"
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <opencv/highgui.h>
namespace fs = boost::filesystem;



ClassificationDataSet::ClassificationDataSet ( std::vector<DataItemPtr> data, std::vector<int> data_class,	std::vector<ClassDescriptorPtr> classes,
												 std::string data_set_path, std::string descriptor_path ) {
		data_ = data;
		data_class_ = data_class;
		classes_ = classes;
		data_set_path_ = data_set_path;
		descriptor_path_ = descriptor_path;
		fs::path dpath ( descriptor_path_ );
		// get the name
				std::string filename = dpath.filename().c_str();
				for ( int i = 0; i < filename.size(); ++ i )
					if ( filename[i] == '.' )
						filename[i] = ' ';

				std::stringstream ss(filename); std::string t1, t2 = filename;
				while ( ss >> t1 ) {
				  name_ = t2;
					t2 = t1;
				}

}
bool ClassificationDataSet::create ( std::string descriptor_path, std::string data_set_path ) {
	data_set_path_ = data_set_path;
	descriptor_path_ = descriptor_path;
	data_.clear();
	classes_.clear();
	data_class_.clear();
	fs::path root ( data_set_path );
	int class_counter = -1, item_counter = -1;
	if ( fs::is_directory ( root ) ) {
		std::cout << "IS DIRECTORY " << root <<  std::endl;
		fs::directory_iterator it = fs::directory_iterator ( root );
		std::cout << "IT init" << std::endl;

		for (; it != fs::directory_iterator(); ++ it ) {
			std::cout << "IN" << std::endl;
			if ( fs::is_directory ( *it ) ) {
				std::cout << "IS 1. DIRECTORY " << std::endl;
				fs::path dir = (*it).path();
				std::cout <<  dir << std::endl;
				ClassDescriptorPtr desc ( new ClassDescriptor() );
				desc->id = ++ class_counter;
				desc->name = dir.filename().string();
				classes_.push_back(desc);
				for ( fs::directory_iterator it2 = fs::directory_iterator ( dir ); it2 != fs::directory_iterator(); ++ it2 ) {
					fs::path file = (*it2).path();
					if ( !fs::is_regular_file ( file ) ) continue;
					// Test extension
					if ( fs::extension(file) != ".jpg" ) continue;
					DataItemPtr item ( new DataItem() );
					item->id = ++ item_counter;
					item->path = file.c_str();
					data_.push_back(item);
					data_class_.push_back(desc->id);
				}
			}
			std::cout << "OUT" << std::endl;

		}
	} else {
		std::cerr << "the input is not a directory" << std::endl;
		return false;
	}
	std::cout << "WORK!!!!! " << std::endl;
	// Writing descriptor
	write(descriptor_path_, data_set_path_);
	return true;
}
bool ClassificationDataSet::read ( std::string descriptor_path ) {
	data_.clear();
	classes_.clear();
	data_class_.clear();

	std::ifstream fin ( descriptor_path.c_str() );
	fin >> name_;
	fin >> data_set_path_;
	int len;

	fin >> len;
	data_.resize(len);
	for ( int i = 0; i < len; ++ i ) {
		data_[i] = DataItemPtr(new DataItem());
		fin >> data_[i]->id >> data_[i]->path;
	}
	fin >> len;
	classes_.resize(len);
	for ( int i = 0; i < len; ++ i ) {
		classes_[i] = ClassDescriptorPtr(new ClassDescriptor());
		fin >> classes_[i]->id >> classes_[i]->name;
	}
	fin >> len;
	data_class_.resize(len);
	for ( int i = 0; i < len; ++ i )
		fin >> data_class_[i];
}
bool ClassificationDataSet::write ( std::string descriptor_path, std::string data_set_path ) {
		std::ofstream fout ( descriptor_path_.c_str() );
		fs::path dpath ( descriptor_path_ );
		// get the name
		std::string filename = dpath.filename().c_str();
		for ( int i = 0; i < filename.size(); ++ i )
			if ( filename[i] == '.' )
				filename[i] = ' ';

		std::stringstream ss(filename); std::string t1, t2 = filename;
		while ( ss >> t1 ) {
		  name_ = t2;
			t2 = t1;
		}
		// writing
	 	fout << name_ << "\n";
	 	fout << data_set_path_ << "\n";
	 	int len = data_.size();
	 	fout << len << "\n";
	 	for ( int i = 0; i < len; ++ i )
	 			fout <<	data_[i]->id << " " << data_[i]->path << "\n";

	 	len = classes_.size();
	 	fout << len << "\n";
		for ( int i = 0; i < len; ++ i )
			fout <<	classes_[i]->id << " " << classes_[i]->name << "\n";

		len = data_class_.size();
		fout << len << "\n";
		for ( int i = 0; i < len; ++ i )
			fout <<	data_class_[i] << "\n";
		fout.flush();
		fout.close();
}
bool ClassificationImageSet::create ( std::string descriptor_path, std::string data_set_path ) {
	ClassificationDataSet::create(descriptor_path, data_set_path);
	//read_images();
}
bool ClassificationImageSet::read ( std::string descriptor_path ) {
	ClassificationDataSet::read(descriptor_path);
	//read_images();
}
bool ClassificationImageSet::write( std::string descriptor_path, std::string data_set_path ) {
	fs::path out_dir ( data_set_path );
	if ( !fs::is_directory ( out_dir ) )
		fs::create_directory ( out_dir );
	out_dir.normalize();
	std::vector<std::string> class_paths ( classes_.size() );
	for ( int i = 0; i < classes_.size(); ++ i ) {
		class_paths[i] = ( out_dir.string() + "/" ) + classes_[i]->name + "/";
		fs::path class_dir ( class_paths[i] );
		fs::create_directory ( class_dir );
	}
	int counter = 0; std::string str_counter;
	for ( int i = 0; i < data_class_.size(); ++ i ) {
		int class_id = data_class_[i];
		std::stringstream ss;
		ss << counter ++;
		ss >> str_counter;
		cv::Mat image = cv::imread( data_[i]->path );
		std::string path = class_paths[class_id] + classes_[class_id]->name + str_counter + ".jpg";
		data_[i]->path = path;
		cv::imwrite(path, image );
	}
	ClassificationDataSet::write( descriptor_path, data_set_path );
}
bool ClassificationImageSet::write() {
	ClassificationImageSet::write( descriptor_path_, data_set_path_ );
}
ClassificationImageSet ClassificationImageSet::clone ( std::string descriptor_path, std::string data_set_path ) {
	ClassificationImageSet output;
	output.set_classes( classes_ );
	output.set_data_class(data_class_);
	output.set_data( data_ );
	output.set_data_set_path ( data_set_path );
	output.set_descriptor_path ( descriptor_path );
	return output;
}

cv::Mat ClassificationImageSet::get_image ( int idx )  const {
	return cv::imread( data_[idx]->path );
}

void ClassificationImageSet::set_image ( const cv::Mat image, int idx ) {
	cv::imwrite ( data_[idx]->path, image );
}
cv::Mat ClassificationImageSet::getImage ( int idx )  const {
	return cv::imread( data_[idx]->path );
}

void ClassificationImageSet::setImage ( const cv::Mat image, int idx ) {
	cv::imwrite ( data_[idx]->path, image );
}



FeatureSet::FeatureSet ( int size , int dims ) {
	resize ( size, dims );
}
FeatureSet::FeatureSet () {
	resize ( 0, 0 );
}
void FeatureSet::resize ( int size, int dim ) {
	size_ = size;
	dim_ = dim;
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
		fout << i;
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
	mat_.resize(size_);
	for ( int i = 0; i < size_; ++ i ) {
		// deprecated
		int id;
		fin >> id;

		mat_[i].resize(dim_);
		for ( int j = 0; j < dim_; ++ j )
			fin >> mat_[i][j];
	}
	cv::Mat ans ( size_, dim_, CV_64F);
	for ( int i = 0; i < size_; ++ i )
		for ( int j = 0; j < dim_; ++ j )
			ans.ptr<double>(i)[j] = mat_[i][j];
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
		for ( int j = 0; j < size_; ++ j ) {
			val = get_value (j, i);
			val = ipcore::normalize ( val, min_, max_, -1, 1 );
			set_value (j, i, val );
		}
	}
}
double FeatureSet::get_value( int i, int j ) const { return mat_[i][j]; }
void FeatureSet::set_value( int i, int j, double val ) { mat_[i][j] = val; }

	/* right way */
void FeatureSet::setValue( int i, int j, double val ) {
	mat_[i][j] = val;
}
double FeatureSet::getValue( int i, int j ) const {
	return mat_[i][j];
}
void FeatureSet::setData ( std::vector<std::vector<double> > data) {
	mat_ = data;
	size_ = mat_.size();
	dim_ = mat_[0].size();
}
std::string FeatureSet::getDataClassName ( int data_id ) {
	return class_names_[ data_class_ids[data_id] ];
}
int FeatureSet::getDataClassId ( int data_id ) {
	return data_class_ids[data_id];
}
void FeatureSet::setDataClassIds( std::vector<int> some ) {
	data_class_ids = some;
}
std::vector<int> FeatureSet::getDataClassIds() {
	return data_class_ids;
}
std::string FeatureSet::getClassName ( int class_id ) {
	return class_names_[ class_id ];
}
std::vector<std::string> FeatureSet::getClassNames() {
	return class_names_;
}
std::vector<double> FeatureSet::getFeatureVector( int idx ) {
	return mat_[idx];
}
void FeatureSet::setClassNames( std::vector<std::string> names ) {
	class_names_ = names;
}
int FeatureSet::size() const { return size_; }
int FeatureSet::dim() const { return dim_; }
int FeatureSet::num_classes() const { return class_names_.size(); }
cv::Mat FeatureSet::cv_mat() const {
	cv::Mat ans ( size_, dim_, CV_64F);
	for ( int i = 0; i < size_; ++ i )
		for ( int j = 0; j < dim_; ++ j )
			ans.ptr<double>(i)[j] = mat_[i][j];
	return ans;
}
// DAOS
void AbstractImageSetDao::setImageSet ( const ImageSet& dset ) {
	images_ = dset;
}
ImageSet AbstractImageSetDao::getImageSet () const {
	return images_;
}

void SimpleImageSetDao::create ( std::string path ) {

}
void SimpleImageSetDao::read ( std::string path ) {

}
void SimpleImageSetDao::update ( std::string path ) {

}
void SimpleImageSetDao::destroy ( std::string path ) {

}

void AbstractFeatureSetDao::setFeatureSet ( const FeatureSet& fset ) {
	features_ = fset;
}
FeatureSet AbstractFeatureSetDao::getFeatureSet ( ) const {
	return features_;
}
void WekaFeatureSetDao::create ( std::string path ) {
	std::ofstream fout ( path.c_str() );
	fout << "@RELATION features" << std::endl;
	for ( int i = 0; i < features_.dim(); ++ i )
		fout << "@ATTRIBUTE " << i << " REAL" << std::endl;
	fout << "@ATTRIBUTE class {";
	for ( int i = 0; i < features_.num_classes(); ++ i )
		fout << (i?",":"") << features_.getClassName(i);
	fout << "}" << std::endl;
	fout << "@DATA" << std::endl;
	for ( int i = 0; i < features_.size(); ++ i ) {
		for ( int j = 0; j < features_.dim(); ++ j )
			fout << features_.getValue(i, j) << ",";
		fout << features_.getClassName( features_.getDataClassId(i) ) << std::endl;
	}
	fout.flush();
	fout.close();
}
void WekaFeatureSetDao::read ( std::string path ) {
	std::ifstream fin ( path.c_str() );
	std::string line, token;
	std::vector < std::string > stream;
	std::vector < std::vector <double> > mat;
	std::vector < int > labels;
	int dim = 0;

	while ( std::getline ( fin, line ) )
		stream.push_back( line );

	for ( int i = 0; i < stream.size(); ++ i ) {
		std::stringstream ss ( stream[i] );
		ss >> token;
		for ( int j = 0; j < token.size(); ++ j )
			token[j] = tolower(token[j]);
		if ( token != "@attribute" ) continue;
		++ dim;
	}
	-- dim;
	bool data = false;
	std::map< std::string, int > mlabels;
	for ( int i = 0; i < stream.size(); ++ i ) {
		if ( !data ) {
			std::stringstream ss2 ( stream[i] );
			ss2 >> token;
			for ( int j = 0; j < token.size(); ++ j )
				token[j] = tolower(token[j]);
			if ( token == "@data" )
				data = true;
			continue;
		}
		line = stream[i];
		for ( int j = 0; j < line.size(); ++ j )
			if ( line[j] == ',' )	line[j] = ' ';
		std::stringstream ss ( line );
		std::vector<double> fv ( dim );
		for ( int j = 0; j < dim; ++ j )
			ss >> fv[j];
		ss >> token;
		mlabels[token] = mlabels.size()-1;
		mat.push_back(fv);
		std::cout << mlabels[token] << std::endl;
		labels.push_back(mlabels[token]);
	}
	features_.setData(mat);
	features_.setDataClassIds(labels);
}
void WekaFeatureSetDao::update ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}
void WekaFeatureSetDao::destroy ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}

void ShogunFeatureSetDao::create ( std::string path ) {
	std::string dpath = path + ".dat";
	std::string lpath = path + ".label";
	std::ofstream foutdata ( dpath.c_str() );
	std::ofstream foutlabels ( lpath.c_str() );
	for ( int i = 0; i < features_.size(); ++ i ) {
		for ( int j = 0; j < features_.dim(); ++ j )
			foutdata << (j?" ":"") << features_.getValue(i, j);
		foutdata << "\n";
	}
	for ( int i = 0; i < features_.size(); ++ i )
		foutlabels << features_.getDataClassId(i) << "\n";
}
void ShogunFeatureSetDao::read ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}
void ShogunFeatureSetDao::update ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}
void ShogunFeatureSetDao::destroy ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}

void AbstractClassifierEvaluatorResultDao::setClassifierEvaluatorResult ( const ClassifierEvaluatorResult & eval ) {
	result_ = eval;
}
ClassifierEvaluatorResult AbstractClassifierEvaluatorResultDao::getClassifierEvaluatorResult() const {
	return result_;
}

SimpleClassifierEvaluatorResultDao::SimpleClassifierEvaluatorResultDao ( ClassifierEvaluatorResult result, ImageSet *image  ) {
	result_ = result;
	image_set_ptr_ = image;
}
void SimpleClassifierEvaluatorResultDao::create ( std::string path ) {

	int num_classes = image_set_ptr_->getClasses().size();
	std::vector<int> data_class_ids = image_set_ptr_->getDataClassIds();
	std::vector<ClassDescriptorPtr> classes = image_set_ptr_->getClasses();
	std::vector<ClassDescriptorPtr> new_classes;
	std::cout << image_set_ptr_->getDataClassIds().size() << std::endl;
	std::cout << image_set_ptr_->getClasses().size() << std::endl;
	for ( int i = 0; i < num_classes; ++ i )
		for ( int j = 0; j < num_classes; ++ j ) {
			std::cout << "ENTRA" << std::endl;
			int id = i*num_classes+j;
			std::vector<int> image_ids = result_.belong_matrix[i][j];
			std::cout << "ohh shit" << std::endl;
			new_classes.push_back( ClassDescriptorPtr ( new ClassDescriptor() ));
			new_classes[id]->id = id;
			new_classes[id]->name += "(REAL_";
			new_classes[id]->name += classes[i]->name;
			new_classes[id]->name += ",OUTCOME_";
			new_classes[id]->name += classes[j]->name;
			new_classes[id]->name += ")";
			std::cout << new_classes[id]->name << std::endl;
			for ( int k = 0; k < image_ids.size(); ++ k)
				data_class_ids[image_ids[k]] = i*num_classes + j;
		}
	std::cout << "called" << std::endl;
	ImageSet new_image_set = image_set_ptr_->clone(path,path +"_desc");
	new_image_set.setClasses(new_classes);
	new_image_set.setDataClass(data_class_ids);
	new_image_set.write();
}
void SimpleClassifierEvaluatorResultDao::read ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}
void SimpleClassifierEvaluatorResultDao::update ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}
void SimpleClassifierEvaluatorResultDao::destroy ( std::string path ) {
	std::cerr << "Not implemented yet" << std::endl;
}
