/*
 * classification_data_set.cpp
 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */

/*
#include "classification_data_set.h"
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
		for ( fs::directory_iterator it = fs::directory_iterator ( root ); it != fs::directory_iterator(); ++ it ) {
			if ( fs::is_directory ( *it ) ) {
				fs::path dir = (*it).path();
				ClassDescriptorPtr desc ( new ClassDescriptor() );
				desc->id = ++ class_counter;
				desc->name = dir.filename().string();
				classes_.push_back(desc);
				for ( fs::directory_iterator it2 = fs::directory_iterator ( dir ); it2 != fs::directory_iterator(); ++ it2 ) {
					fs::path file = (*it2).path();
					if ( !fs::is_regular_file ( file ) ) continue;
					DataItemPtr item ( new DataItem() );
					item->id = ++ item_counter;
					item->path = file.c_str();
					data_.push_back(item);
					data_class_.push_back(desc->id);
				}
			}
		}
	} else
		return false;

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
	output.write( descriptor_path, data_set_path );
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
*/
