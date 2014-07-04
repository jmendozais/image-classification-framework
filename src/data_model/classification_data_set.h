/*
 * classification_data_set.h
 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */

#ifndef CLASSIFICATION_DATA_SET_H_
#define CLASSIFICATION_DATA_SET_H_

#include "../new/data.h"
/*
#include <string>
#include <vector>
#include <opencv/cv.h>
#include <boost/shared_ptr.hpp>
struct DataItem {
	int id;
	std::string path;
};
struct ClassDescriptor {
	int id;
	std::string name;
};
typedef boost::shared_ptr<DataItem> DataItemPtr;
typedef boost::shared_ptr<ClassDescriptor>  ClassDescriptorPtr;

class ClassificationDataSet {
protected:
	std::vector<DataItemPtr> data_;
	std::vector<int> data_class_;
	std::vector<ClassDescriptorPtr> classes_;
	std::string data_set_path_;
	std::string descriptor_path_;
	std::string name_;

	bool write_descriptor();
	std::string generateDataPath ( int idx );
public:
	ClassificationDataSet() {}
	ClassificationDataSet( std::vector<DataItemPtr> data, std::vector<int> data_class,	std::vector<ClassDescriptorPtr> classes,
												 std::string data_set_path, std::string descriptor_path );
	bool create ( std::string descriptor_path, std::string data_set_path );
	bool read ( std::string descriptor_path );
	bool write ( std::string descriptor_path, std::string data_set_path );

	// Accessor & Mutators
	int size() const { return data_.size(); }
	int num_classes() const { return classes_.size(); }
	std::vector<DataItemPtr> data() const { return data_; }
	std::vector<ClassDescriptorPtr> classes() const { return classes_; }
	std::vector<std::string> getClassNames() const {
		std::vector<std::string> resp ( classes_.size() );
		for ( int i = 0; i < classes_.size(); ++ i )
			resp[i] = classes_[i]->name;
		return resp;
	}
	std::vector<int> std_data_classes () const { return data_class_; }
	void set_data_set_path ( std::string p ) { data_set_path_ = p;	}
	void set_descriptor_path ( std::string p ) { descriptor_path_ = p; }
	void set_data ( std::vector<DataItemPtr> data ) { data_ = data; };
	void set_classes ( std::vector<ClassDescriptorPtr> classes ) { classes_ = classes; }
	void set_data_class( std::vector<int> class_ids ) { data_class_ = class_ids; }
	int get_class_id ( int idx ) const { return data_class_[idx]; }
	void set_class_id ( int idx, int class_id ) { data_class_[idx] = class_id; }
	std::vector<int> getDataClassIds() const { return data_class_; }
};

class ClassificationImageSet : public ClassificationDataSet {
private:
	std::vector<cv::Mat> images_;
	bool read_images();
public:
	bool create ( std::string descriptor_path, std::string data_set_path );
	bool read ( std::string descriptor_path );
	bool write ( std::string descriptor_path, std::string data_set_path );
	bool write ();
	ClassificationImageSet clone ( std::string descriptor_path, std::string data_set_path );
	void clone ( std::string descriptor_path, std::string data_set_path, ClassificationImageSet set );
	cv::Mat get_image ( int idx ) const ;
	void set_image ( const cv::Mat image, int idx );
};

*/
#endif /* CLASSIFICATION_DATA_SET_H_ */
