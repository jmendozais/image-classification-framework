/*
 * data.h
 *
 *  Created on: Jan 24, 2013
 *      Author: jmendoza
 */

#ifndef DATA_H_
#define DATA_H_
/****************************************
 * 				 Persistent Objects
 ****************************************/
#include <vector>
#include <string>
#include <opencv/cv.h>
#include <boost/shared_ptr.hpp>
class AbstractParams {};

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

	const std::vector<ClassDescriptorPtr>& getClasses() const {
		return classes_;
	}

	void setClasses(const std::vector<ClassDescriptorPtr>& classes) {
		classes_ = classes;
	}

	const std::vector<DataItemPtr>& getData() const {
		return data_;
	}

	void setData(const std::vector<DataItemPtr>& data) {
		data_ = data;
	}

	const std::vector<int>& getDataClass() const {
		return data_class_;
	}

	void setDataClass(const std::vector<int>& dataClass) {
		data_class_ = dataClass;
	}

	const std::string& getDataSetPath() const {
		return data_set_path_;
	}

	void setDataSetPath(const std::string& dataSetPath) {
		data_set_path_ = dataSetPath;
	}

	const std::string& getDescriptorPath() const {
		return descriptor_path_;
	}

	void setDescriptorPath(const std::string& descriptorPath) {
		descriptor_path_ = descriptorPath;
	}

	const std::string& getName() const {
		return name_;
	}

	void setName(const std::string& name) {
		name_ = name;
	}
};

class ClassificationImageSet : public ClassificationDataSet {
private:
	std::vector<cv::Mat> images_;
	bool read_images();
public:
	/* all this should be cleared */
	bool create ( std::string descriptor_path, std::string data_set_path );
	bool read ( std::string descriptor_path );
	bool write ( std::string descriptor_path, std::string data_set_path );
	bool write ();
	ClassificationImageSet clone ( std::string descriptor_path, std::string data_set_path );
	void clone ( std::string descriptor_path, std::string data_set_path, ClassificationImageSet set );
	cv::Mat get_image ( int idx ) const ;
	void set_image ( const cv::Mat image, int idx );

	/* fixed */
	cv::Mat getImage ( int idx ) const ;
	void setImage ( const cv::Mat image, int idx );
};
typedef ClassificationDataSet DataSet;
typedef ClassificationImageSet ImageSet;
typedef boost::shared_ptr<ImageSet> ImageSetPtr;

class FeatureSet {
private:
	std::string path_;

	std::vector<int> data_class_ids;
	std::vector<std::vector<double> > mat_;
	std::vector<std::string> class_names_;

	int size_;
	int num_classes_;
	int dim_;
public:
	FeatureSet ();
	FeatureSet ( int size, int dim );

	/* to be deleted */
	void write( std::string path );
	void read( std::string path );
	void change_feature_vectors ( cv::Mat fv );
	void set_value( int i, int j, double val );
	double get_value( int i, int j) const;

	/* right way */
	void resize ( int size, int dim );
	void setValue( int i, int j, double val );
	double getValue( int i, int j ) const;

	std::string getDataClassName ( int data_id );
	int getDataClassId ( int data_id );
	std::string getClassName ( int class_id );

	std::vector<std::string> getClassNames();
	std::vector<std::vector<double> > getData();
	std::vector<int> getDataClassIds();
	std::vector<double> getFeatureVector( int idx );

	void setClassNames( std::vector<std::string> );
	void setData ( std::vector<std::vector<double> > );
	void setDataClassIds( std::vector<int> );

	int size() const;
	int dim() const;
	int num_classes() const;
	cv::Mat cv_mat() const;
	std::vector< std::vector<double> > std_mat() const;

	/* methods */
	void normalize ();
	void changeFeatureVectors ( cv::Mat fv );
};
typedef boost::shared_ptr<FeatureSet> FeatureSetPtr;

struct ClassifierEvaluatorResult {
	std::vector<std::vector<int> > confussion_matrix;
	std::vector<std::vector<std::vector<int> > > belong_matrix;
	std::vector<double> accuracy;
	std::vector<double> sensibility;
	std::vector<double> specificity;
	std::vector<double> precision;
	double global_accuracy;
};


/****************************************
 * 					Data Access Objects
 ****************************************/


class AbstractImageSetDao {
public:
	void setImageSet ( const ImageSet& dset );
	ImageSet getImageSet () const ;

	virtual void create ( std::string path ) = 0;
	virtual void read ( std::string path ) = 0;
	virtual void update ( std::string path ) = 0;
	virtual void destroy ( std::string path ) = 0;

protected:
	ImageSet images_;
};

class SimpleImageSetDao {
	virtual void create ( std::string path );
	virtual void read ( std::string path );
	virtual void update ( std::string path );
	virtual void destroy ( std::string path );
};

class AbstractFeatureSetDao {
public:
	void setFeatureSet ( const FeatureSet& fset );
	FeatureSet getFeatureSet ( ) const;

	virtual void create ( std::string path ) = 0;
	virtual void read ( std::string path ) = 0;
	virtual void update ( std::string path ) = 0;
	virtual void destroy ( std::string path ) = 0;

protected:
	FeatureSet features_;
};

class WekaFeatureSetDao : public AbstractFeatureSetDao {
public:
	WekaFeatureSetDao() {}
	~WekaFeatureSetDao() {}
	virtual void create ( std::string path );
	virtual void read ( std::string path );
	virtual void update ( std::string path );
	virtual void destroy ( std::string path );
};

class ShogunFeatureSetDao : public AbstractFeatureSetDao {
public:
	virtual void create ( std::string path );
	virtual void read ( std::string path );
	virtual void update ( std::string path );
	virtual void destroy ( std::string path );
};

class AbstractClassifierEvaluatorResultDao {
public:
	void setClassifierEvaluatorResult ( const ClassifierEvaluatorResult & eval );
	ClassifierEvaluatorResult getClassifierEvaluatorResult() const;

	virtual void create ( std::string path ) = 0;
	virtual void read ( std::string path ) = 0;
	virtual void update ( std::string path ) = 0;
	virtual void destroy ( std::string path ) = 0;

protected:
	ClassifierEvaluatorResult result_;
	ImageSet *image_set_ptr_;
};

class SimpleClassifierEvaluatorResultDao : public AbstractClassifierEvaluatorResultDao {
public:
	SimpleClassifierEvaluatorResultDao ( ClassifierEvaluatorResult result, ImageSet *image );
	virtual void create ( std::string path );
	virtual void read ( std::string path );
	virtual void update ( std::string path );
	virtual void destroy ( std::string path );
};


#endif /* DATA_H_ */
