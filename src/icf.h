/*
 * icf.h
 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */

#ifndef ICF_H_
#define ICF_H_

// Image Classification Framework default headers
#include "classification/classification_model.h"
#include "feature_extraction/gabor_wavelet.h"
#include "feature_extraction/gabor_pca.h"
#include "feature_extraction/cdh.h"
#include "feature_selection/mrMr.h"
#include "feature_selection/abc_selector.h"
#include "classification/knn.h"
#include "validation.h"
#include "image_transformation/image_transform.h"
#include "new/batch.h"
#include "new/data.h"
#include "new/transform.h"
#include "visualization.h"

// SOME FUNCTIONS SHOULD BE REALLOCATED OR CLEANED
void write_feature_vectors_in_libsvm_format ( ClassificationImageSet data, FeatureSet features, std::string output_path ) {
	std::ofstream fout ( output_path.c_str() );
	std::cout << "DATA SIZE " << data.size() << std::endl;
	for ( int i = 0; i < data.size(); ++ i ) {
		fout << data.get_class_id(i);
		for ( int j = 0; j < features.dim(); ++ j ) {
			fout << " " << (j+1) << ":" << features.get_value(i,j);

		}
		fout << std::endl;
	}
	fout.flush();
	fout.close();
}
void read_feature_vectors_in_libsvm_format ( ClassificationImageSet &data, FeatureSet &features, std::string input_path ) {
	std::vector<DataItemPtr> items;
	std::vector<ClassDescriptorPtr> classes;
	std::vector<int> class_ids;
	std::ifstream fin ( input_path.c_str());
	string line;
	int size = 0, class_id;
	std::set<int> class_set;
	std::vector<string> lines;
	while ( getline ( fin, line ) ) {
		lines.push_back(line);
		size ++;
		stringstream ss ( line );
		ss >> class_id;
		class_set.insert( class_id );

		DataItemPtr item ( new DataItem() ); item->id = size;
		items.push_back(item);
		class_ids.push_back(class_id);
	}
	cout << items.size() << " " << class_ids.size() << " " << class_set.size() << std::endl;
	for ( std::set<int>::iterator it = class_set.begin(); it != class_set.end(); ++ it ) {
		ClassDescriptorPtr class_ ( new ClassDescriptor() );
		class_->id = *it;
		class_->name = "unknown" ;
		classes.push_back(class_);
	}
	data.set_classes(classes);
	data.set_data(items);
	data.set_data_class(class_ids);

	int dim = 0;
	std::stringstream ss ( lines[0] );
	std::string tok;
	//while ( ss >> tok ) ++ dim;
	//	-- dim;
	dim = 623;
	features = FeatureSet ( size, dim );
	std::vector<string> class_names;
	for ( int i = 0; i < classes.size(); ++ i)
		class_names.push_back("" + (char)('A' + i));
	features.setClassNames(class_names);
	for ( int i = 0; i < lines.size(); ++ i ) {
		std::stringstream ss2 ( lines[i] );
		ss2 >> tok;
		int idx = 0;
		//std::cout << "LINES " << i << " " << lines[i] << endl;
		while ( ss2 >> tok ) {
			std::string feature = tok.substr( tok.find(':') + 1, tok.size() ).c_str();
			std::stringstream ss3;
			ss3 << feature;
			double value;
			ss3 >> value;
			//std::cout << value << ", ";
			features.set_value( i, idx++, value );
		}
		//std::cout << std::endl;
	}
}
void write_classifier_evaluation_report ( const ClassifierEvaluatorResult &result, const ClassificationDataSet &data, string out_name ) {
// For color
std::vector<int> class_len( data.num_classes() );
for ( int i = 0; i < result.confussion_matrix.size(); ++ i ) {
	for ( int j = 0; j < result.confussion_matrix.size(); ++ j )
		class_len[i] += result.confussion_matrix[i][j];
}

// PRINT IN LATEX FORMAT
string tname = out_name + ".tex";
ofstream out ( tname.c_str() );
out << "\\documentclass[a4paper,10pt]{article}" << endl;
out << "\\usepackage[utf8x]{inputenc}\\usepackage{default}\\usepackage{multirow}\\usepackage[table]{xcolor}" << endl;
out << "\\newcommand{\\mr}[3]{\\multirow{#1}{#2}{#3}}\\newcommand{\\mc}[3]{\\multicolumn{#1}{#2}{#3}}" << endl;
out << "\\begin{document}" << endl;
out << "Exactitud global ( sum(tp)/dataSize ) " << result.global_accuracy << endl;
out << "\\section{Matrix de confusion}" << endl;
out << "\\begin{center}" << endl;

out << "\\begin{tabular}{c|";
int nclases = result.accuracy.size(), val;
for ( int i = 1; i < nclases + 2; ++ i )
	out << "c";
out << "}";
out << " &  & \\mc{" << nclases << "}{c}{Prediccion} \\\\" << endl;
out << "\\hline" << endl;
out << " &  ";
for ( int i = 0; i < nclases; ++ i )
	out << "& " << data.classes()[i]->name.substr(0,2);
out << "\\\\\\cline{" << 3 << "-" << nclases + 2 << "}" << endl;
out << "\\mr{" << nclases << "}{*}{Real}& \\mc{1}{l|}{ " << data.classes()[0]->name << "} ";
for ( int i = 0; i < nclases; ++ i ) {
	val = result.confussion_matrix[0][i];
	double color = abs( 0.5 + 0.5 * (1-val*1.0/class_len[0]) );
	assert ( color > 0 );
	std::cout << val << " " << class_len[0] << std::endl;
	out << "& \\mc{1}{c|}{\\cellcolor[rgb]{ " << color << ", " << color << ",1} " << val << "}";
}
out << "\\\\\\cline{" << 3 << "-" << nclases + 2 << "}" << endl;
for ( int i = 1; i < nclases; ++ i ) {
	out << " & \\mc{1}{l|}{ " << data.classes()[i]->name << "}";
	for ( int j = 0; j < nclases; ++ j ) {
		val = result.confussion_matrix[i][j];
		double color = abs( 0.5 + 0.5 * (1-val*1.0/class_len[i]) );
		out << "& \\mc{1}{c|}{\\cellcolor[rgb]{ " << color << ", " << color << ",1} " << result.confussion_matrix[i][j] << " }";
	}
out << "\\\\\\cline{" << 3 << "-" << nclases + 2 << "}" << endl;
}
out << "\\end{tabular}\\end{center}" << endl;

out << "\\section{Metricas de evaluacion del clasificador}" << endl;
out << "\\begin{center}\\begin{tabular}{";
for ( int i = 0; i < 5; ++ i )
	out << "c";
out << "}\\hline" << endl;
out << fixed << setprecision(4);
out << "Clase & Sensibilidad & Especificidad & Precision & Exactitud \\\\" << endl;
out << "\\hline" << endl;
for ( int i = 0; i < nclases; ++ i ) {
	double sen, spec, pre, acc;
	sen = 0.6 + 0.4 * ( 1 - result.sensibility[i] );
	spec = 0.6 + 0.4 * ( 1 - result.specificity[i] );
	pre = 0.6 + 0.4 * ( 1 - result.precision[i] );
	acc = 0.6 + 0.4 * ( 1 - result.accuracy[i] );

	out << data.classes()[i]->name << " & \\cellcolor[rgb]{"<< sen << ", " << sen << ",1}" << result.sensibility[i] << " & \\cellcolor[rgb]{" << spec << ", " << spec << ", " << "1}" << result.specificity[i] << " & \\cellcolor[rgb]{"<< pre << ", " << pre << ",1}" << result.precision[i] << " & \\cellcolor[rgb]{ 1, "<< acc << ", " << acc << "}"  << result.accuracy[i] << "\\\\" << std::endl;
}
out << "\\hline\\end{tabular}\\end{center} \\end{document}" << endl;
out.flush();
out.close();
// SHOW
string com1 = "pdflatex " + out_name + ".tex";
system ( com1.c_str() );
}
void showImGrid ( std::vector<cv::Mat> images, int rows, int cols ) {
	// switch to double mats
	for ( int i = 0; i < images.size(); ++ i )
		images[i] = cv::Mat_<double>(images[i]);

	// resize
	int max_rows, max_cols;
	max_rows = max_cols = -1;
	int dimx = 1000/cols, dimy = 700/rows;
	for ( int i = 0; i < images.size(); ++ i ) {
		cv::resize(images[i], images[i], cv::Size(images[i].rows, images[i].cols ) );
	}

	for ( int i = 0; i < images.size(); ++ i ) {
		max_rows = max ( max_rows, images[i].rows );
		max_cols = max ( max_cols, images[i].cols );
	}
	cv::Mat resp ( rows * max_rows, cols * max_cols, CV_64FC1 );
	for ( int i = 0; i < rows; ++ i )
		for ( int j = 0; j < cols; ++ j ) {
			if ( i*cols + j >= images.size() ) continue;
			cv::Rect roi ( j*max_cols, i*max_rows, max_cols, max_rows );
			cv::Mat image_roi ( resp, roi );
			ipcore::normalize(images[i*cols+j]).copyTo(image_roi);
		}
	stringstream ss;
	ss << images.size();
	string tok;
	ss >> tok;
	cv::namedWindow("Mosaic" + tok );
	cv::imshow("Mosaic" + tok, resp);
	//cv::waitKey();
}
#endif /* ICF_H_ */
