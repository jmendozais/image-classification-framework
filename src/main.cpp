#include <iostream>
#include "icf.h"

void test_classification_image_set() {
	ClassificationImageSet data_set;
	data_set.create("/home/jmendoza/Escritorio/coffee.data", "/home/jmendoza/projects/cecovasa/data/LotePruebaV1.0/Resultado");
	data_set.read("/home/jmendoza/Escritorio/coffee.data");
	data_set.write("/home/jmendoza/Escritorio/coffee2.data", "/home/jmendoza/Escritorio/coffee2" );
	std::cout << "test_classification_image_set() ... OK" << std::endl;
}
void test_small_data () {
	ClassificationImageSet data_set;
	data_set.create("/home/jmendoza/Escritorio/small.data", "/home/jmendoza/projects/cecovasa/data/small_data_set");
	data_set.read("/home/jmendoza/Escritorio/small.data");
}
void test_feature_set() {
	ClassificationImageSet data_set;
	data_set.read("/home/jmendoza/Escritorio/coffee.data");
	std::cout << "Image set size " << data_set.size() << std::endl;
	FeatureSet features;
	features.read("/home/jmendoza/Escritorio/coffee.feat");
	std::cout << "Feature set size & dim " << features.size() << " " << features.dim() << std::endl;
	features.normalize();
	features.write("/home/jmendoza/Escritorio/coffee_norm.feat");
	write_feature_vectors_in_libsvm_format  ( data_set, features, "/home/jmendoza/Escritorio/coffee_norm.libsvm" );
}
void test_feature_selection() {
	 ClassificationImageSet data_set;
	 data_set.read("/home/jmendoza/Escritorio/coffee.data");
	 FeatureSet features, selected_features;
	 features.read("/home/jmendoza/Escritorio/coffee.feat");

	 MRMR_params params;
	 params.data_set = data_set;
	 params.target_dim = 1;
	 GreedyMRMRSelector selector ( params );

	 selected_features = selector.selectFeatures( features );
	 selected_features.write("/home/jmendoza/Escritorio/coffee_1.feat");
	 write_feature_vectors_in_libsvm_format(data_set, selected_features,
			 	 	 	 	 	 	 	 	 	 	 "/home/jmendoza/Escritorio/coffee_1.libsvm");
}
void test_feature_selection_2() {

	 ClassificationImageSet data_set;
	 FeatureSet features, selected_features;
	 read_feature_vectors_in_libsvm_format ( data_set, features, "/home/jmendoza/Desktop/iris.scale");

	 write_feature_vectors_in_libsvm_format(data_set, features,
			 	 	 	 	 	 	 	 	 	 	 "/home/jmendoza/Desktop/iris.libsvm");
	 MRMR_params params;
	 params.data_set = data_set;
	 params.target_dim = 2;
	 GreedyMRMRSelector selector ( params );

	 selected_features = selector.selectFeatures( features );
	 write_feature_vectors_in_libsvm_format(data_set, selected_features,
			 	 	 	 	 	 	 	 	 	 	 "/home/jmendoza/Desktop/iris2.libsvm");
}
void test_feature_selection_3() {

	 ClassificationImageSet data_set;
	 FeatureSet features, selected_features;
	 read_feature_vectors_in_libsvm_format ( data_set, features, "/home/jmendoza/Desktop/brocado_all");

	 //write_feature_vectors_in_libsvm_format(data_set, features,
			// 	 	 	 	 	 	 	 	 	 	 "/home/jmendoza/Desktop/coffe54.libsm");
	 MRMR_params params;
	 params.data_set = data_set;
	 params.target_dim = 115	;
	 GreedyMRMRSelector selector ( params );

	 selected_features = selector.selectFeatures( features );
	 write_feature_vectors_in_libsvm_format(data_set, selected_features,
			 	 	 	 	 	 	 	 	 	 	 "/home/jmendoza/Desktop/brocado_all_s_115");
}
void test_classification_model() {
	std::string fname = "/Users/mac/research/databases/cecovasa_76/gabor_45";
	ClassificationModel model;

	// Setting gabor wavelets params
	GaborWaveletParams params;
	params.num_orientations = 5;
	params.num_scales = 4;
	params.mask_type = GaborWaveletParams::NO_MASK;
	GaborWavelet gabor_wavelet ( params );
	model.set_feature_extractor( &gabor_wavelet);

	// Classification pipeline
	ClassificationImageSet data_set;
	data_set.read("/Users/mac/research/databases/cecovasa_76/cecovasa.desc_cpp");

	FeatureSet features = model.extractFeatures( data_set );
	features.normalize();
	features.write(fname + ".feat");
	write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".libsvm" );
	features.normalize();
	features.write( fname + ".nfeat");
	write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".nlibsvm" );
	std::cout <<"Test Classification Model... OK" << std::endl;
}

void test_brocados() {
	ClassificationImageSet data;
	//data.create("/home/jmendoza/Desktop/brocado_all.data", "/home/jmendoza/projects/cecovasa/data/lote1.2/Brocados");
	data.read("/home/jmendoza/Desktop/brocado_all.data");
	ClassificationModel model;
	GaborWaveletParams params;
	params.num_orientations = 4;
	params.num_scales = 5;
	params.mask_type = GaborWaveletParams::NO_MASK;
	GaborWavelet extractor ( params );
	model.set_feature_extractor(&extractor);
	FeatureSet features = model.extractFeatures(data);
	features.write("/home/jmendoza/Desktop/brocado_all_45_nomask.feat");
	write_feature_vectors_in_libsvm_format(data,features,"/home/jmendoza/Desktop/brocado_all_45_nomask.libsvm");
}
void test_brocados_2() {
	ClassificationImageSet data_set;
		 data_set.read("/home/jmendoza/Desktop/brocado_all.data");
		 FeatureSet features, selected_features;
		 features.read("/home/jmendoza/Desktop/brocado_all_45_nomask.feat");
		 for ( int i = 80; i < 120; ++ i ) {
			 MRMR_params params;
			 params.data_set = data_set;
			 params.target_dim = i;
			 GreedyMRMRSelector selector ( params );

			 selected_features = selector.selectFeatures( features );
			 stringstream ss;
			 string name = "/home/jmendoza/Desktop/ball.libsvm";
			 ss << name;
			 ss << i;
			 ss >> name;
			 write_feature_vectors_in_libsvm_format(data_set, selected_features, name);
		 }
}

void test_knn () {
	ClassificationImageSet data_set;
	data_set.read("/home/jmendoza/Desktop/brocado_all.data");
	FeatureSet features, selected_features;
	features.read("/home/jmendoza/Desktop/brocado_all_45_nomask.feat");
	std::vector< std::vector<double> > fv = features.std_mat();
	std::vector< int > classes ( data_set.num_classes() );
	KNNClassifier knn;
	knn.train( fv, classes );
	int outcome = knn.predict( fv[350] );
}
void test_all() {
	ClassificationImageSet data;
	data.create("/Users/mac/research/databases/cecovasa_76/cecovasa.desc_cpp", "/Users/mac/research/databases/cecovasa_76");
	/*data.read("/Users/mac/research/projects/coffee1.data");
	ClassificationModel model;
	GaborWaveletParams params;
	params.num_orientations = 4;
	params.num_scales = 3;
	params.mask_type = GaborWaveletParams::NO_MASK;
	GaborWavelet extractor ( params );

	model.set_feature_extractor(&extractor);

	FeatureSet features = model.extractFeatures(data);
	features.write("/home/jmendoza/Desktop/coffee_nomask.feat");
	write_feature_vectors_in_libsvm_format(data, features, "/home/jmendoza/Desktop/coffee_nomask.libsvm");
	*/
}
void test_loocv() {
	ClassificationImageSet data_set;
	data_set.read("/Users/mac/research/databases/cecovasa_76/cecovasa.desc_cpp");
	FeatureSet features, selected_features;
	features.read("/home/jmendoza/Desktop/coffee_mask_gw_5_10.feat");
	std::vector< std::vector<double> > fv = features.std_mat();
	std::vector< int > classes = data_set.std_data_classes();
	KNNClassifier knn(10);
	knn.train( fv, classes );
	KFOLDCV evaluator ( &knn, &data_set, &features, 10);
	ClassifierEvaluatorResult eval = evaluator.evaluate();
	std::cout << "Global Accuracy " << eval.global_accuracy << std::endl;
	write_classifier_evaluation_report( eval, data_set, "/home/jmendoza/Desktop/coffee_mask_gw_5_10" );
}
void cdh() {
	ClassificationImageSet data_set;
	data_set.read("/home/jmendoza/Desktop/coffee.data");
	FeatureSet features;
	features.read("/home/jmendoza/Desktop/coffee_cdh.feat");
	std::vector< std::vector<double> > fv = features.std_mat();
	std::vector< int > classes = data_set.std_data_classes();
	KNNClassifier knn(10);
	//knn.train( fv, classes );
	KFOLDCV evaluator ( &knn, &data_set, &features, 10 );
	ClassifierEvaluatorResult eval = evaluator.evaluate();
	std::cout << "Global Accuracy " << eval.global_accuracy << std::endl;
	write_classifier_evaluation_report( eval, data_set, "/home/jmendoza/Desktop/report_cdh_kfold" );
}
void feat_cdh() {
	ClassificationImageSet data_set;
	data_set.read("/Users/mac/research/databases/cecovasa_76/cecovasa.desc_cpp");
	ClassificationModel model;
	// Setting gabor wavelets params
	CDHFeatureExtractor feature_extractor;
	model.set_feature_extractor( &feature_extractor);
	string fname = "/Users/mac/research/databases/cecovasa_76/cdh_cpp";
	std::cout << fname << std::endl;
	FeatureSet features = model.extractFeatures( data_set );
	features.write( fname + ".feat" );
	write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".libsvm" );
	features.normalize();
	features.write( fname + ".nfeat");
	write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".nlibsvm" );
}
void test_rayner() {
	ClassificationImageSet data_set;
	FeatureSet features, selected_features;
	read_feature_vectors_in_libsvm_format ( data_set, features, "/home/jmendoza/Desktop/features.dat");
	std::vector< std::vector<double> > fv = features.std_mat();
	std::vector< int > classes = data_set.std_data_classes();
	KNNClassifier knn(10);
	//knn.train( fv, classes );
	KFOLDCV evaluator ( &knn, &data_set, &features, 10 );
	ClassifierEvaluatorResult eval = evaluator.evaluate();
	std::cout << "Global Accuracy " << eval.global_accuracy << std::endl;
	//write_classifier_evaluation_report( eval, data_set, "/home/jmendoza/Desktop/report_cdh_kfold" );
}
void test_feature_extractor() {
	// Classification pipeline
	ClassificationImageSet data_set;
	data_set.read("/home/jmendoza/Desktop/coffee.data");
	for ( int i = 3; i < 4; ++ i )
		for ( int j = 5; j < 6; ++ j ) {
			ClassificationModel model;
			// Setting gabor wavelets params
			GaborWaveletParams params;
			params.num_scales = i;
			params.num_orientations = j;
			params.mask_type = GaborWaveletParams::BACKGROUND_REMOVAL_MASK;
			GaborWavelet gabor_wavelet ( params );
			model.set_feature_extractor( &gabor_wavelet);
			stringstream ss;
			ss << i;
			ss << "_";
			ss << j;
			string tok;
			ss >> tok;
			string fname = "/home/jmendoza/Desktop/coffee_mask_gw_" + tok;
			std::cout << fname << std::endl;
			FeatureSet features = model.extractFeatures( data_set );
			features.write( fname + ".feat" );
			write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".libsvm" );
			features.normalize();
			features.write( fname + ".nfeat");
			write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".nlibsvm" );
		}
}
void test_abc() {
	ClassificationImageSet data_set;
	data_set.read("/Users/mac/research/databases/cecovasa_76/cecovasa.desc_cpp");
	FeatureSet features, selected_features;
	read_feature_vectors_in_libsvm_format(data_set, features, "/Users/mac/research/databases/cecovasa_76/cecovasa_fussed_1.libsvm");
	std::vector< std::vector<double> > fv = features.std_mat();
	std::vector< int > classes = data_set.std_data_classes();
	KNNClassifier knn(5);
	knn.train( fv, classes );
	KFOLDCV evaluator ( &knn, &data_set, &features, 5);
	ClassifierFitnessFunction fitness ( &knn, &evaluator, fv, classes );
	//ClassifierEvaluatorResult eval = evaluator.evaluate();
	ABCWrapperSelector selector ( &knn, &evaluator);
	std::vector<std::vector<double> > sfv = selector.selectFeatures(fv,classes);
	//std::cout << "Global Accuracy " << eval.global_accuracy << std::endl;
}
void test_image_transform() {
	ClassificationImageSet data_set, clone_set;
	data_set.read("/home/jmendoza/Desktop/coffee.data");
	clone_set = data_set.clone("/home/jmendoza/Desktop/coffee_clone.data", "/home/jmendoza/Desktop/coffee_clone");
	//clone_set.read("/home/jmendoza/Desktop/coffee_clone.data");
	clone_set.write();
	ClassificationModel model;
	OrientationNormalization transform;
	model.set_image_transform( &transform );
	model.transformImages( data_set, clone_set );
	clone_set.write();
}
/* *************************************************
 *  									EASY
 ***************************************************/
void kfold( std::string dat, std::string feat ) {
	ClassificationImageSet data_set;
	data_set.read(dat);
	FeatureSet features, selected_features;
	features.read(feat);
	features.setClassNames(data_set.getClassNames());
	features.setDataClassIds(data_set.getDataClassIds());

	std::vector< std::vector<double> > fv = features.std_mat();
	std::vector< int > classes = data_set.std_data_classes();
	KNNClassifier knn(21);
	knn.train( fv, classes );
	KFOLDCV evaluator ( &knn, &data_set, &features , 10);
	ClassifierEvaluatorResult eval = evaluator.evaluate();
	//write_classifier_evaluation_report( eval, data_set, feat );
	SimpleClassifierEvaluatorResultDao eval_saver ( eval, &data_set );
	//eval_saver.create( feat + "_kfold" );
	write_classifier_evaluation_report(eval, data_set, feat+"_tex");
}
void test_weka() {
	ClassificationImageSet data_set;
	data_set.read("/home/jmendoza/Desktop/coffee.data");
	FeatureSet f1;
	f1.read("/home/jmendoza/Desktop/coffee_mask_gw_3_5.nfeat");
	f1.setClassNames(data_set.getClassNames());
	f1.setDataClassIds(data_set.getDataClassIds());
	WekaFeatureSetDao dao;
	dao.setFeatureSet(f1);
	dao.create("/home/jmendoza/Desktop/coffee_mask_gw_3_5.arff");
}
void test_gabor( std::string path ) {
	cv::Mat image = cv::imread ( path );
	cv::namedWindow ( "original" );
	cv::imshow ("original", image);
	cvtColor ( image, image, CV_BGR2HSV);
	std::vector<cv::Mat> v_image;
	cv::split ( image, v_image );
	for ( int i = 0; i < v_image.size(); ++ i ) {
		v_image[i] = (cv::Mat) ( cv::Mat_<double>(v_image[i]) );
		v_image[i] = ipcore::normalize(v_image[i]);
	}
	showImGrid ( v_image, 3, 1 );

	GaborWaveletParams params;
	params.num_orientations = 5;
	params.num_scales = 3;
	GaborWavelet wavelets ( params );
	vector<cv::Mat> imgs = wavelets.convolve(  image );
	showImGrid( imgs, 3, 15 );
}
void test_gabor_pca( ) {

	ClassificationImageSet data_set;
	data_set.read("/home/jmendoza/Escritorio/small.data");

	GaborWaveletParams params;
	params.num_orientations = 2;
	params.num_scales = 2;
	params.high_frequency_bound = 0.9;
	params.low_frequency_bound = 0.1;
	GaborWavelet wavelets ( params );

	GaborPCA pca ( wavelets, data_set, 10, 10, 20 );

	ClassificationModel model;
	model.set_feature_extractor(&pca);
	FeatureSet f1 = model.extractFeatures( data_set );
	f1.setClassNames(data_set.getClassNames());
	f1.setDataClassIds(data_set.getDataClassIds());
	f1.normalize();
	WekaFeatureSetDao dao;
	dao.setFeatureSet(f1);
	dao.create("/home/jmendoza/Desktop/gabor_pca.arff");
}
void gabor_FTW ( std::string path ) {
/*
	cv::Mat sample = cv::imread(path);
	cv::resize(sample, sample, cv::Size(1200,1200));
	cv::namedWindow("Original Image");
	cv::imshow("Original Image", sample);
	cv::waitKey();*/
	std::vector<GaborWaveletParams> params ( 3 );
	std::string desc[3] = {"gw_3_5_TUNEADO_entr_mag", "gw_3_5_TUNEADO_entr_real", ""};
	params[0].num_scales = 3;
	params[0].num_orientations = 5;
	params[0].remove_dc = true;
	params[0].high_frequency_bound = 0.125;
	params[0].low_frequency_bound = 0.0175;
	params[0].metrics |= GaborWaveletParams::ENTROPY;
	params[0].response = GaborWaveletParams::MAGNITUDE;

	params[1].num_scales = 3;
	params[1].num_orientations = 5;
	params[1].remove_dc = false;
	params[1].high_frequency_bound = 0.25;
	params[1].low_frequency_bound = 0.0175;
	params[1].metrics |= GaborWaveletParams::ENTROPY;
	params[1].response = GaborWaveletParams::REAL;

	for ( int i = 0; i < 2; ++ i ) {
		ClassificationImageSet data_set;
		data_set.read("/Users/mac/research/databases/cecovasa_76/cecovasa.desc_cpp");
		ClassificationModel model;
		GaborWavelet gabor_wavelet ( params[i] );
		//gabor_wavelet.extractFeatures ( sample );
		model.set_feature_extractor( &gabor_wavelet);
		string fname = "/Users/mac/research/databases/cecovasa_76/gw_entropy_" + desc[i];
		std::cout << fname << std::endl;
		FeatureSet features = model.extractFeatures( data_set );
		features.write( fname + ".feat" );
		write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".libsvm" );
		features.normalize();
		features.write( fname + ".nfeat");
		write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".nlibsvm" );
		// cambio
		/*
		features.setClassNames(data_set.getClassNames());
		features.setDataClassIds(data_set.getDataClassIds());
		WekaFeatureSetDao dao;
		dao.setFeatureSet(features);
		dao.create(fname + ".arff");
		*/
	}
}
void test_gabor_visualizer( string path ) {
	GaborWaveletParams params;
	params.num_scales = 3;
	params.num_orientations = 5;
	params.kernel_side = 60;
	params.remove_dc = true;
	params.high_frequency_bound = 0.3;
	params.low_frequency_bound = 0.1;
	params.metrics |= GaborWaveletParams::ENTROPY;
	params.response = GaborWaveletParams::MAGNITUDE;
	GaborWavelet wavelets ( params );
	//cv::Mat profile = wavelets.getFrequencyProfile();

	GaborWaveletVisualizer visualizer ( wavelets );
	cv::Mat image = cv::imread(path);
	//cv::resize(image, image, cv::Size(512,512));
	cv::namedWindow("Original");
	cv::imshow("Original", image);
	cv::waitKey();
	visualizer.setImage( image );
	visualizer.show();

	while ( true ) {
		//GaborClustererVisualizer clusterer(visualizer.getGaborWavelets());
		//clusterer.clusterize(image, 2);
		cv::waitKey();
	}


}
void show( string window_name, cv::Mat image) {
	cv::Mat temp = image.clone();
	cv::resize(temp, temp, cv::Size(600, 600) );
	cv::imshow(window_name, temp);
}


/// COMPUCORP!!!!!!!!!!!!!!!!!!!!!!!!!!
// canny
cv::Mat canny, closed, gaussian, gray, temp, laplacian, sobel, otsu, bilateral, binary;

void update () {
	cv::threshold(bilateral, otsu, 100, 255, THRESH_OTSU );
	show("otsu", otsu);
	show("closing", closed);
	show("canny", canny);
	show("bilateral", bilateral);
}

int dist;
void track_bilateral(int pos, void *data) {
	if ( dist % 2 == 0 ) return;
	cv::Mat *image = (cv::Mat*) data;
	cv::bilateralFilter ( *image, bilateral, dist, dist*2, dist/2 );
	show("bilateral", bilateral);
	//update();
}
int t1, t2;
void track_canny(int pos, void *data) {
	cv::Mat * image = (cv::Mat*) data;
	cv::Canny(*image, canny, t1,t2);
	show("canny", canny);
}
// closing
int w, h;
void track_closing(int pos, void *data) {
	cv::Mat *image = (cv::Mat *) data;
	Mat element = getStructuringElement( 2, Size( 2*w + 1, 2*h+1 ), Point( w, h ) );
	cv::morphologyEx ( *image, closed, MORPH_CLOSE, element );
	show("closing", closed);
}
int thold1;
void track_thold1(int thold1, void *data) {
	cv::Mat *image = (cv::Mat *) data;
	cv::threshold( *image, binary, (double) thold1, 255, CV_THRESH_BINARY_INV);
	show("thold1", binary);
}
int compucorp() {
	cv::Mat image = cv::imread("/home/jmendoza/projects/compucorp/data/sample.jpg");
	cv::resize(image, image, cv::Size(1000,1000));
	cv::cvtColor(image, image, CV_BGR2GRAY);
	cv::Mat_<double> dimage ( image );
	dimage = ipcore::normalize(dimage,0,255);
	dimage.convertTo( gray, CV_8UC1 );


	dist = 1;

	cv::namedWindow("bilateral");
	cv::createTrackbar("d", "bilateral", &dist, 100, track_bilateral, (void*) &gray );

	/*
	cv::namedWindow("canny");
	t1 = 0, t2 = 0;
	cv::createTrackbar("t1", "canny", &t1, 255, track_canny, (void*) &bilateral );
	cv::createTrackbar("t2", "canny", &t2, 255, track_canny, (void*) &bilateral );
*/
	thold1 = 1;
	cv::namedWindow("thold1");
	cv::createTrackbar("thold1", "thold1", &thold1, 255, track_thold1, (void*) &bilateral );

	cv::namedWindow("closing");
	w = 5, h = 5;
	cv::createTrackbar("w", "closing", &w, 100, track_closing, (void*) &binary);
	cv::createTrackbar("h", "closing", &h, 100, track_closing, (void*) &binary);

	cv::bilateralFilter ( gray, bilateral, dist, dist*2, dist/2 );
	cv::GaussianBlur(gray, gaussian, cv::Size(15,15),2,4);
	cv::Canny(bilateral, canny, t1,t2);
	cv::threshold(bilateral, otsu, 100, 255, THRESH_OTSU );
	cv::waitKey();
}
int mit ( ) {
	std::string data = "/home/jmendoza/projects/image_databases/mit";
	std::string sdesc = data + ".data";
	ClassificationImageSet set; set.create(sdesc, data);

	std::vector<GaborWaveletParams> params ( 3 );
	std::string desc[3] = {"pgabor_mag", "gw_3_5_TUNEADO_entr_real", ""};
	params[0].num_scales = 3;
	params[0].num_orientations = 5;
	params[0].remove_dc = true;
	params[0].high_frequency_bound = 0.25;
	params[0].low_frequency_bound = 0.06;
	params[0].metrics |= GaborWaveletParams::ENTROPY;
	params[0].response = GaborWaveletParams::MAGNITUDE;

	for ( int i = 0; i < 1; ++ i ) {
		ClassificationImageSet data_set;
		data_set.read( sdesc );
		ClassificationModel model;
		PGabor gabor_wavelet ( params[i] );
		model.set_feature_extractor( &gabor_wavelet);
		string fname = data + desc[i];
		std::cout << fname << std::endl;
		FeatureSet features = model.extractFeatures( data_set );
		features.write( fname + ".feat" );
		write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".libsvm" );
		features.normalize();
		features.write( fname + ".nfeat");
		write_feature_vectors_in_libsvm_format  ( data_set, features, fname + ".nlibsvm" );
		// cambio
		features.setClassNames(data_set.getClassNames());
		features.setDataClassIds(data_set.getDataClassIds());
		WekaFeatureSetDao dao;
		dao.setFeatureSet(features);
		dao.create(fname + ".arff");
	}
}
void pgabor() {
	GaborWaveletParams params;
	params.num_scales = 3;
	params.num_orientations = 5;
	params.remove_dc = true;
	params.high_frequency_bound = 0.25;
	params.low_frequency_bound = 0.06;
	params.metrics |= GaborWaveletParams::ENTROPY;
	params.response = GaborWaveletParams::MAGNITUDE;
	PGabor pgabor ( params );
	FeatureExtractorInterface *extractor;
	extractor = &pgabor;
	cv::Mat test = cv::imread("/home/jmendoza/projects/image_databases/caltech256/001.ak47/001_0001.jpg");
	extractor->extractFeatures(test);
}


int main(int argc, char* argv[])
{
	std::string opt = "none";
	if(argc > 2) {
		opt = argv[1];
	}
	//test_small_data();
	//test_classification_model();
	//test_feature_set();
	//test_feature_selection();
	//test_feature_selection_2();
	//test_feature_selection_3();
	//test_brocados_2();
	//test_knn();
	//test_all();
	//feat_cdh();
	//test_loocv();
	//cdh();
	//test_rayner();
	//test_feature_extractor();
	test_abc();
	//test_image_transform();
	//test_weka();
	//kfold ( "/home/jmendoza/Desktop/coffee.data", "/home/jmendoza/Desktop/gw_entropy_gw_3_5_TUNEADO_entr_mag.nfeat" );
	//gabor_FTW("/home/jmendoza/Desktop/LotePruebaV1.3/BS/BrocadoS2888.jpg");
	//gabor_FTW("/home/jmendoza/face.jpg");
	//gabor_FTW("/home/jmendoza/projects/compucorp/data/sample.jpg");
	//test_gabor_visualizer("/Users/mac/research/databases/cecovasa/BrocadoL/BrocadoL1.jpg");
	//compucorp();
	//mit();
	//kfold ( "/home/jmendoza/projects/image_databases/mit.data", "/home/jmendoza/projects/image_databases/mitpgabor_mag.nfeat");
	//test_gabor_pca();
	//pgabor();
	/*
	FeatureSet fs;
	WekaFeatureSetDao weka;
	weka.read("/home/jmendoza/projects/image_databases/mitpgabor_mag.arff");
	fs = weka.getFeatureSet();
	ShogunFeatureSetDao dao;
	dao.setFeatureSet(fs);
	dao.create("/home/jmendoza/Desktop/out");
	*/
	return 0;
}
