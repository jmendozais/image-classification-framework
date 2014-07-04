/*
 * gabor_wavelet.cpp
 *
 *  Created on: Jan 15, 2013
 *      Author: jmendoza
 */

#include "gabor_wavelet.h"

void grid ( std::vector<cv::Mat> images, int rows, int cols ) {
	for ( int i = 0; i < images.size(); ++ i ) {
		cv::resize(images[i], images[i], cv::Size(200,200));
		images[i] = ipcore::normalize(images[i]);
	}
	int max_rows, max_cols;
	max_rows = max_cols = -1;
	for ( int i = 0; i < images.size(); ++ i ) {
		max_rows = max ( max_rows, images[i].rows );
		max_cols = max ( max_cols, images[i].cols );
	}
	std::cout << max_rows << " " << max_cols << std::endl;
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
	cv::waitKey();
}

ManjunathFilterBank::ManjunathFilterBank ( FilterBankParams params ) {
	params_ = params;
}
void ManjunathFilterBank::gabor_filter ( cv::Mat &gabor_real, cv::Mat &gabor_im, int i_scale, int i_orient  ) {
		double a, u0, z, Uvar, Vvar, Xvar, Yvar, X, Y, G, t1, t2, m;
		int x, y;
		double u_h = params_.high_frequency_bound;
		double u_l = params_.low_frequency_bound;
		int num_scales = params_.num_scales;
		int num_orient = params_.num_orientations;
		int side = params_.kernel_side;

		a = pow( u_h / u_l, 1.0/(double)(num_scales-1) );
		u0 = u_h/pow(a, (double) num_scales-i_scale-1);
		Uvar = (a-1.0)*u0/((a+1.0)*sqrt(2.0*log(2.0)));
		z = -2.0*log(2.0)*(Uvar*Uvar)/u0;
		Vvar = tan(PI/(2*num_orient))*(u0+z)/sqrt(2.0*log(2.0)-z*z/(Uvar*Uvar));
		Xvar = 1.0/(2.0*PI*Uvar);
		Yvar = 1.0/(2.0*PI*Vvar);

		t1 = cos(PI/num_orient*i_orient);
		t2 = sin(PI/num_orient*i_orient);

		//fprintf ( stderr, "GF %d %d: Ul %lf, Uh %lf, u0 %lf, Uvar %lf, Vvar %lf, Xvar %lf, Yvar %lf\n", i_scale, i_orient, u_l, u_h, u0, Uvar, Vvar, Xvar, Yvar );
		for (x=0;x<2*side+1;x++) {
			for (y=0;y<2*side+1;y++) {
				X = (double) (x-side)*t1+ (double) (y-side)*t2;
				Y = (double) -(x-side)*t2+ (double) (y-side)*t1;
				G = 1.0/(2.0*PI*Xvar*Yvar)*pow(a, (double) num_scales-i_scale-1)*exp(-0.5*((X*X)/(Xvar*Xvar)+(Y*Y)/(Yvar*Yvar)));
				gabor_real.ptr<double>(x)[y] = G*cos(2.0*PI*u0*X);
				gabor_im.ptr<double>(x)[y] = G*sin(2.0*PI*u0*X);
			}
		}

		if ( params_.remove_dc ) {
			m = 0;
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					m += gabor_real.ptr<double>(x)[y];
			m /= pow((double) 2.0*side+1, 2.0);
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					gabor_real.ptr<double>(x)[y] -= m;

			m = 0;
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					m += gabor_im.ptr<double>(x)[y];
			m /= pow((double) 2.0*side+1, 2.0);
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					gabor_im.ptr<double>(x)[y] -= m;
		}
	}
std::vector<cv::Mat> ManjunathFilterBank::getFilterBank ( FilterBankParams params ) {
	params_ = params;
#define p_ params_
	std::vector<cv::Mat> filter_bank ( p_.num_scales * p_.num_orientations * 2 );
	for ( int i = 0; i < p_.num_scales * p_.num_orientations; ++ i ) {
		filter_bank[2*i] = cv::Mat ( 2*p_.kernel_side+1, 2*p_.kernel_side+1, CV_64FC1 );
		filter_bank[2*i+1] = cv::Mat ( 2*p_.kernel_side+1, 2*p_.kernel_side+1, CV_64FC1 );
		gabor_filter ( filter_bank[2*i], filter_bank[2*i+1], i/p_.num_orientations, i%p_.num_orientations );
	}
	return filter_bank;
}
	// Generate a filter bank with the current parameters
	void GaborWavelet::update() {
		filter_bank = filter_bank_builder_ptr_->getFilterBank(params_);
		/*
		filter_bank  = std::vector<cv::Mat> ( params_.num_scales * params_.num_orientations * 2 );
		for ( int i = 0; i < params_.num_scales * params_.num_orientations; ++ i ) {
			filter_bank[2*i] = cv::Mat ( 2*params_.kernel_side+1, 2*params_.kernel_side+1, CV_64FC1 );
			filter_bank[2*i+1] = cv::Mat ( 2*params_.kernel_side+1, 2*params_.kernel_side+1, CV_64FC1 );
			gabor_filter ( filter_bank[2*i], filter_bank[2*i+1], i/params_.num_orientations, i%params_.num_orientations );
		}
		*/
	}
	// Generate a gabor filter
	/*
	void GaborWavelet::gabor_filter ( cv::Mat &gabor_real, cv::Mat &gabor_im, int i_scale, int i_orient  ) {
		double a, u0, z, Uvar, Vvar, Xvar, Yvar, X, Y, G, t1, t2, m;
		int x, y;


		double u_h = params_.high_frequency_bound;
		double u_l = params_.low_frequency_bound;
		int num_scales = params_.num_scales;
		int num_orient = params_.num_orientations;
		int side = params_.kernel_side;

		a = pow( u_h / u_l, 1.0/(double)(num_scales-1) );
		u0 = u_h/pow(a, (double) num_scales-i_scale-1);
		Uvar = (a-1.0)*u0/((a+1.0)*sqrt(2.0*log(2.0)));
		z = -2.0*log(2.0)*(Uvar*Uvar)/u0;
		Vvar = tan(PI/(2*num_orient))*(u0+z)/sqrt(2.0*log(2.0)-z*z/(Uvar*Uvar));
		Xvar = 1.0/(2.0*PI*Uvar);
		Yvar = 1.0/(2.0*PI*Vvar);

		t1 = cos(PI/num_orient*i_orient);
		t2 = sin(PI/num_orient*i_orient);

		//fprintf ( stderr, "GF %d %d: Ul %lf, Uh %lf, u0 %lf, Uvar %lf, Vvar %lf, Xvar %lf, Yvar %lf\n", i_scale, i_orient, u_l, u_h, u0, Uvar, Vvar, Xvar, Yvar );
		for (x=0;x<2*side+1;x++) {
			for (y=0;y<2*side+1;y++) {
				X = (double) (x-side)*t1+ (double) (y-side)*t2;
				Y = (double) -(x-side)*t2+ (double) (y-side)*t1;
				G = 1.0/(2.0*PI*Xvar*Yvar)*pow(a, (double) num_scales-i_scale-1)*exp(-0.5*((X*X)/(Xvar*Xvar)+(Y*Y)/(Yvar*Yvar)));
				gabor_real.ptr<double>(x)[y] = G*cos(2.0*PI*u0*X);
				gabor_im.ptr<double>(x)[y] = G*sin(2.0*PI*u0*X);
			}
		}

		if ( params_.remove_dc ) {

			m = 0;
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					m += gabor_real.ptr<double>(x)[y];
			m /= pow((double) 2.0*side+1, 2.0);
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					gabor_real.ptr<double>(x)[y] -= m;

			m = 0;
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					m += gabor_im.ptr<double>(x)[y];
			m /= pow((double) 2.0*side+1, 2.0);
			for (x=0;x<2*side+1;x++)
				for (y=0;y<2*side+1;y++)
					gabor_im.ptr<double>(x)[y] -= m;
		}
	}
	*/
	GaborWavelet::GaborWavelet ( int sc, int ors, double u_l, double u_h, int side, bool remove_dc = false ) {
		params_.num_scales = sc;
		params_.num_orientations = ors;
		params_.low_frequency_bound = u_l;
		params_.high_frequency_bound = u_h;
		params_.kernel_side = side;
		params_.remove_dc = false;
		filter_bank_builder_ptr_ = FilterBankBuilderPtr(new ManjunathFilterBank(params_));
		update();
	}
	GaborWavelet::GaborWavelet ( GaborWaveletParams params ) {
			mask_type = params.mask_type;
			params_ = params;
			filter_bank_builder_ptr_ = FilterBankBuilderPtr(new ManjunathFilterBank(params_));
			update();
	}
	GaborWavelet::GaborWavelet () {
		GaborWaveletParams params;
		params_ = params;
		mask_type = params.mask_type;
		filter_bank_builder_ptr_ = FilterBankBuilderPtr(new ManjunathFilterBank(params_));
		update();
	}
	bool GaborWavelet::persist ( std::string filename ) {

	}
	bool GaborWavelet::load ( std::string filename ) {

	}
	std::vector<cv::Mat> GaborWavelet::convolveOneChannel ( const cv::Mat &image ) const {
		std::vector<cv::Mat> resp;
		Mat feature_vector ( 2 * params_.num_orientations * params_.num_scales, 1, CV_64F );
		int n = getOptimalDFTSize ( image.rows + 2*params_.kernel_side );
		int m = getOptimalDFTSize ( image.cols + 2*params_.kernel_side );
		Mat image_padded, I;
		// Padding
		//copyMakeBorder ( (Mat)(Mat_<double>(image)), image_padded, side, n - image.rows - side, side, m - image.cols - side, BORDER_REFLECT );
		copyMakeBorder ( (Mat)(Mat_<double>(image)), image_padded, params_.kernel_side, n - image.rows - params_.kernel_side, params_.kernel_side, m - image.cols - params_.kernel_side, BORDER_REFLECT );
		// DFT
		dft ( image_padded, I, DFT_COMPLEX_OUTPUT );
		// Por ahora solo trabajamos con la parte real del filtro de gabor :S
		for ( int i = 0; i < filter_bank.size(); i += 2 ) {
			Mat g_padded, result, R, G;

			Mat g_complex;
			Mat gp_arr[2];

			// Padding
			//copyMakeBorder ( filter_bank[i], g_padded, 0, n - filter_bank[i].rows, 0, m - filter_bank[i].cols, BORDER_CONSTANT );

			//copyMakeBorder ( filter_bank[i], gp_arr[0], side, n - filter_bank[i].rows - side, side, m - filter_bank[i].cols - side, BORDER_CONSTANT );
			//copyMakeBorder ( filter_bank[i+1], gp_arr[1], side, n - filter_bank[i+1].rows - side, side, m - filter_bank[i+1].cols - side, BORDER_CONSTANT );
			copyMakeBorder ( filter_bank[i], gp_arr[0], 0, n - filter_bank[i].rows, 0, m - filter_bank[i].cols, BORDER_CONSTANT );
			copyMakeBorder ( filter_bank[i+1], gp_arr[1], 0, n - filter_bank[i+1].rows, 0, m - filter_bank[i+1].cols, BORDER_CONSTANT );


			merge ( gp_arr, 2, g_complex );
			dft ( g_complex, G, DFT_COMPLEX_OUTPUT );

			// use linear filter for kernels width side <= 11, and DFT for kernels with side > 11
			mulSpectrums ( I, G, R, 1);
			dft ( R, result, DFT_INVERSE + DFT_SCALE + DFT_COMPLEX_OUTPUT );

			//filter2D ( image_padded , result, -1, filter_bank[i] );
			Mat rarr[2];
			split ( result, rarr );
			Mat cresult1 ( rarr[0], Rect ( 2*params_.kernel_side, 2*params_.kernel_side, image.cols, image.rows ) );
			Mat cresult2 ( rarr[1], Rect ( 2*params_.kernel_side, 2*params_.kernel_side, image.cols, image.rows ) );
			/*
			printf ( "TO SPLIT\n" );
			cv::namedWindow ( "filter" );
			cv::imshow ( "filter", ipcore::normalize(cresult));
			cv::namedWindow ( "padded" );
			cv::imshow ( "padded", ipcore::normalize(gp_arr[0]) );
			cv::waitKey();
			*/
			//Mat cropped ( result, Rect ( 2*side, 2*side, image.cols, image.rows ) );
			resp.push_back ( cresult1 );
			resp.push_back ( cresult2 );
		}
		return resp;
	}
	std::vector<cv::Mat> GaborWavelet::convolve( const cv::Mat &image ) const {
		std::vector<cv::Mat> v_image;
		cv::split ( image, v_image );
		std::vector<cv::Mat> ans;
		for ( int i = 0; i < v_image.size(); ++ i ) {
			std::vector<cv::Mat> curconv = convolveOneChannel(v_image[i]);
			ans.insert( ans.begin(), curconv.begin(), curconv.end() );
		}
		return ans;
	}
	std::vector<cv::Mat> GaborWavelet::getMagnitudes ( const std::vector<cv::Mat> &real_imag ) const {
		std::vector<cv::Mat> ans ( real_imag.size()/2 );
		for ( int i = 0; i < real_imag.size(); i += 2 ) {
			cv::Mat mag ( real_imag[i].rows, real_imag[i].cols, CV_64FC1 );
			for ( int x = 0; x < real_imag[i].rows; ++ x )
				for ( int y = 0; y < real_imag[i].cols; ++ y ) {
					double a = real_imag[i].ptr<double>(x)[y];
					double b = real_imag[i+1].ptr<double>(x)[y];
					mag.ptr<double>(x)[y] = sqrt( a*a + b*b );
				}
			ans[i/2] = mag;
		}
		return ans;
	}
	std::vector<cv::Mat> GaborWavelet::getMagnitudes ( const cv::Mat &image ) const {
		std::vector<cv::Mat> convs = convolve ( image );
		return GaborWavelet::getMagnitudes(convs);
	}
	const int GaborWavelet::getNumMetrics() const {
		int num_metrics = 0;
		for ( int i = 0; i < 31; ++ i )
			if ( params_.metrics & (1<<i) )
				++ num_metrics;
		return num_metrics;
	}
	cv::Mat GaborWavelet::measure( std::vector<cv::Mat> inputs ) {
		int num_metrics = GaborWavelet::getNumMetrics();
		cv::Mat resp ( num_metrics * inputs.size(), 1, CV_64FC1 );
		int idx = 0;
		for ( int i = 0; i < inputs.size(); ++ i ) {
			double imean, sd, entropy;
			if ( mask_type != GaborWaveletParams::NO_MASK ) {
				imean = ipcore::mean ( inputs[i], mask_ );
				sd = ipcore::standard_deviation ( inputs[i], mask_, imean );
				entropy = ipcore::entropy ( inputs[i], mask_ );
			} else {
				imean = ipcore::mean ( inputs[i] );
				sd = ipcore::standard_deviation ( inputs[i], imean );
				entropy = ipcore::entropy ( inputs[i] );
			}
			if ( params_.metrics & GaborWaveletParams::MEAN )
				resp.ptr<double>( idx++ )[0] = imean;
			if ( params_.metrics & GaborWaveletParams::STD )
				resp.ptr<double>( idx++ )[0] = sd;
			if ( params_.metrics & GaborWaveletParams::ENTROPY )
				resp.ptr<double>( idx++ )[0] = entropy;
		}
		return resp;
	}
	/*
	cv::Mat GaborWavelet::extract_features_one_channel ( const cv::Mat & image ) {
		std::vector<cv::Mat> temp, resp;
    temp = GaborWavelet::getMagnitudes(image);

    if ( params_.response & GaborWaveletParams::MAGNITUDE )
            resp = temp;
    temp = GaborWavelet::convolveOneChannel(image);
    if ( params_.response & GaborWaveletParams::REAL )
            for ( int i = 0; i < temp.size(); i += 2 )
                    resp.push_back ( temp[i] );
    if ( params_.response & GaborWaveletParams::IMAGINARY )
            for ( int i = 0; i+1 < temp.size(); i += 2 )
                    resp.push_back ( temp[i+1] );
    //grid(resp,3,5);
		int num_metrics = 0;
		for ( int i = 0; i < 30; ++ i )
			if ( (1<<i) & params_.metrics )
				++ num_metrics;
		Mat feature_vector ( num_metrics * resp.size(), 1, CV_64F );

		int idx = 0;
		for ( int i = 0; i < resp.size(); ++i ) {
			double imean, sd, entropy;
			if ( mask_type != GaborWaveletParams::NO_MASK ) {
				imean = ipcore::mean ( resp[i], mask_ );
				sd = ipcore::standard_deviation ( resp[i], mask_, imean );
				entropy = ipcore::entropy ( resp[i], mask_ );
			} else {
				imean = ipcore::mean ( resp[i] );
				sd = ipcore::standard_deviation ( resp[i], imean );
				entropy = ipcore::entropy ( resp[i] );
			}
			if ( params_.metrics & GaborWaveletParams::MEAN )
				feature_vector.ptr<double>( idx++ )[0] = imean;
			if ( params_.metrics & GaborWaveletParams::STD )
				feature_vector.ptr<double>( idx++ )[0] = sd;
			if ( params_.metrics & GaborWaveletParams::ENTROPY )
				feature_vector.ptr<double>( idx++ )[0] = entropy;
		}
		return feature_vector;
	}
	*/
	void GaborWavelet::setMaskType( int type ) {
		mask_type = type;
		params_.mask_type = type;
	}
 std::vector<cv::Mat> GaborWavelet::getResponses ( const cv::Mat &image ) const {
	 std::vector<cv::Mat> responses;
	 std::vector<cv::Mat> real_imag = convolve ( image );
	 if ( params_.response & GaborWaveletParams::REAL )
		 for ( int i = 0; i < real_imag.size(); i += 2 )
			 responses.push_back ( real_imag[i] );
	 if ( params_.response & GaborWaveletParams::IMAGINARY )
		 for ( int i = 0; i+1 < real_imag.size(); i += 2 )
			 responses.push_back ( real_imag[i+1] );
	 if ( params_.response & GaborWaveletParams::MAGNITUDE ) {
		 std::vector<cv::Mat> mags = GaborWavelet::getMagnitudes(real_imag);
	   for ( int i = 0; i < mags.size(); ++ i )
	  	 responses.push_back(mags[i]);
	 }
	 return responses;
 }
 cv::Mat GaborWavelet::extractFeatures ( cv::Mat &image ) {
	 	if ( mask_type == GaborWaveletParams::CUSTOM_MASK );
	 	else if ( mask_type != GaborWaveletParams::NO_MASK )
	 		GaborWavelet::generate_mask ( image, mask_, mask_type );
	 	std::vector<cv::Mat> responses = GaborWavelet::getResponses(image);
	 	return GaborWavelet::measure(responses);
	}


	void GaborWavelet::generate_mask ( const cv::Mat& input, cv::Mat& mask, int mask_type ) {
		if ( mask_type == GaborWaveletParams::BACKGROUND_REMOVAL_MASK );

		// otsu thresholding
		Mat gray_input, u_input, temp;
		cvtColor( input, gray_input, CV_BGR2GRAY );
		gray_input.convertTo( u_input, CV_8UC1 );
		cv::threshold(u_input, temp, 100, 255, THRESH_OTSU );

		// invert scalars
		for ( int i = 0; i < temp.rows; ++ i )
			for ( int j = 0; j < temp.cols; ++ j )
				if ( temp.ptr<uchar>(i)[j] == 0 )
					temp.ptr<uchar>(i)[j] = 255;
				else temp.ptr<uchar>(i)[j] = 0;

		// closing holes ( max rad of holes 5px )
		int side = 5;
		cv::copyMakeBorder(temp, temp, side + 1, side + 1, side + 1, side + 1, BORDER_CONSTANT );
		Mat element = getStructuringElement( 2, Size( 2*side + 1, 2*side+1 ), Point( side, side ) );
		cv::morphologyEx ( temp, temp, MORPH_CLOSE, element );
		Mat cropped ( temp, Rect ( side + 1, side + 1, input.cols, input.rows ) );
		mask = cropped;
}
const cv::Mat& GaborWavelet::getFrequencyProfile() const {
	int rows = filter_bank[0].rows;
	int cols = filter_bank[0].cols;
	cv::Mat ans ( rows, cols, CV_64FC1 );
	for ( int i = 0; i < ans.rows; ++ i )
		for ( int j = 0; j < ans.cols; ++ j )
			ans.ptr<double>(i)[j] = 0;
	for ( int i = 0; i < filter_bank.size(); i += 2 ) {

		std::vector<cv::Mat> g_complex (2);
		g_complex[0] = filter_bank[i];
		g_complex[1] = filter_bank[i+1];

		cv::Mat g, G;
		cv::merge(g_complex, g);

		cv::dft( g, G );
		cv::Mat G_arr[2], mag;
		cv::split( G, G_arr );
		cv::magnitude ( G_arr[0], G_arr[1], mag );

		double k = -1e100, K;
		for ( int x = 0; x < mag.rows; ++ x )
			for ( int y = 0; y < mag.cols; ++ y )
				k = std::max( mag.ptr<double>(x)[y], k );

		for ( int x = 0; x < mag.rows; ++ x )
			for ( int y = 0; y < mag.cols; ++ y )
				if ( abs(mag.ptr<double>(x)[y] - k/2) < k/32 )
					mag.ptr<double>(x)[y] = 1;

		G = mag + cv::Scalar::all(1);
		cv::log(G, G);
		G = ipcore::normalize(G);
		ans = ans + G;
	}
	int midx = ans.rows/2;
	int midy = ans.cols/2;
	cv::Mat tl ( ans, cv::Rect(0, 0, midy, midx) );
	cv::Mat tr ( ans, cv::Rect(midy, 0, midy, midx) );
	cv::Mat bl ( ans, cv::Rect(0, midx, midy, midx) );
	cv::Mat br ( ans, cv::Rect(midy, midx, midy, midx) );

	cv::Mat temp;
	tl.copyTo(temp);
	br.copyTo(tl);
	temp.copyTo(br);
	tr.copyTo(temp);
	bl.copyTo(tr);
	temp.copyTo(bl);

	cv::namedWindow("Frequencies");
	cv::imshow("Frequencies", ipcore::normalize(ans));
	cv::waitKey();
	return ans;
}

const int GaborWavelet::getNumScales() const {
	return params_.num_scales;
}
const int GaborWavelet::getNumOrientations() const {
	return params_.num_orientations;
}
const double GaborWavelet::getHighFrequencyBound() const {
	return params_.high_frequency_bound;
}
const double GaborWavelet::getLowFrequencyBound() const {
	return params_.low_frequency_bound;
}
double GaborWavelet::setHighFrequencyBound(double value) {
	params_.high_frequency_bound = value;
	update();
}
double GaborWavelet::setLowFrequencyBound(double value) {
	params_.low_frequency_bound = value;
	update();
}
void GaborWavelet::setFilterBank(std::vector<cv::Mat> filters) {
	filter_bank = filters;
}
std::vector<cv::Mat> GaborWavelet::getFilterBank() {
	return filter_bank;
}
void GaborWavelet::setMask( const cv::Mat &mask ) {
	mask_ = mask;
}
void GaborWavelet::setFilterBankBuilderPtr( FilterBankBuilderPtr ptr) {
	filter_bank_builder_ptr_ = ptr;
}
FilterBankBuilderPtr GaborWavelet::getFilterBankBuilderPtr() {
	return filter_bank_builder_ptr_;
}

PGabor::PGabor ( GaborWaveletParams params ) : GaborWavelet(params) {
	deep = 3;
}
cv::Mat PGabor::extractFeatures ( cv::Mat &image ) {
	int rows = image.rows;
	int cols = image.cols;
	std::vector<cv::Mat> responses = GaborWavelet::getResponses(image);
	int counter = 0;
	for ( int i = 0; i < deep; ++ i )
		counter += (1<<i)*(1<<i);
	int num_metrics = GaborWavelet::getNumMetrics();
	cv::Mat fv ( counter * num_metrics * responses.size(), 1, CV_64FC1 );

	GaborWavelet::setMaskType(GaborWaveletParams::CUSTOM_MASK);
	int idx = 0;
	for ( int i = 0; i < deep; ++ i ) {
		int srows = rows/(1<<i);
		int scols = cols/(1<<i);
		for ( int j = 0; j < (1<<i); ++ j )
			for ( int k = 0; k < (1<<i); ++ k ) {
				std::vector<double> temp;
				cv::Mat mask ( rows, cols, CV_8UC1 );
				mask = cv::Scalar(0);
				int currows = min ( (j+1)*srows, rows);
				int curcols = min ( (k+1)*scols, cols);
				cv::Mat roi = mask( cv::Rect( k*scols, j*srows, curcols - k*scols, currows - j*srows ) );
				roi = cv::Scalar(255);
				GaborWavelet::setMask(mask);
				cv::Mat fv_temp = GaborWavelet::measure(responses);
				for ( int i = 0; i < fv_temp.rows; ++ i )
					fv.ptr<double>(idx++)[0] = fv_temp.ptr<double>(i)[0];
			}
	}
	return fv;
}
