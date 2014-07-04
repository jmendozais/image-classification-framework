/*
 * mrMr.cpp
 *
 *  Created on: Jan 16, 2013
 *      Author: jmendoza
 */

#include "mrMr.h"

#include <iostream>
#include <vector>
#include "../ipcore/util.h"

GreedyMinRedundanceMaxRelevanceSelector::GreedyMinRedundanceMaxRelevanceSelector ( MRMR_params params ) {
		params_ = params;
	}
FeatureSet GreedyMinRedundanceMaxRelevanceSelector::selectFeatures ( const FeatureSet& features ) {
		cv::Mat feature_set = features.cv_mat();
		cv::Mat classes ( features.size(), 1, CV_64FC1 );
		for ( int i = 0; i < features.size(); ++ i  )
			classes.ptr<double>(i)[0] = (double) params_.data_set.get_class_id(i);

		// incremental search algorithm
		int N = feature_set.rows;
		int M = feature_set.cols;
		assert ( params_.target_dim <= M );
		std::vector < bool > S ( M );
		cv::Mat corr = ipcore::correlation ( feature_set );
		for ( int m = 0; m < params_.target_dim; ++ m ) {
			int imax = 0;
			double vmax = -1e100;
			for ( int i = 0; i < M; ++ i ) {
				// belong
				if ( S[i] ) continue;
				// Calculating mrMr
				double mrMr = 0;
				cv::Mat x_i = feature_set ( cv::Range::all(), cv::Range( i, i + 1 ) );
				// s_m-1
				if ( m ) {
					for ( int j = 0; j < M; ++ j )
						// belong ?
						if ( S[j] ) {
							cv::Mat x_j = feature_set ( cv::Range::all(), cv::Range( j, j + 1 ) );
							mrMr += abs ( corr.ptr<double>(i)[j] );
						}
					mrMr /= ( m*m );
				}
				// difference or quotient ? quotient :D
				mrMr = ipcore::f_test ( x_i , classes ) / mrMr;
				// update max D
				if ( vmax < mrMr ) {
					vmax = mrMr;
					imax = i;
				}
			}
			S[imax] = true;
		}
		cv::Mat resp =  cv::Mat( N, params_.target_dim, CV_64F ).clone();
		for ( int i = 0; i < N; ++ i ) {
			int idx = 0;
			for ( int j = 0; j < M; ++ j )
				if ( S[j] )
					resp.ptr<double>(i)[idx++] = feature_set.ptr<double>(i)[j];
		}
		FeatureSet ans = features;
		ans.change_feature_vectors(resp);
		return ans;
}
