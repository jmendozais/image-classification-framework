/*
 * abc_selector.cpp
 *
 *  Created on: Jan 19, 2013
 *      Author: jmendoza
 */

#include "abc_selector.h"
std::vector< std::vector< double > > ABCWrapperSelector::selectFeatures ( const std::vector<std::vector<double> > &data, const std::vector<int> &responses ) {
	ClassifierFitnessFunction fitness ( classifier_, evaluator_, data, responses );
	int dim = data[0].size();
	BinaryABCOptimizacion optimizer ( dim, &fitness );
	optimizer.ejecutar();
	std::cout << "END SELECTION" << std::endl;
}

