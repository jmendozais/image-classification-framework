/*
 * base.h
 *
 *  Created on: Jan 16, 2013
 *      Author: jmendoza
 */

#ifndef BASE_H_
#define BASE_H_
class Task {
	void perform();
};
template<typename I, typename O>
class Batch : Task {
	Batch ( const I& input, O& output );
};


#endif /* BASE_H_ */
