/*
 * AlshaFastqFileParser.h
 *
 *  Created on: 30-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHAFASTQFILEPARSER_H_
#define ALSHAFASTQFILEPARSER_H_

#include "AlshaFileParser.h"

class AlshaFastqFileParser : public AlshaFileParser
{
public:
	AlshaFastqFileParser(const char* fileName);
	~AlshaFastqFileParser();

	char* getNextSeq(int* seqLength);
	int getNextSeq(char* sequence);
private:
	char buffer[10000];
};

#endif /* ALSHAFASTQFILEPARSER_H_ */
