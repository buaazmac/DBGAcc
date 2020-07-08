/*
 * AlshaFastaFileParser.h
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHAFASTAFILEPARSER_H_
#define ALSHAFASTAFILEPARSER_H_

#include "AlshaFileParser.h"

class AlshaFastaFileParser : public AlshaFileParser
{
public:
	AlshaFastaFileParser(const char* fileName);
	~AlshaFastaFileParser();

	int getNextSeq(char* sequence);
	char* getNextSeq(int* length);

private:
	char buffer[10000];
};
#endif /* ALSHAFASTAFILEPARSER_H_ */
