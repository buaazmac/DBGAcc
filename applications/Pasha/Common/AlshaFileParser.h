/*
 * AlshaFileParser.h
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHAFILEPARSER_H_
#define ALSHAFILEPARSER_H_

#include "../zlib/zlib.h"
#include "AlshaTypes.h"
#include "AlshaUtils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

class AlshaFileParser
{
public:
	static AlshaFileParser* getParser(const char* fileName, int type);

	virtual ~AlshaFileParser();

	virtual int getNextSeq(char* sequence);
	virtual char* getNextSeq(int* seqLength);
	NumType getNumSeqs(){return numSeqs;}
	NumType getRealNumSeqs(){return realNumSeqs;}
	void  	reset();
protected:
	//the output file is opened in append mode
	AlshaFileParser();
	AlshaFileParser(const char* fileName);

	//protected member functions
	bool encoding(char* str);

	NumType numSeqs;
	NumType realNumSeqs;
	
	gzFile file;
	bool first;
};
#endif /* ALSHAFILEPARSER_H_ */
