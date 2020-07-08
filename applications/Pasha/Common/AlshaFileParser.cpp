/*
 * AlshaFileParser.cpp
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#include "AlshaFileParser.h"
#include "AlshaUtils.h"
#include <iostream>
using namespace std;
#include "AlshaFastaFileParser.h"
#include "AlshaFastqFileParser.h"

AlshaFileParser* AlshaFileParser::getParser(const char* fileName, int fileType)
{
	AlshaFileParser* parser;
	switch(fileType){
	case ALSHA_FILE_FORMAT_FASTA:
 		parser = new AlshaFastaFileParser(fileName);
		break;
	case ALSHA_FILE_FORMAT_FASTQ:
		parser = new AlshaFastqFileParser(fileName);
		break;
	default:
		AlshaUtils::exitProgram("Unknown file format");
   	}
	return parser;
}
AlshaFileParser::AlshaFileParser()
{
	//initilaize varaibles
	this->numSeqs = 0;
	this->realNumSeqs = 0;

	this->file = NULL;
	this->first = true;
}
AlshaFileParser::AlshaFileParser(const char* fileName)
{
	//initilaize varaibles
	this->numSeqs = 0;
	this->realNumSeqs = 0;
	this->file = NULL;
	this->first = true;
	//
	if(fileName == NULL){
		cout << "Input file is not specified" << endl;
		AlshaUtils::exitProgram();
	}
	//open the input file wiht reading mode
	this->file = gzopen(fileName, "rb");
	if(this->file == NULL){
		cout << "Failed to open input file: " << fileName << endl;
		AlshaUtils::exitProgram();
	}else{
		gzseek(this->file, 0, SEEK_SET);
		cout <<"Reading file " << fileName << endl;
	}
}
AlshaFileParser::~AlshaFileParser()
{
	if(this->file != NULL){
		gzclose(this->file);
		this->file = NULL;
	}
}
int AlshaFileParser::getNextSeq(char* sequence)
{
	AlshaUtils::exitProgram("getNextSeq is not supported for this file format");
	return 0;
}
char* AlshaFileParser::getNextSeq(int* seqLength)
{
	return 0;
}
void AlshaFileParser::reset()
{
	if(this->file){
		gzseek(this->file, 0, SEEK_SET);
	}
	this->numSeqs = 0;
	this->realNumSeqs = 0;
	this->first = true;
}
bool AlshaFileParser::encoding(char* str)
{
	for( ; *str != '\0'; str++){
		switch(*str){
		case 'A':
		case 'a':
			*str = ADENINE;
			break;
		case 'C':
		case 'c':
			*str = CYTOSINE;
			break;
		case 'G':
		case 'g':
			*str = GUANINE;
			break;
		case 'T':
		case 't':
			*str = THYMINE;
			break;
		default:
			*str = ADENINE;
		}
	}
	return true;
}
