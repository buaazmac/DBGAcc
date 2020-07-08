/*
 * AlshaFastaFileParser.cpp
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#include "AlshaFastaFileParser.h"

AlshaFastaFileParser::AlshaFastaFileParser(const char* fileName) : AlshaFileParser(fileName)
{
}
int AlshaFastaFileParser::getNextSeq(char* buffer)
{
	int seqLength = 0;
	const int maxline = 10000;
	if(first){
		//check the file format
		int c = gzgetc(this->file);
		if(c != '>'){
			AlshaUtils::exitProgram("NOT in FASTA format");
		}
		gzungetc(c, this->file);
		first = false;
	}
	//get the sequence name
	if(gzgets(this->file, buffer, maxline) && buffer[0] != '\n' && buffer[0] != '\r'){
		//increase the number of sequences
		this->numSeqs++;
	}else{
		//reaching the end of the file
		cout << "Found " << this->numSeqs << " sequences" << endl;
		return 0;
	}
	//get the sequence body
	if(gzgets(this->file, buffer, maxline) == NULL){
		AlshaUtils::exitProgram("incomplete file\n");
	}
	//converting the strings
	int i;
	for(i = strlen(buffer) - 1; i >= 0 && (buffer[i] == '\n' || buffer[i] == '\r'); --i){
		buffer[i] = '\0';
	}
	seqLength = i + 1;
	this->encoding(buffer);
	
	return seqLength;
}

char* AlshaFastaFileParser::getNextSeq(int* seqLength)
{
	const int maxline = 10000;
	if(first){
		//check the file format
		int c = gzgetc(this->file);
		if(c != '>'){
			AlshaUtils::exitProgram("NOT in FASTA format");
		}
		gzungetc(c, this->file);
		first = false;
	}
	//get the sequence name
	if(gzgets(this->file, buffer, maxline) && buffer[0] != '\n' && buffer[0] != '\r'){
		//increase the number of sequences
		this->numSeqs++;
	}else{
		//reaching the end of the file
		cout << "Found " << this->numSeqs << " sequences" << endl;
		*seqLength = 0;
		return 0;
	}
	//get the sequence body
	if(gzgets(this->file, buffer, maxline) == NULL){
		AlshaUtils::exitProgram("incomplete file\n");
	}
	//converting the strings
	int i;
	for(i = strlen(buffer) - 1; i >= 0 && (buffer[i] == '\n' || buffer[i] == '\r'); --i){
		buffer[i] = '\0';
	}
	*seqLength = i + 1;

	this->encoding(buffer);

	return buffer;
}
AlshaFastaFileParser::~AlshaFastaFileParser()
{
}
