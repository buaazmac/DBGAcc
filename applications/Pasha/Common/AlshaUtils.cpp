/*
 * AlshaUtils.cpp
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#include <stdarg.h>
#include <stdio.h>
#include <errno.h>

#include "AlshaTypes.h"
#include "AlshaUtils.h"

void* AlshaUtils::memAlloc(NumType size)
{
	void* addr = malloc(size);
	if(addr == NULL){
		memAllocFailed(size, "memAlloc");
	}
	return addr;
}
void* AlshaUtils::memCalloc(NumType nmemb, int size)
{
	void* addr = calloc(nmemb, size);
	if(addr == NULL){
		memAllocFailed(size, "memCalloc");
	}
	return addr;
}

void* AlshaUtils::memRealloc(void* addr, NumType size)
{
	void* newAddr = realloc(addr, size);
	if(newAddr == NULL){
		memAllocFailed(size, "memRealloc");
	}
	return newAddr;
}

void AlshaUtils::memFree(void* addr)
{
	if(addr != NULL){
		free(addr);
	}
}
void AlshaUtils::memVerify(void* addr)
{
	if(addr == NULL){
		AlshaUtils::exitProgram("Invalid object");
	}
}
void AlshaUtils::memAllocFailed(NumType size, const char* function)
{
	cout<< function << " failed to allocate memory of size " << size << endl;

	//quit the program
	exitProgram("Program exits due to memory allocation failure");
}

void AlshaUtils::exitProgram(const char* msg)
{
	cout<< "*****************************************" << endl;
	if(msg != NULL && strlen(msg) > 0){
		cout<< "OOPS: " << msg << endl;
	}
	cout<< "*****************************************" << endl;

#ifdef ALSHA_MPI
	//exiting the program
	MPI_Finalize();
#endif

	exit(1);
}
//test the existence of the file
bool AlshaUtils::testFile(const char* fileName)
{
	if(fileName == NULL){
		return false;
	}
	FILE* fp = fopen(fileName, "rb");
	if(fp == NULL){
		return false;
	}
	fclose(fp);

	return true;
}
const char* AlshaUtils::getProgramName()
{
	return ALSHA_NAME;
}
const char* AlshaUtils::getProgramVersion()
{
	return ALSHA_VERSION;
}
const char* AlshaUtils::getProgramAuthor()
{
	return ALSHA_AUTHOR;
}

