/*
 * AlshaUtils.h
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHAUTILS_H_
#define ALSHAUTILS_H_

#include "AlshaTypes.h"
class AlshaUtils
{
public:
	//memory operations
	static void* 	memAlloc(NumType size);
	static void*	memCalloc(NumType nmemb, int size);
	static void* 	memRealloc(void* addr, NumType size);
	static void 	memFree(void* addr);
	static void		memVerify(void* addr);

	//exiting programs
	static void 	exitProgram(const char* msg = "");

	//file operations
	static bool	testFile(const char* fileName);

	//program name and version
	static const char* getProgramName();
	static const char* getProgramVersion();
	static const char* getProgramAuthor();

private:
	static void	memAllocFailed(NumType size, const char* function);
};

#define ALSHA_NAME 		"PASHA"
#define ALSHA_VERSION 	"1.0.10"
#define ALSHA_AUTHOR 	"LIU YONGCHAO (liuy0039@ntu.edu.sg)"

#endif /* ALSHAUTILS_H_ */
