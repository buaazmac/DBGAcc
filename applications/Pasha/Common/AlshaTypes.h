/*
 * AlahsTypes.h
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALAHSTYPES_H_
#define ALAHSTYPES_H_

#include <algorithm>
#include <map>
#include <ext/hash_set>
#include <set>
#include <list>
#include <deque>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#ifdef ALSHA_MPI
#include <mpi.h>
#endif

using namespace std;
using namespace __gnu_cxx;
#ifdef _WIN32
typedef unsigned __int64		uint64_t;
typedef __int64					int64_t;
typedef unsigned int			uint32_t;
typedef signed int				int32_t;
typedef unsigned short			uint16_t;
typedef signed short			int16_t;
typedef unsigned char			uint8_t;
typedef signed char				int8_t;
#else
#include <stdint.h>
#endif

typedef int64_t 	IDtype;
typedef int64_t 	NumType;

typedef uint64_t		AlshaKmerKey;
typedef uint64_t		AlshaKmerData;
#define AlshaMap		map
#define AlshaSet		hash_set
#define AlshaQueue		list


#include "tbb/tbb_allocator.h"
#include "tbb/scalable_allocator.h"
#include "tbb/cache_aligned_allocator.h"
#define AlshaAllocator tbb::tbb_allocator

#include "AlshaDefs.h"
#endif /* ALAHSTYPES_H_ */
