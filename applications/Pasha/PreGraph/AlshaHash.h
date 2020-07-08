#ifndef _ALSHA_HASH_H
#define _ALSHA_HASH_H

#include "AlshaUtils.h"

class AlshaHash
{
public:
	static unsigned int hashCRC(const char* buf, int len);
private:
	static unsigned int crc_table[256];
	static unsigned int crc32_v(const char *buf, int len);
};
#endif

