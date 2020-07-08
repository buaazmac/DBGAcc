/*
 * AlshaKmer.h
 *
 *  Created on: 31-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHAKMER_H_
#define ALSHAKMER_H_

#include "AlshaTypes.h"
#include "AlshaUtils.h"
#include "AlshaParams.h"

typedef uint8_t 	Nucleotide;
typedef struct tagAlshaKmerProps
{
	//linkage information in two directions
	uint8_t linkages [ALSHA_KMER_DIR_NUM];
	//multiplicity of the kmer
	uint16_t multiplicity;
	//k-mer node status
	uint8_t status;
}__attribute__((packed)) AlshaKmerProps;

class AlshaKmerUtils
{
public:
	static void printKmer(AlshaKmerData kmerData)
	{
		for(int i = AlshaParams::KMERLENGTH - 1, j = 0; i >= 0; --i, ++j){
			Nucleotide nucleotide = (kmerData >> (i << 1)) & 3;
			cout << (char)decoding(nucleotide);
		}
		cout << endl;
	}
	static void printKmerProps(AlshaKmerProps& props)
	{
		cout <<"DIR: " << hex << (uint32_t)props.linkages[ALSHA_FORWARD_DIR] << " " << hex << (uint32_t)props.linkages[ALSHA_REVERSE_DIR]
			<< " Status: " << hex << (uint32_t)props.status << endl;
		cout << (int)props.multiplicity << endl;
	}
	static string getKmerString(AlshaKmerData kmerData)
	{
 		string kmer;
		kmer.reserve(AlshaParams::KMERLENGTH);
		for(int i = AlshaParams::KMERLENGTH - 1, j = 0; i >= 0; --i, ++j){
			Nucleotide nucleotide = (kmerData >> (i << 1)) & 3;
			kmer.push_back((char)decoding(nucleotide));
		}
		return kmer;
	}
	static string getKmerString(AlshaKmerData kmerData, int low, int high)
	{
 		string kmer;
		kmer.reserve(AlshaParams::KMERLENGTH);
		for(int i = high, j = 0; i >= low; --i, ++j){
			Nucleotide nucleotide = (kmerData >> (i << 1)) & 3;
			kmer.push_back((char)decoding(nucleotide));
		}
		return kmer;
	}

	static int qsortKmerCmp(const void* a, const void* b)
	{
		AlshaKmerData* p = (AlshaKmerData*)a;
		AlshaKmerData* q = (AlshaKmerData*)b;

		AlshaKmerData aData = *p;
		AlshaKmerData bData = *q;
		if(aData > bData){
			return 1;
		}else if(aData < bData){
			return -1;
		}
		return 0;
	}
	static int compareKmer(AlshaKmerData aData, AlshaKmerData bData)
	{
		if(aData > bData){
			return 1;
		}else if(aData < bData){
			return -1;
		}
		return 0;
	}
	static AlshaKmerData reverseComplementKmer(AlshaKmerData kmer)
	{
		uint8_t baseByte;
		uint8_t byteNum;
		AlshaKmerData rcKmer;

		rcKmer = 0;
		byteNum = (AlshaParams::KMERLENGTH + 3) / 4;
		//for the lowest seven bytes
		for(int i = 0; i < byteNum; ++i){
			rcKmer <<= 8;

			baseByte = byteRC[kmer & 0xFF];
			kmer >>= 8;
			
	#ifdef ALSHA_COLOR
			rcKmer += (uint64_t)baseByte;
	#else
			baseByte = ~baseByte;
			rcKmer += (uint64_t)(baseByte);
	#endif
		}
		//shift to the right
		rcKmer >>= (byteNum * 8 - AlshaParams::KMERLENGTH * 2);

		return rcKmer;
	}
	static Nucleotide shiftLeft(AlshaKmerData& data)
	{
		Nucleotide saved;

		saved = (data >> AlshaParams::KMERLENGTHSHIFT) & 3;

		data <<= 2;
		data &= AlshaParams::KMERLENGHTFILTER;

		return saved;
	}
	static Nucleotide shiftRight(AlshaKmerData& data)
	{
		Nucleotide nucleotide;

		nucleotide = (Nucleotide)(data & 3);

		data >>= 2;

		return nucleotide;
	}
	static Nucleotide popNucleotide(AlshaKmerData& data)
	{
		Nucleotide nucleotide;

		//get the popped nucleotide
		nucleotide = (Nucleotide)(data & 3);

		//shift to the right by one nucleotide
		data >>= 2;

		return nucleotide;
	}
	static Nucleotide pushNucleotideRet(AlshaKmerData& data, Nucleotide nucleotide)
	{
		Nucleotide saved;
		//save the left-most nucleotide
		saved = (data >> AlshaParams::KMERLENGTHSHIFT) & 3;
		//shift to the left by one nucleotide
		data <<= 2;

		//put the nucleotide in the 0~1 bits
		data |= (uint64_t) nucleotide;
		data &= AlshaParams::KMERLENGHTFILTER;

		return saved;
	}
	static void pushNucleotide(AlshaKmerData& data, Nucleotide nucleotide)
	{
		//shift to the left by one nucleotide
		data <<= 2;

		//put the nucleotide in the 0~1 bits
		data += (uint64_t)nucleotide;
		data &= AlshaParams::KMERLENGHTFILTER;
	}
	static Nucleotide reversePushNucleotideRet(AlshaKmerData& data, Nucleotide nucleotide)
	{
		Nucleotide saved;
		//save the right-most nucleotide
		saved = data & 3;

		//shift to the right by one nucleotide
		data >>= 2;

		data |= (((uint64_t)nucleotide) << AlshaParams::KMERLENGTHSHIFT);

		return saved;
	}
	static void reversePushNucleotide(AlshaKmerData& data, Nucleotide nucleotide)
	{
		//shift to the right by one nucleotide
		data >>= 2;

		data |= (((uint64_t) nucleotide) << AlshaParams::KMERLENGTHSHIFT);
	}
	static Nucleotide getNucleotide(AlshaKmerData& data, unsigned int pos)
	{
		return ((data >> (pos << 1)) & 3);
	}
	static void putNucleotide(AlshaKmerData& data, unsigned int pos, Nucleotide nucleotide)
	{
		data |= ((uint64_t)nucleotide) << (pos << 1);
	}
	static Nucleotide decoding(Nucleotide nucleotide)
	{
		switch(nucleotide){
		case ADENINE: 	return 'A';
		case CYTOSINE: 	return 'C';
		case GUANINE: 	return 'G';
		case THYMINE: 	return 'T';
		}
		return 'N';
	}
	static Nucleotide encoding(Nucleotide nucleotide)
	{
		switch(nucleotide){
		case 'A':
		case 'a': 	return ADENINE;
		case 'C':
		case 'c':	return CYTOSINE;
		case 'G':
		case 'g':	return GUANINE;
		case 'T':
		case 't':	return THYMINE;
		}
		return ADENINE;
	}
	static void reverseProps(AlshaKmerProps& props)
	{
		uint8_t forward = props.linkages[ALSHA_FORWARD_DIR];
		props.linkages[ALSHA_FORWARD_DIR] = complements[props.linkages[ALSHA_REVERSE_DIR]];
		props.linkages[ALSHA_REVERSE_DIR] = complements[forward];
	}
	static int getTransmitSize(AlshaKmerProps& props)
	{
		return sizeof(props.linkages) + sizeof(props.multiplicity) + sizeof(props.status);
	}
	static unsigned int serialize(AlshaKmerProps& props, uint8_t* buffer)
	{
		unsigned int offset = 0;

		//save linages
		buffer[offset++] = props.linkages[ALSHA_FORWARD_DIR];
		buffer[offset++] = props.linkages[ALSHA_REVERSE_DIR];
		
		//save multiplicity
		uint16_t data = props.multiplicity;
   		for(int i = 0; i < sizeof(props.multiplicity); i++){
       		buffer[offset++] = (uint8_t)(data & 0xFF);
       		data >>= 8;
   		}

		//save status
		buffer[offset++] = props.status;

		return offset;
	}
	static unsigned int unserialize(AlshaKmerProps& props, const uint8_t* buffer)
	{
		unsigned int offset = 0;

		//save linages
		props.linkages[ALSHA_FORWARD_DIR] = buffer[offset++];
		props.linkages[ALSHA_REVERSE_DIR] = buffer[offset++];
		
		//save multiplicity
  		props.multiplicity = 0;
   		for (int i = 0; i < sizeof(props.multiplicity); i++) {
   			uint16_t data = buffer[offset++];
       		props.multiplicity |= data << (i << 3);
   		}

		//save status
		props.status = buffer[offset++];

		return offset;

	}
	static void addMultiplicity(AlshaKmerProps& props)
	{
		if(props.multiplicity < 65535){
			props.multiplicity++;
		}
	}
	static unsigned int getMultiplicity(AlshaKmerProps& props)
	{
		return props.multiplicity;
	}
	static void setLinkage(AlshaKmerProps& props, uint8_t dir, Nucleotide nucleotide)
	{
		if(nucleotide >= ALSHA_NUCLEOTIDE_NUM){
			assert(nucleotide == ALSHA_INVALID_NUCLEOTIDE);
			return;
		}
		uint8_t mask = 1 << nucleotide;
		props.linkages[dir] |= mask;
	}
	static void clearLinkage(AlshaKmerProps& props, uint8_t dir, Nucleotide nucleotide)
	{
		assert(nucleotide < ALSHA_NUCLEOTIDE_NUM);
		uint8_t mask = ~(1 << nucleotide);

		props.linkages[dir] &= mask;
	}
	static void clearLinkage(AlshaKmerProps& props, uint8_t dir)
	{
		props.linkages[dir] = 0;
	}

	static bool getLinkage(AlshaKmerProps& props, uint8_t dir, Nucleotide nucleotide)
	{
		assert(nucleotide < ALSHA_NUCLEOTIDE_NUM);
		return (props.linkages[dir] >> nucleotide) & 1;
	}
	static int checkLinkage(AlshaKmerProps& props)
	{
		uint8_t forward = bitsCountTable[props.linkages[ALSHA_FORWARD_DIR]]; 
		uint8_t reverse = bitsCountTable[props.linkages[ALSHA_REVERSE_DIR]];

		if(forward == 0 && reverse == 0){
			return ALSHA_KMER_ISOLATED;
		}else if(forward == 0){
			return ALSHA_KMER_FORWARD_TIP;
		}else if(reverse == 0){
			return ALSHA_KMER_REVERSE_TIP;
		}
		return ALSHA_KMER_CONTINUOUS;
	}
	static int checkLinkage(AlshaKmerProps& props, uint8_t& toDir)
	{
		uint8_t forward = bitsCountTable[props.linkages[ALSHA_FORWARD_DIR]]; 
		uint8_t reverse = bitsCountTable[props.linkages[ALSHA_REVERSE_DIR]];

		if(forward == 0 && reverse == 0){
			toDir = ALSHA_INVALID_DIR;
			return ALSHA_KMER_ISOLATED;
		}else if(forward == 0){
			toDir = ALSHA_REVERSE_DIR;
			return ALSHA_KMER_FORWARD_TIP;
		}else if(reverse == 0){
			toDir = ALSHA_FORWARD_DIR;
			return ALSHA_KMER_REVERSE_TIP;
		}
		toDir = ALSHA_INVALID_DIR;
		return ALSHA_KMER_CONTINUOUS;
	}
	static int checkLinkage(AlshaKmerProps& props, uint8_t& forward, uint8_t& reverse, uint8_t& toDir)
	{	
		forward = bitsCountTable[props.linkages[ALSHA_FORWARD_DIR]]; 
		reverse = bitsCountTable[props.linkages[ALSHA_REVERSE_DIR]];

		if(forward == 0 && reverse == 0){
			toDir = ALSHA_INVALID_DIR;
			return ALSHA_KMER_ISOLATED;
		}else if(forward == 0){
			toDir = ALSHA_REVERSE_DIR;
			return ALSHA_KMER_FORWARD_TIP;
		}else if(reverse == 0){
			toDir = ALSHA_FORWARD_DIR;
			return ALSHA_KMER_REVERSE_TIP;
		}
		toDir = ALSHA_INVALID_DIR;
		return ALSHA_KMER_CONTINUOUS;
	}

	static Nucleotide getFirstLinkage(AlshaKmerProps& props, uint8_t dir)
	{
		for (int i = 0; i < ALSHA_NUCLEOTIDE_NUM; ++i){
			uint8_t mask = 1 << i;
			if((props.linkages[dir] & mask) != 0) {
				return i;
			}
		}
		return ALSHA_INVALID_NUCLEOTIDE;
	}
	static Nucleotide getNextLinkage(AlshaKmerProps& props, uint8_t dir, Nucleotide nucleotide)
	{
		for (int i = nucleotide + 1; i < ALSHA_NUCLEOTIDE_NUM; ++i) {
			uint8_t mask = 1 << i;
			if((props.linkages[dir] & mask) != 0){
				return i;
			}
		}
		return ALSHA_INVALID_NUCLEOTIDE;
	}
	static inline uint8_t getLinkageCount(AlshaKmerProps& props, uint8_t dir)
	{
		return bitsCountTable[props.linkages[dir]];
	}

	static inline bool isDead(AlshaKmerProps& props)
	{
		return (props.status & (1 << ALSHA_STATUS_DEAD_OFF));
	}
	static bool isMarked(AlshaKmerProps& props, uint8_t dir)
	{
		if(dir == ALSHA_FORWARD_DIR){
			return (props.status & (1 << ALSHA_STATUS_MARK_FORWARD_OFF));
		}
		return (props.status & (1 << ALSHA_STATUS_MARK_REVERSE_OFF));
	}
	static inline void setStatus(AlshaKmerProps& props, uint8_t kmerStatusOff)
	{
		props.status |= 1 << kmerStatusOff;
	}
	static void mark(AlshaKmerProps& props, uint8_t dir)
	{
		if(dir == ALSHA_FORWARD_DIR){
			props.status |= 1 << ALSHA_STATUS_MARK_FORWARD_OFF;
		}else{
			props.status |= 1 << ALSHA_STATUS_MARK_REVERSE_OFF;
		}
	}
	static uint8_t getStatus(AlshaKmerProps& props, uint8_t kmerStatusOff)
	{
		return (props.status >> kmerStatusOff) & 1;
	}
private:
	static const uint8_t bitsCountTable[16];
	static const uint8_t complements[16];
	static const uint8_t byteRC[256];
};

typedef pair<AlshaKmerData, AlshaKmerProps> AlshaKmer;

#endif /* ALSHAKMER_H_ */

