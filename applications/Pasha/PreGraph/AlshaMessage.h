/*
 * AlshaMessage.h
 *
 *  Created on: 05-Apr-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHAMESSAGE_H_
#define ALSHAMESSAGE_H_

#include "AlshaTypes.h"
#include "AlshaUtils.h"
#include "AlshaKmer.h"


class AlshaMessage
{
public:
	AlshaMessage(int type);
	virtual ~AlshaMessage(){};

	static	uint8_t		readMessageType(const uint8_t* buffer);
	inline	int			getType(){ return type;}
	virtual	uint32_t	getTransmitSize() = 0;
	virtual uint32_t 	serialize(uint8_t* buffer) = 0;
	virtual uint32_t 	unserialize(const uint8_t* buffer) = 0;

protected:
	int8_t type;
	static AlshaAllocator<char> allocator;
	static AlshaAllocator<uint8_t > uint8Allocator;
	static AlshaAllocator<AlshaKmerData> kmerDataAllocator;
};

class AlshaMessageControl : public AlshaMessage
{
public:
	AlshaMessageControl() : AlshaMessage(ALSHA_MSG_CONTROL_COMMAND){}
	AlshaMessageControl(uint32_t cmd);
	AlshaMessageControl(uint32_t cmd, uint64_t descriptor);

	uint32_t 	getTransmitSize();
	uint32_t	serialize(uint8_t* buffer);
	uint32_t	unserialize(const uint8_t* buffer);

	inline uint32_t	getCommand(){return command;}
	inline uint64_t	getDescriptor(){return descriptor;}

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

private:
	uint32_t	command;
	uint64_t	descriptor;
};
class AlshaMessageKmerInsertion: public AlshaMessage
{
public:
	AlshaMessageKmerInsertion(vector<string>* batch);
	~AlshaMessageKmerInsertion();

	uint32_t getTransmitSize() {return 0;}
	uint32_t serialize(uint8_t* buffer){return 0;}
	uint32_t unserialize(const uint8_t* buffer) {return 0;}

	inline vector<string>* getBatch(){return batch;}

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

private:
	vector<string>* batch;
};

class AlshaMessageKmerLinkage : public AlshaMessage
{
public:
	AlshaMessageKmerLinkage(int procId = -1);
	virtual ~AlshaMessageKmerLinkage();

	uint32_t getTransmitSize();
	uint32_t serialize(uint8_t* buffer);
	uint32_t unserialize(const uint8_t* buffer);

	void	insert(AlshaKmerData kmerData, uint8_t dir, Nucleotide base);
	void	sort();
	void	clear();
	bool	ready(int maxElements = 1);

	inline int	getProcId(){return procId;}
	inline vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& getKmers(){return kmers;}
	inline vector<uint8_t, AlshaAllocator<uint8_t> >& getDirs(){return dirs;}
	inline vector<Nucleotide, AlshaAllocator<Nucleotide> >& getBases(){return bases;}

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

protected:
	int procId;
	vector <AlshaKmerData, AlshaAllocator<AlshaKmerData> > kmers;
	vector <uint8_t, AlshaAllocator<uint8_t> > dirs;
	vector <Nucleotide, AlshaAllocator<Nucleotide> > bases;

	void	reserve();
};
class AlshaMessageKmerLinkageRemoval : public AlshaMessage
{
public:
	AlshaMessageKmerLinkageRemoval(int procId = -1);
	~AlshaMessageKmerLinkageRemoval();

	uint32_t	getTransmitSize();
	uint32_t 	serialize(uint8_t* buffer);
	uint32_t 	unserialize(const uint8_t* buffer);

	void	insertLinkage(AlshaKmerData kmer, uint8_t dir, Nucleotide base);
	void	clear();
	bool	ready(int maxElements = 1);

	inline int getProcId(){return procId;}
	inline vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& getLinkageKmers(){return linkageKmers;}
	inline vector<uint8_t, AlshaAllocator<uint8_t> >& getLinkageDirs(){return linkageDirs;}
	inline vector<Nucleotide, AlshaAllocator<Nucleotide> >& getLinkageBases(){return linkageBases;}

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

private:
	int procId;
	//the three vectors are for k-mer linkages
	vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> > linkageKmers;
	vector<uint8_t, AlshaAllocator<uint8_t> > linkageDirs;
	vector<Nucleotide, AlshaAllocator<Nucleotide> > linkageBases;

	void	reserve();

};

class AlshaMessageKmerRemoval: public AlshaMessage
{
public:
	AlshaMessageKmerRemoval(int procId = -1);
	~AlshaMessageKmerRemoval();

	uint32_t	getTransmitSize();
	uint32_t 	serialize(uint8_t* buffer);
	uint32_t 	unserialize(const uint8_t* buffer);

	void	insertKmerKey(AlshaKmerData kmerKey);
	void	insertLinkage(AlshaKmerData kmer, uint8_t dir, Nucleotide base);
	void	clear();
	bool	ready(int maxElements = 1);

	inline int getProcId(){return procId;}
	inline vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& getKmerKeys(){return kmerKeys;};
	inline vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& getLinkageKmers(){return linkageKmers;}
	inline vector<uint8_t, AlshaAllocator<uint8_t> >& getLinkageDirs(){return linkageDirs;}
	inline vector<Nucleotide, AlshaAllocator<Nucleotide> >& getLinkageBases(){return linkageBases;}

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

private:
	int procId;
	//first: the key of the k-mer; second: it is reverse complement or not
	vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> > kmerKeys;	//the keys of kmers to be removed 
	
	//the three vectors are for k-mer linkages
	vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> > linkageKmers;
	vector<uint8_t, AlshaAllocator<uint8_t> > linkageDirs;
	vector<Nucleotide, AlshaAllocator<Nucleotide> > linkageBases;

	void	reserve();
};
class AlshaMessageReadBatch : public AlshaMessage
{
public:
	AlshaMessageReadBatch(vector< pair<char*, int> >* batch);
	~AlshaMessageReadBatch();
  uint32_t  getTransmitSize(){return 0;}
  uint32_t  serialize(uint8_t* buffer){return 0;}
  uint32_t  unserialize(const uint8_t* buffer){return 0;}

	vector< pair<char*, int> >* getBatch(){return batch;}
private:
	vector< pair<char*, int>  >* batch;
};
class AlshaMessageString : public AlshaMessage
{
public:
	AlshaMessageString(string& str, NumType stringIndex);
	uint32_t 	getTransmitSize();
	uint32_t 	serialize(uint8_t* buffer);
	uint32_t	unserialize(const uint8_t* buffer);

	inline NumType getStringIndex(){return strIndex;}
	inline string& getString(){return str;}

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

private:
	NumType strIndex;
	string str;
};
class AlshaMessageKmerPropsRequest : public AlshaMessage
{
public:
	AlshaMessageKmerPropsRequest(AlshaKmerData kmerKey, bool rc, int procId);

    uint32_t    getTransmitSize();
    uint32_t    serialize(uint8_t* buffer);
    uint32_t    unserialize(const uint8_t* buffer);

	inline AlshaKmerData getKmerKey(){return kmerKey;}
	inline uint8_t isRC(){return rc;}
	inline int getProcId(){return procId;};

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

private:
	int procId;
	uint8_t rc;
	AlshaKmerData kmerKey;
};

class AlshaMessageKmerPropsResponse : public AlshaMessage
{
public:
	AlshaMessageKmerPropsResponse();
	AlshaMessageKmerPropsResponse(AlshaKmerProps& props);

	uint32_t    getTransmitSize();
 	uint32_t    serialize(uint8_t* buffer);
 	uint32_t    unserialize(const uint8_t* buffer);

	inline AlshaKmerProps& getKmerProps(){return props;}

	void* operator new(size_t size){
		return allocator.allocate(size);
	}
	void  operator delete(void* ptr, size_t size){
		allocator.deallocate((char*)ptr, size);
	}

private:
	AlshaKmerProps props;
};

#if 0
class AlshaMessageRemoveLinkage : public AlshaMessage
{
public:
	AlshaMessageRemoveLinkage() : AlshaMessage(ALSHA_MSG_REMOVE_LINKAGE) {}
	AlshaMessageRemoveLinkage(AlshaKmerData kmerKey, AlshaKmerData kmerData, uint8_t dir, Nucleotide linkage);

	uint32_t	getTransmitSize();
	uint32_t 	serialize(uint8_t* buffer);
	uint32_t 	unserialize(const uint8_t* buffer);

	inline AlshaKmerData 	getKmerKey(){return kmerKey;}
	inline AlshaKmerData 	getKmerData(){return kmerData;}
	inline uint8_t			getDir(){return dir;}
	inline Nucleotide 		getLinkage(){return linkage;}
private:
	AlshaKmerData kmerKey;
	AlshaKmerData kmerData;
	uint8_t dir;
	Nucleotide linkage;
};
/*********************************************************
 * Message pair to obtain the kmer properties; one is sending
 * message, and the other is receiving message
 * *******************************************************/
#endif

#endif /* ALSHAMESSAGE_H_ */
