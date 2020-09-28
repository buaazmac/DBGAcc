/*
 * AlshaMessage.cpp
 *
 *  Created on: 05-Apr-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#include "AlshaMessage.h"
#include "AlshaComms.h"
#include "AlshaParams.h"

AlshaAllocator<char> AlshaMessage::allocator;
AlshaAllocator<uint8_t> AlshaMessage::uint8Allocator;
AlshaAllocator<AlshaKmerData> AlshaMessage::kmerDataAllocator;
AlshaMessage::AlshaMessage(int type)
{
	this->type = type;
}
uint32_t AlshaMessage::getTransmitSize()
{
	//type is unsigned char type
	return sizeof(type);
}
uint8_t AlshaMessage::readMessageType(const uint8_t* buffer)
{
	return (*buffer);
}
/***********************************************
 *	ALSHA_MSG_CONTROL_COMMAND
 **********************************************/
AlshaMessageControl::AlshaMessageControl(uint32_t cmd) : AlshaMessage(ALSHA_MSG_CONTROL_COMMAND)
{
	this->command = cmd;
	this->descriptor = 0;
}
AlshaMessageControl::AlshaMessageControl(uint32_t cmd, uint64_t descriptor) : AlshaMessage(ALSHA_MSG_CONTROL_COMMAND)
{
    this->command = cmd;
    this->descriptor = descriptor;
}
uint32_t AlshaMessageControl::getTransmitSize()
{
	return AlshaMessage::getTransmitSize() + sizeof(this->command) + sizeof(descriptor);
}
uint32_t AlshaMessageControl::serialize(uint8_t* buffer)
{
	uint32_t offset = 0;

	//save message type
	buffer[offset++] = this->type;

	//save the command
	uint32_t data = this->command;
	for(uint32_t i = 0; i < sizeof(this->command); i++){
		buffer[offset++] = (uint8_t)(data & 0xFF);
		data >>= 8;
	}
	//save the descriptor
	uint64_t data1= this->descriptor;
    for(uint32_t i = 0; i < sizeof(this->descriptor); i++){
        buffer[offset++] = (uint8_t)(data1 & 0xFF);
        data1 >>= 8;
    }

	return offset;
}
uint32_t AlshaMessageControl::unserialize(const uint8_t* buffer)
{
	uint32_t offset = 0;

	//skip the type
	offset++;

	//the commmand
	this->command = 0;
	for (uint32_t i = 0; i < sizeof(this->command); i++) {
		uint32_t data = buffer[offset++];
		this->command |= data << (i << 3);
	}
	//the descriptor
	this->descriptor = 0;
	for (uint32_t i = 0; i < sizeof(this->descriptor); i++) {
		uint64_t data = buffer[offset++];
		this->descriptor |= data << (i << 3);
	}

	return offset;
}
AlshaMessageKmerInsertion::AlshaMessageKmerInsertion(vector<string>* batch) : AlshaMessage(ALSHA_MSG_KMER_INSERTION)
{
	this->batch = batch;
}
AlshaMessageKmerInsertion::~AlshaMessageKmerInsertion()
{
	if(batch){
		delete batch;
	}
}
AlshaMessageKmerLinkage::AlshaMessageKmerLinkage(int procId): AlshaMessage(ALSHA_MSG_KMER_LINKAGE),
				dirs(uint8Allocator), bases(uint8Allocator), kmers(kmerDataAllocator)
{
	this->procId = procId;
	reserve();
}
AlshaMessageKmerLinkage::~AlshaMessageKmerLinkage()
{
	clear();
}
uint32_t AlshaMessageKmerLinkage::getTransmitSize()
{
	/*the number of kmers + kmerNum * (kmerData + dir + base)*/
	return AlshaMessage::getTransmitSize() + kmers.size() * (sizeof(AlshaKmerData)
					+ sizeof(uint8_t) + sizeof(Nucleotide)) + sizeof(uint32_t);
}
uint32_t AlshaMessageKmerLinkage::serialize(uint8_t* buffer)
{
	uint32_t offset = 0;

	//save message type
	buffer[offset++] = this->type;

	//transfer the number of kmers
	uint32_t size = kmers.size();
	for(int i = 0; i < sizeof(size); ++i){
		buffer[offset++] = (uint8_t)(size & 0xFF);
		size >>= 8;
	}

	//transfer the kmer data
	for(size_t i = 0; i < kmers.size(); ++i){

		//save the data of kmers
		AlshaKmerData kmerData = kmers[i];
		for(int j = 0; j < sizeof(kmerData); ++j){
			buffer[offset++] = (uint8_t)(kmerData & 0xFF);
			kmerData >>= 8;
		}

		//save the dirs
		buffer[offset++] = dirs[i];

		//save the bases
		buffer[offset++] = bases[i];
	}

	return offset;
}
uint32_t AlshaMessageKmerLinkage::unserialize(const uint8_t* buffer)
{
	uint32_t offset = 0;

	//skip the type
	offset++;

	//get the number of elements
	uint32_t size = 0;
	for(int i = 0; i < sizeof(size); ++i){
		uint32_t data = buffer[offset++];
		size |= data << (i << 3);
	}
	
	//resize the vectors
	kmers.resize(size);
	dirs.resize(size);
	bases.resize(size);

	//get the kmer data
	for(int i = 0; i < size; ++i){

		//get the kmer
		AlshaKmerData kmerData = 0;
		for(int j = 0; j < sizeof(kmerData); ++j){
			uint64_t data = buffer[offset++];
			kmerData |= data << (j << 3);
		}

		kmers[i] = kmerData;

		//get the dir
		dirs[i] = buffer[offset++];

		//get the base
		bases[i] = buffer[offset++];
	}

	return offset;
}

void AlshaMessageKmerLinkage::insert(AlshaKmerData kmerData, uint8_t dir, Nucleotide base)
{
	kmers.push_back(kmerData);
	dirs.push_back(dir);
	bases.push_back(base);
}
void AlshaMessageKmerLinkage::sort()
{
#if 0
	typedef pair<uint8_t, Nucleotide> AlshaKmerLinkage;
	multimap <AlshaKmerData, AlshaKmerLinkage> linkageMap;

	for(size_t i = 0; i < kmers.size(); ++i){
		AlshaKmerData kmerData = kmers[i].x;
		kmerData = (kmerData << 32) | kmers[i].y;
		linkageMap.insert(make_pair(kmerData, make_pair(dirs[i], bases[i])));
	}
	clear();

	//output the results
	multimap <AlshaKmerData, AlshaKmerLinkage>::iterator iter;
	for(iter = linkageMap.begin(); iter != linkageMap.end(); ++iter){
		CudaKmerData kmerData = {iter->first >> 32, iter->first & 0x0FFFFFFFF};
		kmers.push_back(kmerData);
		dirs.push_back(iter->second.first);
		bases.push_back(iter->second.second);
	}
	linkageMap.erase(linkageMap.begin(), linkageMap.end());
	linkageMap.clear();
#endif
}
void AlshaMessageKmerLinkage::clear()
{
	kmers.clear();
	dirs.clear();
	bases.clear();
}
void AlshaMessageKmerLinkage::reserve()
{
	//reserve 
	int maxElements = AlshaParams::getLinkageBatchSize() * 2;
	kmers.reserve(maxElements);
	dirs.reserve(maxElements);
	bases.reserve(maxElements);
}
bool AlshaMessageKmerLinkage::ready(int maxElements)
{
	if(kmers.size() > 0 && kmers.size() >= maxElements){
		return true;
	}
	return false;
}

AlshaMessageReadBatch::AlshaMessageReadBatch(vector< pair<char*, int> >*batch) : AlshaMessage(ALSHA_MSG_READ_BATCH)
{
	this->batch = batch;
}
AlshaMessageReadBatch::~AlshaMessageReadBatch()
{
	for(size_t i = 0; i < batch->size(); i++){
		AlshaUtils::memFree((*batch)[i].first);
	}
	delete batch;
}
AlshaMessageString::AlshaMessageString(string& str, NumType strIndex) : AlshaMessage(ALSHA_MSG_STRING)
{
	this->str = str;
	this->strIndex = strIndex;
}
uint32_t AlshaMessageString::getTransmitSize()
{
	return AlshaMessage::getTransmitSize();
}
uint32_t AlshaMessageString::serialize(uint8_t* buffer)
{
	return 0;
}
uint32_t AlshaMessageString::unserialize(const uint8_t* buffer)
{
	return 0;
}
/***********************************************
 *	ALSHA_MSG_REMOVE_KMER
 **********************************************/
AlshaMessageKmerRemoval::AlshaMessageKmerRemoval(int procId) : AlshaMessage(ALSHA_MSG_KMER_REMOVAL),
			linkageDirs(uint8Allocator), linkageBases(uint8Allocator), linkageKmers(kmerDataAllocator)
{
	this->procId = procId;
	reserve();
}
AlshaMessageKmerRemoval::~AlshaMessageKmerRemoval()
{
	clear();
}
uint32_t AlshaMessageKmerRemoval::getTransmitSize()
{
	/*#kmers + #kmers * sizeof(AlshaKmerData) + #linkages + #linkages *
 		(sizeof(AlshaKmerData) + sizeof(uint8_t) + sizeof(Nucleotide)
	*/
	return AlshaMessage::getTransmitSize() + sizeof(uint32_t)
		+ kmerKeys.size() * sizeof(AlshaKmerData) + sizeof(uint32_t) +
		linkageKmers.size() * (sizeof(AlshaKmerData) + sizeof(uint8_t) + sizeof(Nucleotide));
}
uint32_t AlshaMessageKmerRemoval::serialize(uint8_t* buffer)
{
	uint32_t offset = 0;

	//save message type
	buffer[offset++] = this->type;

	//transfer the number of kmer keys
	uint32_t size = kmerKeys.size();
	for(int i = 0; i < sizeof(size); ++i){
		buffer[offset++] = (uint8_t)(size & 0xFF);
		size >>= 8;
	}

	//transfer the kmer data
	for(size_t i = 0; i < kmerKeys.size(); ++i){
		//save the data of kmer keys
		AlshaKmerData kmerData = kmerKeys[i];
		for(int j = 0; j < sizeof(kmerData); ++j){
			buffer[offset++] = (uint8_t)(kmerData & 0xFF);
			kmerData >>= 8;
		}
	}

	//transfer the number of linkages
	size = linkageKmers.size();
	for(int i = 0; i < sizeof(size); ++i){
		buffer[offset++] = (uint8_t)(size & 0xFF);
		size >>= 8;
	}
	
	//transfer the linkage data
	for(size_t i = 0; i < linkageKmers.size(); ++i){
		//save the data of kmers
		AlshaKmerData kmerData = linkageKmers[i];
		for(int j = 0; j < sizeof(kmerData); ++j){
			buffer[offset++] = (uint8_t)(kmerData & 0xFF);
			kmerData >>= 8;
		}
		//save the dirs
		buffer[offset++] = linkageDirs[i];

		//save the bases
		buffer[offset++] = linkageBases[i];
	}

	return offset;
}
uint32_t AlshaMessageKmerRemoval::unserialize(const uint8_t* buffer)
{
	uint32_t offset = 0;

	//skip the type
	offset++;

	//get the number of kmers
	uint32_t size = 0;
	for(int i = 0; i < sizeof(size); ++i){
		uint32_t data = buffer[offset++];
		size |= data << (i << 3);
	}
	
	//resize the kmer vector
	kmerKeys.resize(size);

	//get the kmer data
	for(int i = 0; i < size; ++i){

		//get the kmer
		AlshaKmerData kmerData = 0;
		for(int j = 0; j < sizeof(kmerData); ++j){
			uint64_t data = buffer[offset++];
			kmerData |= data << (j << 3);
		}
		kmerKeys[i] = kmerData;
	}
	
	//get the number of linkages
	size = 0;
	for(int i = 0; i < sizeof(size); ++i){
		uint32_t data = buffer[offset++];
		size |= data << (i << 3);
	}

	//resize the linkage vectors
	linkageKmers.resize(size);
	linkageDirs.resize(size);
	linkageBases.resize(size);

	//get the linkage data
	for(int i = 0; i < size; ++i){
		//get the kmer
		AlshaKmerData kmerData = 0;
		for(int j = 0; j < sizeof(kmerData); ++j){
			uint64_t data = buffer[offset++];
			kmerData |= data << (j << 3);
		}
		linkageKmers[i] = kmerData;

		//get the dir
		linkageDirs[i] = buffer[offset++];

		//get the base
		linkageBases[i] = buffer[offset++];
	}

	return offset;
}
void AlshaMessageKmerRemoval::insertKmerKey(AlshaKmerData kmerKey)
{
	kmerKeys.push_back(kmerKey);
}
void AlshaMessageKmerRemoval::insertLinkage(AlshaKmerData kmerKey, uint8_t dir, Nucleotide base)
{
	linkageKmers.push_back(kmerKey);
	linkageDirs.push_back(dir);
	linkageBases.push_back(base);
}
void AlshaMessageKmerRemoval::clear()
{
	kmerKeys.clear();
	linkageKmers.clear();
	linkageDirs.clear();
	linkageBases.clear();
}
void AlshaMessageKmerRemoval::reserve()
{
    int maxElements = AlshaParams::getLinkageBatchSize() * 2;

	kmerKeys.reserve(maxElements);
	linkageKmers.reserve(maxElements);
	linkageDirs.reserve(maxElements);
	linkageBases.reserve(maxElements);
}
bool AlshaMessageKmerRemoval::ready(int maxElements)
{
	if(kmerKeys.size() >= maxElements || linkageKmers.size() >= maxElements){
		return true;
	}
	return false;
}
/***********************************************
 *	ALSHA_MSG_KMER_LINKAGE_REMOVAL
 **********************************************/
AlshaMessageKmerLinkageRemoval::AlshaMessageKmerLinkageRemoval(int procId) : AlshaMessage(ALSHA_MSG_KMER_LINKAGE_REMOVAL),
		linkageDirs(uint8Allocator), linkageBases(uint8Allocator), linkageKmers(kmerDataAllocator)
{
	this->procId = procId;
	reserve();
}
AlshaMessageKmerLinkageRemoval::~AlshaMessageKmerLinkageRemoval()
{
	clear();
}
uint32_t AlshaMessageKmerLinkageRemoval::getTransmitSize()
{
	/*#kmers + #linkages + #linkages *
 		(sizeof(AlshaKmerData) + sizeof(uint8_t) + sizeof(Nucleotide)
	*/
	return AlshaMessage::getTransmitSize() + sizeof(uint32_t) +
		linkageKmers.size() * (sizeof(AlshaKmerData) + sizeof(uint8_t) + sizeof(Nucleotide));
}
uint32_t AlshaMessageKmerLinkageRemoval::serialize(uint8_t* buffer)
{
	uint32_t offset = 0;

	//save message type
	buffer[offset++] = this->type;

	//transfer the number of linkages
	uint32_t size = linkageKmers.size();
	for(int i = 0; i < sizeof(size); ++i){
		buffer[offset++] = (uint8_t)(size & 0xFF);
		size >>= 8;
	}
	
	//transfer the linkage data
	for(size_t i = 0; i < linkageKmers.size(); ++i){
		//save the data of kmers
		AlshaKmerData kmerData = linkageKmers[i];
		for(int j = 0; j < sizeof(kmerData); ++j){
			buffer[offset++] = (uint8_t)(kmerData & 0xFF);
			kmerData >>= 8;
		}
		//save the dirs
		buffer[offset++] = linkageDirs[i];

		//save the bases
		buffer[offset++] = linkageBases[i];
	}

	return offset;
}
uint32_t AlshaMessageKmerLinkageRemoval::unserialize(const uint8_t* buffer)
{
	uint32_t offset = 0;

	//skip the type
	offset++;

	//get the number of linkages
	uint32_t size = 0;
	for(int i = 0; i < sizeof(size); ++i){
		uint32_t data = buffer[offset++];
		size |= data << (i << 3);
	}

	//resize the linkage vectors
	linkageKmers.resize(size);
	linkageDirs.resize(size);
	linkageBases.resize(size);

	//get the linkage data
	for(int i = 0; i < size; ++i){
		//get the kmer
		AlshaKmerData kmerData = 0;
		for(int j = 0; j < sizeof(kmerData); ++j){
			uint64_t data = buffer[offset++];
			kmerData |= data << (j << 3);
		}
		linkageKmers[i] = kmerData;

		//get the dir
		linkageDirs[i] = buffer[offset++];

		//get the base
		linkageBases[i] = buffer[offset++];
	}

	return offset;
}
void AlshaMessageKmerLinkageRemoval::insertLinkage(AlshaKmerData kmerKey, uint8_t dir, Nucleotide base)
{
	linkageKmers.push_back(kmerKey);
	linkageDirs.push_back(dir);
	linkageBases.push_back(base);
}
void AlshaMessageKmerLinkageRemoval::clear()
{
	linkageKmers.clear();
	linkageDirs.clear();
	linkageBases.clear();
}
void AlshaMessageKmerLinkageRemoval::reserve()
{
    int maxElements = AlshaParams::getLinkageBatchSize() * 2;

	linkageKmers.reserve(maxElements);
	linkageDirs.reserve(maxElements);
	linkageBases.reserve(maxElements);
}
bool AlshaMessageKmerLinkageRemoval::ready(int maxElements)
{
	if(linkageKmers.size() >= maxElements){
		return true;
	}
	return false;
}
/*******************************************************
	KMER PROPERTY REQUEST
*******************************************************/
AlshaMessageKmerPropsRequest::AlshaMessageKmerPropsRequest(AlshaKmerData kmerKey, bool rc, int procId)
					: AlshaMessage(ALSHA_MSG_KMER_PROPS_REQUEST)
{
	this->kmerKey = kmerKey;
	this->rc = rc;
	this->procId = procId;
}
uint32_t AlshaMessageKmerPropsRequest::getTransmitSize()
{
	return AlshaMessage::getTransmitSize() + sizeof(kmerKey) + sizeof(rc) + sizeof(procId);
}
uint32_t AlshaMessageKmerPropsRequest::serialize(uint8_t* buffer)
{
	unsigned int offset = 0;

	//save message type
	buffer[offset++] = this->type;

	//save the data of kmer keys
	AlshaKmerData kmerData = this->kmerKey;
	for(int i = 0; i < sizeof(kmerData); ++i){
		buffer[offset++] = (uint8_t)(kmerData & 0xFF);
		kmerData >>= 8;
	}

	//save the rc
	buffer[offset++] = this->rc;

	//save the procId
	unsigned int data = this->procId;
	for(int i = 0; i < sizeof(data); ++i){
		buffer[offset++] = (uint8_t)(data & 0xFF);
		data >>= 8;
	}
	return offset;
}
uint32_t AlshaMessageKmerPropsRequest::unserialize(const uint8_t* buffer)
{
	uint32_t offset = 0;

	//skip the type
	offset++;

	//get the kmer
	kmerKey = 0;
	for(int i = 0; i < sizeof(kmerKey); ++i){
		uint64_t data = buffer[offset++];
		kmerKey |= data << (i << 3);
	}

	//get the rc
	rc = buffer[offset++];

	//get the procId
	procId = 0;
	for(int i = 0; i < sizeof(procId); ++i){
		unsigned int data = buffer[offset++];
		procId |= data << (i << 3);
	}

	return offset;
}
/***********************************************
 *	ALSHA_MSG_KMER_PROPS_RESPONSE
 **********************************************/

AlshaMessageKmerPropsResponse::AlshaMessageKmerPropsResponse()
					: AlshaMessage(ALSHA_MSG_KMER_PROPS_RESPONSE)
{
}
AlshaMessageKmerPropsResponse::AlshaMessageKmerPropsResponse(AlshaKmerProps& props)
					: AlshaMessage(ALSHA_MSG_KMER_PROPS_RESPONSE)
{
	this->props = props;
}
uint32_t AlshaMessageKmerPropsResponse::getTransmitSize()
{
	return AlshaMessage::getTransmitSize() + AlshaKmerUtils::getTransmitSize(props);
}
uint32_t AlshaMessageKmerPropsResponse::serialize(uint8_t* buffer)
{
	uint32_t offset = 0;

	//save message type
	buffer[offset++] = this->type;

	//save the kmer props
	offset += AlshaKmerUtils::serialize(props, buffer + offset);

	return offset;
}
uint32_t AlshaMessageKmerPropsResponse::unserialize(const uint8_t* buffer)
{
	uint32_t offset = 0;

	//skip the type
	offset++;

	//extract the kmer props
	offset += AlshaKmerUtils::unserialize(props, buffer + offset);

	return offset;
}

#if 0
/***********************************************
 *	ALSHA_MSG_REMOVE_LINKAGE
 **********************************************/
AlshaMessageRemoveLinkage::AlshaMessageRemoveLinkage(AlshaKmerData kmerKey, AlshaKmerData kmerData, uint8_t dir, uint8_t linkage)
					: AlshaMessage(ALSHA_MSG_REMOVE_LINKAGE)
{
	this->kmerKey = kmerKey;
	this->kmerData = kmerData;
	this->dir = dir;
	this->linkage = linkage;
}
uint32_t AlshaMessageRemoveLinkage::getTransmitSize()
{
	return AlshaMessage::getTransmitSize() + sizeof(this->kmerKey) + sizeof(this->kmerData)
									+ sizeof(this->dir) + sizeof(this->linkage);
}
uint32_t AlshaMessageRemoveLinkage::serialize(uint8_t* buffer)
{
	uint32_t offset = 0;

	//save message type
	buffer[offset++] = this->type;

	//save the kmer key
	AlshaKmerData key = this->kmerKey;
	for(uint32_t i = 0; i < sizeof(this->kmerKey); ++i){
		buffer[offset++] = (uint8_t)(key & 0xFF);
		key >>= 8;
	}
	//save the kmer data
	AlshaKmerData data = this->kmerData;
	for(uint32_t i = 0; i < sizeof(this->kmerData); ++i){
		buffer[offset++] = (uint8_t)(data & 0xFF);
		data >>= 8;
	}
	//save the dir and linkage
	buffer[offset++] = this->dir;
	buffer[offset++] = this->linkage;

	return offset;
}
uint32_t AlshaMessageRemoveLinkage::unserialize(const uint8_t* buffer)
{
	uint32_t offset = 0;

	//skip the type
	offset++;

	//extract the kmer key
	this->kmerKey = 0;
	for(uint32_t i = 0; i < sizeof(this->kmerKey); ++i){
		uint64_t data = (unsigned char)buffer[offset++];
		this->kmerKey |= data << (i << 3);
	}

	//extract the kmer data
	this->kmerData = 0;
	for(uint32_t i = 0; i < sizeof(this->kmerData); ++i){
		uint64_t data = (unsigned char) buffer[offset++];
		this->kmerData |= data << (i << 3);
	}

	//extract the dir and linkage
	this->dir = buffer[offset++];
	this->linkage = buffer[offset++];

	return offset;
}

#endif
