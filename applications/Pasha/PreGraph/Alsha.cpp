/*
 * Alsha.cpp
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#include "Alsha.h"
#include "AlshaUtils.h"
#include "AlshaHash.h"
#include <mpi.h>


#define KMER_BLOCK_SIZE			512
Alsha::Alsha(AlshaParams* params)
{
	this->params = params;

	//initialize the sorted vector information
	vecKmerNum = 0;
	vecKmerData = NULL;
	vecKmerProps = NULL;
	vecKmerAccelerator = NULL;

	taskThread = NULL;
	pthread_mutex_init(&threadMutex, NULL);

	//the number of kmers eroded
	numErodedKmers = 0;
	numSeqs = 0;
	preGraphNodeIndex = 0;

	//initialize comms
	AlshaComms::init();
}
Alsha::~Alsha()
{
	//auxiliary thread
	if(taskThread){
		delete taskThread;
	}

	if(this->params){
		delete this->params;
	}

	//release k-mer nodes
	releaseKmerNodes();

	pthread_mutex_destroy(&threadMutex);

	AlshaComms::destroy();
}
void Alsha::createThread()
{
	//auxiliary thread
	taskThread = new AlshaThread(this);
}
void Alsha::destroyThread()
{
	if(taskThread){
		delete taskThread;
		taskThread = NULL;
	}
}
void Alsha::lock()
{
	pthread_mutex_lock(&threadMutex);
}
void Alsha::unlock()
{
	pthread_mutex_unlock(&threadMutex);
}
int Alsha::calcDestinationProc(AlshaKmerData kmerData)
{
	unsigned int factor = 19;
	unsigned int numBytes = (AlshaParams::KMERLENGTH + 3) / 4;
	
	unsigned int sum = 0;
	for(int i = 0; i < numBytes; i++){
		sum = sum * factor + (kmerData & 0xFF);
		kmerData >>= 8;
	}
	return sum % AlshaParams::numProcs;
}
int Alsha::calcDestinationProc(unsigned int hashValue)
{
	unsigned int factor = 19;
	unsigned int numBytes = (AlshaParams::KMERLENGTH + 3) / 4;
	
	unsigned int sum = 0;
	for(int i = 0; i < numBytes; i++){
		sum = sum * factor + (hashValue & 0xFF);
		hashValue >>= 8;
	}
	return sum % AlshaParams::numProcs;
}

NumType Alsha::getKmerNum(AlshaKmerSet* kmerSet)
{
	return kmerSet->size();
}
AlshaKmerData Alsha::getKmerKey(AlshaKmerData kmer, bool& rc)
{
	//compute the reverse complement of this kmer
	AlshaKmerData rcKmer = AlshaKmerUtils::reverseComplementKmer(kmer);

	//if(AlshaKmerUtils::compareKmer(kmer, rcKmer) > 0){
	if(kmer > rcKmer){
		rc = true;
		return rcKmer;
	}
	
	rc = false;
	return kmer;
}

NumType Alsha::findKmer(AlshaKmerData kmerData)
{
	NumType left;
	NumType right;

	if(vecKmerAccelerator){
		NumType keyIndex = kmerData >> vecKmerAcceleratorShift;
		left = vecKmerAccelerator[keyIndex];
		right = vecKmerAccelerator[keyIndex + 1] - 1;
	}else{
		left = 0;
		right = vecKmerNum - 1;
	}
	
	while(left <= right){
		NumType median = (left + right) / 2;
		AlshaKmerData value = vecKmerData[median];
		if(value < kmerData){
			left = median + 1;
		}else if(value > kmerData){
			right = median - 1;
		}else{
			return median;
		}
	}

	return -1;
}

//remove a kmer path that is invoked by the task thread, so it does not allow communication functions
void Alsha::removePath(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
			AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
			AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
			AlshaQueue <NumType, AlshaAllocator<NumType> >& kmerIndexPath)
{
	//remove the k-mers and the linkages connected to them
	while(!kmerPath.empty()){
		AlshaKmer& kmer = kmerPath.back();
		AlshaKmerData& kmerKey = kmerKeyPath.back();
		NumType& kmerIndex = kmerIndexPath.back();
		
		removeKmer(msg2Proc, kmer, kmerKey, kmerIndex);

		kmerPath.pop_back();
		kmerKeyPath.pop_back();
		kmerIndexPath.pop_back();
	}

	//send the messages to processes
	for(size_t procId = 0; procId < msg2Proc.size(); ++procId){
		if(procId == AlshaParams::procId){
			kmerRemovalMsgHandler(*msg2Proc[procId]);
		}else{
			//send this message to the master thread and then send it to the destination process
			taskThread->sendResponse(msg2Proc[procId]);
			msg2Proc[procId] = new AlshaMessageKmerRemoval(procId);
		}
	}
	for(size_t i = 0; i < msg2Proc.size(); ++i){
		msg2Proc[i]->clear();
	}
	//fprintf(stderr, "function %s line: %d\n", __FUNCTION__, __LINE__);
}
void Alsha::removePathDirect(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
				AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath, 
				AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
				AlshaQueue <NumType, AlshaAllocator<NumType> >& kmerIndexPath)
{
	//remove the k-mers and the linkages connected to them
	AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >::iterator pathIter = kmerPath.begin();
	AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >::iterator keyIter = kmerKeyPath.begin();
	AlshaQueue<NumType, AlshaAllocator<NumType> >::iterator indexIter = kmerIndexPath.begin();
	for( ; pathIter != kmerPath.end(); ++pathIter){
		removeKmer(msg2Proc, *pathIter, *keyIter, *indexIter);
		++keyIter;
		++indexIter;
	}

	//send the messages to processes
	for(size_t procId = 0; procId < msg2Proc.size(); ++procId){
		if(procId == AlshaParams::procId){
			kmerRemovalMsgHandler(*msg2Proc[procId]);
		}else{
			//send this message to the process
			AlshaComms::sendMessage(procId, msg2Proc[procId]);
		}
	}
	for(size_t i = 0; i < msg2Proc.size(); ++i){
		msg2Proc[i]->clear();
	}
	//fprintf(stderr, "function %s line: %d\n", __FUNCTION__, __LINE__);
}

//remove a kmer, which does not have communications,and can be invoked by any threads
void Alsha::removeKmer(NumType kmerIndex)
{
	assert(kmerIndex >= 0 && kmerIndex < vecKmerNum);

	//set the status of this k-mer to DEAD
	lock();
	if(!AlshaKmerUtils::isDead(vecKmerProps[kmerIndex])){
		AlshaKmerUtils::setStatus(vecKmerProps[kmerIndex], ALSHA_STATUS_DEAD_OFF);
		++numErodedKmers;
	}
	unlock();

}

//remove a kmer and its linkages
void Alsha::removeKmer(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
				AlshaKmer& kmer, AlshaKmerData kmerKey, NumType kmerIndex)
{	
	int procId;
	Nucleotide base;
	AlshaKmerData kmerData;

	//determing the destination of the kmer to be removed
	if(kmerIndex >= 0){
		//this kmer belongs to the process itself;
		kmerKey = vecKmerData[kmerIndex];
		procId = AlshaParams::procId;
	}else{
		procId = calcDestinationProc(kmerKey);
	}
	msg2Proc[procId]->insertKmerKey(kmerKey);

	//determine the linkages connected to the kmer
	for(uint8_t dir = ALSHA_FORWARD_DIR; dir <= ALSHA_REVERSE_DIR; ++dir){
		//get the linkage nucleotide
		base = (dir == ALSHA_FORWARD_DIR) ? AlshaKmerUtils::getNucleotide(kmer.first, AlshaParams::KMERLENGTH - 1):
					AlshaKmerUtils::getNucleotide(kmer.first, 0);
		//iterate each linkage
		for(Nucleotide nuc = AlshaKmerUtils::getFirstLinkage(kmer.second, dir); nuc != ALSHA_INVALID_NUCLEOTIDE; 
					nuc = AlshaKmerUtils::getNextLinkage(kmer.second, dir, nuc)){

			kmerData = kmer.first;	//get a copy of the kmer data
			//for a new k-mer
			if(dir == ALSHA_FORWARD_DIR){
				AlshaKmerUtils::pushNucleotide(kmerData, nuc);
			}else{
				AlshaKmerUtils::reversePushNucleotide(kmerData, nuc);
			}
			bool rc;
			kmerKey = getKmerKey(kmerData, rc);
			procId = calcDestinationProc(kmerKey);
			
			//insert this linkage
			if(rc){
				msg2Proc[procId]->insertLinkage(kmerKey, dir, toReverseNucleotide(base));
			}else{
				msg2Proc[procId]->insertLinkage(kmerKey, toReverseDir(dir), base);
			}
		}
	}
}
void Alsha::removeKmerLinkage(vector<AlshaMessageKmerLinkageRemoval*, AlshaAllocator<AlshaMessageKmerLinkageRemoval*> >& msg2Proc,
							AlshaKmer& kmer, NumType kmerIndex, uint8_t dir)
{
	bool rc;
	int procId;
	Nucleotide base;
	AlshaKmerData kmerData, kmerKey;

	//get the linkage nucleotide
	base = (dir == ALSHA_FORWARD_DIR) ? AlshaKmerUtils::getNucleotide(kmer.first, AlshaParams::KMERLENGTH - 1):
				AlshaKmerUtils::getNucleotide(kmer.first, 0);
	//iterate each linkage
	for(Nucleotide nuc = AlshaKmerUtils::getFirstLinkage(kmer.second, dir); nuc != ALSHA_INVALID_NUCLEOTIDE; 
				nuc = AlshaKmerUtils::getNextLinkage(kmer.second, dir, nuc)){

		kmerData = kmer.first;	//get a copy of the kmer data
		//for a new k-mer
		if(dir == ALSHA_FORWARD_DIR){
			AlshaKmerUtils::pushNucleotide(kmerData, nuc);
		}else{
			AlshaKmerUtils::reversePushNucleotide(kmerData, nuc);
		}
		bool rc;
		kmerKey = getKmerKey(kmerData, rc);
		procId = calcDestinationProc(kmerKey);
		if(procId == AlshaParams::procId){
			int kmerIndex = findKmer(kmerKey);
			assert(kmerIndex >= 0);
			lock();
			if(rc){
				AlshaKmerUtils::clearLinkage(vecKmerProps[kmerIndex], dir, toReverseNucleotide(base));
			}else{
				AlshaKmerUtils::clearLinkage(vecKmerProps[kmerIndex], toReverseDir(dir), base);
			}
			unlock();
		}else{
			//insert this linkage
			if(rc){
				msg2Proc[procId]->insertLinkage(kmerKey, dir, toReverseNucleotide(base));
			}else{
				msg2Proc[procId]->insertLinkage(kmerKey, toReverseDir(dir), base);
			}
		}
	}
	
	//clear the linkages in the direction
	lock();
	AlshaKmerUtils::clearLinkage(vecKmerProps[kmerIndex], dir);
	unlock();
}
//removing a local k-mer and its linkages to the dir direction
void Alsha::removeLocalKmer(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
					AlshaKmer& kmer, NumType kmerIndex, uint8_t dir)
{

	int procId;
	Nucleotide base;
	AlshaKmerData kmerData, kmerKey;

	assert(kmerIndex >= 0);

	//this kmer belongs to the process itself;
	kmerKey = vecKmerData[kmerIndex];
	removeKmer(kmerIndex);

	if(dir == ALSHA_INVALID_DIR){
		return;
	}

	//get the linkage nucleotide
	base = (dir == ALSHA_FORWARD_DIR) ? AlshaKmerUtils::getNucleotide(kmer.first, AlshaParams::KMERLENGTH - 1):
					AlshaKmerUtils::getNucleotide(kmer.first, 0);

	//iterate each linkage
	for(Nucleotide nuc = AlshaKmerUtils::getFirstLinkage(kmer.second, dir); nuc != ALSHA_INVALID_NUCLEOTIDE; 
					nuc = AlshaKmerUtils::getNextLinkage(kmer.second, dir, nuc)){

		kmerData = kmer.first;	//get a copy of the kmer data
		//for a new k-mer
		if(dir == ALSHA_FORWARD_DIR){
			AlshaKmerUtils::pushNucleotide(kmerData, nuc);
		}else{
			AlshaKmerUtils::reversePushNucleotide(kmerData, nuc);
		}
		bool rc;
		kmerKey = getKmerKey(kmerData, rc);
		procId = calcDestinationProc(kmerKey);
		if(procId == AlshaParams::procId){
			kmerIndex = findKmer(kmerKey);
			assert(kmerIndex >= 0);
			lock();
			if(rc){
				AlshaKmerUtils::clearLinkage(vecKmerProps[kmerIndex], dir, toReverseNucleotide(base));
			}else{
				AlshaKmerUtils::clearLinkage(vecKmerProps[kmerIndex], toReverseDir(dir), base);
			}
			unlock();
		}else{
			//insert this linkage
			if(rc){
				msg2Proc[procId]->insertLinkage(kmerKey, dir, toReverseNucleotide(base));
			}else{
				msg2Proc[procId]->insertLinkage(kmerKey, toReverseDir(dir), base);
			}
		}
	}
}
//this function is used by the main function for communication purpose
void Alsha::getRemoteKmerProps(AlshaKmerProps& props, AlshaMessageKmerPropsRequest* request)
{
	AlshaKmerData kmerData = request->getKmerKey();
	int procId = request->getProcId();

	//send the request message to remote process
	AlshaComms::sendCommand(procId, ALSHA_CONTROL_KMER_PROPS_REQUEST, kmerData);
}
//This function is used by the task thread
void Alsha::getRemoteKmerProps(AlshaKmerProps& props, AlshaKmerData kmerData, bool rc, int procId)
{
	//printf("%s: send..\n", __FUNCTION__);
	//send the request to the main thread
	taskThread->sendResponse(new AlshaMessageKmerPropsRequest(kmerData, rc, procId));
	
	//wait for the data from the main thread
	AlshaMessage* msg = taskThread->recvData();
	//printf("%s: recv..\n", __FUNCTION__);

	assert(msg->getType() == ALSHA_MSG_KMER_PROPS_RESPONSE);

	AlshaMessageKmerPropsResponse* response = dynamic_cast<AlshaMessageKmerPropsResponse*>(msg);
	
	props = response->getKmerProps();

	//if it is the reverse complement
	if(rc){
		AlshaKmerUtils::reverseProps(props);
	}

	delete msg;
}
void Alsha::getRemoteKmerPropsDirect(AlshaKmerProps& props, AlshaKmerData kmerData, bool rc, int procId)
{
	uint8_t* buffer;
	unsigned int msgSize;

	assert(procId != AlshaParams::procId);

	//send the request message to remote process
	AlshaComms::sendCommand(procId, ALSHA_CONTROL_KMER_PROPS_REQUEST, kmerData);

	//receive the results
	MPI_Status status;
	
	//check the existence of message
	while(AlshaComms::checkMessage() == false);
	
	//receiving the message
	AlshaComms::recvMessage(buffer, msgSize, status);

	uint8_t type = AlshaMessage::readMessageType(buffer);
	if(type != ALSHA_MSG_KMER_PROPS_RESPONSE){
		cout <<"wroing message type (should be ALSHA_MSG_KMER_PROPS_RESPONSE): " << (int)type 
			<< " in process " << AlshaParams::procId << endl;
		AlshaUtils::exitProgram();
	}

	AlshaMessageKmerPropsResponse response;
	response.unserialize(buffer);

	props = response.getKmerProps();
	
	//if it is the reverse complement
	if(rc){
		AlshaKmerUtils::reverseProps(props);
	}
	AlshaUtils::memFree(buffer);
}

void Alsha::kmerPropsRequestHandler(int procId, AlshaKmerData kmerData)
{
	NumType kmerIndex = findKmer(kmerData);

	if(kmerIndex >= 0){
		lock();
		AlshaMessageKmerPropsResponse response(vecKmerProps[kmerIndex]);
		unlock();

		AlshaComms::sendMessage(procId, &response);
	}else{
		cout<< "Failed to find k-mer: " << kmerData << " by " << AlshaParams::procId << endl;
		AlshaKmerUtils::printKmer(kmerData);
		AlshaUtils::exitProgram();
	}
}
//the kmer linkage message handler
void Alsha::kmerLinkageMsgHandler(AlshaMessageKmerLinkage& linkage)
{
	//build the linkages of the k-mers in this message
	vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmers = linkage.getKmers();
	vector<uint8_t, AlshaAllocator<uint8_t> >& dirs = linkage.getDirs();
	vector<Nucleotide, AlshaAllocator<Nucleotide> >& bases = linkage.getBases();
	for(size_t i = 0; i < kmers.size(); ++i){
		NumType kmerIndex = findKmer(kmers[i]);
		if(kmerIndex >= 0){
			//it is possible that this kmer does not exist
			lock();
			AlshaKmerUtils::setLinkage(vecKmerProps[kmerIndex], dirs[i], bases[i]);
			unlock();
		}else{
			cout << "Failed to find k-mer: " << kmers[i] << "by " << AlshaParams::procId << endl;
			AlshaKmerUtils::printKmer(kmers[i]);
			AlshaUtils::exitProgram();
		}
	}
}
//the kmer removal message handler
void Alsha::kmerRemovalMsgHandler(AlshaMessageKmerRemoval& msg)
{
	vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeys = msg.getKmerKeys();
	vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& linkageKmers = msg.getLinkageKmers();
	vector<uint8_t, AlshaAllocator<uint8_t> >& linkageDirs = msg.getLinkageDirs();
	vector<Nucleotide, AlshaAllocator<Nucleotide> >& linkageBases = msg.getLinkageBases();

	//remove the k-mers
	for(size_t i = 0; i < kmerKeys.size(); ++i){
		NumType kmerIndex = findKmer(kmerKeys[i]);
		
		if(kmerIndex < 0){
			cout<< "Failed to remove kmer: " << kmerKeys[i] << endl;
			AlshaKmerUtils::printKmer(kmerKeys[i]);
			AlshaUtils::exitProgram();
		}
		removeKmer(kmerIndex);
	}

	//remove the linkages
	for(size_t i = 0; i < linkageKmers.size(); ++i){
		NumType kmerIndex = findKmer(linkageKmers[i]);

		if(kmerIndex < 0){
			cout<< "Failed to find linkage kmer: " << linkageKmers[i] << endl;
			AlshaKmerUtils::printKmer(linkageKmers[i]);
			AlshaUtils::exitProgram();
		}

		lock();
		AlshaKmerUtils::clearLinkage(vecKmerProps[kmerIndex], linkageDirs[i], linkageBases[i]);
		unlock();
	}
}
//the kmer linkage removal message handler
void Alsha::kmerLinkageRemovalMsgHandler(AlshaMessageKmerLinkageRemoval& msg)
{
	vector<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& linkageKmers = msg.getLinkageKmers();
	vector<uint8_t, AlshaAllocator<uint8_t> >& linkageDirs = msg.getLinkageDirs();
	vector<Nucleotide, AlshaAllocator<Nucleotide> >& linkageBases = msg.getLinkageBases();

	//remove the linkages
	for(size_t i = 0; i < linkageKmers.size(); ++i){
		NumType kmerIndex = findKmer(linkageKmers[i]);

		if(kmerIndex < 0){
			cout<< "Failed to find linkage kmer: " << linkageKmers[i] << endl;
			AlshaKmerUtils::printKmer(linkageKmers[i]);
			AlshaUtils::exitProgram();
		}

		lock();
		AlshaKmerUtils::clearLinkage(vecKmerProps[kmerIndex], linkageDirs[i], linkageBases[i]);
		unlock();
	}
}

/****************************************************
 * distribute k-mers over all processes
 ****************************************************/
static int cmpKmer(const void* a, const void* b)
{
	AlshaKmerData aData = *((AlshaKmerData*)a);
	AlshaKmerData bData = *((AlshaKmerData*)b);
	if(aData > bData){
		return 1;
	}else if(aData < bData){
		return -1;
	}
	return 0;
}
void Alsha::distributeKmers()
{
	FILE* file;

	//create the kmer set
#ifdef HAVE_GOOGLE_SPARSE_HASH_SET
	AlshaKmerHashSet* kmerSet = new AlshaKmerHashSet(200000000);
	kmerSet->max_load_factor(0.432);
	//set the deletete key
	kmerSet->set_deleted_key((uint64_t)-1);
#else
	AlshaKmerHashSet* kmerSet = new AlshaKmerHashSet;
#endif

	AlshaUtils::memVerify(kmerSet);

	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"Start distributing k-mer among MPI processes" << endl;
	}

	vector< pair<char*, int> >* batch;
	AlshaMessageReadBatch* batchMsg;
  //send request
	taskThread->sendRequest(new AlshaMessageControl(ALSHA_CONTROL_READ_FILE));

	//get the error positions of a batch of reads
	while(1){
	 	AlshaMessage* response = taskThread->recvResponse();
 		uint8_t type = response->getType();
		if(type == ALSHA_MSG_CONTROL_COMMAND){
 			AlshaMessageControl* control = dynamic_cast<AlshaMessageControl*>(response);
 			if(control->getCommand() == ALSHA_CONTROL_READ_FILE_DONE){
   			delete control;
    		break;
    		}else{
      	  cerr << "Unrecognized message command " << control->getCommand() << "in filtering out function" << endl;
       	 AlshaUtils::exitProgram();
      	}
    }else if (type == ALSHA_MSG_READ_BATCH){
     	batchMsg = dynamic_cast<AlshaMessageReadBatch*>(response);
      batch = batchMsg->getBatch();

			size_t readIdx;
			for(readIdx = 0; readIdx < batch->size(); readIdx++){
				const char* sequence = (*batch)[readIdx].first;
				int readLength = (*batch)[readIdx].second;

       	AlshaKmerData kmer = 0;
       	AlshaKmerData antiKmer = 0;
       	Nucleotide nucleotide;

        	//initialize the first k-mer
        	for(int i = 0; i < AlshaParams::KMERLENGTH - 1; ++i){
            	nucleotide = sequence[i];
            	AlshaKmerUtils::pushNucleotide(kmer, nucleotide);
           		AlshaKmerUtils::reversePushNucleotide(antiKmer, toReverseNucleotide(nucleotide));
        	}

 	       for(int i = AlshaParams::KMERLENGTH - 1; i < readLength; ++i){
    	   		nucleotide = sequence[i];

        	 	//get a k-mer and its reverse complement
         		AlshaKmerUtils::pushNucleotide(kmer, nucleotide);
	          AlshaKmerUtils::reversePushNucleotide(antiKmer, toReverseNucleotide(nucleotide));

						AlshaKmerData kmerData = kmer;
        	 	//if(AlshaKmerUtils::compareKmer(kmer, antiKmer) > 0){
						if(kmer > antiKmer){
							kmerData = antiKmer;
						}
						if (calcDestinationProc(kmerData) == AlshaParams::procId){
							kmerSet->insert(kmerData);
            }
        	}
				}
			this->numSeqs += batch->size();
			if(this->numSeqs % 20000000 == 0){
				cout << " Processed " << this->numSeqs << " reads by " << AlshaParams::procId << endl;
			}
		}
		delete batchMsg;
	}

	/********************************************
 	Each process stores its kmers to a k-mer file, 
	and release the memory occupied by the k-mer set
 	********************************************/
	char kmerFile[256];
	sprintf(kmerFile, "%s/pasha-Kmers-%d", AlshaParams::getPathName().c_str(), AlshaParams::procId);
	
	file = fopen(kmerFile, "wb");
	if(!file){
		cout<< "Failed to open file: " << kmerFile << endl;
		AlshaUtils::exitProgram();
	}
	//write the number of kmers
	vecKmerNum = kmerSet->size();
	if(fwrite(&vecKmerNum, sizeof(NumType), 1, file) != 1){
			AlshaUtils::exitProgram("fwrite failed");
	}
	//iterate each k-mer
	vecKmerNum = 0;	//record the number of kmers
	int kmerBlockSize = 0;
	AlshaKmerData kmerBlock[KMER_BLOCK_SIZE];
	for(AlshaKmerHashSet::const_iterator iter = kmerSet->begin(); iter != kmerSet->end(); ++iter){
		if(kmerBlockSize == KMER_BLOCK_SIZE){
			if(fwrite(kmerBlock, sizeof(AlshaKmerData), kmerBlockSize, file) != kmerBlockSize){
				AlshaUtils::exitProgram("fwrite failed");
			}
			kmerBlockSize = 0;
		}
		kmerBlock[kmerBlockSize++] = *iter;
		vecKmerNum++;
	}
	if(kmerBlockSize > 0){
  	if(fwrite(kmerBlock, sizeof(AlshaKmerData), kmerBlockSize, file) != kmerBlockSize){
     	AlshaUtils::exitProgram("fwrite failed...");
   	}
	}
	fclose(file);

	//calculate the total number of kmers
	NumType totalKmerNum = 0;
	NumType kmerNum = kmerSet->size();
	MPI_Reduce(&kmerNum, &totalKmerNum, 1, 	MPI_LONG_LONG, MPI_SUM, ALSHA_MASTER_RANK, MPI_COMM_WORLD);
	if(AlshaParams::procId == ALSHA_MASTER_RANK){

		//writes the reads information to a file
		string infoFileName = params->getReadsFileName();
		infoFileName.append(".info");
		file = fopen(infoFileName.c_str(), "wb");
		if(!file){
			cout <<"Failed to create file: " << infoFileName << endl;
			AlshaUtils::exitProgram();
		}
		fprintf(file, "%lld\t%lld\t%d\t%d\n", this->numSeqs, totalKmerNum, AlshaParams::KMERLENGTH, AlshaParams::numProcs);
		fclose(file);

		cout <<"Total number of sequences: " << this->numSeqs << endl;
		cout<< "Total number of k-mers: " << totalKmerNum << endl;
	}

	//delete the kmer set
	kmerSet->erase(kmerSet->begin(), kmerSet->end());
	kmerSet->clear();
	delete kmerSet;

	//synchronize all processes
	MPI_Barrier(MPI_COMM_WORLD);
}
void Alsha::buildSortedKmerVector()
{
	int numProcs;
	int kmerLength;
	FILE* file;
	char buffer[1024];
	//read the number of sequences
	string infoFileName = params->getReadsFileName();
	infoFileName.append(".info");
	file = fopen(infoFileName.c_str(), "rb");
	if(!file){
		cout <<"Failed to open file: " << infoFileName << endl;
		AlshaUtils::exitProgram();
	}
	if(fgets(buffer, 1024, file) == NULL){
		cerr <<"Failed to read file " << infoFileName << endl;
		AlshaUtils::exitProgram();
	}
	fclose(file);

	sscanf(buffer, "%lld\t%*lld\t%d\t%d\n", &this->numSeqs, &kmerLength, &numProcs);
	//check the number of processors
	if(numProcs != AlshaParams::numProcs){
		cout << "Different number of processes " << AlshaParams::numProcs << " != " << numProcs <<" (old value)" << endl;
		AlshaUtils::exitProgram();
	}
	//set the kmer size
	AlshaParams::setKmerSize(kmerLength);
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout << "Number of sequences: " << this->numSeqs << endl;
		cout << "Kmer size: " << AlshaParams::KMERLENGTH << endl;
	}
 	/********************************************
	Build a sorted vector of k-mers from the k-mer file
 	********************************************/
	char kmerFile[256];
  sprintf(kmerFile, "%s/pasha-Kmers-%d", AlshaParams::getPathName().c_str(), AlshaParams::procId);

	//load all k-mers from the file to the sorted vector
	file = fopen(kmerFile, "rb");
	if(!file){
		cout<< "Failed to open file: " << kmerFile << endl;
		AlshaUtils::exitProgram();
	}
	rewind(file);

	//read the total number of kmers
	if(fread(&vecKmerNum, sizeof(NumType), 1, file) != 1){
		AlshaUtils::exitProgram("fread failed");
	}
	
	NumType totalKmerNum = 0;
	MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&vecKmerNum, &totalKmerNum, 1,  MPI_LONG_LONG, MPI_SUM, ALSHA_MASTER_RANK, MPI_COMM_WORLD);
  if(AlshaParams::procId == ALSHA_MASTER_RANK){
    cout<< "Total number of k-mer nodes: " << totalKmerNum << endl;
	}

	//build a sorted vector for all k-mers
	vecKmerData = (AlshaKmerData*)AlshaUtils::memAlloc((vecKmerNum + 1) * sizeof(AlshaKmerData));
	vecKmerProps = (AlshaKmerProps*)AlshaUtils::memAlloc((vecKmerNum + 1) * sizeof(AlshaKmerProps));
	memset(vecKmerProps, 0, (vecKmerNum + 1)* sizeof(AlshaKmerProps));
	
	//read all k-mers from the file
	NumType kmerIndex= 0;
  int kmerBlockSize = 0;
  AlshaKmerData kmerBlock[KMER_BLOCK_SIZE];
	do{
		kmerBlockSize = fread(kmerBlock, sizeof(AlshaKmerData), KMER_BLOCK_SIZE, file);
		for(int i = 0; i < kmerBlockSize; i++){
			vecKmerData[kmerIndex] = kmerBlock[i];
			++kmerIndex;
		}
	}while(kmerBlockSize == KMER_BLOCK_SIZE);

	if(kmerIndex != vecKmerNum){
		cerr << "Incorrect number of kmers are read from the file " << kmerIndex << " / " << vecKmerNum << endl;
		AlshaUtils::exitProgram("");
	}
	fclose(file);

	//sorting the vector
	qsort(vecKmerData, vecKmerNum, sizeof(AlshaKmerData),  cmpKmer);

	//build kmer vector accelerator
	buildKmerVectorAccelerator();

	//synchronize all processes
	MPI_Barrier(MPI_COMM_WORLD);
}

/****************************************************
 * build linkages for all nodes (k-mers)
 ****************************************************/
void Alsha::buildLinkages()
{
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"Start building k-mer linkages from the input short reads" << endl;
	}

	//building the linkages from the reads
	buildLinkagesCore(AlshaParams::getReadsFileName().c_str());

	//synchronize all processes
	MPI_Barrier(MPI_COMM_WORLD);

}
void Alsha::buildLinkagesCore(const char* readsFileName)
{
	NumType numOfReads = 0;
	vector< pair<char*, int> >* batch;
	AlshaMessageReadBatch* batchMsg;
  //send request
	taskThread->sendRequest(new AlshaMessageControl(ALSHA_CONTROL_READ_FILE));

	//get the error positions of a batch of reads
	while(1){
	 	AlshaMessage* response = taskThread->recvResponse();
 		uint8_t type = response->getType();
		if(type == ALSHA_MSG_CONTROL_COMMAND){
 			AlshaMessageControl* control = dynamic_cast<AlshaMessageControl*>(response);
 			if(control->getCommand() == ALSHA_CONTROL_READ_FILE_DONE){
   			delete control;
    		break;
    		}else{
      	  cerr << "Unrecognized message command " << control->getCommand() << "in filtering out function" << endl;
       	 AlshaUtils::exitProgram();
      	}
    }else if (type == ALSHA_MSG_READ_BATCH){
     	batchMsg = dynamic_cast<AlshaMessageReadBatch*>(response);
      batch = batchMsg->getBatch();

    	size_t readIdx;
      for(readIdx = 0; readIdx < batch->size(); readIdx++){
        const char* sequence = (*batch)[readIdx].first;
        int readLength = (*batch)[readIdx].second;

     	  AlshaKmerData kmer = 0;
       	AlshaKmerData antiKmer = 0;
        Nucleotide nucleotide, leftBase, rightBase;

        //initialize the first k-mer
   	    for(int i = 0; i < AlshaParams::KMERLENGTH - 1; ++i){
     	    nucleotide = sequence[i];
         	AlshaKmerUtils::pushNucleotide(kmer, nucleotide);
         	AlshaKmerUtils::reversePushNucleotide(antiKmer, toReverseNucleotide(nucleotide));
   	    }
				int leftBaseIndex = 0;
				int rightBaseIndex = AlshaParams::KMERLENGTH;
       	for(int i = AlshaParams::KMERLENGTH - 1; i < readLength; ++i){
         	nucleotide = sequence[i];

					//get the left and right bases
					leftBase = (i >= AlshaParams::KMERLENGTH) ? sequence[leftBaseIndex++] : ALSHA_INVALID_NUCLEOTIDE;
					rightBase = ( i + 1 < readLength) ? sequence[rightBaseIndex++] : ALSHA_INVALID_NUCLEOTIDE;

    	   	//get a k-mer and its reverse complement
        	AlshaKmerUtils::pushNucleotide(kmer, nucleotide);
         	AlshaKmerUtils::reversePushNucleotide(antiKmer, toReverseNucleotide(nucleotide));
	
					bool rc = false;
					AlshaKmerData kmerData = kmer;
       		//if(AlshaKmerUtils::compareKmer(kmer, antiKmer) > 0){
					if(kmer > antiKmer){	//using the reverse complement as the canonical kmer
						rc = true;
						kmerData = antiKmer;
					}
					if (calcDestinationProc(kmerData) == AlshaParams::procId){
						NumType kmerIndex = findKmer(kmerData);	
						if(kmerIndex < 0){
							cout <<"Failed to find kmer: " << kmerData << endl;
							AlshaKmerUtils::printKmer(kmerData);
							AlshaUtils::exitProgram();
						}	
						if(!rc){
							AlshaKmerUtils::setLinkage(vecKmerProps[kmerIndex], ALSHA_FORWARD_DIR, rightBase);
							AlshaKmerUtils::setLinkage(vecKmerProps[kmerIndex], ALSHA_REVERSE_DIR, leftBase);
						}else{
							if(rightBase != ALSHA_INVALID_NUCLEOTIDE){
								AlshaKmerUtils::setLinkage(vecKmerProps[kmerIndex], ALSHA_REVERSE_DIR, toReverseNucleotide(rightBase));
							}	
							if(leftBase != ALSHA_INVALID_NUCLEOTIDE){
								AlshaKmerUtils::setLinkage(vecKmerProps[kmerIndex], ALSHA_FORWARD_DIR, toReverseNucleotide(leftBase));
							}
						}
						AlshaKmerUtils::addMultiplicity(vecKmerProps[kmerIndex]);
          }
       	}
			}
			numOfReads += batch->size();
			if(numOfReads % 20000000 == 0){
				cout << " Processed " << numOfReads << " reads by " << AlshaParams::procId << endl;
			}
		}
		delete batchMsg;
	}
}

/***********************************************
	Estimate coverage from k-mer multiplicities
************************************************/
void Alsha::estimateCoverage()
{
	//estimate the coverage
	if(AlshaParams::numProcs > 1){
		vector<uint64_t> hVec = buildHistogram();
		vector<uint64_t> hVecSum (hVec.size());
		MPI_Allreduce(&hVec[0], &hVecSum[0], hVec.size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
		
		setCoverage(AlshaHistogram(hVecSum));
	}else{
		setCoverage(buildHistogram());
	}

	//synchronize all processes
	MPI_Barrier(MPI_COMM_WORLD);
}
/*******************************************
	Remove low-multiplicity tip k-mers
*******************************************/
void Alsha::erodeKmers()
{
	int erodeRounds;
	NumType	erodeKmers = 0;
	NumType nodeErodedKmers = 0;

	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"Start remvoing low-multiplicity ( < " << AlshaParams::ERODEMULTI << " ) tip k-mers" << endl;
	}

	for(erodeRounds = 0; ; ++erodeRounds){
		numErodedKmers = 0;
		if(AlshaParams::numProcs > 1){
			//send message to start the tip clipping
			taskThread->sendRequest(new AlshaMessageControl(ALSHA_CONTROL_KMER_ERODE_START));

			//initailize the number of eroded k-mers in this iteration by itself
			int numProcs = AlshaParams::numProcs;
			while(numProcs > 0){
				erodeKmersLoop(numProcs);
			}
			//send the completion command to all other processes
			if(AlshaParams::procId == ALSHA_MASTER_RANK){
				AlshaComms::bcastCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_KMER_ERODE_DONE);
			}
			//syncrhonize all the processes
			MPI_Barrier(MPI_COMM_WORLD);

 			//check the completion of all message in flight
   			while(!AlshaComms::checkCompletion()){
     			erodeKmersLoop(numProcs);
 			}
       	}else{
			//directly invoke the core function
			erodeKmersCore();
		}
		NumType sumErodedKmers;
		MPI_Allreduce(&numErodedKmers, &sumErodedKmers, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		if(sumErodedKmers == 0){
			break;
		}
		if(AlshaParams::procId == ALSHA_MASTER_RANK){
			cout << sumErodedKmers <<" kmers are eroded in round " << erodeRounds + 1 << endl;
		}
		erodeKmers += sumErodedKmers;
		nodeErodedKmers += numErodedKmers;
		if(nodeErodedKmers >= vecKmerNum / 2){
			cerr << "resizing the kmer vector (size: " << vecKmerNum - nodeErodedKmers << " / " << vecKmerNum << " )" << endl;
			resizeKmerVector();
			nodeErodedKmers = 0;
		}
		//synchronization
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"Totally removed " << erodeKmers << " kmers in " << erodeRounds << " rounds" << endl;
	}

	//synchronize all threads
	MPI_Barrier(MPI_COMM_WORLD);

}
void Alsha::erodeKmersLoop(int& numProcs)
{
	//waiting for messages from procId
	int flag;
	unsigned int msgSize;
	uint8_t* msgBuffer;
	MPI_Status status;
	
	flag = taskThread->getResponseNum();
	if(flag > 0){
		AlshaMessage* response = taskThread->recvResponse();
		uint8_t type = response->getType();
		if(type == ALSHA_MSG_CONTROL_COMMAND){
			AlshaMessageControl* control = dynamic_cast<AlshaMessageControl*>(response);
			switch(control->getCommand()){
			case ALSHA_CONTROL_KMER_ERODE_DONE:
				//send message to the master process
				if(AlshaParams::procId != ALSHA_MASTER_RANK){
					AlshaComms::sendCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_KMER_ERODE_DONE);
				}else{
					--numProcs;
				}
				break;
			default:
				fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
				cout<<"unknown message: " << type << endl;
				AlshaUtils::exitProgram();
			}
		}else if(type == ALSHA_MSG_KMER_REMOVAL){
			AlshaMessageKmerRemoval* removal = dynamic_cast<AlshaMessageKmerRemoval*>(response);
		
			if(removal->getProcId() != AlshaParams::procId){
				AlshaComms::sendMessage(removal->getProcId(), removal);
			}else{
				kmerRemovalMsgHandler(*removal);
			}
		}else{
			fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
			cout<<"unknown message: " << type << endl;
			AlshaUtils::exitProgram();
		}	
	
		delete response;
	}
	//check the message
	if(AlshaComms::checkMessage() == false){
		return;
	}
	//receive the message
	AlshaComms::recvMessage(msgBuffer, msgSize, status);

	//get the message type
	int msgType = AlshaMessage::readMessageType(msgBuffer);
	if(msgType == ALSHA_MSG_KMER_REMOVAL){
		AlshaMessageKmerRemoval removal;
				
		//unserialize the message
		removal.unserialize(msgBuffer);
		
		//handle this message
		kmerRemovalMsgHandler(removal);
	}else if(msgType == ALSHA_MSG_CONTROL_COMMAND){
		AlshaMessageControl control;
				
		//unserialize the message
		control.unserialize(msgBuffer);
				
		switch(control.getCommand()){
		case ALSHA_CONTROL_KMER_ERODE_DONE:
			if(AlshaParams::procId == ALSHA_MASTER_RANK){
				--numProcs;
			}else{
				numProcs  = 0;
			}
			break;
		default:
			fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
			cout<< "Unmatched command (" << AlshaParams::procId << "): " << control.getCommand() << endl;
			AlshaUtils::exitProgram();
			break;
		}
	}else{
		fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
		cout<< "Unknown message type: " << msgType << endl;
		AlshaUtils::exitProgram();
	}

	//release msgBuffer
	AlshaUtils::memFree(msgBuffer);

}
void Alsha::erodeKmersCore()
{
	int procId;

	const unsigned int maxElements = AlshaParams::getLinkageBatchSize();
	vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> > msg2Proc;
	msg2Proc.resize(AlshaParams::numProcs);
	for(procId = 0; procId < msg2Proc.size(); ++procId){
		msg2Proc[procId] = new AlshaMessageKmerRemoval(procId);
	}

	//iterate each k-mer to construct linages between k-mer graph nodes
	int iterations = 0;
	AlshaKmer kmer;
	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
		lock();
		kmer = make_pair(vecKmerData[kmerIndex], vecKmerProps[kmerIndex]);
		unlock();

		if(AlshaKmerUtils::isDead(kmer.second)){
			continue;
		}
		//check the linkage of this kmer
		uint8_t toDir;
		int linkage = AlshaKmerUtils::checkLinkage(kmer.second, toDir);
		if(linkage == ALSHA_KMER_CONTINUOUS){
			continue;
		}
        if(linkage == ALSHA_KMER_ISOLATED
			|| (AlshaKmerUtils::getMultiplicity(kmer.second) < AlshaParams::ERODEMULTI
            && AlshaKmerUtils::getLinkageCount(kmer.second, toDir) < 2)){

			//remove the kmer and the linkage
			removeLocalKmer(msg2Proc, kmer, kmerIndex, toDir);
			
			iterations++;
		}
		//check if this message can be sent out
		if(iterations >= maxElements){
			for(size_t procId = 0; procId < msg2Proc.size(); ++procId){
				if(msg2Proc[procId]->ready(maxElements) == true){
					taskThread->sendResponse(msg2Proc[procId]);
					msg2Proc[procId] = new AlshaMessageKmerRemoval(procId);
				}
			}
			iterations = 0;
		}
	}

	//flush all the message vecotrs
	for(size_t procId = 0; procId < msg2Proc.size(); ++procId){
		if(msg2Proc[procId]->ready() == true){
			taskThread->sendResponse(msg2Proc[procId]);
			msg2Proc[procId] = NULL;
		}
	}
	msg2Proc.clear();
}

/****************************************************
 * clip short and low-coverage dead ends
 ****************************************************/
void Alsha::clipTips()
{
	int tipCutoff = AlshaParams::getTipCutoff();
	int clipTipsRounds;
	NumType	clippedKmers = 0;
	NumType	nodeClippedKmers = 0;


	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"Start clip short and low-coverage dead ends" << endl;
	}
	
	for(clipTipsRounds = 0; ; ++clipTipsRounds){
		numErodedKmers = 0;
		//send message to start the tip clipping
		taskThread->sendRequest(new AlshaMessageControl(ALSHA_CONTROL_CLIP_TIPS_START, tipCutoff));

		//initailize the number of eroded k-mers in this iteration by itself
		int numProcs = AlshaParams::numProcs;
		while(numProcs > 0){
			clipTipsLoop(numProcs);
		}
		//send the completion command to all other processes
		if(AlshaParams::procId == ALSHA_MASTER_RANK){
			AlshaComms::bcastCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_CLIP_TIPS_DONE);
		}
		//syncrhonize all the processes
		MPI_Barrier(MPI_COMM_WORLD);

   	//check the completion of all message in flight
  	while(!AlshaComms::checkCompletion()){
			clipTipsLoop(numProcs);
		}

		NumType sumErodedKmers;
		MPI_Allreduce(&numErodedKmers, &sumErodedKmers, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		if(sumErodedKmers == 0){
			break;
		}
		if(AlshaParams::procId == ALSHA_MASTER_RANK){
			cout << sumErodedKmers <<" kmers are eroded in round " << clipTipsRounds + 1 << endl;
		}
		clippedKmers += sumErodedKmers;
		nodeClippedKmers += numErodedKmers;
    if(nodeClippedKmers >= vecKmerNum / 2){
      cerr << "resizing the kmer vector (size: " << vecKmerNum - nodeClippedKmers << " / " << vecKmerNum << " )" << endl;
      resizeKmerVector();
      nodeClippedKmers = 0;
    }
		//synchronization
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"Totally removed " << clippedKmers << " kmers in " << clipTipsRounds << " rounds" << endl;
		cout << "resizing the k-mer vector" << endl;
	}

	//resize the kmer vector
	resizeKmerVector();

	//synchronize all threads
	MPI_Barrier(MPI_COMM_WORLD);
}
void Alsha::clipTipsLoop(int& numProcs)
{
	//waiting for messages from procId
	int flag;
	unsigned int msgSize;
	uint8_t* msgBuffer;
	MPI_Status status;
	
	//check the message queue from the auxiliary thread
	flag = taskThread->getResponseNum();
	if(flag > 0){
		AlshaMessage* response = taskThread->recvResponse();
		uint8_t type = response->getType();
		if(type == ALSHA_MSG_CONTROL_COMMAND){
			AlshaMessageControl* control = dynamic_cast<AlshaMessageControl*>(response);
			switch(control->getCommand()){
			case ALSHA_CONTROL_CLIP_TIPS_DONE:
				//send message to the master process
				if(AlshaParams::procId != ALSHA_MASTER_RANK){
					AlshaComms::sendCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_CLIP_TIPS_DONE);
				}else{
					--numProcs;
				}
				break;
			}
		}else if(type == ALSHA_MSG_KMER_REMOVAL){
			AlshaMessageKmerRemoval* removal = dynamic_cast<AlshaMessageKmerRemoval*>(response);
			
			assert(AlshaParams::procId != removal->getProcId());
			AlshaComms::sendMessage(removal->getProcId(), removal);

		}else if(type == ALSHA_MSG_KMER_PROPS_REQUEST){
			AlshaKmerProps props;
			AlshaMessageKmerPropsRequest* request = dynamic_cast<AlshaMessageKmerPropsRequest*>(response);
			//get the props
			getRemoteKmerProps(props, request);
		}else{
			fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
			cout<<"unknown message: " << type << endl;
			AlshaUtils::exitProgram();
		}	
		
		delete response;
	}
	//check the existence of message
	if(AlshaComms::checkMessage() == false){
		return;
	}
	//receive the message
	AlshaComms::recvMessage(msgBuffer, msgSize, status);

	//get the message type
	int msgType = AlshaMessage::readMessageType(msgBuffer);
	if(msgType == ALSHA_MSG_KMER_REMOVAL){
		AlshaMessageKmerRemoval removal;
		
		//unserialize the message
		unsigned int offset = removal.unserialize(msgBuffer);
		assert(offset == msgSize);

		//do it by itself
		kmerRemovalMsgHandler(removal);

	}else if(msgType == ALSHA_MSG_KMER_PROPS_RESPONSE){
		AlshaMessageKmerPropsResponse* response = new AlshaMessageKmerPropsResponse;
		
		//unserialize the message
		unsigned int offset = response->unserialize(msgBuffer);
		assert(offset == msgSize);

		//send the kmer props to the task thread
		taskThread->sendData(response);

	}else if(msgType == ALSHA_MSG_CONTROL_COMMAND){
		AlshaMessageControl control;
	
		//unserialize the message
		unsigned int offset = control.unserialize(msgBuffer);
		if(offset != msgSize){
			cout<<"offset: " << offset << " msgSize: " << msgSize << endl;
			AlshaUtils::exitProgram();
		}
		
		switch(control.getCommand()){
		case ALSHA_CONTROL_CLIP_TIPS_DONE:
			//indicating the completion
			if(AlshaParams::procId == ALSHA_MASTER_RANK){
				--numProcs;
			}else{
				numProcs = 0;
			}
			break;

		case ALSHA_CONTROL_KMER_PROPS_REQUEST:

			//deal with this request
			kmerPropsRequestHandler(status.MPI_SOURCE, control.getDescriptor());

			break;
		default:
			fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
			cout<< "Unmatched command (" << AlshaParams::procId << "): " << control.getCommand() << endl;
			AlshaUtils::exitProgram();
			break;
		}
	}else{
		fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
		cout<< "Unknown message type: " << msgType << endl;
		AlshaUtils::exitProgram();
	}
	AlshaUtils::memFree(msgBuffer);
}
void Alsha::clipTipsCore(int tipCutoff)
{
	AlshaKmer kmer;
	AlshaKmerSet kmerSet;
	AlshaQueue <AlshaKmer, AlshaAllocator<AlshaKmer> > kmerPath;
	AlshaQueue <AlshaKmerData, AlshaAllocator<AlshaKmerData> > kmerKeyPath;
	AlshaQueue <NumType, AlshaAllocator<NumType> > kmerIndexPath;

	//message for kmer and linkage removal
    vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> > msg2Proc;
    msg2Proc.resize(AlshaParams::numProcs);
    for(int procId = 0; procId < msg2Proc.size(); ++procId){
        msg2Proc[procId] = new AlshaMessageKmerRemoval(procId);
    }

	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
		lock();
		kmer = make_pair(vecKmerData[kmerIndex], vecKmerProps[kmerIndex]);
		unlock();

		//check the validity of this k-mer
		if(AlshaKmerUtils::isDead(kmer.second)){
			continue;
		}

		//check the linkage of this k-mer
		uint8_t toDir;
		int linkage = AlshaKmerUtils::checkLinkage(kmer.second, toDir);
		if(linkage == ALSHA_KMER_CONTINUOUS){
			continue;
		}
		
		//if this k-mer is isolated, remove it
		if(linkage == ALSHA_KMER_ISOLATED){
			removeKmer(kmerIndex);
			continue;
		}

		//add the current k-mer to the k-mer node path
		kmerSet.insert(kmer.first);
		kmerPath.push_back(kmer);
		kmerKeyPath.push_back((AlshaKmerData)-1); //for local k-mers, no need of this key
		kmerIndexPath.push_back(kmerIndex);

		//check the multiplicity of this node
		if(isEligibleTip(kmerSet, kmerPath, kmerKeyPath, kmerIndexPath, tipCutoff, toDir)){
			//remove all k-mer nodes in the k-mer path
			removePath(msg2Proc, kmerPath, kmerKeyPath, kmerIndexPath);
		}

		//clear the kmer path
		kmerSet.erase(kmerSet.begin(), kmerSet.end());
		kmerSet.clear();
		kmerPath.erase(kmerPath.begin(), kmerPath.end());
		kmerPath.clear();
		kmerKeyPath.erase(kmerKeyPath.begin(), kmerKeyPath.end());
		kmerKeyPath.clear();
		kmerIndexPath.erase(kmerIndexPath.begin(), kmerIndexPath.end());
		kmerIndexPath.clear();
	}

	//release the msg2Proc vector
	for(size_t i = 0; i < msg2Proc.size(); ++i){
		if(msg2Proc[i]){
			delete msg2Proc[i];
			msg2Proc[i] = NULL;
		}
	}
	msg2Proc.clear();
}

bool Alsha::isEligibleTip(AlshaKmerSet& kmerSet, AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
				AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
				AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath, int tipCutoff, uint8_t dir)
{
	bool rc;
	int procId;
	NumType kmerIndex;
	AlshaKmer* kmer;
	int totalLength = AlshaParams::KMERLENGTH;

	uint8_t reverseDir = toReverseDir(dir);

	kmer = &kmerPath.front();
	while(AlshaKmerUtils::getLinkageCount(kmer->second, dir) < 2
			&& AlshaKmerUtils::getLinkageCount(kmer->second, reverseDir) < 2){

		Nucleotide base = AlshaKmerUtils::getFirstLinkage(kmer->second, dir);
		//reaching the end of one linear path
		if(base == ALSHA_INVALID_NUCLEOTIDE){
			return true;
		}

		//get the k-mer linked to from the current k-mer
		AlshaKmer toKmer;
		toKmer.first = kmer->first;
		if(dir == ALSHA_FORWARD_DIR){
			AlshaKmerUtils::pushNucleotide(toKmer.first, base);
		}else{
			AlshaKmerUtils::reversePushNucleotide(toKmer.first, base);
		}
		//get the k-mer information
		AlshaKmerData kmerKey = getKmerKey(toKmer.first, rc);
		procId = calcDestinationProc(kmerKey);

		if(procId == AlshaParams::procId){
			kmerIndex = findKmer(kmerKey);
			assert(kmerIndex >= 0);

			lock();
			toKmer.second = vecKmerProps[kmerIndex];
			unlock();
			if(rc){
				AlshaKmerUtils::reverseProps(toKmer.second);
			}
		}else{
			kmerIndex = -1;
			getRemoteKmerProps(toKmer.second, kmerKey, rc, procId);
		}
		//AlshaKmerUtils::printKmer(toKmer.first);
		//AlshaKmerUtils::printKmer(kmerKey);
	
		//check the validity of this kmer
		if(AlshaKmerUtils::isDead(toKmer.second)){
			if(AlshaParams::numProcs == 1){
				cout << "reaching a dead k-mer, and it is impossible" << endl;
				AlshaUtils::exitProgram();
			}
			return true;
		}
		//loop detecting
		pair<AlshaKmerSet::iterator, bool> ret;
		ret = kmerSet.insert(toKmer.first);
		if(ret.second == false){
			//if a small loop is detected, remove this small repeats
			return true;
		}
		//push this k-mer to the k-mer path
		kmerPath.push_back(toKmer);
		kmerIndexPath.push_back(kmerIndex);
		kmerKeyPath.push_back(kmerKey);

		//increase the total length of this tip
		++totalLength;
		//cout<<"totalLength: " << totalLength << endl;
		if(totalLength >= tipCutoff){
			return false;
		}
		
		//exending the tip
		kmer = &kmerPath.back();
	}
	/*****************************
		--A--B--C
	E--F
		--M--N
	For this case, do not remove this tip
	******************************/
	if(AlshaKmerUtils::getLinkageCount(kmer->second, reverseDir) < 2){
		return false;
	}
	//if it is a singleton
	if(AlshaKmerUtils::getMultiplicity(kmer->second) == 1){
		return true;
	}

	//excluding the last one
	kmerPath.pop_back();	
	kmerIndexPath.pop_back();
	kmerKeyPath.pop_back();

	if(kmerPath.size() == 0){
		return true;
	}
	unsigned int kmerCounts = 0;
	for(AlshaQueue<AlshaKmer>::iterator iter = kmerPath.begin(); iter != kmerPath.end(); ++iter){
		kmerCounts += AlshaKmerUtils::getMultiplicity(iter->second);
	}
	double coverage = kmerCounts / kmerPath.size();
	if(coverage >= AlshaParams::COVERAGE){
		return false;
	}
	return true;
}
/****************************************************
	mark branching
****************************************************/
void Alsha::splitBranching()
{

	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"splitting linear paths" << endl;
	}

	//mark branching
	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
		if(AlshaKmerUtils::isDead(vecKmerProps[kmerIndex])){
			continue;
		}

		//check each direction
		for(uint8_t dir = ALSHA_FORWARD_DIR; dir <= ALSHA_REVERSE_DIR; ++dir){
			if(AlshaKmerUtils::getLinkageCount(vecKmerProps[kmerIndex], dir) > 1){
				AlshaKmerUtils::mark(vecKmerProps[kmerIndex], dir);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(AlshaParams::numProcs > 1){
		//send message to split the branches
		taskThread->sendRequest(new AlshaMessageControl(ALSHA_CONTROL_SPLIT_BRANCHING_START));

		int numProcs = AlshaParams::numProcs;
		while(numProcs > 0){
			splitBranchingLoop(numProcs);
		}
		//send the complete command to all the other processes
		if(AlshaParams::procId == ALSHA_MASTER_RANK){
			AlshaComms::bcastCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_SPLIT_BRANCHING_DONE);
		}
		//synchronize all processes
		MPI_Barrier(MPI_COMM_WORLD);
	
		//check the completion of all message in flight
	   	while(!AlshaComms::checkCompletion()){
   			splitBranchingLoop(numProcs);
   		}
	}else{
		//directly invoke the core function
		splitBranchingCore();
	}

	//synchronize all processes
	MPI_Barrier(MPI_COMM_WORLD);
}
void Alsha::splitBranchingLoop(int& numProcs)
{
	//waiting for messages from procId
	int flag;
	unsigned int msgSize;
	uint8_t* msgBuffer;
	MPI_Status status;
	
	//check the message from other processes
	if(AlshaComms::checkMessage() == false){
		//check the message queue from the auxiliary thread
		flag = taskThread->getResponseNum();
		if(flag > 0){
			AlshaMessage* response = taskThread->recvResponse();
			uint8_t type = response->getType();
			if(type == ALSHA_MSG_CONTROL_COMMAND){
				AlshaMessageControl* control = dynamic_cast<AlshaMessageControl*>(response);
				switch(control->getCommand()){
				case ALSHA_CONTROL_SPLIT_BRANCHING_DONE:
					//send message to the master process
					if(AlshaParams::procId != ALSHA_MASTER_RANK){
						AlshaComms::sendCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_SPLIT_BRANCHING_DONE);
					}else{
						--numProcs;
					}
					break;
				default:
					fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
					cout<<"unknown message: " << type << endl;
					AlshaUtils::exitProgram();
				}
			}else if(type == ALSHA_MSG_KMER_LINKAGE_REMOVAL){
				AlshaMessageKmerLinkageRemoval* linkage = dynamic_cast<AlshaMessageKmerLinkageRemoval*>(response);
			
				if(AlshaParams::procId != linkage->getProcId()){
					AlshaComms::sendMessage(linkage->getProcId(), linkage);
				}else{
					kmerLinkageRemovalMsgHandler(*linkage);
				}
			}else{
				fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
				cout<<"unknown message: " << type << endl;
				AlshaUtils::exitProgram();
			}	
	
			delete response;
		}
		return;
	}

	//receive the message
	AlshaComms::recvMessage(msgBuffer, msgSize, status);

	//get the message type
	int msgType = AlshaMessage::readMessageType(msgBuffer);
	if(msgType == ALSHA_MSG_KMER_LINKAGE_REMOVAL){
		AlshaMessageKmerLinkageRemoval linkage;
				
		//unserialize the message
		linkage.unserialize(msgBuffer);
		
		//handle this message
		kmerLinkageRemovalMsgHandler(linkage);
	}else if(msgType == ALSHA_MSG_CONTROL_COMMAND){
		AlshaMessageControl control;
				
		//unserialize the message
		control.unserialize(msgBuffer);
				
		switch(control.getCommand()){
		case ALSHA_CONTROL_SPLIT_BRANCHING_DONE:
			if(AlshaParams::procId == ALSHA_MASTER_RANK){
				--numProcs;
			}else{
				numProcs = 0;
			}
			break;
		default:
			fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
			cout<< "Unmatched command (" << AlshaParams::procId << "): " << control.getCommand() << endl;
			AlshaUtils::exitProgram();
			break;
		}
	}else{
		fprintf(stderr ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
		cout<< "Unknown message type: " << msgType << endl;
		AlshaUtils::exitProgram();
	}

	//release msgBuffer
	AlshaUtils::memFree(msgBuffer);
}
void Alsha::splitBranchingCore()
{
	int procId;
	//remove linkages
    const unsigned int maxElements = AlshaParams::getLinkageBatchSize();
    vector<AlshaMessageKmerLinkageRemoval*, AlshaAllocator<AlshaMessageKmerLinkageRemoval*> > msg2Proc;
    msg2Proc.resize(AlshaParams::numProcs);
    for(procId = 0; procId < msg2Proc.size(); ++procId){
        msg2Proc[procId] = new AlshaMessageKmerLinkageRemoval(procId);
    }

	int iterations = 0;
	AlshaKmer kmer;
	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
		
		lock();
		kmer = make_pair(vecKmerData[kmerIndex], vecKmerProps[kmerIndex]);
		unlock();

		if(AlshaKmerUtils::isDead(kmer.second)){
			continue;
		}

		//check each direction
		for(uint8_t dir = ALSHA_FORWARD_DIR; dir <= ALSHA_REVERSE_DIR; ++dir){
			if(AlshaKmerUtils::isMarked(kmer.second, dir)){
				removeKmerLinkage(msg2Proc, kmer, kmerIndex, dir);
				iterations++;
			}
		}
 		//check if this message can be sent out
		if(iterations >= maxElements){
           	for(procId = 0; procId < msg2Proc.size(); ++procId){
               	if(msg2Proc[procId]->ready(maxElements) == true){
                   	taskThread->sendResponse(msg2Proc[procId]);
                   	msg2Proc[procId] = new AlshaMessageKmerLinkageRemoval(procId);
            	}
			}
			iterations = 0;
		}
	}
 	//check if this message can be sent out
	for(procId = 0; procId < msg2Proc.size(); ++procId){
    	if(msg2Proc[procId]->ready() == true){
           	taskThread->sendResponse(msg2Proc[procId]);
          	msg2Proc[procId] = NULL;
      	}
 	}
	msg2Proc.clear();
}
/****************************************************
 * creating a pre-graph
 ****************************************************/
void Alsha::createPreGraph()
{
	int rounds;
	NumType totalPreGraphNodes;
	NumType totalErodedKmers;
	NumType lastTotalLinearPaths;
	string preGraphFileName = AlshaParams::getPreGraphFileName();

	//split branching
	splitBranching();
	
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout <<"Start the generation of linear paths" << endl;
	}

	//generate the file name
	char* fileName = new char [preGraphFileName.length() + 32];
	sprintf(fileName, "%s-%d", preGraphFileName.c_str(), AlshaParams::procId);

	//open the file for pregraph
	mpiFile = fopen(fileName, "wb");
	if(!mpiFile){
		cout << "Failed to open file " << fileName << endl;
		AlshaUtils::exitProgram();
	}

	//
	preGraphNodeIndex = 0;

	//parallel generation of linear paths
	for(rounds = 0; ; ++rounds){
		
		//zero the number of removed kmers
		numErodedKmers = 0;

		//send start create pregraph command to the task thread
		taskThread->sendRequest(new AlshaMessageControl(ALSHA_CONTROL_CREATE_PREGRAPH_START));
		//
		int numProcs = AlshaParams::numProcs;
		while(numProcs > 0){
			createPreGraphLoop(numProcs);
		}
		//send the completion command to all other processes
		if(AlshaParams::procId == ALSHA_MASTER_RANK){
			AlshaComms::bcastCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_CREATE_PREGRAPH_DONE);
		}
		//syncrhonize all the processes
		MPI_Barrier(MPI_COMM_WORLD);
	
  	//check the completion of all message in flight
  	while(!AlshaComms::checkCompletion()){
   		createPreGraphLoop(numProcs);
 		}

		//reduction to get the total number of kmers that have been removed
		MPI_Allreduce(&numErodedKmers, &totalErodedKmers, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		if(totalErodedKmers == 0){
			break;
		}
	}
	//
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&preGraphNodeIndex, &totalPreGraphNodes, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout << "Generate " << totalPreGraphNodes << " linear paths in parallel in " << rounds << " runs" << endl;
		cout << "Start the following sequential generation of linear paths" << endl;
	}

	//sequential generation of linear paths
	for(int procId = 0; procId < AlshaParams::numProcs; ++procId){
		if(procId == AlshaParams::procId){
			//create the pre-graph from the its own k-mer nodes
			sequentialCreatePreGraphCore(mpiFile);

			//send DONE message to all the other processes
			AlshaComms::bcastCommand(AlshaParams::procId, ALSHA_CONTROL_CREATE_PREGRAPH_DONE);
		}else{
			int numProcs = 1;
			while(numProcs > 0){
				createPreGraphLoop(numProcs);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

	  	//check the completion of all message in flight
 		while(!AlshaComms::checkCompletion()){
			int numProcs;
      		createPreGraphLoop(numProcs);
   		}
	}
	//close the file
	fclose(mpiFile);
	delete [] fileName;

	//reduce to get the total number of linear paths
	MPI_Allreduce(&preGraphNodeIndex, &totalPreGraphNodes, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout << "Totally Generate: " << totalPreGraphNodes << " linear paths" << endl;
	}
	//write the total number of graph nodes to the index file
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		string fileName = preGraphFileName;
		fileName.append(".info");
		FILE* file = fopen(fileName.c_str(),"wb");
		if(!file){
			cout << "Failed to create file: " << fileName << endl;
			AlshaUtils::exitProgram();
		}
		int double_strand = 1;
		fprintf(file, "%lld\t%lld\t%d\t%d", totalPreGraphNodes, this->getNumSeqs(), AlshaParams::KMERLENGTH, double_strand);
		fclose(file);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

}
void Alsha::createPreGraphLoop(int& numProcs)
{
	//waiting for messages from procId
	int flag;
	unsigned int msgSize;
	uint8_t* msgBuffer;
	MPI_Status status;
	
	//check the message queue from the auxiliary thread
	flag = taskThread->getResponseNum();
	if(flag > 0){
		AlshaMessage* response = taskThread->recvResponse();
		uint8_t type = response->getType();
		if(type == ALSHA_MSG_CONTROL_COMMAND){
			AlshaMessageControl* control = dynamic_cast<AlshaMessageControl*>(response);
			switch(control->getCommand()){
			case ALSHA_CONTROL_CREATE_PREGRAPH_DONE:
				//send message to the master process
				if(AlshaParams::procId != ALSHA_MASTER_RANK){
					AlshaComms::sendCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_CREATE_PREGRAPH_DONE, preGraphNodeIndex);
				}else{
					--numProcs;
				}
				break;
			}
		}else if(type == ALSHA_MSG_KMER_REMOVAL){
			AlshaMessageKmerRemoval* removal = dynamic_cast<AlshaMessageKmerRemoval*>(response);
			
			assert(AlshaParams::procId != removal->getProcId());
			AlshaComms::sendMessage(removal->getProcId(), removal);

		}else if(type == ALSHA_MSG_KMER_PROPS_REQUEST){
			AlshaKmerProps props;
			AlshaMessageKmerPropsRequest* request = dynamic_cast<AlshaMessageKmerPropsRequest*>(response);
			//get the props
			getRemoteKmerProps(props, request);
		/*}else if(type == ALSHA_MSG_STRING){
			AlshaMessageString* str = dynamic_cast<AlshaMessageString*>(response);
			//write the string to the file
			writeString(mpiFile, str->getString(), str->getStringIndex());*/
		}else{
			fprintf(stdout ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
			cout<<"unknown message: " << type  << endl;
			AlshaUtils::exitProgram();
		}	
		
		delete response;
	}
	//check the existence of message
	if(AlshaComms::checkMessage() == false){
		return;
	}
	//receive the message
	AlshaComms::recvMessage(msgBuffer, msgSize, status);

	//get the message type
	int msgType = AlshaMessage::readMessageType(msgBuffer);
	if(msgType == ALSHA_MSG_KMER_REMOVAL){
		AlshaMessageKmerRemoval removal;
				
		//unserialize the message
		unsigned int offset = removal.unserialize(msgBuffer);
		assert(offset == msgSize);
		//do it by itself
		kmerRemovalMsgHandler(removal);
	
	}else if(msgType == ALSHA_MSG_KMER_PROPS_RESPONSE){
		AlshaMessageKmerPropsResponse* response = new AlshaMessageKmerPropsResponse;
		
		//unserialize the message
		unsigned int offset = response->unserialize(msgBuffer);
		assert(offset == msgSize);

		//send the kmer props to the task thread
		taskThread->sendData(response);

	}else if(msgType == ALSHA_MSG_CONTROL_COMMAND){
		AlshaMessageControl control;
			
		//unserialize the message
		control.unserialize(msgBuffer);
		
		switch(control.getCommand()){
		case ALSHA_CONTROL_CREATE_PREGRAPH_START:
			if(AlshaParams::procId != ALSHA_MASTER_RANK){
				//perform the sequential routine
				sequentialCreatePreGraphCore(mpiFile);
				/*//send the completion message
				AlshaComms::sendCommand(ALSHA_MASTER_RANK, ALSHA_CONTROL_CREATE_PREGRAPH_DONE, preGraphNodeIndex);
				*/
				AlshaComms::bcastCommand(AlshaParams::procId, ALSHA_CONTROL_CREATE_PREGRAPH_DONE, preGraphNodeIndex);
			}
			break;
		case ALSHA_CONTROL_CREATE_PREGRAPH_DONE:
			if(AlshaParams::procId == ALSHA_MASTER_RANK){
				--numProcs;
			}else{
				numProcs = 0;
			}
			break;
		case ALSHA_CONTROL_KMER_PROPS_REQUEST:
			kmerPropsRequestHandler(status.MPI_SOURCE, control.getDescriptor());
			break;
		default:
			fprintf(stdout ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
			cout<< "Unmatched command (" << AlshaParams::procId << "): " << control.getCommand() << endl;
			AlshaUtils::exitProgram();
			break;
		}
	}else{
		fprintf(stdout ,"FUNCTION:%s LINE: %d\n", __FUNCTION__, __LINE__);
		cout<< "Unknown message type: " << msgType << "from process: " << status.MPI_SOURCE << endl;
		AlshaUtils::exitProgram();
	}

	//release msgBuffer
	AlshaUtils::memFree(msgBuffer);
}
void Alsha::parallelCreatePreGraphCore()
{
	AlshaKmer kmer;
  AlshaKmerSet kmerSet;
	AlshaQueue <AlshaKmer, AlshaAllocator<AlshaKmer> > kmerPath;
	AlshaQueue <AlshaKmerData, AlshaAllocator<AlshaKmerData> > kmerKeyPath;
	AlshaQueue <NumType, AlshaAllocator<NumType> > kmerIndexPath;

	//message for kmer and linkage removal
	vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> > msg2Proc;
 	msg2Proc.resize(AlshaParams::numProcs);
	for(int procId = 0; procId < msg2Proc.size(); ++procId){
 		msg2Proc[procId] = new AlshaMessageKmerRemoval(procId);
 	}

	
	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
	
		lock();
		kmer = make_pair(vecKmerData[kmerIndex], vecKmerProps[kmerIndex]);
		unlock();

		//check the validity of this k-mer
		if(AlshaKmerUtils::isDead(kmer.second)){
			continue;
		}

		//if this k-mer is isolated, remove it
		uint8_t toDir;
		int linkage = AlshaKmerUtils::checkLinkage(kmer.second, toDir);
		if(linkage == ALSHA_KMER_CONTINUOUS){
			continue;
		}
		
		//add the current k-mer to the k-mer node path
		kmerSet.insert(kmer.first);
		kmerPath.push_back(kmer);
		kmerKeyPath.push_back((AlshaKmerData)-1); //for local k-mers, no need of this key
		kmerIndexPath.push_back(kmerIndex);
		
		if(linkage == ALSHA_KMER_ISOLATED){
			//output the linear path and remove the k-mers on this path
			outputLinearPath(mpiFile, msg2Proc, kmerPath, true);
			//remove the k-mers
			removeKmer(kmerIndex);
		}else if(parallelFormLinearPath(kmerSet, kmerPath, kmerKeyPath, kmerIndexPath, toDir) == true){
			//output the linear path and remove the k-mers on this path
			outputLinearPath(mpiFile, msg2Proc, kmerPath, true);
			
			//remove the k-mers
			removePath(msg2Proc, kmerPath, kmerKeyPath, kmerIndexPath);
		}
		//clear the kmer path
		kmerSet.erase(kmerSet.begin(), kmerSet.end());
		kmerSet.clear();
		kmerPath.erase(kmerPath.begin(), kmerPath.end());
		kmerPath.clear();
		kmerKeyPath.erase(kmerKeyPath.begin(), kmerKeyPath.end());
		kmerKeyPath.clear();
		kmerIndexPath.erase(kmerIndexPath.begin(), kmerIndexPath.end());
		kmerIndexPath.clear();
	}

	//release the msg2Proc vector
	for(size_t i = 0; i < msg2Proc.size(); ++i){
		if(msg2Proc[i]){
			delete msg2Proc[i];
			msg2Proc[i] = NULL;
		}
	}
	msg2Proc.clear();

}
void Alsha::sequentialCreatePreGraphCore(AlshaFile file)
{
	AlshaKmer kmer;
  AlshaKmerSet kmerSet;
	AlshaQueue <AlshaKmer, AlshaAllocator<AlshaKmer> > kmerPath;
	AlshaQueue <AlshaKmerData, AlshaAllocator<AlshaKmerData> > kmerKeyPath;
	AlshaQueue <NumType, AlshaAllocator<NumType> > kmerIndexPath;

	//message for kmer and linkage removal
    vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> > msg2Proc;
    msg2Proc.resize(AlshaParams::numProcs);
    for(int procId = 0; procId < msg2Proc.size(); ++procId){
        msg2Proc[procId] = new AlshaMessageKmerRemoval(procId);
    }

	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
	
		lock();
		kmer = make_pair(vecKmerData[kmerIndex], vecKmerProps[kmerIndex]);
		unlock();

		//check the validity of this k-mer
		if(AlshaKmerUtils::isDead(kmer.second)){
			continue;
		}

		//if this k-mer is isolated, remove it
		int linkage = AlshaKmerUtils::checkLinkage(kmer.second);
			
		//add the current k-mer to the k-mer node path
		kmerSet.insert(kmer.first);
		kmerPath.push_back(kmer);
		kmerKeyPath.push_back((AlshaKmerData)-1); //for local k-mers, no need of this key
		kmerIndexPath.push_back(kmerIndex);
		
		if(linkage == ALSHA_KMER_ISOLATED){
			//output the linear path and remove the k-mers on this path
			outputLinearPath(file, msg2Proc, kmerPath, false);
			//remove the k-mers
			removeKmer(kmerIndex);
		}else if(sequentialFormLinearPath(kmerSet, kmerPath, kmerKeyPath, kmerIndexPath) == true){
			//output the linear path and remove the k-mers on this path
			outputLinearPath(file, msg2Proc, kmerPath, false);
			
			//remove the k-mers
			removePathDirect(msg2Proc, kmerPath, kmerKeyPath, kmerIndexPath);
		}

		//clear the kmer path
		kmerSet.erase(kmerSet.begin(), kmerSet.end());
		kmerSet.clear();
		kmerPath.erase(kmerPath.begin(), kmerPath.end());
		kmerPath.clear();
		kmerKeyPath.erase(kmerKeyPath.begin(), kmerKeyPath.end());
		kmerKeyPath.clear();
		kmerIndexPath.erase(kmerIndexPath.begin(), kmerIndexPath.end());
		kmerIndexPath.clear();
	}

	//release the msg2Proc vector
	for(size_t i = 0; i < msg2Proc.size(); ++i){
		if(msg2Proc[i]){
			delete msg2Proc[i];
			msg2Proc[i] = NULL;
		}
	}
	msg2Proc.clear();
}
bool Alsha::parallelFormLinearPath(AlshaKmerSet& kmerSet, AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
					AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
					AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath, uint8_t dir)
{
	bool ret = true;
	bool rc;
	int procId;
	NumType kmerIndex;
	AlshaKmer* kmer;

	uint8_t reverseDir = toReverseDir(dir);

	//in the forward direction
	kmer = &kmerPath.front();
	//if this kmer has a single linkage in this direction
	while(AlshaKmerUtils::getLinkageCount(kmer->second, dir) == 1){
		Nucleotide base = AlshaKmerUtils::getFirstLinkage(kmer->second, dir);
		if(base == ALSHA_INVALID_NUCLEOTIDE){
			cout<< "Inconsistent k-mer states" << endl;
			AlshaUtils::exitProgram();
		}
		//get the k-mer linked to from this k-mer
		AlshaKmer toKmer;
		toKmer.first = kmer->first;
		if(dir == ALSHA_FORWARD_DIR){
			AlshaKmerUtils::pushNucleotide(toKmer.first, base);
		}else{
			AlshaKmerUtils::reversePushNucleotide(toKmer.first, base);
		}
		//get the k-mer information
		AlshaKmerData kmerKey = getKmerKey(toKmer.first, rc);
		procId = calcDestinationProc(kmerKey);
		if(procId == AlshaParams::procId){
			kmerIndex = findKmer(kmerKey);
			assert(kmerIndex >= 0);

			lock();
			toKmer.second = vecKmerProps[kmerIndex];
			unlock();
			if(rc){
				AlshaKmerUtils::reverseProps(toKmer.second);
			}
		}else{
			kmerIndex = -1;
			/*AlshaKmerUtils::printKmer(toKmer.first);
			AlshaKmerUtils::printKmer(kmerKey);*/
			getRemoteKmerProps(toKmer.second, kmerKey, rc, procId);
		}
		
		//check the validity of this kmer
		if(AlshaKmerUtils::isDead(toKmer.second)){
			if(AlshaParams::numProcs == 1){
				cout << "reaching a dead k-mer, and it is impossible" << endl;
				AlshaUtils::exitProgram();
			}

			ret = false;
			break;
		}
#if 0
		AlshaKmerUtils::printKmer(kmer->first);
		AlshaKmerUtils::printKmerProps(kmer->second);
		AlshaKmerUtils::printKmer(toKmer.first);
		AlshaKmerUtils::printKmerProps(toKmer.second);
		//check the reverse linkage of the k-mer toKmer
		cout<< (int)AlshaKmerUtils::getLinkageCount(kmer->second, dir) << " "
			<< (int)AlshaKmerUtils::getLinkageCount(toKmer.second, reverseDir) <<endl;
#endif
		if(AlshaKmerUtils::getLinkageCount(toKmer.second, reverseDir) == 1){
			//check whether this k-mer is an end point
			if(AlshaKmerUtils::getLinkageCount(toKmer.second, dir) == 0 
				&& procId > AlshaParams::procId){
				ret = false;
				break;
			}
			//loop detecting
			pair<AlshaKmerSet::iterator, bool> res;
			res = kmerSet.insert(kmerKey);
			if(res.second == false){
				//a loop is detected
				break;
			}
			//push this k-mer to the k-mer path
			if(dir == ALSHA_FORWARD_DIR){
				kmerPath.push_back(toKmer);
				kmerIndexPath.push_back(kmerIndex);
				kmerKeyPath.push_back(kmerKey);
			}else{
				//push this k-mer to the k-mer path
				kmerPath.push_front(toKmer);
				kmerIndexPath.push_front(kmerIndex);
				kmerKeyPath.push_front(kmerKey);
			}
		}else{
			break;
		}
		//do the next loop
		if(dir == ALSHA_FORWARD_DIR){
			kmer = &kmerPath.back();
		}else{
			kmer = &kmerPath.front();
		}
	}

	return ret;
}
bool Alsha::sequentialFormLinearPath(AlshaKmerSet& kmerSet, AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
						AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
						AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath)
{
	bool rc;
	int procId;
	NumType kmerIndex;
	AlshaKmer* kmer;

	bool ret = true;
	uint8_t dir, reverseDir;

	//in the forward direction
	kmer = &kmerPath.front();
	dir = ALSHA_FORWARD_DIR;
	reverseDir = toReverseDir(dir);

	//if this kmer has a single linkage in this direction
	while(AlshaKmerUtils::getLinkageCount(kmer->second, dir) == 1){
		Nucleotide base = AlshaKmerUtils::getFirstLinkage(kmer->second, dir);
		if(base == ALSHA_INVALID_NUCLEOTIDE){
			cout<< "Inconsistent k-mer states" << endl;
			AlshaUtils::exitProgram();
		}
		//get the k-mer linked to from this k-mer
		AlshaKmer toKmer;
		toKmer.first = kmer->first;
		if(dir == ALSHA_FORWARD_DIR){
			AlshaKmerUtils::pushNucleotide(toKmer.first, base);
		}else{
			AlshaKmerUtils::reversePushNucleotide(toKmer.first, base);
		}
		//get the k-mer information
		AlshaKmerData kmerKey = getKmerKey(toKmer.first, rc);
		procId = calcDestinationProc(kmerKey);
		if(procId == AlshaParams::procId){
			kmerIndex = findKmer(kmerKey);
			assert(kmerIndex >= 0);

			lock();
			toKmer.second = vecKmerProps[kmerIndex];
			unlock();
			if(rc){
				AlshaKmerUtils::reverseProps(toKmer.second);
			}
		}else{
			kmerIndex = -1;
			/*AlshaKmerUtils::printKmer(toKmer.first);
			AlshaKmerUtils::printKmer(kmerKey);*/
			getRemoteKmerPropsDirect(toKmer.second, kmerKey, rc, procId);
		}
		
		//check the validity of this kmer
		if(AlshaKmerUtils::isDead(toKmer.second)){
			cout << "reaching a dead k-mer, and it is impossible" << endl;
			AlshaUtils::exitProgram();
			break;
		}
		/*AlshaKmerUtils::printKmerProps(kmer->second);
		AlshaKmerUtils::printKmerProps(toKmer.second);
		//check the reverse linkage of the k-mer toKmer
		cout<< (int)AlshaKmerUtils::getLinkageCount(kmer->second, dir) << " "
			<< (int)AlshaKmerUtils::getLinkageCount(toKmer.second, reverseDir) <<endl;*/
		if(AlshaKmerUtils::getLinkageCount(toKmer.second, reverseDir) == 1){
			//a single linkage, concaten the two k-mers together
			//loop detecting
			pair<AlshaKmerSet::iterator, bool> res;
			res = kmerSet.insert(kmerKey);
			if(res.second == false){
				//a loop is detected
				break;
			}
			//push this k-mer to the k-mer path
			kmerPath.push_back(toKmer);
			kmerIndexPath.push_back(kmerIndex);
			kmerKeyPath.push_back(kmerKey);
		}else{
			break;
		}
		//do the next loop
		kmer = &kmerPath.back();
	}

	//in the reverse direction
	kmer = &kmerPath.front();
	dir = ALSHA_REVERSE_DIR;
	reverseDir = toReverseDir(dir);

	//if this kmer has a single linkage in this direction
	while(AlshaKmerUtils::getLinkageCount(kmer->second, dir) == 1){
		Nucleotide base = AlshaKmerUtils::getFirstLinkage(kmer->second, dir);
		if(base == ALSHA_INVALID_NUCLEOTIDE){
			cout<< "Inconsistent k-mer states" << endl;
			AlshaUtils::exitProgram();
		}

		//get the k-mer linked to from this k-mer
		AlshaKmer toKmer;
		toKmer.first = kmer->first;
		if(dir == ALSHA_FORWARD_DIR){
			AlshaKmerUtils::pushNucleotide(toKmer.first, base);
		}else{
			AlshaKmerUtils::reversePushNucleotide(toKmer.first, base);
		}
		//get the k-mer information
		AlshaKmerData kmerKey = getKmerKey(toKmer.first, rc);
		procId = calcDestinationProc(kmerKey);
		if(procId == AlshaParams::procId){
			kmerIndex = findKmer(kmerKey);
			assert(kmerIndex >= 0);

			lock();
			toKmer.second = vecKmerProps[kmerIndex];
			unlock();
			if(rc){
				AlshaKmerUtils::reverseProps(toKmer.second);
			}
		}else{
			kmerIndex = -1;
			/*AlshaKmerUtils::printKmer(toKmer.first);
			AlshaKmerUtils::printKmer(kmerKey);*/
			getRemoteKmerPropsDirect(toKmer.second, kmerKey, rc, procId);
		}
		
		//check the validity of this kmer
		if(AlshaKmerUtils::isDead(toKmer.second)){
            cout << "reaching a dead k-mer, and it is impossible" << endl;
			AlshaUtils::exitProgram();
			break;
		}

		//check the reverse linkage of the k-mer toKmer
		if(AlshaKmerUtils::getLinkageCount(toKmer.second, reverseDir) == 1){
			//a single linkage or zero, concaten the two k-mers together
			//loop detecting
			pair<AlshaKmerSet::iterator, bool> res;
			res = kmerSet.insert(kmerKey);
			if(res.second == false){
				//a loop is detected
				break;
			}
			//push this k-mer to the k-mer path
			kmerPath.push_front(toKmer);
			kmerIndexPath.push_front(kmerIndex);
			kmerKeyPath.push_front(kmerKey);
		}else{
			break;
		}
		//do the next loop
		kmer = &kmerPath.front();
	}
	return ret;
}

void Alsha::outputLinearPath(AlshaFile file, vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
					AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath, bool parallel)
{
    string contig;
    contig.reserve(AlshaParams::KMERLENGTH + kmerPath.size());;

	unsigned int kmerCounts = 0;
	//generate the string
	AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >::iterator iter = kmerPath.begin();
	
	contig = AlshaKmerUtils::getKmerString(iter->first);
	kmerCounts += AlshaKmerUtils::getMultiplicity(iter->second);
	
    for(++iter; iter != kmerPath.end(); ++iter){
		Nucleotide base = AlshaKmerUtils::getNucleotide(iter->first, 0);
        char nucleotide = (char)AlshaKmerUtils::decoding(base);
        contig.append(1, nucleotide);
		kmerCounts += AlshaKmerUtils::getMultiplicity(iter->second);
    }
   	double coverage = kmerCounts / kmerPath.size();
    if(coverage < AlshaParams::COVERAGE){
        return;
    }

#if 0
	if(parallel){
    	//send this linear path to the master thread
	    taskThread->sendResponse(new AlshaMessageString(contig, preGraphNodeIndex));
	}else{
	    //write the linear path to the file
		writeString(file, contig, preGraphNodeIndex);
	}
#else
	writeString(file, contig, preGraphNodeIndex);
#endif

	preGraphNodeIndex++;
	if(preGraphNodeIndex % 5000 == 0){
		cout << "Node index: " << preGraphNodeIndex << "(process " << AlshaParams::procId << ")" << endl;
	}
}
void Alsha::writeString(AlshaFile file, string& contig, NumType preGraphNodeIndex)
{
	static AlshaAllocator<char> allocator;

	int bufferSize = contig.length() + 256;
	//allocate memory
	char* buffer = (char*)allocator.allocate(bufferSize);
	
 	sprintf(buffer, "NODE\t%d\t%lld\t%lld\n", AlshaParams::procId, preGraphNodeIndex, contig.length() - AlshaParams::KMERLENGTH + 1);
	sprintf(buffer + strlen(buffer), "%s\n", contig.c_str());
	
	int bufferLength = strlen(buffer);
	if(fwrite(buffer, sizeof(char), bufferLength, file) != bufferLength){
		cout <<"file write failed for process " << AlshaParams::procId << endl;
		AlshaUtils::exitProgram();
	}

	//deallocate memory
	allocator.deallocate(buffer, bufferSize);
}

void Alsha::releaseKmerNodes()
{
	if(vecKmerData){
		AlshaUtils::memFree(vecKmerData);
		vecKmerData = NULL;
	}
	if(vecKmerProps){
		AlshaUtils::memFree(vecKmerProps);
		vecKmerProps = NULL;
	}
	if(vecKmerAccelerator){
		AlshaUtils::memFree(vecKmerAccelerator);
		vecKmerAccelerator = NULL;
	}
}
AlshaHistogram Alsha::buildHistogram()
{
	AlshaHistogram h;
	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
		if(AlshaKmerUtils::isDead(vecKmerProps[kmerIndex])){
			continue;
		}
		h.insert(AlshaKmerUtils::getMultiplicity(vecKmerProps[kmerIndex]));
	}
	return h;
}

/** Calculate a k-mer coverage threshold from the given k-mer coverage
 * histogram. */
float Alsha::calculateCoverageThreshold(const AlshaHistogram& h)
{
	float cov = h.firstLocalMinimum();
	for (unsigned iteration = 0; iteration < 100; iteration++) {
		AlshaHistogram trimmed = h.trimLow((unsigned int)roundf(cov));

		uint64_t median = trimmed.median();
		float cov1 = sqrt(median);
		if (cov1 == cov) {
			return cov;
		}
		cov = cov1;
	}
	if (AlshaParams::procId == ALSHA_MASTER_RANK){
		cerr << "warning: coverage threshold did not converge" << endl;
	}
	return 0;
}
/** Set the coverage-related parameters e and c from the given k-mer
 * coverage histogram. */
void Alsha::setCoverage(const AlshaHistogram& h)
{
	float minCov = calculateCoverageThreshold(h);
	if (minCov < 2){
		minCov = 2;
	}
 	if (AlshaParams::ERODEMULTI < 0) {
       	AlshaParams::ERODEMULTI = (int)(floorf(minCov * 0.5));
				if(AlshaParams::ERODEMULTI < 2){
					AlshaParams::ERODEMULTI = 2;
				}
        if (AlshaParams::procId == ALSHA_MASTER_RANK){
            cout << "Set the lowest multiplicity to "<< AlshaParams::ERODEMULTI << endl;
		}
    }

	if (AlshaParams::COVERAGE < 0) {
		AlshaParams::COVERAGE = minCov;
		if (AlshaParams::procId == ALSHA_MASTER_RANK){
			cout << "Set coverage threshold to " << AlshaParams::COVERAGE << endl;
		}
	}
}
void Alsha::resizeKmerVector()
{

	FILE* file;
	char kmerFileName[256];
	sprintf(kmerFileName, "%s/Kmers-%d", AlshaParams::getPathName().c_str(), AlshaParams::procId);

	file = fopen(kmerFileName, "wb");
	if(!file){
		cerr << "Failed to open file " << kmerFileName << endl;
		AlshaUtils::exitProgram();
	}

	lock();
	//find the number of kmers alive
	NumType kmerNum = 0;
	int kmerBlockSize = 0;
	AlshaKmerData kmerBlock[KMER_BLOCK_SIZE];
	AlshaKmerProps kmerPropsBlock[KMER_BLOCK_SIZE];
	for(NumType kmerIndex = 0; kmerIndex < vecKmerNum; ++kmerIndex){
		if(!AlshaKmerUtils::isDead(vecKmerProps[kmerIndex])){
			++kmerNum;
			if(kmerBlockSize == KMER_BLOCK_SIZE){
				//write the kmer to the file
				if(fwrite(kmerBlock, sizeof(AlshaKmerData), kmerBlockSize, file) != kmerBlockSize){
					AlshaUtils::exitProgram("kmer writting failed");
				}
				if(fwrite(kmerPropsBlock, sizeof(AlshaKmerProps), kmerBlockSize, file) != kmerBlockSize){
					AlshaUtils::exitProgram("kmer props writting failed");
				}
				kmerBlockSize = 0;
			}
			kmerBlock[kmerBlockSize] = vecKmerData[kmerIndex];
			kmerPropsBlock[kmerBlockSize] = vecKmerProps[kmerIndex];
			kmerBlockSize++;
		}
	}
	if(kmerBlockSize > 0){
		//write the kmer to the file
		if(fwrite(kmerBlock, sizeof(AlshaKmerData), kmerBlockSize, file) != kmerBlockSize){
			AlshaUtils::exitProgram("kmer writting failed");
		}
		if(fwrite(kmerPropsBlock, sizeof(AlshaKmerProps), kmerBlockSize, file) != kmerBlockSize){
			AlshaUtils::exitProgram("kmer props writting failed");
		}
	}

	//close the file
	fclose(file);

	//release memory
	AlshaUtils::memFree(vecKmerData);
 	AlshaUtils::memFree(vecKmerProps);
	
	//if it is an emplty array
	if(kmerNum == 0){
		vecKmerNum = 0;	
		vecKmerData = NULL;
		vecKmerProps = NULL;
		
		unlock();
		return;
	}

	//reallocate memory
	vecKmerNum = kmerNum;
	vecKmerData = (AlshaKmerData*)AlshaUtils::memAlloc((vecKmerNum + 1) * sizeof(AlshaKmerData));
	vecKmerProps = (AlshaKmerProps*)AlshaUtils::memAlloc((vecKmerNum + 1) * sizeof(AlshaKmerProps));

	//open file kmer file
	file = fopen(kmerFileName, "rb");
	if(!file){
		unlock();
		cerr << "Failed to open file for read " << kmerFileName << endl;
		AlshaUtils::exitProgram();
	}
	rewind(file);

	//read the file
	NumType index = 0;
	NumType numKmerBlocks = vecKmerNum / KMER_BLOCK_SIZE;
	for(int iter = 0; iter < numKmerBlocks; iter++){
		kmerBlockSize = fread(kmerBlock, sizeof(AlshaKmerData), KMER_BLOCK_SIZE, file);

		NumType current = index;
		for(int i = 0; i < kmerBlockSize; i++){
			vecKmerData[current] = kmerBlock[i];
			current++;
		}
		if(fread(kmerPropsBlock, sizeof(AlshaKmerProps), KMER_BLOCK_SIZE, file) != kmerBlockSize){
			unlock();
			AlshaUtils::exitProgram("resize() fread failed");
		}
		current = index;
		for(int i = 0; i < kmerBlockSize; i++){
			vecKmerProps[current] = kmerPropsBlock[i];
			current++;
		}
		//increase the index
		index += kmerBlockSize;
	}
	kmerBlockSize = vecKmerNum % KMER_BLOCK_SIZE;
	if(kmerBlockSize > 0){
		if(fread(kmerBlock, sizeof(AlshaKmerData), kmerBlockSize, file) != kmerBlockSize){
			unlock();
			AlshaUtils::exitProgram("resize() fread failed 1");
		}
		
		NumType current = index;
    for(int i = 0; i < kmerBlockSize; i++){
      vecKmerData[current] = kmerBlock[i];
      current++;
    }
		if(fread(kmerPropsBlock, sizeof(AlshaKmerProps), kmerBlockSize, file) != kmerBlockSize){
      unlock();
      AlshaUtils::exitProgram("resize() fread failed");
    }
    current = index;
    for(int i = 0; i < kmerBlockSize; i++){
      vecKmerProps[current] = kmerPropsBlock[i];
      current++;
    }
		//increase the index
		index += kmerBlockSize;
	}
	//close the file
	fclose(file);

	if(index != vecKmerNum){
		cerr << "Error occurrend for file " << kmerFileName << "(" << index << "/" << vecKmerNum << ")" << endl;
		AlshaUtils::exitProgram();
	}
	
	//build kmer vector accelerator
	buildKmerVectorAccelerator();

	unlock();
}

void Alsha::buildKmerVectorAccelerator()
{
	NumType* ptr;
	NumType current, last;

	if(vecKmerAccelerator == NULL){
		vecKmerAcceleratorBits = 24;
		vecKmerAcceleratorShift = 2 * AlshaParams::KMERLENGTH - vecKmerAcceleratorBits;
		vecKmerAccelerator = (NumType*)AlshaUtils::memAlloc(((1 << vecKmerAcceleratorBits) + 1) * sizeof(NumType));
	}
	ptr = vecKmerAccelerator;
	*ptr = 0;
	last = 0;
	for(NumType index = 0; index < vecKmerNum; index++){
		current = vecKmerData[index] >> vecKmerAcceleratorShift;
		while(last < current){
			last++;
			ptr++;
			*ptr = index;
		}
	}
	while(last < (NumType)(1 << vecKmerAcceleratorBits)){
		last++;
		ptr++;
		*ptr = vecKmerNum;
	}
}
