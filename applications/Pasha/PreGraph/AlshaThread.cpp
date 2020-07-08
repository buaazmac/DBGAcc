/*
 * AlshaThread.cpp
 *
 *  Created on: 03-May-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#include "AlshaThread.h"
#include "Alsha.h"

#define MEMORY_CONTROL
#define USE_COND_WAIT
unsigned int AlshaThread::LIST_SIZE_THRESHOLD = 100;

AlshaThread::AlshaThread(Alsha* alsha)
{
	this->alsha = alsha;
  this->isResponseOverSize = false;

	pthread_mutex_init(&this->requestMutex, NULL);
	pthread_mutex_init(&this->requestCondMutex, NULL);
	pthread_mutex_lock(&this->requestCondMutex);
	pthread_mutex_init(&this->responseMutex, NULL);
	pthread_mutex_init(&this->responseCondMutex, NULL);	
	pthread_mutex_lock(&this->responseCondMutex);

	pthread_mutex_init(&this->dataMutex, NULL);
	pthread_mutex_lock(&this->dataMutex);
	pthread_mutex_init(&this->dataCondMutex, NULL);
	pthread_mutex_lock(&this->dataCondMutex);
	pthread_mutex_init(&this->overSizeCondMutex, NULL);
	pthread_mutex_lock(&this->overSizeCondMutex);
	pthread_barrier_init(&this->threadBarrier, NULL, 2);

	this->requestList.clear();
	this->responseList.clear();
	this->dataList.clear();

	//create a thread
	pthread_create(&this->threadIdx, NULL, threadFunc, this);
}

AlshaThread::~AlshaThread()
{
	//destroy the thread
	this->sendRequest(new AlshaMessageControl(ALSHA_CONTROL_THREAD_EXIT));
	this->waitForBarrier();

	//destroy the mutexes
	pthread_mutex_destroy(&this->requestMutex);
	pthread_mutex_destroy(&this->requestCondMutex);

	pthread_mutex_destroy(&this->responseMutex);
	pthread_mutex_destroy(&this->responseCondMutex);

    pthread_mutex_destroy(&this->dataMutex);
    pthread_mutex_destroy(&this->dataCondMutex);

	pthread_mutex_destroy(&this->overSizeCondMutex);

	pthread_barrier_destroy(&this->threadBarrier);
}
void* AlshaThread::threadFunc(void* arg)
{
	AlshaThread* thread = static_cast<AlshaThread*>(arg);
	Alsha* alsha = thread->getAlsha();

	//waiting for request from the master thread
	bool done = false;
	while(!done){
		//get a request
		AlshaMessage* request = thread->recvRequest();
		uint8_t type = request->getType();
		if(type == ALSHA_MSG_CONTROL_COMMAND){
			AlshaMessageControl* control = dynamic_cast<AlshaMessageControl*>(request);
			switch(control->getCommand()){
			case ALSHA_CONTROL_READ_FILE:
				thread->readFile();
				break;
			case ALSHA_CONTROL_KMER_ERODE_START:
				thread->erodeKmers();
				break;
			case ALSHA_CONTROL_CLIP_TIPS_START:
				thread->clipTips(control->getDescriptor());
				break;
			case ALSHA_CONTROL_SPLIT_BRANCHING_START:
				thread->splitBranchingCore();
				break;
			case ALSHA_CONTROL_CREATE_PREGRAPH_START:
				thread->createPreGraph();
				break;
			case ALSHA_CONTROL_THREAD_EXIT:
				done = true;
				break;
			case ALSHA_CONTROL_THREAD_SYNC:
				thread->waitForBarrier();
				break;
			default:
				cout <<"Unrecoginized task request in the auxiliary thread" << endl;
				AlshaUtils::exitProgram();
				break;
			}
		}else if(type == ALSHA_MSG_KMER_LINKAGE){

			AlshaMessageKmerLinkage* linkage = dynamic_cast<AlshaMessageKmerLinkage*>(request);

			thread->buildKmerLinkage(linkage);
		}else if(type == ALSHA_MSG_KMER_REMOVAL){

			//fprintf(stderr, "function %s line: %d\n", __FUNCTION__, __LINE__);
			AlshaMessageKmerRemoval* msg = dynamic_cast<AlshaMessageKmerRemoval*>(request);
			thread->kmerRemoval(msg);

			//fprintf(stderr, "function %s line: %d\n", __FUNCTION__, __LINE__);
		}else{
			cout <<"Unexpected message type: " << type << "in the auxiliary thread" << endl;
			AlshaUtils::exitProgram();
		}

		//release the request data
		delete request;
	}

	//wait for the thread barrier
	thread->waitForBarrier();

	return 0;
}
void AlshaThread::readFile()
{
	AlshaParams* params = alsha->getParams();
	char* sequence;
	int readLength;
	unsigned int batchSize = 1000;
	
	//read a batch of reads
	vector< pair<char*, int> >* batch = NULL;

	NumType numSequences = 0;
	AlshaFileParser* parser;
	vector< pair<AlshaPairedFile, int> > inFileNames = params->getInputFiles();
	for(size_t fileIndex = 0; fileIndex < 2 * inFileNames.size(); fileIndex++){
		AlshaPairedFile fileName = inFileNames[fileIndex / 2].first;
		int fileType = inFileNames[fileIndex / 2].second;

		if(fileIndex % 2 == 0){
			if(fileName.first.length() == 0){
				continue;
			}
			parser = AlshaFileParser::getParser(fileName.first.c_str(), fileType);
		}else{
			if(fileName.second.length() == 0){
				continue;
			}
			parser = AlshaFileParser::getParser(fileName.second.c_str(), fileType);
		}
		
		//get sequence
		while((sequence = parser->getNextSeq(&readLength))){
			//allocate a new batch of reads
			if(batch == NULL){
					batch = new vector< pair<char*, int> >;
			}
    	if(readLength < AlshaParams::KMERLENGTH){
    		//cout<< "Two short read (skipped): " << line << endl;
      	continue;
     	}
#if 0
			//DEBUGGING
			if(numSequences >= 3000000000){
				break;
			}
#endif
			numSequences ++;
			//save this read
			char* buffer = (char*)AlshaUtils::memAlloc(readLength);
			memcpy(buffer, sequence, readLength);
     	batch->push_back(make_pair(buffer, readLength));

			//check the batch size
			if(batch->size() == batchSize){
				//send the batch of reads to the master thread
				sendResponse(new AlshaMessageReadBatch(batch));
				//
				batch = NULL;
			}
		}
    delete parser;
	}
	if(batch  && batch->size() > 0){
		//send the batch of reads to the master thread
		sendResponse(new AlshaMessageReadBatch(batch));
		//
		batch = NULL;
	}
	//send the DONE message
	sendResponse(new AlshaMessageControl(ALSHA_CONTROL_READ_FILE_DONE));

}
void AlshaThread::erodeKmers()
{
	alsha->erodeKmersCore();
	
	//send the number of eroded kmers
	sendResponse(new AlshaMessageControl(ALSHA_CONTROL_KMER_ERODE_DONE));
}
void AlshaThread::clipTips(int tipCutoff)
{
	alsha->clipTipsCore(tipCutoff);

	//send the number of eroded kmers
	sendResponse(new AlshaMessageControl(ALSHA_CONTROL_CLIP_TIPS_DONE));
}
void AlshaThread::splitBranchingCore()
{
	alsha->splitBranchingCore();

	//send the number of eroded kmers
	sendResponse(new AlshaMessageControl(ALSHA_CONTROL_SPLIT_BRANCHING_DONE));
}
void AlshaThread::createPreGraph()
{
	alsha->parallelCreatePreGraphCore();

  //send the number of eroded kmers
	sendResponse(new AlshaMessageControl(ALSHA_CONTROL_CREATE_PREGRAPH_DONE));
}
//building linkages
void AlshaThread::buildKmerLinkage(AlshaMessageKmerLinkage* linkage)
{
	alsha->kmerLinkageMsgHandler(*linkage);
}
//kmer removal
void AlshaThread::kmerRemoval(AlshaMessageKmerRemoval* msg)
{
	alsha->kmerRemovalMsgHandler(*msg);
}
//communication functions
void AlshaThread::waitForRequestCond()
{
	pthread_mutex_lock(&this->requestCondMutex);
}
void AlshaThread::triggerRequestCond()
{
	pthread_mutex_unlock(&this->requestCondMutex);
}
void AlshaThread::waitForResponseCond()
{
	pthread_mutex_lock(&this->responseCondMutex);
}
void AlshaThread::triggerResponseCond()
{
	pthread_mutex_unlock(&this->responseCondMutex);
}
void AlshaThread::waitForDataCond()
{
    pthread_mutex_lock(&this->dataCondMutex);
}
void AlshaThread::triggerDataCond()
{
  pthread_mutex_unlock(&this->dataCondMutex);
}
void AlshaThread::waitForOverSizeCond()
{
	pthread_mutex_lock(&this->overSizeCondMutex);
}
void AlshaThread::triggerOverSizeCond()
{
  pthread_mutex_unlock(&this->overSizeCondMutex);
}
void AlshaThread::waitForBarrier()
{
	pthread_barrier_wait(&this->threadBarrier);
}
unsigned int AlshaThread::getRequestNum()
{
	unsigned int size;

	pthread_mutex_lock(&this->requestMutex);

	size = this->requestList.size();

	pthread_mutex_unlock(&this->requestMutex);

	return size;
}
void AlshaThread::sendRequest(AlshaMessage* request)
{
	pthread_mutex_lock(&this->requestMutex);

	if(this->requestList.size() == 0){
		this->requestList.push_front(request);
#ifdef USE_COND_WAIT
		this->triggerRequestCond();
#endif
	}else{
		this->requestList.push_front(request);
	}

	pthread_mutex_unlock(&this->requestMutex);
}
AlshaMessage* AlshaThread::recvRequest()
{
	AlshaMessage* request;

#ifdef USE_COND_WAIT
	if(getRequestNum() == 0){
		waitForRequestCond();
	}
#else
	while(getRequestNum() == 0);
#endif

	pthread_mutex_lock(&this->requestMutex);
	if(this->requestList.size() == 0){
		cout <<"The number of requests cannot be zero" << endl;
		AlshaUtils::exitProgram();
	}

	request = this->requestList.back();
	this->requestList.pop_back();

	pthread_mutex_unlock(&this->requestMutex);

	return request;
}
unsigned int AlshaThread::getResponseNum()
{
	unsigned int size;

	pthread_mutex_lock(&this->responseMutex);

	size = this->responseList.size();

#ifdef MEMORY_CONTROL
#ifdef USE_COND_WAIT
  //inform the auxiliary thread
  //if(this->responseList.size() <= AlshaThread::LIST_SIZE_THRESHOLD){
  if(this->isResponseOverSize){
    if(size <= AlshaThread::LIST_SIZE_THRESHOLD / 2){
      this->isResponseOverSize = false;
			this->triggerOverSizeCond();
    }
  }
#endif
#endif

	pthread_mutex_unlock(&this->responseMutex);

	return size;
}
void AlshaThread::sendResponse(AlshaMessage* response)
{
#ifdef MEMORY_CONTROL
#ifdef USE_COND_WAIT
  if(getResponseNum() >= AlshaThread::LIST_SIZE_THRESHOLD){
    pthread_mutex_lock(&this->responseMutex);
    this->isResponseOverSize = true;    //set the oversize mark
    pthread_mutex_unlock(&this->responseMutex);
    this->waitForOverSizeCond();
  }
#else
	while(getResponseNum() >= AlshaThread::LIST_SIZE_THRESHOLD);
#endif
#endif

	pthread_mutex_lock(&this->responseMutex);

	if(this->responseList.size() == 0){
		this->responseList.push_front(response);
#ifdef USE_COND_WAIT
		this->triggerResponseCond();
#endif
	}else{
		this->responseList.push_front(response);
	}

	pthread_mutex_unlock(&this->responseMutex);
}
AlshaMessage* AlshaThread::recvResponse()
{
	AlshaMessage* response;

#ifdef USE_COND_WAIT
	while(getResponseNum() == 0){
		waitForResponseCond();
	}
#else
	while(getResponseNum() == 0);
#endif

	pthread_mutex_lock(&this->responseMutex);
	if(this->responseList.size() == 0){
    cout <<"The number of responses cannot be zero" << endl;
    AlshaUtils::exitProgram();
	}
	response = this->responseList.back();
	this->responseList.pop_back();
#ifdef MEMORY_CONTROL
  //inform the auxiliary thread
#ifdef USE_COND_WAIT
  //if(this->responseList.size() <= AlshaThread::LIST_SIZE_THRESHOLD){
  if(this->isResponseOverSize){
    if(this->responseList.size() <= AlshaThread::LIST_SIZE_THRESHOLD / 2){
      this->isResponseOverSize = false;
			this->triggerOverSizeCond();
    }
  }
#endif
#endif
	pthread_mutex_unlock(&this->responseMutex);

	return response;
}

void AlshaThread::sendData(AlshaMessage* data)
{
	//save the data to the data list
	this->dataList.push_front(data);

	//unlock the mutex
	pthread_mutex_unlock(&this->dataMutex);
}
AlshaMessage* AlshaThread::recvData()
{
	AlshaMessage* data;

	//wait for message
	pthread_mutex_lock(&this->dataMutex);

	data = this->dataList.back();
	this->dataList.pop_back();

	return data;
}



