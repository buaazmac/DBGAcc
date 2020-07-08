/*
 * AlshaThread.h
 *
 *  Created on: 03-May-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHATHREAD_H_
#define ALSHATHREAD_H_

#include "AlshaTypes.h"
#include "AlshaUtils.h"
#include "AlshaMessage.h"

class Alsha;
class AlshaThread
{
public:
    AlshaThread(Alsha* alsha);
    ~AlshaThread();

	static void		setListOverSize(unsigned int size){LIST_SIZE_THRESHOLD = size;}
    Alsha*          getAlsha(){return alsha;}
    void			sendRequest(AlshaMessage* request);
    AlshaMessage*	recvRequest();
    unsigned int 	getRequestNum();

    void            sendResponse(AlshaMessage* response);
    AlshaMessage*   recvResponse();
    unsigned int 	getResponseNum();

		//blocking functions
		void			sendData(AlshaMessage* data);
		AlshaMessage*	recvData();

    void	waitForRequestCond();
    void	waitForResponseCond();
		void	waitForOverSizeCond();
		void	waitForDataCond();
	 	void	triggerRequestCond();
    void	triggerResponseCond();
		void	triggerOverSizeCond();
		void	triggerDataCond();

    void	waitForBarrier();
	
	void readFile();
	void	erodeKmers();
	void	clipTips(int tipCutoff);
	void	splitBranchingCore();
	void	createPreGraph();
	void	kmerRemoval(AlshaMessageKmerRemoval* msg);
	void	buildKmerLinkage(AlshaMessageKmerLinkage* linkage);
private:
	static void* 	threadFunc(void*);
	static unsigned int LIST_SIZE_THRESHOLD;

    Alsha* alsha;
	//barrier for synchronization
    pthread_barrier_t   threadBarrier;

	//for task responses passing from worker thread to master thread
	list<AlshaMessage*, AlshaAllocator<AlshaMessage*> > requestList;

	//for task requests passing from master thread to worker thread
	list<AlshaMessage*, AlshaAllocator<AlshaMessage*> > responseList;
	bool isResponseOverSize;

	//for data passing from master threads to worker thread
	list<AlshaMessage*, AlshaAllocator<AlshaMessage*> > dataList;

	//internal mutex
    pthread_mutex_t requestMutex;
    pthread_mutex_t requestCondMutex;

	pthread_mutex_t responseMutex;
    pthread_mutex_t responseCondMutex;
 	
	pthread_mutex_t	dataMutex;
	pthread_mutex_t	dataCondMutex;
    pthread_mutex_t overSizeCondMutex;


	//thread id
	pthread_t threadIdx;
};


#endif

