/*
 * Alsha.h
 *
 *  Created on: 29-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHA_H_
#define ALSHA_H_

#include <mpi.h>
#include <omp.h>
#include <iostream>
using namespace std;

#include "AlshaTypes.h"
#include "AlshaUtils.h"
#include "AlshaParams.h"
#include "AlshaFileParser.h"
#include "AlshaKmer.h"
#include "AlshaMessage.h"
#include "AlshaComms.h"
#include "AlshaThread.h"
#include "AlshaHistogram.h"

#ifdef HAVE_GOOGLE_SPARSE_HASH_SET
#include <google/sparse_hash_set>
#include <google/dense_hash_set>
using google::sparse_hash_set;
using google::dense_hash_set;

#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
#define final(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}

struct sparseKmerHashFunc
{
  size_t operator()(const AlshaKmerData& data) const
	{
		unsigned int a, b, c;
		unsigned int numBytes = (AlshaParams::KMERLENGTH + 3 ) >> 2;
		
  	a = b = c = 0xdeadbeef + numBytes + 131;
		
		a += data & 0xFFFFFFFF;
		b += (data >> 32) & 0xFFFFFFFF;

		final(a, b, c);
		
		return c;
	}
};

typedef sparse_hash_set<AlshaKmerData, sparseKmerHashFunc> AlshaKmerHashSet;
#else
typedef AlshaSet<AlshaKmerData> AlshaKmerHashSet;
#endif
typedef AlshaSet<AlshaKmerData> AlshaKmerSet;

typedef FILE*	AlshaFile;
class Alsha
{
public:
	Alsha(AlshaParams* params);
	~Alsha();

	static int	getProcID() {return AlshaParams::procId;}
	static int	getNumProcs() {return AlshaParams::numProcs;}
	static int	calcDestinationProc(AlshaKmerData kmerData);
	static int	calcDestinationProc(unsigned int hashValue);
	AlshaParams* getParams(){return params;}
	
	/********************************
 	kernel functions
 	********************************/
	NumType	getNumSeqs(){return numSeqs;}
	//create thread
	void createThread();
	void destroyThread();
	//distribute k-mers among all the processes
	void 	distributeKmers();
	//build sorted kmer vector
	void buildSortedKmerVector();
	//build linkages between k-mer nodes
	void 	buildLinkages();
	//estimate coverage
	void	estimateCoverage();
	//remove low-multiplicity kmer
	void 	erodeKmers();
	//clip short and error-prone tips
	void 	clipTips();
	void 	splitBranching();	
	//generate a pre-graph
	void 	createPreGraph();
	/****************************
 		Utility functions
 	****************************/
	//release kmer nodes
	void	releaseKmerNodes();
	//release graph nodes
	void	releaseGraphNodes();
	//the kmer linkage message handler;
	void	kmerLinkageMsgHandler(AlshaMessageKmerLinkage& linkage);
	//the kmer removal message handler
	void	kmerRemovalMsgHandler(AlshaMessageKmerRemoval& msg);
	//the kmer linkage removal message handler
	void 	kmerLinkageRemovalMsgHandler(AlshaMessageKmerLinkageRemoval& msg);
protected:
	//varaibles
	AlshaParams* params;
	//sorted vector of k-mers
	AlshaKmerData* vecKmerData;
	AlshaKmerProps* vecKmerProps;
	NumType vecKmerNum;
	NumType* vecKmerAccelerator;
	uint32_t vecKmerAcceleratorBits;
	uint32_t vecKmerAcceleratorShift;

	//auxiliary thread
	AlshaThread* taskThread;
	pthread_mutex_t threadMutex;
	void	lock();
	void	unlock();

	//number of reads
	NumType numSeqs;
	//number of kmers eroded
	NumType numErodedKmers;

	//index of preGraph node
	AlshaFile mpiFile;
	NumType preGraphNodeIndex;

	/***********************************************
 	*	private functions
 	***********************************************/
	NumType	getKmerNum(AlshaKmerSet* kmerSet);	//get the number of k-mers in the collection
	AlshaKmerData getKmerKey(AlshaKmerData kmer, bool& rc);

	//find a k-mer from the k-mer sorted vector
	NumType	findKmer(AlshaKmerData kmerData);
	
	//remove a k-mer and all linkages connected to it
	void	removeKmer(NumType kmerIndex);	//used to remove isolated kmers in the sorted vector of this process
	void 	removePath(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
								AlshaQueue<AlshaKmer,  AlshaAllocator<AlshaKmer> >& kmerPath,
								AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
        				AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath);
	void 	removePathDirect(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
								AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
								AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
								AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath);

	void	removeKmer(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
								AlshaKmer& kmer, AlshaKmerData kmerKey, NumType kmerIndex);
	void 	removeKmerLinkage(vector<AlshaMessageKmerLinkageRemoval*, AlshaAllocator<AlshaMessageKmerLinkageRemoval*> >& msg2Proc,
								AlshaKmer& kmer, NumType kmerIndex, uint8_t dir);
	void	removeLocalKmer(vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
								AlshaKmer& kmer, NumType kmerIndex, uint8_t dir);
	
	//get the remote k-mer property
	void	getRemoteKmerProps(AlshaKmerProps& props, AlshaMessageKmerPropsRequest* request);
	void	getRemoteKmerProps(AlshaKmerProps& props, AlshaKmerData kmerKey, bool rc, int procId);
	void	getRemoteKmerPropsDirect(AlshaKmerProps& props, AlshaKmerData kmerKey, bool rc, int procId);
	void	kmerPropsRequestHandler(int procId, AlshaKmerData kmerData);

	//the core function of k-mer node linkage buliding
	void 	buildLinkagesCore(const char* readsFileName);
	void	buildLinkagesLoop(int& numProcs);
	
	//erode low-multiplicity kmers
	void	erodeKmersLoop(int& numProcs);
	void	erodeKmersCore();
	//the core function of clipping tips
	void	clipTipsCore(int tipCutoff);
	void	clipTipsLoop(int& numProcs);
	bool	isEligibleTip(AlshaKmerSet& kmerSet, AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
					AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
					AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath, int tipCutoff, uint8_t dir);

	//split branching
	void	splitBranchingLoop(int& numProcs);
	void	splitBranchingCore();

	//resize the kmer vector
	void	resizeKmerVector();
	//build kmer vector accelerator
	void	buildKmerVectorAccelerator();

	//the core function of creating pre-graph
	void	createPreGraphLoop(int& numProcs);
	void	writeString(AlshaFile file, string& str, NumType index);
	void	parallelCreatePreGraphCore();
	void	sequentialCreatePreGraphCore(AlshaFile file);
	void	outputLinearPath(AlshaFile file,
					vector<AlshaMessageKmerRemoval*, AlshaAllocator<AlshaMessageKmerRemoval*> >& msg2Proc,
					AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath, bool parallel);

	bool	parallelFormLinearPath(AlshaKmerSet& kmerSet, AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
					AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
					AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath, uint8_t dir);
	bool	sequentialFormLinearPath(AlshaKmerSet& kmerSet, AlshaQueue<AlshaKmer, AlshaAllocator<AlshaKmer> >& kmerPath,
					AlshaQueue<AlshaKmerData, AlshaAllocator<AlshaKmerData> >& kmerKeyPath,
					AlshaQueue<NumType, AlshaAllocator<NumType> >& kmerIndexPath);

	//build coverage histogram
	AlshaHistogram buildHistogram();
	static float calculateCoverageThreshold(const AlshaHistogram& h);
	void setCoverage(const AlshaHistogram& h);

	friend class AlshaThread;
};
#endif /* ALSHA_H_ */
