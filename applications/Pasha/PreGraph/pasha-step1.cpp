#include <mpi.h>
#include "Alsha.h"
#include <omp.h>

int main(int argc, char* argv[])
{

    	//============ KmerGen start

	double timeStart, timeCurr, timeInit;

	//initialize MPI environment
	MPI_Init(&argc, &argv);
	//get the rank of this process
	MPI_Comm_rank(MPI_COMM_WORLD, &AlshaParams::procId);
	//get the total number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &AlshaParams::numProcs);

	//parsing parameters;
	AlshaParams* params = new AlshaParams(argc, argv);

	//generate entity for master and slave processes
	Alsha* alsha = new Alsha(params);
	
	//create thread
	alsha->createThread();

	timeStart = MPI_Wtime();
	timeInit = timeStart;	//save the starting time

	//distribute all k-mers over all processes
	alsha->distributeKmers();

	timeCurr = MPI_Wtime();
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout << "Total runtime of distributing kmers is " << timeCurr - timeInit << " seconds" << endl;
	}
	timeStart = timeCurr;

    	//============ KmerGen stop, PreGraph start


    	//build sorted kmer vector
	alsha->buildSortedKmerVector();
	timeCurr = MPI_Wtime();
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout<< "Build sorted kmer vector takes " << timeCurr - timeStart << " seconds" << endl;
	}
	timeStart = timeCurr;

	//start the asynchronous receiving
	AlshaComms::startup();

	//build the linkages of k-mer nodes
	alsha->buildLinkages();
	timeCurr = MPI_Wtime();
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout<< "Building k-mer linkages takes " << timeCurr - timeStart << " seconds" << endl;
	}
	timeStart = timeCurr;


    	//============ PreGraph end, common part that end the threads and processes start


	//release resources
	delete alsha;

	//synchronize all processes
	MPI_Barrier(MPI_COMM_WORLD);
	timeCurr = MPI_Wtime();
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout << "Total runtime of step1 is " << timeCurr - timeInit << " seconds" << endl;
	}

	//finalize MPI environment
	MPI_Finalize();

	return 0;
}

