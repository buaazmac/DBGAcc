#include <mpi.h>
#include "Alsha.h"
#include <omp.h>

int main(int argc, char* argv[])
{
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

	//estimate coverage
	alsha->estimateCoverage();

	//optional feature 
	if(AlshaParams::isKmerErodible()){
		alsha->erodeKmers();
		timeCurr = MPI_Wtime();
		if(AlshaParams::procId == ALSHA_MASTER_RANK){
			cout<< "Remove low-multiplicity tip k-mers takes " << timeCurr - timeStart << " seconds" << endl;
		}
		timeStart = timeCurr;
	}
	//clip tips
	alsha->clipTips();
	timeCurr = MPI_Wtime();
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout<< "Clip short and low-coverage dead ends takes " << timeCurr - timeStart << " seconds" << endl;
	}
	timeStart = timeCurr;

	//delete intermediate k-mer file
	char kmerFile[1024];
  sprintf(kmerFile, "%s/Kmers-%d", AlshaParams::getPathName().c_str(), AlshaParams::procId);
	unlink(kmerFile);

	//create pre-graph
	alsha->createPreGraph();
	timeCurr = MPI_Wtime();
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout<< "Creating linear paths takes " << timeCurr - timeStart << " seconds" << endl;
	}
	timeStart = timeCurr;

	//release k-mer nodes
	alsha->releaseKmerNodes();

	//release resources
	delete alsha;

	//concatenate the pregraph file
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
 		string preGraphFileName = AlshaParams::getPreGraphFileName();
		char* fileName = new char [preGraphFileName.length() + 32];
		char* command = new char [preGraphFileName.length() * 2 + 256];

		//remove the original file
		unlink(preGraphFileName.c_str());

		//form the new file
		for(int procId = 0; procId < AlshaParams::numProcs; procId++){
			//get the file name;
			sprintf(fileName, "%s-%d", preGraphFileName.c_str(), procId);
			//form the command
			sprintf(command, "cat %s >> %s", fileName, preGraphFileName.c_str());
			//execute the command
			system(command);

			//remove the pre-graph file of each process
			unlink(fileName);
		}
		
		delete [] fileName;
		delete [] command;
	}

	//synchronize all processes
	MPI_Barrier(MPI_COMM_WORLD);
	timeCurr = MPI_Wtime();
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		cout << "Total runtime of building pre-graph is " << timeCurr - timeInit << " seconds" << endl;
	}

	//finalize MPI environment
	MPI_Finalize();

	return 0;
}
