/*
 * AlshaParams.cpp
 *
 *  Created on: 31-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */
#include <sys/time.h>
#include "AlshaTypes.h"
#include "AlshaUtils.h"
#include "AlshaParams.h"
#include <omp.h>

#define DEFAULT_KMERLENGTH 	21			//kmer length must be odd number between 1 and 31
//static member variables initialization
unsigned int AlshaParams::KMERLENGTH = DEFAULT_KMERLENGTH;
uint64_t AlshaParams::KMERLENGTHSHIFT = (uint64_t)(AlshaParams::KMERLENGTH - 1) * 2;
uint64_t AlshaParams::KMERLENGHTFILTER = ~(((uint64_t)-1) << (AlshaParams::KMERLENGTH * 2));
double AlshaParams::COVERAGE = -1.0f;
double AlshaParams::EXPCOVERAGE = -1.0f;
int AlshaParams::ERODEMULTI = -1;
int AlshaParams::CATEGORIES = 1;
bool AlshaParams::doScaffolding = true;

#define DEFAULT_READS_FILENAME				"pasha-Sequences"
#define DEFAULT_PREGRAPH_FILENAME			"pasha-PreGraph"
#define DEFAULT_CONTIG_FILENAME				"pasha-Contig"
#define DEFAULT_SCAFFOLD_FILENAME			"pasha-Scaffold"
#define DEFAULT_SCAFFOLD_GRAPH_FILENAME		"plsha-ScaffoldGraph"
#define DEFAULT_LOG_FILENAME				"pasha-Log"
#define DEFAULT_MAPPED_READS_FILENAME		"pasha-MappedReads"
string AlshaParams::pathName="";
string AlshaParams::readsFileName = DEFAULT_READS_FILENAME;
string AlshaParams::preGraphFileName = DEFAULT_PREGRAPH_FILENAME;
string AlshaParams::contigFileName = DEFAULT_CONTIG_FILENAME;
string AlshaParams::scaffoldFileName = DEFAULT_SCAFFOLD_FILENAME;
string AlshaParams::scaffoldGraphFileName = DEFAULT_SCAFFOLD_GRAPH_FILENAME;
string AlshaParams::logFileName = DEFAULT_LOG_FILENAME;
string AlshaParams::mappedReadsFileName = DEFAULT_MAPPED_READS_FILENAME;

//static member variables for MPI processes
int AlshaParams::procId = 0;
int AlshaParams::numProcs = 1;
int AlshaParams::numThreads = 1;

//cutoff length of tips
int AlshaParams::tipCutoff = 2 * AlshaParams::KMERLENGTH;	//2kmer by default

//linkage batch size
unsigned int AlshaParams::linkageBatchSize = 256;

//paired
bool AlshaParams::paired = false;
//erodible kmers
bool AlshaParams::kmerErodible = true;

void AlshaParams::printUsage()
{
#if defined(ALSHA_KMERGEN)
	printKmerGenUsage();
#elif defined(ALSHA_PREGRAPH)
	printPreGraphUsage();
#elif defined(ALSHA_GRAPH)
	printGraphUsage();
#endif
}
void AlshaParams::printKmerGenUsage()
{
	cout << "PASHA is parallelized short read assembler for large genomes using de Bruijn graphs" << endl;
	cout << "Usage:" << endl;
	cout << "\tpasha-kmergen directory [options]" << endl;
	cout << endl << "\tdirectory:\t\t working directory name" << endl;
	cout << "Options:" << endl;
	cout << "\t-fasta <string>(input reads file in FASTA format)" << endl;
	cout << "\t-fastq <string>(input reads file in FASTQ format)" << endl;
	cout << "\t-k <integer> (the kmer size (odd number between 1 and 31), default value: " << DEFAULT_KMERLENGTH
		 << " )" << endl;
	cout << "\t-help or -? (print out the usage information)" << endl;
	cout << "\t-version (print out the version)" << endl;
}
void AlshaParams::printPreGraphUsage()
{
	cout << "PASHA is parallelized short read assembler for large genomes using de Bruijn graphs" << endl;
	cout << "Usage:" << endl;
	cout << "\tpasha-pregraph directory [options]" << endl;
	cout << endl << "\tdirectory:\t\t working directory name" << endl;
	cout << "Options:" << endl;
	cout << "\t-fasta <string>(input reads file in FASTA format)" << endl;
	cout << "\t-fastq <string>(input reads file in FASTQ format)" << endl;
	cout << "\t-help or -? (print out the usage information)" << endl;
	cout << "\t-version (print out the version)" << endl;
}
void AlshaParams::printGraphUsage()
{
	cout << "PASHA is parallelized short read assembler for large genomes using de Bruijn graphs" << endl;
	cout << "Usage:" << endl;
	cout << "\tpasha-graph directory [options]" << endl;
	cout << endl << "\tdirectory:\t\t working directory name" << endl;
	cout << "Options:" << endl;
	cout << "\t-fasta <string> (input reads file in FASTA format)" << endl;
	cout << "\t-fastaPaired <string string> (input paired-end reads file in FASTA format)" << endl;
	cout << "\t-fastaPairedFile <string> (input interleaved paired-end reads file in FASTA format)" << endl;
	cout << "\t-fastq <string> (input reads file in FASTQ format)" << endl;
	cout << "\t-fastqPaired <string string> (input paired-end reads file in FASTQ format)" << endl;
	cout << "\t-fastqPairedFile <string> (input interleaved paired-end reads file in FASTQ format)" << endl;
	cout << "\t-numthreads <integer> (the number of threads for OpenMP (default value: " << AlshaParams::numThreads << " )"<< endl;
	cout << "\t-help or -? (print out the usage information)" << endl;
	cout << "\t-version (print out the version)" << endl;
}

AlshaParams::AlshaParams(int argc, char* argv[])
{
	static bool first = true;

	int index = 1;
	if(argc < 3){
		if(argc < 2){
				printUsage();
				exit(1);
		}
		if(!strcmp(argv[index], "-version")){
			cout <<"Version " << ALSHA_VERSION << endl;
			exit(1);
		//help
		}else if(!strcmp(argv[index], "-help") || !strcmp(argv[index], "-?")){
			printUsage();
			exit(1);
		}else{
			printUsage();
			exit(1);
		}
	}
	//generate the absolute path
	if(argv[index][0] == '-'){
		printUsage();
		exit(1);
	}

	//get the absolute path
	pathName = argv[index];
	readsFileName = argv[index] + (string)"/" + readsFileName;
	preGraphFileName = argv[index] + (string)"/" + preGraphFileName;
	contigFileName = argv[index] + (string)"/" + contigFileName;
	scaffoldFileName = argv[index] + (string)"/" + scaffoldFileName;
	scaffoldGraphFileName = argv[index] + (string)"/" + scaffoldGraphFileName;
	logFileName = argv[index] + (string)"/" + logFileName;
	mappedReadsFileName = argv[index] + (string)"/" + mappedReadsFileName;

	//test the directory
	if(AlshaParams::procId == ALSHA_MASTER_RANK){
		DIR * dir;
		dir = opendir(argv[index]);
		if(dir == NULL){
			mkdir(argv[index], 0777);
		}else{

#if	defined(ALSHA_PREGRAPH) || defined(ALSHA_KMERGEN)
   	unlink(preGraphFileName.c_str());
   	unlink(contigFileName.c_str());
   	unlink(scaffoldFileName.c_str());
   	unlink(scaffoldGraphFileName.c_str());
   	unlink(logFileName.c_str());
   	unlink(mappedReadsFileName.c_str());
#endif

		closedir(dir);
		}
	}
#ifdef ALSHA_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	for(++index; index < argc; ){
		if(argv[index][0] == '-'){
			//FASTA file
			if(!strcmp(argv[index], "-fasta")){
				if(first){
					first = false;
					setPaired(false);
				}else if(getPaired()){
					cout << "Cannot specify single-end and paired-end reads together" << endl;
					AlshaUtils::exitProgram();
				}
				++index;
				bool yes = false;
				if(index < argc && argv[index][0] != '-'){
					AlshaPairedFile file = make_pair(string(argv[index]), string(""));
					inFileNames.push_back(pair<AlshaPairedFile, int>(file, ALSHA_FILE_FORMAT_FASTA));
					++index;
					yes = true;
				}
				if(!yes){
					cout <<"please specify a file name" << endl;
					printUsage();
					exit(0);
				}
			//FASTQ file
			}else if(!strcmp(argv[index], "-fastq")){
                if(first){
                    first = false;
                    setPaired(false);
                }else if(getPaired()){
                    cout << "Cannot specify single-end and paired-end reads together" << endl;
                    AlshaUtils::exitProgram();
                }

				++index;
				bool yes = false;
				if(index < argc && argv[index][0] != '-'){
					AlshaPairedFile file = make_pair(string(argv[index]), string(""));
					inFileNames.push_back(pair<AlshaPairedFile, int>(file, ALSHA_FILE_FORMAT_FASTQ));
					++index;

					yes = true;
				}
       	if(!yes){
          	cout <<"please specify a file name" << endl;
                    printUsage();
                    exit(0);
                }

#if defined(ALSHA_GRAPH)
			//FASTA paired file
			}else if(!strcmp(argv[index], "-fastaPaired")){
                if(first){
                    first = false;
                    setPaired(true);
                }else if(!getPaired()){
                    cout << "Cannot specify single-end and paired-end reads together" << endl;
                    AlshaUtils::exitProgram();
                }
				++index;
				bool yes = false;
				if(index + 1 < argc && argv[index][0] != '-' && argv[index + 1][0] != '-'){
					AlshaPairedFile file = make_pair(string(argv[index]), string(argv[index + 1]));
					inFileNames.push_back(pair<AlshaPairedFile, int>(file, ALSHA_FILE_FORMAT_FASTA));
					index += 2;

					yes = true;
				}
				if(!yes){
					cout <<"please specify two file names" << endl;
					printUsage();
					exit(0);
				}
			//FASTQ paired file
			}else if(!strcmp(argv[index], "-fastqPaired")){
                if(first){
                    first = false;
                    setPaired(true);
                }else if(!getPaired()){
                    cout << "Cannot specify single-end and paired-end reads together" << endl;
                    AlshaUtils::exitProgram();
                }
				++index;
				bool yes = false;
				if(index + 1 < argc && argv[index][0] != '-' && argv[index + 1][0] != '-'){
					AlshaPairedFile file = make_pair(string(argv[index]), string(argv[index + 1]));
					inFileNames.push_back(pair<AlshaPairedFile, int>(file, ALSHA_FILE_FORMAT_FASTQ));
					index += 2;

					yes = true;
				}
                if(!yes){
                    cout <<"please specify two file names" << endl;
                    printUsage();
                    exit(0);
                }
			//FASTA paired reads in a single file
			}else if(!strcmp(argv[index], "-fastaPairedFile")){
                if(first){
                    first = false;
                    setPaired(true);
                }else if(!getPaired()){
                    cout << "Cannot specify single-end and paired-end reads together" << endl;
                    AlshaUtils::exitProgram();
                }
				++index;
				bool yes = false;
				if(index < argc && argv[index][0] != '-'){
					AlshaPairedFile file = make_pair(string(argv[index]), "");
					inFileNames.push_back(pair<AlshaPairedFile, int>(file, ALSHA_FILE_FORMAT_FASTA));
					index++;
					yes = true;
				}
				if(!yes){
					cout <<"please specify a file name" << endl;
					printUsage();
					exit(0);
				}
			//FASTQ paired file in a single file
			}else if(!strcmp(argv[index], "-fastqPairedFile")){
                if(first){
                    first = false;
                    setPaired(true);
                }else if(!getPaired()){
                    cout << "Cannot specify single-end and paired-end reads together" << endl;
                    AlshaUtils::exitProgram();
                }
				++index;
				bool yes = false;
				if(index < argc && argv[index][0] != '-'){
					AlshaPairedFile file = make_pair(string(argv[index]), "");
					inFileNames.push_back(pair<AlshaPairedFile, int>(file, ALSHA_FILE_FORMAT_FASTQ));
					index++;

					yes = true;
				}
                if(!yes){
                    cout <<"please specify a file name" << endl;
                    printUsage();
                    exit(0);
                }
			//number of threads
			}else if(!strcmp(argv[index], "-numthreads")){
				++index;
				sscanf(argv[index], "%d", &numThreads);
				if(numThreads < 1){
					numThreads = 1;
				}
				++index;

#elif defined(ALSHA_KMERGEN)
			//kmer size
			}else if(!strcmp(argv[index], "-k")){
				++index;
				if(index >= argc){
					printUsage();
					exit(0);
				}
				int kmerSize;
				sscanf(argv[index], "%d", &kmerSize);
				if(kmerSize < 0){
					cout << "Kmer size cannot be negative integer and will use default value: " << DEFAULT_KMERLENGTH << endl;
					kmerSize = DEFAULT_KMERLENGTH;
				}else if(kmerSize > 31){
					cout << "Kmer size cannot be over 31 and will be clamped to 31" << endl;
					kmerSize = 31;
				}
				if((kmerSize & 1) == 0){
					kmerSize --;
					cout <<"Kmer size must be odd integer, and thus use " << kmerSize << " instead" << endl;
				}
				index++;
				setKmerSize(kmerSize);
#endif
			//version
			}else if(!strcmp(argv[index], "-version")){
				cout <<"Version " << ALSHA_VERSION << endl;
				exit(1);
			//help
			}else if(!strcmp(argv[index], "-help") || !strcmp(argv[index], "-?")){
				printUsage();
				exit(1);
			}else{
				cout <<"Unknown option: " << argv[index] << endl;
				AlshaUtils::exitProgram();
			}
		}else{
			cout <<"Unknown option: " << argv[index] << endl;
			AlshaUtils::exitProgram();
		}
	}

	//record the number of categories
	CATEGORIES = inFileNames.size();
#ifdef ALSHA_GRAPH
	//if not doing scaffolding
	if(!(isPaired() && toDoScaffolding())){
		CATEGORIES = 1;
	}
	insertLengths.resize(CATEGORIES);
	for(size_t i = 0; i < insertLengths.size(); i++){
		insertLengths[i] = 0;
	}
	cout << "Number of libraries: " << CATEGORIES << endl;
#endif
}
AlshaParams::~AlshaParams()
{
}
void AlshaParams::setKmerSize(int kmerSize)
{
	AlshaParams::KMERLENGTH = (unsigned int )kmerSize;
	AlshaParams::KMERLENGTHSHIFT = (uint64_t)(AlshaParams::KMERLENGTH - 1) * 2;
	AlshaParams::KMERLENGHTFILTER = ~(((uint64_t)-1) << (AlshaParams::KMERLENGTH * 2));
}
void AlshaParams::setTipCutoff(int cutoff)
{
	tipCutoff = cutoff;
}
void AlshaParams::getSysTime(double *dtime)
{
  struct timeval tv;

  gettimeofday (&tv, NULL);

  *dtime = (double) tv.tv_sec;
  *dtime = *dtime + (double) (tv.tv_usec) / 1000000.0;
}

