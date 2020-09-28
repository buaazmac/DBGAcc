/*
 * AlshaParams.h
 *
 *  Created on: 31-Mar-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHAPARAMS_H_
#define ALSHAPARAMS_H_

#include "AlshaTypes.h"
#include "AlshaUtils.h"

typedef pair<string, string> AlshaPairedFile;
class AlshaParams
{
public:
	AlshaParams(int argc, char* argv[]);
	~AlshaParams();

	/********************************
 	* public functions
 	* ******************************/
	//static public member functions
	static void		getSysTime(double* dtime);
	static void 	printUsage();
	static void printKmerGenUsage();
	static void 	printPreGraphUsage();
	static void 	printGraphUsage();
	
	//get the intermediate reads file name
	static string& getPathName(){return pathName;}
	static string& getReadsFileName(){return readsFileName;}
	static string& getPreGraphFileName(){return preGraphFileName;}
	static string& getContigFileName(){return contigFileName;}
	static string& getScaffoldFileName(){return scaffoldFileName;}
	static string& getScaffoldGraphFileName(){return scaffoldGraphFileName;}
	static string& getLogFileName(){return logFileName;}
	static string& getMappedReadsFileName(){return mappedReadsFileName;}

	//get the input file list
	vector< pair<AlshaPairedFile, int> >& getInputFiles(){return inFileNames;}
	int	getInsertLength(int cat){return insertLengths[cat];}
	//get the intput files for each category
	pair<AlshaPairedFile, int>& getInputFiles(int cat){assert(cat >= 0 && cat <CATEGORIES); return inFileNames[cat];}
	
	//MPI process information
	static int	getProcID() {return procId;}
	static int	getNumProcs() {return numProcs;}
	static int 	getNumThreads(){return numThreads;}

	//the batch size for k-mer linkage building
	static unsigned int getLinkageBatchSize(){return linkageBatchSize;}
	//the cutoff length of tips
	static int getTipCutoff(){return tipCutoff;}
	static double getCoverageCutoff(){return COVERAGE;}
	static void setCoverageCutoff(double cutoff){COVERAGE = cutoff;}
	static double getExpectedCoverage(){return EXPCOVERAGE;}
	static void	setExpectedCoverage(double cov){EXPCOVERAGE = cov;}

	//number of categories
	static int getNumCategories(){return CATEGORIES;}
	static void setNumCategories(int cat){CATEGORIES = cat;}
	//set the k-mer size
	static void	setKmerSize(int kmerSize);
	//is paired
	static bool isPaired(){return paired;}
	//do scaffolding
	static bool toDoScaffolding(){return doScaffolding;}
	//erode low-multiplicity k-mers
	static bool isKmerErodible(){return kmerErodible;}
	/********************************
 	* public variables
 	* ******************************/

	static int			procId;	//the rank of this process
	static int			numProcs;	//the number of processes
	static int 			numThreads;

	//static public varables for k-mer
	static unsigned int	KMERLENGTH;
	static uint64_t		KMERLENGTHSHIFT;
	static uint64_t		KMERLENGHTFILTER;
	static const unsigned int MAX_KMERLENGTH = 32;
	static double 		COVERAGE;
	static double		EXPCOVERAGE;
	static int	ERODEMULTI;
private:
	//set the cutoff length of tips
	static void	setTipCutoff(int cutoff);	
	static void setPaired(bool value){paired = value;}
	static bool getPaired(){return paired;}
	static void setDoScaffolding(bool value){doScaffolding = value;}
	static void setKmerErodible(bool value){kmerErodible = value;}

	//input file lists
	vector< pair<AlshaPairedFile, int> > inFileNames;
	vector<int> insertLengths;

	//static private member variable
	/*static char readsFileName[256];
	static char preGraphFileName[256];*/
	static string readsFileName;
	static string preGraphFileName;
	static string contigFileName;
	static string scaffoldFileName;
	static string scaffoldGraphFileName;
	static string logFileName;
	static string mappedReadsFileName;
	static string pathName;

	//the cutoff length for tips
	static int tipCutoff;
	static bool paired;
	static bool doScaffolding;
	static bool kmerErodible;
	static unsigned int linkageBatchSize;
	static int CATEGORIES;

};

#endif /* ALSHAPARAMS_H_ */
