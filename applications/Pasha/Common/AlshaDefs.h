/*
 * AlshaDefs.h
 *
 *  Created on: 07-Apr-2010
 *      Author: Yongchao Liu, School of Computer Engineering,
 *				Nanyang Technology University.
 *		Emails: liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 */

#ifndef ALSHADEFS_H_
#define ALSHADEFS_H_

//data transmission message
enum
{
	//for control commands
	ALSHA_MSG_CONTROL_COMMAND = 100,
	ALSHA_MSG_READ_BATCH,
	//insert kmers
	ALSHA_MSG_KMER_INSERTION,
	//build linkage for the k-mers
	ALSHA_MSG_KMER_LINKAGE,
	/*remove a kmer itself*/
	ALSHA_MSG_KMER_REMOVAL,
	/*removing linkages of a kmer*/
	ALSHA_MSG_KMER_LINKAGE_REMOVAL,
	/*the message contains a string*/
	ALSHA_MSG_STRING,
	/*message pair for kmer properties request*/
	ALSHA_MSG_KMER_PROPS_REQUEST,
	ALSHA_MSG_KMER_PROPS_RESPONSE,
};
//control commands carried through the control message
enum
{
	//indicating the startup and completiong of linking kmers
	ALSHA_CONTROL_KMER_LINKAGE_STARTUP = 128,
	ALSHA_CONTROL_KMER_LINKAGE_REBUILD,
	ALSHA_CONTROL_KMER_LINKAGE_DONE,

	//
	ALSHA_CONTROL_READ_FILE,
	ALSHA_CONTROL_READ_FILE_DONE,
	//
	ALSHA_CONTROL_LOAD_KMER_VECTOR,
	
	//for the erosion of tip kmers
	ALSHA_CONTROL_KMER_ERODE_START,
	ALSHA_CONTROL_KMER_ERODE_DONE,

	//splitting branching
	ALSHA_CONTROL_SPLIT_BRANCHING_START,
	ALSHA_CONTROL_SPLIT_BRANCHING_DONE,
	//request for the k-mer property
	ALSHA_CONTROL_KMER_PROPS_REQUEST,

	//for the operations on tip clipping
	ALSHA_CONTROL_CLIP_TIPS_START,
	ALSHA_CONTROL_CLIP_TIPS_DONE,
	
	//for the creation of pregraph
	ALSHA_CONTROL_CREATE_PREGRAPH_START,
	ALSHA_CONTROL_CREATE_PREGRAPH_DONE,

	//echo information
	ALSHA_CONTROL_ECHO,
	//exit the auxiliary thread
	ALSHA_CONTROL_THREAD_EXIT,
	//synchronize the auxiliary thread
	ALSHA_CONTROL_THREAD_SYNC
};

//sending model determining whether sending message immediately or enqueueing the message
enum
{
	ALSHA_SEND_IMMEDIATE,
	ALSHA_SEND_QUEUED
};
//tags types
enum
{
	ALSHA_TAG_NONE,
	ALSHA_TAG_DATA_TRANSMIT,
};

enum
{
	ALSHA_KMER_ISOLATED,
	ALSHA_KMER_FORWARD_TIP,
	ALSHA_KMER_REVERSE_TIP,
	ALSHA_KMER_CONTINUOUS,
};

//input file formats supported
enum
{
    ALSHA_FILE_FORMAT_FASTA,
    ALSHA_FILE_FORMAT_FASTQ,
    ALSHA_FILE_FORMAT_FASTQGZ
};

#define ADENINE     0
#define CYTOSINE    1
#define GUANINE     2
#define THYMINE     3

//#define ALSHA_COLOR	

#define ALSHA_FORWARD_DIR		0		//forward edge (A-->B), where the (k-1) bp suffix of A is identical to the (k-1) bp prefix of B.
#define ALSHA_REVERSE_DIR		1		//reverse edge (B<--A), where the (k-1) bp prefix of A is identical to the (k-1) bp suffixe of B.
#define ALSHA_KMER_DIR_NUM		2		//total number of directions
#define ALSHA_INVALID_DIR		ALSHA_KMER_DIR_NUM

//k-mer status
#define ALSHA_STATUS_DEAD_OFF			0
#define ALSHA_STATUS_MARK_FORWARD_OFF	1
#define ALSHA_STATUS_MARK_REVERSE_OFF	2

#define ALSHA_NUCLEOTIDE_NUM		4
#define ALSHA_INVALID_NUCLEOTIDE	ALSHA_NUCLEOTIDE_NUM

#ifdef ALSHA_COLOR
#define	toReverseNucleotide(nuc)	(nuc)
#else
#define	toReverseNucleotide(nuc)	(3 - (nuc))
#endif

#define toReverseDir(dir)			(1 - (dir))

#define ALSHA_MASTER_RANK		0

#endif /* ALSHADEFS_H_ */
