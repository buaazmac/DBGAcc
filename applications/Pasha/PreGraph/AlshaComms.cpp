#include "AlshaComms.h"
#include "AlshaParams.h"
#include <mpi.h>

MPI_Request AlshaComms::msgRequest = MPI_REQUEST_NULL;
uint8_t* AlshaComms::RX_BUFFER = NULL;
NumType AlshaComms::txPackets = 0;
NumType AlshaComms::rxPackets = 0;

void AlshaComms::init()
{
	RX_BUFFER = (uint8_t*)AlshaUtils::memAlloc(RX_BUFFER_SIZE);
}
void AlshaComms::destroy()
{
	AlshaUtils::memFree(RX_BUFFER);
}
void AlshaComms::startup()
{
    assert(msgRequest == MPI_REQUEST_NULL);
    MPI_Irecv(RX_BUFFER, RX_BUFFER_SIZE, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgRequest);
}
void AlshaComms::sendCommand(int procId, uint32_t command, uint64_t descriptor, int tag)
{
	if(procId < 0) return;

	AlshaMessageControl msg(command, descriptor);

	mpiSend(procId, &msg, tag);
}
void AlshaComms::bcastCommand(int root, uint32_t command, uint64_t descriptor, int tag)
{
	if(root < 0 || AlshaParams::numProcs == 1) return;
	
	AlshaMessageControl msg(command, descriptor);

	//broadcast the command
	mpiBcast(root, &msg, tag);
}
void AlshaComms::sendMessage(int procId, AlshaMessage* msg, int tag)
{
	if(procId < 0) return;

	AlshaUtils::memVerify(msg);

	mpiSend(procId, msg, tag);
}

void AlshaComms::mpiSend(int procId, AlshaMessage* msg, int tag)
{
	if(procId < 0) return;

	unsigned int totalSize = msg->getTransmitSize();
	
	uint8_t* buffer = (uint8_t*)AlshaUtils::memAlloc(totalSize);
	
	//serialize the message
	unsigned int offset = msg->serialize(buffer);
	assert(offset == totalSize);

	//send the message to the destination process
	MPI_Send(buffer, totalSize, MPI_BYTE, procId, tag, MPI_COMM_WORLD);
	AlshaUtils::memFree(buffer);

	txPackets++;
}
void AlshaComms::mpiBcast(int root, AlshaMessage* msg, int tag)
{
	if(root < 0 || AlshaParams::numProcs == 1) return;

	unsigned int totalSize = msg->getTransmitSize();
	
	uint8_t* buffer = (uint8_t*)AlshaUtils::memAlloc(totalSize);
	
	//serialize the message
	unsigned int offset = msg->serialize(buffer);
	assert(offset == totalSize);

	//broadcasting the message
	int reqIndex = 0;
	MPI_Request* requests = (MPI_Request*)AlshaUtils::memAlloc(sizeof(MPI_Request) * AlshaParams::numProcs);
	for(int procId = 0; procId < AlshaParams::numProcs; ++procId){
		if(procId == root) continue;
		
		MPI_Send_init(buffer, totalSize, MPI_BYTE, procId, tag, MPI_COMM_WORLD, &requests[reqIndex]);
		++reqIndex;
	}
	MPI_Startall(reqIndex, requests);
	MPI_Waitall(reqIndex, requests, MPI_STATUSES_IGNORE);

	txPackets += reqIndex;

	AlshaUtils::memFree(requests);
	AlshaUtils::memFree(buffer);
}
bool AlshaComms::checkCompletion()
{
    //perform reduction
    NumType diff = txPackets - rxPackets;
    NumType sumDiff;
    MPI_Allreduce(&diff, &sumDiff, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if(sumDiff == 0){
        txPackets = rxPackets = 0;
	}
    //cout << txPackets << " " << rxPackets << " " << sumDiff << endl;

    return (sumDiff == 0);
}

bool AlshaComms::checkMessage()
{
	int flag;
	MPI_Status status;
	
	MPI_Request_get_status(msgRequest, &flag, &status);
	if(!flag){
		MPI_Request_get_status(msgRequest, &flag, &status);
	}
	return (flag == 1) ? true : false;
}
bool AlshaComms::recvMessage(uint8_t*& buffer, unsigned int& bufferSize, MPI_Status& status)
{
	int flag;

	MPI_Test(&msgRequest, &flag, &status);
	assert(flag == 1);
	
	int size;
	MPI_Get_count(&status, MPI_BYTE, &size);

	bufferSize = size;
	buffer = (uint8_t*)AlshaUtils::memAlloc(bufferSize);
	memcpy(buffer, RX_BUFFER, bufferSize);

 	assert(msgRequest == MPI_REQUEST_NULL);
    MPI_Irecv(RX_BUFFER, RX_BUFFER_SIZE, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgRequest);
	
	rxPackets++;

	return true;
}
