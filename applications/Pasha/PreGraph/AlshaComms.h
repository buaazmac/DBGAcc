#ifndef ALSHA_COMMS_H_
#define ALSHA_COMMS_H_

#include "AlshaTypes.h"
#include "AlshaUtils.h"
#include "AlshaMessage.h"
#include <mpi.h>

class AlshaComms
{
public:
	static void init();
	static void destroy();
	static void startup();
	static void sendCommand(int procId, uint32_t command, uint64_t descriptor = 0, int tag = ALSHA_TAG_DATA_TRANSMIT);
	static void bcastCommand(int root, uint32_t command, uint64_t descriptor = 0, int tag = ALSHA_TAG_DATA_TRANSMIT);
	static void sendMessage(int procId, AlshaMessage* msg, int tag = ALSHA_TAG_DATA_TRANSMIT); 

	static bool checkCompletion();
	static bool checkMessage();
	static bool recvMessage(uint8_t*& buffer, unsigned int& bufferSize, MPI_Status& status);
private:
	static MPI_Request msgRequest;
	static const unsigned int RX_BUFFER_SIZE = 16 * 1024 * 1024;
	static uint8_t* RX_BUFFER;
	static void mpiSend(int procId, AlshaMessage* msg, int tag);
	static void mpiBcast(int root, AlshaMessage* msg, int tag);

    static NumType txPackets;
    static NumType rxPackets;
};
#endif

