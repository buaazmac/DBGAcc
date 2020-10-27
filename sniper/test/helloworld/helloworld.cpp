// Your First C++ Program

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


#include "../../include/sim_api.h"


int main() {
	
	SimRoiStart();
	std::cout << "Hello World!"<<std::endl;
	SimRoiEnd();
	return 0;
}
