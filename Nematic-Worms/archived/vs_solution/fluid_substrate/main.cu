
//#include "NWSim.h"
//#include "NWHalloweenSimulation.h"
#include "NWSimulation.h"

int main(int argc, char *argv[])
{
	NWSimulation mySim(argc, argv);
	mySim.Run();
	return 0;
}