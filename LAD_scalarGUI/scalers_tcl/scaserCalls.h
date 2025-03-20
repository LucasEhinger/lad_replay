#include "scaserRpc.h"

CLIENT *scaserOpen(char *host);
int scaserClose(CLIENT *clnt);
int scaserClear(CLIENT *clnt);
int scaserGetInfo(CLIENT *clnt, int *slots, int *channels, int *rip,
		  int *clock);
int scaserReadScalers(CLIENT *clnt, int first, int count, int *values,
		      int *overflows, int *events);
int scaserReadEorScalers(CLIENT *clnt, int first, int count, int *values,
		      int *overflows, int *events);

