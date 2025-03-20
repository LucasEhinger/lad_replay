/* Routines to read Hall A scalers */

#include <sys/ioctl.h>
#include <sys/types.h>
#include "time.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include "vxWscaler.h"


/* Actually open connection and read scalers.  Return number of channels */
/* Number might be hardwired per host */
struct request myRequest;
struct request bobreply;
struct sockaddr_in serverSockAddr;    /*  socket addr of server */
static char  *defMsg = "Bob wants some data";
int sFd; /* Socket file descriptor */

genscaOpen(char *serverhost) {
  int nRead, nRead1;
  int clearall, checkend;
  clearall = 0;
  checkend = 0;

  serverSockAddr.sin_family = PF_INET;
  serverSockAddr.sin_port = htons (SERVER_PORT_NUM);
  /* Need to do DNS lookup, hardwire for now */
  if((sFd = socket(PF_INET, SOCK_STREAM, 0))==ERROR) {
    return(-1);
  }
  serverSockAddr.sin_addr.s_addr = inet_addr ("129.57.192.5");
  if(connect (sFd, (struct sockaddr *) &serverSockAddr, sizeof(serverSockAddr)) == ERROR) {
    fprintf(stderr,"Error opening socket\n");
    close (sFd);
    return (-2); /* Return a negative channel count */
  }
  myRequest.reply = 1;
  myRequest.clearflag = htonl(clearall);
  myRequest.checkend = htonl(checkend);
  strcpy(myRequest.message,defMsg);
  if(write(sFd, (char *)&myRequest, sizeof(myRequest)) == ERROR) {
    printf(stderr,"Error sendign write\n");
    close (sFd);
    return (-3); /* Return a negative channel count */
  }
  nRead = read(sFd, (char *)&bobreply, sizeof(bobreply));
  if(nRead < 0) {
    close(sFd);
    return(-4);
  }
  /* Keep reading if we didn't get the full expected packet */
  /* Note: should this ever happen?  Seems to be this will break if server
     wants to return a different size buffer than we are expecting. */
  while (nRead < sizeof(bobreply)) {
    nRead1 = read (sFd, ((char *) &bobreply)+nRead,
		   sizeof(bobreply)-nRead);int clearall;
    int checkend;

    /*if (debug) printf("reading some more \n");*/
    if(nRead1 < 0) {
      close(sFd);
      return(-5);
    }
    nRead += nRead1;
  }
  return NUMBLOCKS*16;
}  

  


genscaCopy(int first, int count, int *scaler_list)
{
  int i;
  for(i=0;i<count;i++) {
    scaler_list[i] = ntohl(bobreply.ibuf[first+i]);
  }
}

genscaClose()
{
  close (sFd);
}
