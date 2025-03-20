/*
 * Scalers --
 *
 *    This file contains a tcl interface to a Hall A style scaler server
 *    as implemented for GeN.
 *
 */
#include "tcl.h"
/*#include "../../RPC/scaserCalls.h"*/
#include <stdio.h>
#include <stdlib.h>

static int GenscalerObjCmd _ANSI_ARGS_((ClientData dummy,
	    Tcl_Interp *interp, int objc,  Tcl_Obj * CONST   objv[]));

static int GenscalerObjCmd(ClientData dummy, Tcl_Interp *interp, int objc,
		       Tcl_Obj * CONST objv[])
{

  Tcl_Obj *listObjv;
  int result;

  if (objc < 3) {
    Tcl_WrongNumArgs(interp, 1, objv, "{get} host [args]");
    return TCL_ERROR;
  }

# if 0
  if (strcmp(Tcl_GetString(objv[1]), "clear") == 0) {
    int result;
    clnt = scaserOpen(Tcl_GetString(objv[2]));	/* For now open and close each call */
    if(!clnt) {
      /*      interp->result = "failed to open scaler server";*/
      return TCL_ERROR;
    }
    if(scaserClear(clnt)==1) {
      result = 1;
    } else {
      /* Clear was not successful, probably because run in progress */
      result = 0;
    }
    scaserClose(clnt);
    Tcl_SetObjResult(interp, Tcl_NewIntObj(result));
    return TCL_OK;
  } else
#endif
    if ((strcmp(Tcl_GetString(objv[1]), "get") == 0) ||
	     (strcmp(Tcl_GetString(objv[1]), "eor") == 0) ||
	     (strcmp(Tcl_GetString(objv[1]), "eor_overflow") == 0)) {
    int first=0;
    int count=0;
    int actualcount,i;
    int *scaler_list;
    int channels;

    channels = genscaOpen(Tcl_GetString(objv[2]));	/* For now open and close each call */
    if(channels <= 0) {
      /*      interp->result = "failed to read from scaler server";*/
      return TCL_ERROR;
    }
    if(objc > 3) {
      first = atol(Tcl_GetString(objv[3]));
      if(objc > 4) {
	count = atol(Tcl_GetString(objv[4]));
      }
    }
    if(count <= 0) {
      count = channels - first;
    }
    scaler_list = (int *) malloc(count*sizeof(int));
    actualcount = genscaCopy(first, count, scaler_list);
    genscaClose(); /* Close the network connection */
    listObjv = Tcl_NewListObj(0, (Tcl_Obj **) NULL);
    /* Copy scaler array into TCL list */
    for(i=0; i<actualcount; i++) {
      Tcl_ListObjAppendElement(interp, listObjv,
			       Tcl_NewIntObj(scaler_list[i]));
    }
    free(scaler_list);
    Tcl_SetObjResult(interp, listObjv);
    return TCL_OK;
  }
  TCL_ERROR;
}

int Genscaler_Init(Tcl_Interp *interp)
{
  int code;

  if (Tcl_InitStubs(interp, TCL_VERSION, 1) == NULL) {
    return TCL_ERROR;
  }
  code = Tcl_PkgProvide(interp, "Genscaler", "1.0");
  if (code != TCL_OK) {
    return code;
  }
  Tcl_CreateObjCommand(interp, "scaler", GenscalerObjCmd,
		     (ClientData) 0, (Tcl_CmdDeleteProc *) NULL);
  return TCL_OK;
}
#if 0
  } else if (strcmp(Tcl_GetString(objv[1], "get") == 0) { /* Return a list of variables */

    if(!clnt) {
      interp->result = "failed to open scaler server";
      return TCL_ERROR;
    }
    
#endif

