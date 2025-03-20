/* scaler2.h          R.Michaels
   July 1998 port to Radstone PPC

   Executible statement definitions, to put after declarations. 
   Contains:
      VME offset addresses vmeoff[] 
      Types of scalers scalertype[] = Types of scalers --
           1151   =  LeCroy model 1151 
           7200   =  Struck model 7200
           560    =  CAEN model V560           
      The correspondence between VME addresses and scaler
      type is specified here.

   Sept 2003.  Removed the 2nd redundant helicity scaler.
               Added 2 units of SIS3801, remove a 1151.
               vxWorks 5.4: 0xe0000000 -> 0xfa000000
   Mod. Oct 2003
               Use sysBusToLocalAdrs.

*/

int ix, res;
unsigned long laddr;
static int addr_assign=0;

/* Address starts for Motorola 2400 */
vmeoff[0]=0xab1000;
vmeoff[1]=0xab1000;
vmeoff[2]=0xab2000;
vmeoff[3]=0xab3000;
vmeoff[4]=0xab4000;
vmeoff[5]=0xab5000;
vmeoff[6]=0xab6000;

/* VME header is a traditional header used by decoding, and
   is independent of the cpu offset */
vmeheader[0]=0xb0d00000;
vmeheader[1]=0xb0d10000;
vmeheader[2]=0xb0d20000;
vmeheader[3]=0xb0d30000;
vmeheader[4]=0xb0d40000;
vmeheader[5]=0xb0d50000;
vmeheader[6]=0xb0d60000;
vmeheader[7]=0xb0d70000;
vmeheader[8]=0xb0d80000;
vmeheader[9]=0xb0d90000;

scalertype[0]=380101; /* 1st scaler, plus helicity */
scalertype[1]=380102; /* 1st scaler, minus helicity */
scalertype[2]=3800;
scalertype[3]=3800; 
scalertype[4]=3800;
scalertype[5]=3800;
scalertype[6]=3800;








