#!/usr/bin/env python3

SCREENMAP = "GRINCH.screen"
NPERCOLUMN = 43
NVETROCS = 4

def getindex(slot,chan):
    '''Convert vetroc module # and chan into a fake
    slot and channel assuming fake 32 channel modules and
    accounting for the 3 extra words and the end of each of the 128
    channels per vetroc.'''
    index = slot*(128+3)+chan
    return(index)

scalermap=open(SCREENMAP,"w")

for ivetroc in range(NVETROCS):
    for ichan in range(128):
        index=getindex(ivetroc,ichan)
        #column = (ivetroc*128+ichan)//NPERCOLUMN + 1
        column = ivetroc*3 + ichan//NPERCOLUMN + 1
        #row = (ivetroc*128+ichan)%NPERCOLUMN
        row = ichan%NPERCOLUMN
        print(f"{column} {row} {ivetroc}/{ichan} {index}",file=scalermap)
