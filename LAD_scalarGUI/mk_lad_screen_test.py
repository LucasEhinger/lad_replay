#!/usr/bin/env python3

SCREENMAP = "LAD_test.screen"
NPLANE = 5
BARPPLANE = 11

def getindex(plane,bar,pmt,src):#pmt: top or buttom: 0 or 1; src:fADC or discriminator 0 or 1
    '''Convert vetroc module # and chan into a fake
    slot and channel assuming fake 32 channel modules and
    accounting for the 3 extra words and the end of each of the 128
    channels per vetroc.'''
    index = ((plane*(BARPPLANE)+bar)*2+pmt)*2+src
    return(index)

scalermap=open(SCREENMAP,"w")
walls = ['00','01','10','11','20']
for iplane in range(5):
    for ibar in range(BARPPLANE):
        for ipmt in range(2):
            for isrc in range(2):
                index=getindex(iplane,ibar,ipmt,isrc)
                wall = walls[iplane]
                print(f"{iplane} {'W'}{wall} {ibar}{'T' if ipmt == 0 else 'B'}{'F' if isrc == 0 else 'D'} {index}", file=scalermap)
