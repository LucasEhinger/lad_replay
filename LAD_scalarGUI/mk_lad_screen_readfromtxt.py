#!/usr/bin/env python3

import os

input_file_name = 'crate_config.txt'
output_file_name = 'LAD_readfromtxt.screen'
wall = 0
layer = 0
bar = 0
if not os.path.exists(input_file_name):
    print(f"Error: The file '{input_file_name}' does not exist in the current directory.")
else:
    with open(input_file_name, 'r') as file:
        lines = file.readlines()
    with open(output_file_name, 'w') as output_file:
        for line in lines:
            # Clean the line by removing leading/trailing whitespace and skipping empty lines
            cleaned_line = line.strip()
            if not cleaned_line:
                continue
            parts = cleaned_line.split()
            if parts[0]=="Wall":
                wall = int(parts[1])
            if parts[0]=="Layer":
                layer = int(parts[1])
            if parts[0].isdigit():
                bar = int(parts[0])
            if len(parts)>6 and parts[2]=="U":
                output_line1 = str(wall*2+layer)+" "+"W"+str(wall)+str(layer)+" "+str(bar)+parts[2]+"F"+" "+str((int(parts[4])-13)*16+8*32+int(parts[5]))
                output_line2 = str(wall*2+layer)+" "+"W"+str(wall)+str(layer)+" "+str(bar)+parts[2]+"D"+" "+str((int(parts[7])-3)*32+int(parts[8])-1)
                print(output_line1, output_line2)
                output_file.write(output_line1+'\n')
                output_file.write(output_line2+'\n')
            if len(parts)>6 and parts[1]=="D":
                output_line1 = str(wall*2+layer)+" "+"W"+str(wall)+str(layer)+" "+str(bar)+parts[1]+"F"+" "+str((int(parts[3])-13)*16+8*32+int(parts[4]))
                output_line2 = str(wall*2+layer)+" "+"W"+str(wall)+str(layer)+" "+str(bar)+parts[1]+"D"+" "+str((int(parts[6])-3)*32+int(parts[7])-1)
                print(output_line1, output_line2)
                output_file.write(output_line1+'\n')
                output_file.write(output_line2+'\n')


    print(f"Output has been written to '{output_file_name}'.")
