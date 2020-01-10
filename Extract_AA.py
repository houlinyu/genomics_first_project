#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 09:06:29 2020

@author: houlin.yu
"""

import os
os.getcwd()
os.chdir('/Users/houlin.yu/Desktop')


def extractAA(X):

    import re
    keywords1 = ['start', 'protein sequence', ']', '[A-Z]{98}'] 
    pattern1 = re.compile('|'.join(keywords1))

#keywords2 = ['protein sequence', ']', '[A-Z]{98}']
#pattern2 = re.compile('|'.join(keywords2))
     
    import os
    inputFilepath = X
    filename_w_ext = os.path.basename(inputFilepath)
    filename, file_extension = os.path.splitext(filename_w_ext)
    

    gff_file = open(X, 'r')
    faa_file = open(''.join([filename,'.faa']), 'w')
    for line in gff_file:
 #   if pattern2.search(line):
 #       faa_file.rstrip(line)
         if '# ' in line:
             if pattern1.search(line):
                 faa_file.write(line)
    gff_file.close()
    faa_file.close()
             
    faa_file0 = open(''.join([filename,'.faa']), 'r')
    faa_file2 = open(''.join([filename,'2.faa']), 'w')
 
    stripped_lines = [l.lstrip(l[0:2]) for l in faa_file0.readlines()] 
    faa_file2.write("".join(stripped_lines)) 

    faa_file0.close()
    faa_file2.close()



    faa_file2 = open(''.join([filename,'2.faa']), 'r')
    faa_file3 = open(''.join([filename,'3.faa']), 'w')
    
    prefix = ''.join([filename,'_'])
    seprefix = ''.join(['\n>',prefix])


    for line in faa_file2:
        if 'protein ' in line:
            [line] = [line.lstrip(line[0:20])]
        elif ']' in line:
            [line] = [line.replace(']','')]        
        elif 'start' in line:
            [line] = [line.replace('start gene ',prefix)]

        if ''.join([filename,'_']) not in line:
            [line] = [line.rstrip(line[-1])]
                 
        if ''.join([filename,'_']) in line:
            [line] = [line.replace(prefix,seprefix)]
                     
        faa_file3.write(line)
     
    faa_file2.close()
    faa_file3.close()
















