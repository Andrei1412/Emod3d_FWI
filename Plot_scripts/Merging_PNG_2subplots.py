import os
from glob import glob
from time import time
import numpy as np
import matplotlib.pyplot as plt
from shutil import rmtree
from fpdf import FPDF
import math

init_time=time()

#directory = 'TestDirectory'
#version1 = '20p5p1'
#version2 = '20p5p8'
#version3 = '19p5p8_emp_All_trim'
##Concatenate predictor variable plots
##simtype = ['sim']
##variabletype = ['Mw']
##variabletype = ['Mw','CD','FM','Rrup','Vs30']
#variabletype = ['Mw','CD','rrup','vs30','H250','H1250']
#biastype = ['Biased']
#for sim in simtype:
#for variable in variabletype:
#    for bias in biastype:
#        loc1 = "/".join([os.getcwd(),directory,'Version_v' + version1,'Plots','Predictor_Variables',variable])
#loc1 = "/home/andrei/Documents/NZVM_paper/NZVM_2010_Oct_revision/figures"
#loc1 = "/home/andrei/Documents/NZVM_paper/DONNA_Updated_NZVM_2020/figures"
loc1 = "/home/andrei/Documents/NZVM_paper/NZVM_2020_Nov_revision/figures"
#loc1 = "/home/andrei/Documents/NZVM_paper/NZVM_2020_Nov_revision/figures/Geology"

pdf = FPDF()
pdf.add_page()

#pdf.image("/".join([loc1, 'data_misfit_m09_2021.png']), 15, 5, 80, 80)
#pdf.image("/".join([loc1, 'data_mean_m09.png']), 95, 5, 80, 80)

pdf.image("/".join([loc1, 'data_misfit_sub1_12iter.png']), 15, 5, 80, 80)
pdf.image("/".join([loc1, 'data_misfit_sub2_12iter.png']), 95, 5, 80, 80)

#pdf.image("/".join([loc1, 'data_misfit_sub2.png']), 15, 5, 80, 80)
#pdf.image("/".join([loc1, 'data_mean_m12.png']), 95, 5, 80, 80)

#pdf.image("/".join([loc1, 'GIS_Marlborough.png']), 15, 5, 80, 80)
#pdf.image("/".join([loc1, 'GIS_Marlborough.png']), 15, 5, 75, 75)
#pdf.image("/".join([loc1, 'Raypath_Marlborough.png']), 100, 5, 80, 80)

#pdf.set_font('Arial', 'B', 16)
pdf.set_font('Arial', '', 10)
#pdf.set_font('normal', '', 12)
pdf.cell(1)
pdf.cell(1, -5, '(a)',0)
#pdf.multi_cell(3, -4, '(a)',0)
pdf.cell(80)
#pdf.cell(90)
#pdf.cell(1, -5, '(b)',0)
pdf.cell(1, -5, '(b)',0)
#pdf.output("/".join([os.getcwd(),'PNG_2subplots_m09.pdf']),'F')
pdf.output("/".join([os.getcwd(),'PNG_2subplots_m09_new.pdf']),'F')
#pdf.output("/".join([os.getcwd(),'PNG_2subplots_m15.pdf']),'F')

#final_time = time()
#print("Done in {:10.1f} secs".format(final_time-init_time))