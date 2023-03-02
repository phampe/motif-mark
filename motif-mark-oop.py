#!/usr/bin/env python

'''
This code is meant to go through a fasta file and a motif file and output out exon location that are based off of capitalized letters as well as motif location on each nucleotide sequence
'''
import argparse
import re
import cairo
import bioinfo
import math

def get_args():
    parser=argparse.ArgumentParser(description="This code will parse through a give SAM file and output a SAM file where there are only unique PCR reads")
    parser.add_argument("-f", help= "FASTA file to input" )
    parser.add_argument("-m", help= "Motif file" )
    return parser.parse_args()

args=get_args()

fasta_file = (args.f) ### working
motif_file = (args.m) ### working

#########################################################################################################################
### Creating our classes

class sequence:
    '''
    This class sequence is taking in the class and keeping the objects with 
    '''
    def __init__(self, nucleotide_sequence , header):

        self.nucleotide_sequence = nucleotide_sequence
        self.header = header

    def exon(self):
        exon = re.search("[A-Z]+", self.nucleotide_sequence)
        start_stop_position = exon.span() ### python is zero base so you have to add one to both values
        return(start_stop_position)


class motif:
    '''
    This will take in the motif and adjust the bases as need if there are any that are ambiguous 
    This will also run through and find the motif and save it so we can refer to it later on
    '''
    def __init__(self, motif_sequence):

        self.motif_sequence = motif_sequence
        self.change_motif_bases()

    def change_motif_bases(self):
        ###replace functions needs to be output to a new variable.
        replaced_base_motif = self.motif_sequence.upper()
        replaced_base_motif = replaced_base_motif.replace("A","[Aa]")
        replaced_base_motif = replaced_base_motif.replace("C","[Cc]")
        replaced_base_motif = replaced_base_motif.replace("G","[Gg]")
        replaced_base_motif = replaced_base_motif.replace("T","[Tt]")
        replaced_base_motif = replaced_base_motif.replace("U","[UuTt]")
        replaced_base_motif = replaced_base_motif.replace("W","[ATatUu]")
        replaced_base_motif = replaced_base_motif.replace("S","[CGcg]")
        replaced_base_motif = replaced_base_motif.replace("M","[ACac]")
        replaced_base_motif = replaced_base_motif.replace("K","[GTgtUu]")
        replaced_base_motif = replaced_base_motif.replace("R","[AGag]")
        replaced_base_motif = replaced_base_motif.replace("Y","[CTctUu]")
        replaced_base_motif = replaced_base_motif.replace("B","[CGTcgtUu]")
        replaced_base_motif = replaced_base_motif.replace("D","[AGTagtUu]")
        replaced_base_motif = replaced_base_motif.replace("H","[ACTactUu]")
        replaced_base_motif = replaced_base_motif.replace("V","[ACGacg]")
        replaced_base_motif = replaced_base_motif.replace("N","[ACGTacgtUu]")
        replaced_base_motif = "(?=(" + replaced_base_motif + "))" 

        ### look ahead will pull look ahead to see if the same sequence is there and if it is get it as well without skipping over

        self.replaced_base_motif = replaced_base_motif

        ### you want to run ^^^ when trying to convert the motif into possible bases (replaced_base_motif)

    def motif_finder(self, sequence):

        motif_spot = re.finditer(self.replaced_base_motif , sequence.nucleotide_sequence)
    
        matches_motif_list = []
        for match in motif_spot:

            matches = match.span() 
            matches = (matches[0], matches[1]+len(self.motif_sequence))
            matches_motif_list.append(matches)

        return(matches_motif_list)

################################################################################################################
### putting my fasta file into a dictionary to use later. Only getting the gene name and sequence. 
fasta_gene_name_nucleotide_dict = {}
sequence_object_list = []
with open(fasta_file,"r") as fa:

    for line in fa:

        if line[0]==">":

            line = line.split(" ")
            gene_name = line[0].strip(">")
            nucleotide = ""

        else:

            nucleotide += line.strip("\n")
        insert = {gene_name:nucleotide}
        fasta_gene_name_nucleotide_dict.update(insert)

for gene in fasta_gene_name_nucleotide_dict:
    sequence_object_list.append(sequence(fasta_gene_name_nucleotide_dict[gene],gene))
    ### we are putting our now nucleotide and gene into our object from our sequence class.

################################################################################################################
### motif list
### we want this to be a list of objects

motif_list = []
with open (motif_file,"r") as motif_f:
        for line in motif_f:
            line = line.strip("\n")
            motif_list.append(motif(line))

#################################################################################################################
### color dictionary for genes. There will only be 4 different colors set here
### total color is 255

color_list = [(0.2, 0.23, 0.9), (0.9, 0.1, 0.1), (0.4, 0.9, 0.4), (0.7, 0.59, 0), (180, 150, 0)]
### blue, red, green, mud brown, yellow

#################################################################################################################
### our pycairo
height_of_figure = int(len(sequence_object_list))
### height is the total number of sequences 

nucleotide_length = []
for i in fasta_gene_name_nucleotide_dict:
    longest = len(fasta_gene_name_nucleotide_dict[i])
    nucleotide_length.append(longest)

width_of_figure = max(nucleotide_length)
### finding the max nucleotide length and having it be the width of the image

width,height = ((100)+width_of_figure),((150)*height_of_figure)

surface = cairo.PDFSurface(str(fasta_file) + ".png", width, height)

image = cairo.Context(surface)

start_y = 100
start_x = 50
key_start = 10
test = 0

for nucleotide in sequence_object_list:
    ### this section will draw the length of the nucleotide as well as the location of the exons

    a = nucleotide.exon() ### gives you exon location on the nucleotide
    # print(nucleotide.header) ### gives you the header for checking

    image.set_source_rgb(0, 0, 0)
    image.move_to(start_x, start_y - 25)
    image.show_text(nucleotide.header)

    image.set_line_width(2)
    image.set_source_rgb(0, 0, 0)
    image.move_to(start_x,start_y)
    image.line_to(start_x + len(nucleotide.nucleotide_sequence),start_y)
    image.stroke()

    image.set_line_width(40)
    image.set_source_rgb(0, 0, 0)
    image.move_to(a[0]+start_x,start_y)
    image.line_to(a[1]+start_x,start_y)
    image.stroke()

    for i, motifs in enumerate(motif_list) :
        ### Here we are going through our motifs and setting up for the next section

        motif_check = (motifs.motif_sequence) ## print motifs and it works (goes through each of the motifs)
        # print(motif_check)
        # print(i)
        replaced_bases_check = motifs.replaced_base_motif ### since we are already running the code inside the function we can just call it and it works
        # print(replaced_bases_check)
        motif_location = (motifs.motif_finder(nucleotide)) ### this going to output the motif location
        # print(motif_location)
        
        if test < len(motif_list):
            ### This will draw our key with a set limit so that we don't keep writing it out

            image.set_source_rgb(color_list[i][0], color_list[i][1], color_list[i][2])
            image.move_to(key_start + 50, 20)
            image.show_text(motif_check)

            image.set_line_width(20)
            image.move_to(key_start + 20 , 20)
            image.line_to(key_start + 40, 20)
            image.stroke()
            key_start += 100
            test+=1

        for location in motif_location:
            ### Here we are drawing the motif at their given locaiton on the nucleotide.

            r = color_list[i][0]
            b = color_list[i][1]
            g = color_list[i][2]
            image.set_line_width(30)
            image.set_source_rgb(r, b, g)
            image.move_to(location[0]+start_x,start_y)
            image.line_to(location[1]+start_x,start_y)
            image.stroke()

    start_y += 125
