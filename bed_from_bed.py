# !/usr/bin/python2
# 1 March 2018
# A. Vershinina together with https://github.com/AlexKnyshov

# Goal: having a known CONST length of regions of interest (say 1kb) and a step (say 30kb), create a set of coordinates (in BED format) using another .bed file as a guide.
# This particular version implies that BED coordinats tell us where the regions should be. Example: if .bed file:
# chr1	4	35000	feature1
# chr1	40000	43000	feature2
# chr1	100000	220000	feature3
# Region of interest is 1kb and step is 30kb, the output:
# [[4, 1004], [31005, 32005], [100000, 101000], [131001, 132001], [162002, 163002], [193003, 194003]]
#######
# NB! This script works only if step > region. 
# I do no not know how it beahaves if intervals between features are < region length.
#######
# Important details about BED coordinates
# * chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.  
# * chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
# b - bed file
# a - array
# l - length
# f - feature
# r - region
# s - step
# c - coordinate

import pybedtools
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="bed file with guide coordinates")
parser.add_argument("-r", help="length of target regions", type=int)
parser.add_argument("-s", help="min step between regions", type=int)
parser.add_argument("-o", help="output filename")
args = parser.parse_args()
print 'Input file:', args.f
print 'Creating regions of length', args.r, 'bp.'
print 'Minimum step between the regions', args.s, 'bp.'
print 'Saving output as', args.o



b = pybedtools.BedTool(args.f)

r = args.r
s = args.s + r + 1 # desired step + region + 1


def extract_region(start, end, array):
    
    """ 
    Extract coordinates of a region knowing two coordinates. Output - an array with coordinates. 
    """
    a_temp = [start, end]
    array.append(a_temp)
    return array
    
def jump_to_next_interval(bedfile, iterator):
    """
     Update while-loop iterator based on a condition. In this case condition is checking where is the nearest usable feature. Output - increment for an iterator
    """ 
    dists = []
    for feature in bedfile:
        dist = feature.start - iterator
        if dist > 0:
            dists.append(dist)
    m = min(dists)
    return m       

def ordered_set(in_list): 
    """
    what are the chromosomes?
    stolen from https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists
    Takes a list of chromosomes, reduces it to a an ordered set where each chromosome is repeated only one. 
    Order is the same as the order of chroms in the input file.
    """
    out_list = []
    added = set()
    for val in in_list:
        if not val in added:
            out_list.append(val)
            added.add(val)
    return out_list

def read_chroms(f):
    """
    reads bed file creates an ordered list of chromosomes.
    Creates inpit list for the main loop below.
    """

    infile=open(f,"r")
    lines=infile.readlines()
    chroms=[]
    for x in lines:
        chroms.append(x.split('\t')[0]) # bed file should be tab delimited, thus \t is used.
    infile.close()
    chrom_list = ordered_set(chroms)
    return chrom_list

def get_regions(bed, chr):
    """
    this is the main function of this program.
    It takes coordinates from a bed file, and creates a set of new coordinates, using input bed file as a guide.
    """
    features = []
    COORDs=[]
    for f in bed:
        if f.chrom == chr:
            features.append(f)
    f_l = len(features)
    start_1st_f=features[0].start
    end_last_f=features[f_l-1].end
    i = start_1st_f 
    while i < end_last_f:
        start_c = i
        end_c = start_c + r # !NB! BED coordinates do not include the last bp (0 _ 100 is 0-99). Stop coordinate in the output will be 1 bp longer than the desired region length. Let's say you want it to be 3 bp long, starting from 4. Output will be 4-7, i.e. 4,5,6, and 7.   
        ok = False
        for f in features:
            #print f.chrom
            if start_c >= f.start and end_c <= f.end:
                ok = True
                break
        if ok:
            COORDs = extract_region(start_c, end_c, COORDs)
            i = i + s
        else:
            m = jump_to_next_interval(features, i)
            i = i + m
    return COORDs

def flatten(l):
    """
    unlist list of lists
    stolen from stackoverflow
    """
    flat_list = []
    for sublist in l:
        for item in sublist:
            flat_list.append(item)
    return flat_list

def get_regions_by_chr(chrom_list, b): 
    """
    run a loop, iterating over each chromosome in a list looking like
    ['chr1','chr2'];
    b is a bed file read by pybedtools
    """

    final_regions = []
    for C in chrom_list:
        r1 = get_regions(b, C)
        n_regions = len(r1)
        key = n_regions*[C] #  chromosome is the key
        r2 = zip(key, r1) # connect chromosome with its corresponding region, create a tuple looking like ('chr1', [12, 1012])
        temp_regions = [[i[0]]+list(i[1]) for i in r2] # convert tuple to a list ['chr1', 12, 1012]
        final_regions.append(temp_regions)
    return final_regions

CHROMs = read_chroms(args.f)

final_regions = get_regions_by_chr(CHROMs, b)
flat_regions = flatten(final_regions)

with open(args.o, 'w') as outfile:
    np.savetxt(outfile, flat_regions, fmt="%s", delimiter="\t")

print len(b), 'coordinates processed from the guide bed file.'

l=[]
for f in flat_regions:
    l.append(f[0])

unique = set(l)
print 'Count of neutral regions created for each chromosome:'
for item in unique:
    print item, l.count(item)
print 'Total', len(flat_regions), 'putative neutral regions found.'