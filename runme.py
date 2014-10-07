#
# A script to find functional loci in cancer genes, using
# reconstructed ancestors of kinase genes and a library
# of mutations obeserved in many real cancer strains.
#
# Victor Hanson-Smith
# 2014
# victorhansonsmith@gmail.com
#

import os, sys, re

from argParser import *
ap = ArgParser(sys.argv)

msa_path = "../2014-Feb25/msaprobs/cmgc.msaprobs.phylip"

seed_cbiopath = {"Homo.sapiens.CDK1":"cBioPortal/cdkl.txt",
                 "Homo.sapiens.CDK2":"cBioPortal/cdk2.txt",
                 "Homo.sapiens.CDK4":"cBioPortal/cdk4.txt",
                 "Homo.sapiens.CDK6":"cBioPortal/cdk6.txt",
                 "Homo.sapiens.MAPK1":"cBioPortal/mapk1.txt",
                 "Homo.sapiens.MAPK3":"cBioPortal/mapk3.txt",
                 "Homo.sapiens.MAP2K1":"cBioPortal/map2k1.txt",
                 "Homo.sapiens.MAP2k2":"cBioPortal/map2k2.txt"}


branch_dirs = [] # a list of directories, one directory for each branch, containing mutation information
for d in os.listdir("../2014-Feb25"):
    if d.startswith("Anc") and d.__contains__("to"):
        branch_dirs.append(d)
               

def build_site_map(mpath):
    """Returns (ref2seed, seed2ref)"""
    
    """msapath is the filepath to a PHYLIP aligned file."""
    if False == os.path.exists(mpath):
        print "\n. Error: I can't find your MSA file at", msapth
        exit()
    
    seed2ref = {} # key = taxa name, value = hash; key = site, value = site
    ref2seed = {}

    fin = open(mpath, "r")
    line1 = fin.readline()
    nsites = int( line1.split()[1] )
    for l in fin.xreadlines():
        tokens = l.split()
        this_taxa = tokens[0]
        this_seq = tokens[1]
        seed2ref[this_taxa] = {}
        ref2seed[this_taxa] = {}
        jj = -1
        for ii in range(0, nsites):
            c = this_seq[ii]
            if c == "-":
                continue
            else:
                jj += 1
                ref2seed[this_taxa][ii] = jj
                seed2ref[this_taxa][jj] = ii
    return (ref2seed, seed2ref)


def read_cbio_file(cpath):
    if False == os.path.exists(cpath):
        print "\n. Error, I cannot find your cBioPortal data file at", cpath
        exit()
    
    mu_data = {}
        
    fin = open(cpath, "r")
    for l in fin.xreadlines():
        if l.__len__() <= 1:
            continue
        tokens = l.split()
        if tokens.__len__() <= 8:
            continue

        mu_name = tokens[0]

        type_ii = None
        for tt in range(0, tokens.__len__()):
            if tokens[tt] == "Missense" or tokens[tt] == "Nonsense":
                type_ii = tt
        if type_ii == None:
            #print tokens
            continue


        next = 1
        study = ""
        mu_ii = None
        if tokens[type_ii-1] == "3D":
            study = " ".join( tokens[next:type_ii-2] )
            mu_ii = type_ii-2
        else:
            study = " ".join( tokens[next:type_ii-1] )
            mu_ii = type_ii - 1

        mu = tokens[mu_ii]
        fromstate = mu[0]
        tostate = mu[ mu.__len__()-1 ]
        #print fromstate, tostate, l
        site = None
        if mu.__len__() == 3:
            site = int( mu[1] )
        else:
            site = int( mu[1:mu.__len__()-2] )
        type = tokens[type_ii]
        next = type_ii + 1
        copy = tokens[next]
        next += 1
        if tokens[next] != "U" and tokens[next] != "V":
            next += 2
        else:
            next += 1

        #print l
        #print tokens
        #print next, fromstate, tostate, site, type, copy
        #print tokens[next]


        if tokens.__len__() - next == 3:        
            next += 1 # skip mutation assessor

        allele_freq = tokens[next]

        next += 1
        n_in_sample = tokens[next]
        mu_data[mu_name] = (study, fromstate, tostate, site, type, allele_freq, n_in_sample)
    return mu_data


def parse_df(dpath):
    pass # continue here!

(ref2seed, seed2ref) = build_site_map(msa_path)


#
#
#
for seed in seed_cbiopath:
    if seed not in ref2seed:
        print "\n. I can't find the taxa", seed, "in your MSA. I'm skipping it!"
        continue

    if False == os.path.exists( seed_cbiopath[seed] ):
        print "\n. Error, I can't find your cBioPortal data directory at", seed_cbiopath[seed]
        continue
    
    mu_data = read_cbio_file( seed_cbiopath[seed] )
    if mu_data.__len__() <= 0:
        print "\n. I didn't find any data in", seed_cbiopath[seed]
    
    """Open an Excel file and write data as we gather it."""
    fout = open(seed + ".xls", "w")
    for muname in mu_data:
        data = mu_data[muname]
        line = muname + "\t"
        line += data[1].__str__() + "\t" + data[2].__str__() + "\t" + data[3].__str__()
        refsite = seed2ref[seed][ data[3] ]
        #print line, refsite


    fout.close()


    #print seed
    #for x in mu_data:
    #    print x, mu_data[x]

        

exit()
    





exit()
print ref2seed[seed_taxa]
print seed2ref[seed_taxa]

#for branch in branch_dirs():




