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

seed_cbiopath = {"Homo.sapiens.CDK1":"cBioPortal/cdk1.txt",
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
        branch_dirs.append("../2014-Feb25/" + d)
               

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
            site = int( mu[1:mu.__len__()-1] )
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
    if False == os.path.exists(dpath):
        print "\n. Error, I cannot find your Df path at", dpath
        exit()
    
    key_model = "msaprobs.PROTCATLG"

    refsite_df = {}
    refsite_fromstate = {}
    refsite_frompp = {}
    refsite_tostate = {}
    refsite_topp = {}

    fin = open(dpath, "r")
    last_df = None
    last_refsite = None
    keyc = 0
    for l in fin.xreadlines():
        if l.startswith("-->"):
            tokens = l.split()
            last_df = float(tokens[3])
            last_ref = None
            last_refsite = None
        elif l.startswith( key_model ):
            #print key_model, l
            tokens = l.split()
            last_refsite = int( tokens[2] )
            refsite_df[last_refsite] = last_df
            keyc = 2
        elif keyc == 2:
            #print 2, l
            keyc = 1
            tokens = l.split()
            refsite_fromstate[ last_refsite ] = tokens[3]
            refsite_frompp[ last_refsite ] = tokens[4]
        elif keyc == 1:
            #print 1, l
            keyc = 0
            tokens = l.split()
            refsite_tostate[ last_refsite ] = tokens[3]
            refsite_topp[ last_refsite ] = tokens[4]
            #print "174:", refsite_topp[ last_refsite ], refsite_frompp[ last_refsite ]
    return (refsite_df, refsite_fromstate, refsite_frompp, refsite_tostate, refsite_topp)
            

(ref2seed, seed2ref) = build_site_map(msa_path)

#
#
#
branch_data = {}
for dir in branch_dirs:
    x = parse_df( dir + "/Df.details.txt")
    branch_data[dir] = x

branches_sorted = branch_data.keys()
branches_sorted.sort()


#
# test:
#
#print branch_data[ branches_sorted[0] ][0]
#print branch_data[ branches_sorted[0] ][1]
#print branch_data[ branches_sorted[0] ][2]
#print branch_data[ branches_sorted[0] ][3]
#print branch_data[ branches_sorted[0] ][4]
#exit()


#
#
#
for seed in seed_cbiopath:
    if seed not in ref2seed:
        print "\n. I can't find the taxa", seed, "in your MSA. I'm skipping it!"
        continue
    
    mu_data = read_cbio_file( seed_cbiopath[seed] )
    if mu_data.__len__() <= 0:
        print "\n. I didn't find any data in", seed_cbiopath[seed]
    
    """Open an Excel file and write data as we gather it."""
    fout = open(seed + ".xls", "w")
    header = "Mutation\tfrom\tto\tsite in cBioPortal\tsite in msaprobs\t"
    for branch in branches_sorted:
        btok = branch.split("/")
        branch_short = btok[ btok.__len__()-1 ]
        header += branch_short.__str__() + "\t\t"
    fout.write(header + "\n")
    

    for muname in mu_data:
        data = mu_data[muname]
        line = muname + "\t"
        """site, to state, from state"""
        line += data[1].__str__() + "\t" + data[2].__str__() + "\t" + data[3].__str__() + "\t"
        
        refsite = seed2ref[seed][ data[3]-1 ]
        #print "229:", refsite

        line += (refsite+1).__str__() + "\t"

        for branch in branches_sorted:
            if refsite+1 in branch_data[branch][0]:
                this_df = branch_data[branch][0][refsite+1]
                line += "%.3f"%this_df + "\t"
                from_state = branch_data[branch][1][refsite+1]
                from_pp = branch_data[branch][2][refsite+1]
                to_state = branch_data[branch][3][refsite+1]
                to_pp = branch_data[branch][4][refsite+1]
                line += from_state + "(" + from_pp + ") -> " + to_state + "(" + to_pp + ")\t"
                
            else:
                #q = branch_data[branch][0].keys()
                #q.sort()
                #print "236:", refsite, q 
                #exit()
                line += "NA\tNA\t"

            #(refsite_df, refsite_fromstate, refsite_frompp, refsite_tostate, refsite_topp)

        line += "\n"
        fout.write(line)

    fout.close()


    #print seed
    #for x in mu_data:
    #    print x, mu_data[x]

        

exit()
    





exit()
print ref2seed[seed_taxa]
print seed2ref[seed_taxa]

#for branch in branch_dirs():




