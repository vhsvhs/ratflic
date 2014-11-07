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
    
    """mpath is the filepath to a PHYLIP aligned file."""
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
    
    mu_data = {} # key = site
        
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


def collapse_mu_data(mu_data):
    site_mu_data = {} # key = site, value = mu_data hashes corresponding to that site
    for muname in mu_data:
        site = mu_data[muname][3]
        if site not in site_mu_data:
            site_mu_data[site] = {}
        site_mu_data[site][muname] = mu_data[muname]
    return site_mu_data

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


def get_branch_data(branchname):
    """x = (refsite_df, refsite_fromstate, refsite_frompp, refsite_tostate, refsite_topp)"""
    x = parse_df( branchname + "/Df.details.txt")
    return x

def get_branch_musites(branchname, ppthresh = 0.6):
    """x = (refsite_df, refsite_fromstate, refsite_frompp, refsite_tostate, refsite_topp)"""
    x = parse_df( branchname + "/Df.details.txt")
    branch_musites[branchname] = []
    for site in x[0]:
        if x[1][site] != x[3][site]:
            if float( x[2][site] ) > 0.6 or float( x[4][site] ) > 0.6:
                branch_musites[branchname].append( site )
    return branch_musites

#
# branch_data
#
branch_data = {}
branch_musites = {}
for branchname in branch_dirs:
    branch_data[branchname] = get_branch_data(branchname)
    branch_musites = get_branch_musites(branchname)


branches_sorted = branch_data.keys()
branches_sorted.sort()

for branch in branch_musites:
    print branch, branch_musites[branch].__len__()

#
# Excel styles:
#
from xlwt import Workbook, easyxf, Borders
hit_style1 = easyxf('pattern: pattern solid, fore_colour orange;')
hit_style2 = easyxf('pattern: pattern solid, fore_colour yellow;')
hit_style3 = easyxf('pattern: pattern solid, fore_colour yellow;')

stripe_style1 = easyxf('pattern: pattern solid, fore_colour white;')
stripe_style2 = easyxf('pattern: pattern solid, fore_colour white;')
header_style = easyxf('font: bold True;')
wrap_style = easyxf('align: wrap 1;')
center_style = easyxf('alignment: horizontal center;')

#
#
#
book = Workbook()

for seed in seed_cbiopath:
    if seed not in ref2seed:
        print "\n. I can't find the taxa", seed, "in your MSA. I'm skipping it!"
        continue
    
    mu_data = read_cbio_file( seed_cbiopath[seed] )
    if mu_data.__len__() <= 0:
        print "\n. I didn't find any data in", seed_cbiopath[seed]
    
    site_mu_data = collapse_mu_data( mu_data )
    
    branch_counthits = {}
    for branch in branches_sorted:
        branch_counthits[branch] = 0

    branch_countvalidsites = {} # count the number of cBio sites that actually existed on this branch
    for branch in branch_countvalidsites:
        branch_countvalidsites[branch] = 0

    """Open an Excel file and write data as we gather it."""
    sheet1 = book.add_sheet(seed)
    sheet1.write(0,0,"Site (cBioPortal)", header_style)
    sheet1.col(0).width = 3700
    sheet1.write(0,1,"Observed mutations (cBioPortal)", header_style)
    sheet1.col(1).width = 7000
    sheet1.write(0,2,"site (MSAProbs)", header_style)
    sheet1.col(2).width = 3800

    col = 3
    for branch in branches_sorted:
        btok = branch.split("/")
        branch_short = btok[ btok.__len__()-1 ]
        sheet1.write(0,col,branch_short, header_style)
        sheet1.col(col).width = 6000
        col += 1
    
    sheet1.write(0,col,"N hits", header_style)

    row = 1
    col = None
    """One site per row."""
    for site in site_mu_data:
        count_mu_branches = 0 # count the number of branches that have a mutation at this site.
        found_mu_for_row = False
        found_convergent_mu_for_row = False

        sheet1.write(row,0,site,center_style)

        fromaas = []
        toaas = []

        mutations = []
        for muname in site_mu_data[site]:
            data = mu_data[muname]
            fromaa = data[1]
            toaa = data[2]
            if fromaa not in fromaas:
                fromaas.append( fromaa )
            if toaa not in toaas:
                toaas.append( toaa )
            mutations.append( fromaa + "->" + toaa)
        sheet1.write(row,1, ", ".join( mutations ), wrap_style)

        """site, to state, from state"""
        refsite = seed2ref[seed][ site_mu_data[site][ site_mu_data[site].keys()[0] ][3] - 1]
        sheet1.write(row, 2, refsite+1, center_style)

        col = 2
        branch_count = 0
        for branch in branches_sorted:
            branch_count += 1
            if refsite+1 in branch_data[branch][0]:
                branch_countvalidsites[branch] += 1
                
                this_df = branch_data[branch][0][refsite+1]
                from_state = branch_data[branch][1][refsite+1]
                from_pp = branch_data[branch][2][refsite+1]
                to_state = branch_data[branch][3][refsite+1]
                to_pp = branch_data[branch][4][refsite+1]
                
                if from_state == "-":
                    from_pp = "n/a"
                    from_state = "-"
                if to_state == "-":
                    to_pp = "n/a"
                    to_state = "-"

                """Reversion?"""
                if to_state in fromaas and from_state != to_state:
                    found_convergent_mu_for_row = True
                    count_mu_branches += 1
                    branch_counthits[branch] += 1
                    col += 1
                    sheet1.write(row, col, from_state + "(" + from_pp + ") -> " + to_state + "(" + to_pp + ")", hit_style1)
                elif (from_state != to_state) and to_pp != "n/a" and ( float(to_pp) > 0.6 ):
                    found_mu_for_row = True
                    count_mu_branches += 1
                    branch_counthits[branch] += 1
                    col += 1
                    sheet1.write(row, col, from_state + "(" + from_pp + ") -> " + to_state + "(" + to_pp + ")", hit_style2)
                else:
                    st = stripe_style2
                    if branch_count%2 == 0:
                        st = stripe_style1

                    col += 1
                    sheet1.write(row, col, from_state + "(" + from_pp + ") -> " + to_state + "(" + to_pp + ")", st)
 
            else:
                st = stripe_style2
                if branch_count%2 == 0:
                    st = stripe_style1
                col += 1
                sheet1.write(row,col, "NA", st)

        if count_mu_branches > 0:
            sheet1.write(row,col+1, count_mu_branches, center_style)
        row += 1
    
    row += 1
    col = 2
    sheet.write(row, col, "hits:")
    for branch in branches_sorted:
        col += 1
        sheet1.write(row,col, branch_counthits[branch], st)

    row += 1
    col = 2
    sheet.write(row, col, "hit max possible:")
    for branch in branches_sorted:
        col += 1
        sheet1.write(row,col, branch_countvalidsites[branch], st)
        
    row += 1
    col = 2
    sheet.write(row, col, "Count mu on branch:")
    for branch in branches_sorted:
        col += 1
        sheet1.write(row,col, branch_musites[branch].__len__(), st)

    row += 1
    col = 2
    sheet.write(row, col, "Total sites:")
    for branch in branches_sorted:
        col += 1
        sheet1.write(row,col, branch_data[branch][0].keys().__len__(), st)

    book.save("ASR-cBio.11-7-2014.xls")

exit()




