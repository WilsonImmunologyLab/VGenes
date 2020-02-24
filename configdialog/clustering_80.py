import sys
import os
import csv
import sqlite3

import numpy as np
import scipy as sp
import scipy.cluster

import multiprocessing as mp
import itertools

import Levenshtein

def pdist(X):
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            dm[k] = Levenshtein.distance(X[i], X[j])
            k += 1
    return dm

def cluster_seqs(seqs,cutoff,linkage='single'):
    if len(seqs) == 0:
        return (np.array([]),{})

    #checks if there is only 1 unique seq
    if len(seqs) == 1:
        T = np.array([1]*len(seqs))
        return T

    #compute distance matrix
    Y = pdist(seqs)
    
    #compute linkage
    Z = sp.cluster.hierarchy.linkage(Y,method=linkage)

    # determine the clusters at level cutoff
    T = sp.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')

    return T

#get list of subgroups for each pool of clonal assignments
def get_subgroups(c, subject):
    query = "SELECT subgroup, count(*) FROM " + subject + " GROUP BY subgroup ORDER BY count(*);"
    results = c.execute(query).fetchall()
    subgroup = [x[0].encode('ascii', 'ignore') for x in results]
    return subgroup

#get sequence, celltype and CDR3 len info for clustering and post-clustering analysis
def get_subgroup_seqs(c, subgroup):
    query = 'SELECT Sequence, cell_type, CDR3_len, CDR3_NT, "V.REGION.Nb.of.AA.changes" FROM ' + subject + " WHERE subgroup = '" + subgroup + "';"
    results = c.execute(query).fetchall()
    print results
    seqs = [x[0] for x in results]
    celltype = [x[1] for x in results]
    cdr3_len = results[0][2]
    cdr3_nt = [x[3] for x in results]
    mut = [x[4] for x in results]
    return [seqs, celltype, cdr3_len, cdr3_nt, mut]

#group sequences into clones with max edit distance of the CDR3 length
def clones(data):
    max_edit = int(data[2])*.2
    print int(max_edit)
    results = cluster_seqs(data[3], int(max_edit))
    t = [int(x) for x in results]
    return t

#format data to write to csv
def format_data(subgroup_list, data, results):
    if len(data) != len(results):
        return []
    rv = []
    for i in range(len(data)):
        subgroup = subgroup_list[i]
        seqs = data[i][0]
        celltype = data[i][1]
        cdr3_len = data[i][2]
        cdr3_nt = data[i][3]
        mut = data[i][4]
        clone_assignments = results[i]
        for j in range(len(seqs)):
            if len(seqs) != len(clone_assignments):
                print("not correct order!!")
                return []
            row = [subgroup, cdr3_len, seqs[j], celltype[j], cdr3_nt[j], clone_assignments[j], mut[j]]
            rv.append(row)
    
    return rv
        
    
def main(db, subject, outfile):

    connection = sqlite3.connect(db)
    c = connection.cursor()
    connection.text_factory = str
    
    print "getting data to analyze"
    subgroup_list = get_subgroups(c, subject)
    
    data = []
    
    for subgroup in subgroup_list:
        x = get_subgroup_seqs(c, subgroup)
        data.append(x)
    
    pool = mp.Pool(processes=4)

    print "assigning clones"
    results = pool.map(clones, data)

    rv = format_data(subgroup_list, data, results)

    out = open(outfile, 'wb')
    csv_out = csv.writer(out)
    csv_out.writerows(rv)

    connection.close()
    return rv
    
if __name__ == "__main__":
    db = '/Users/denise/Documents/RepSeq2/IMGT_parsed2.sqlite'
    subject = 'IMGT_011'
    outfile = '/Users/denise/Documents/RepSeq2/RepSeq/vdj/test.csv'

    test = main(db, subject, outfile)

