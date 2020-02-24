#!/usr/local/bin/python
# filename: clonify.py


###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################


import os
import time
import math
import argparse
import itertools
from multiprocessing import Pool, cpu_count

from pymongo import MongoClient

import numpy as np

import fastcluster as fc
from scipy.cluster.hierarchy import fcluster

from Bio import pairwise2
from Levenshtein import distance



parser = argparse.ArgumentParser("")
parser.add_argument('-d', '--db', dest='db', required=True,
					help="The MongoDB database to be queried. Required.")
parser.add_argument('-c', '--coll', dest='coll', default=None,
					help="The MongoDB collection to be queried. If ommitted, all collections in the database will be processed.")
parser.add_argument('-i', '--ip', dest='ip', default='localhost',
					help="The IP address of the MongoDB server. Defaults to 'localhost'.")
parser.add_argument('-p', '--port', dest='port', default=27017, type=int,
					help="The port used to connect to the MongoDB server. Defaults to '27017'.")
parser.add_argument('-o', '--out', dest='output', default='',
					help="Directory for the output files. Files will be named '<collection>_clones.txt'. \
					Failing to provide an output directory will result in no output files being written.")
parser.add_argument('-s', '--split_by', dest='split_by', default='none', choices=['none', 'fam', 'gene'],
					help="Define how to split the input sequence pool. Choices are 'fam' and 'gene'. \
					Splitting by 'gene' may cause a small number of false negatives due to IgBLAST germline misassignment, \
					but can provide large reductions in processing time and memory use. Default is 'fam'.")
parser.add_argument('-t', '--threads', dest='threads', type=int, default=None,
					help="Number of threads to use. Default is max available.")
parser.add_argument('-x', '--dist', dest='distance_cutoff', default=0.35, type=float,
					help="The cutoff adjusted edit distance (aED) for segregating sequences into clonal families. \
					Defaults to 0.35.")
parser.add_argument('-z', '--no_split', dest='no_split', action='store_true', default=False,
					help="Use to process all sequences using a single thread. Can reduce errors with very small sequence sets.")
parser.add_argument('-n', '--nt', dest='is_aa', action='store_false', default=True,
					help="Uses nucleotide CDR3 sequences instead of amino acid sequences (default).")
parser.add_argument('-u', '--no_update', dest='update', action='store_false', default=True,
					help="Use to skip updating the MongoDB database with clonality info.")
parser.add_argument('-k', '--chunksize', dest='chunksize', type=int, default=500,
					help="Splits the input_seqs into roughly <chunksize>-sized pieces and builds distance submatrices \
					for each pair of chunks separately. Typically, this shouldn't be changed.")
args = parser.parse_args()


class Seq(object):
	'''
	Contains genetic characteristics for a single sequence.
	Input:
	data = a MongoDB result (dict-like) containing the following fields:
				[seq_id, v_gene, j_gene, <junc_query>, var_muts_nt]
	   	   where <junc_query> is the sequence of the nucleotide or AA junction.
	junc_query = either 'junc_aa' or 'junc_nt' for nucleotide or AA junctions, respectively.
	'''
	def __init__(self, data, junc_query):
		self.id = data['seq_id']
		self.v_fam = data['v_gene']['fam']
		self.v_gene = data['v_gene']['gene']
		self.v_all = data['v_gene']['all']
		self.j_gene = data['j_gene']['gene']
		self.j_all = data['j_gene']['all']
		self.junc = data[junc_query]
		self.junc_len = len(self.junc)
		self.muts = []
		if 'var_muts_nt' in data.keys():
			self.muts = ['{}{}'.format(d['loc'], d['mut']) for d in data['var_muts_nt']]

	def v_gene_string(self):
		return 'v{0}-{1}'.format(self.v_fam, self.v_gene)

	def v_fam_string(self):
		return 'v{0}'.format(self.v_fam)




# -----------------------------
#           CLONIFY
# -----------------------------



def clonify((ichunk, jchunk)):
	results = []
	for i in ichunk:
		results.append(get_scores(i, jchunk))
	return results


def get_scores(i, jchunk):
	results = []
	for j in jchunk:
		if i.id == j.id:
			results.append(0.0)
			continue
		LD = get_LD(i, j)
		vPenalty = vCompare(i, j)
		jPenalty = jCompare(i, j)
		lenPenalty = math.fabs(i.junc_len - j.junc_len) * 2
		editLength = min(i.junc_len, j.junc_len)
		mutBonus = sharedMuts(i, j)
		if mutBonus > (LD + vPenalty + jPenalty):
			mutBonus = (LD + vPenalty + jPenalty - 0.001)  # distance values can't be negative
		results.append((LD + vPenalty + jPenalty + lenPenalty - mutBonus) / editLength)
	return results


def get_LD(i, j):
	'''
	Calculate sequence distance between a pair of Seq objects
	'''
	# pairwise2 is used to force 'gapless' distance when sequence pair is of the same length
	if i.junc_len == j.junc_len:
		identity = pairwise2.align.globalms(i.junc, j.junc, 1, 0, -50, -50,
											score_only=True,
											one_alignment_only=True)
		return i.junc_len - identity
	# Levenshtein distance is used for sequence pairs of different lengths
	else:
		return distance(i.junc, j.junc)


def vCompare(i, j):
	'Calculate penalty for mismatches in Variable segment.'
	if i.v_gene != j.v_gene:
		return 10
	return 0


def jCompare(i, j):
	'Calculate penalty for mismatches in Joining segment.'
	if i.j_gene != j.j_gene:
		return 8
	return 0


def sharedMuts(i, j):
	'Calculate bonus for shared mutations.'
	if i.id == j.id:
		return 0.0
	bonus = 0.0
	for mut in i.muts:
		if mut == '':
			continue
		if mut in j.muts:
			bonus += 0.35
	return bonus


# -----------------------------
#             MONGO
# -----------------------------



def get_seqs(database, collection):
	conn = MongoClient(args.ip, args.port)
	db = conn[database]
	c = db[collection]
	if args.is_aa:
		junc_query = 'junc_aa'
	else:
		junc_query = 'junc_nt'
	results = c.find({'chain': 'heavy'},
					 {'_id': 0, 'seq_id': 1, 'v_gene': 1, 'j_gene': 1, junc_query: 1, 'var_muts_nt': 1})
	output = []
	for r in results:
		if r is not None:
			if junc_query not in r.keys():
				continue
			output.append(Seq(r, junc_query))
	return output


def get_collections():
	if args.coll:
		return [args.coll, ]
	conn = MongoClient(args.ip, args.port)
	db = conn[args.db]
	subjects = db.collection_names(include_system_collections=False)
	return sorted(subjects)


def update_db(database, collection, clusters):
	conn = MongoClient(args.ip, args.port)
	db = conn[database]
	c = db[collection]
	clust_sizes = []
	clust_count = 0
	for clust_id in clusters:
		clust_size = len(clusters[clust_id])
		seq_ids = [s.id for s in clusters[clust_id]]
		if clust_size > 1:
			if args.update:
				c.update({'seq_id': {'$in': seq_ids}},
						 {'$set': {'clonify': {'id': clust_id, 'size': clust_size}}},
						 multi=True)
			clust_count += 1
			clust_sizes.append(clust_size)
	clustered_seqs = sum(clust_sizes)
	avg_clust_size = float(clustered_seqs) / clust_count
	max_clust_size = max(clust_sizes)
	return [clust_count, clustered_seqs, avg_clust_size, max_clust_size]




# -----------------------------
#         CLUSTERING
# -----------------------------



def count_cpus():
	if args.threads:
		return args.threads
	return cpu_count()


def split_by_gene(seqs):
	'''
	Splits sequences by variable gene.
	Returns a dict with var genes as keys, seqs as values.
	'''
	split = {}
	for seq in seqs:
		if seq.v_gene_string() in split:
			split[seq.v_gene_string()].append(seq)
		else:
			split[seq.v_gene_string()] = [seq, ]
	return split


def split_by_fam(seqs):
	'''
	Splits sequences by variable family.
	Returns a dict with var fams as keys, seqs as values.
	'''
	split = {}
	for seq in seqs:
		if seq.v_fam_string() in split:
			split[seq.v_fam_string()].append(seq)
		else:
			split[seq.v_fam_string()] = [seq, ]
	return split


def get_chunksize(input):
	'''
	Calculates an appropriate sequence 'chunk' size, based on the number of
	input sequences and the number of available processors.
	'''
	if args.no_split:
		return len(input)
	if len(input) < args.chunksize:
		return len(input)
	# To cover a rare case where taking math.ceil() will result in too few chunks of sequence data.
	# Overly simple example:
	# 10 seqs, 6 CPUs will result in creation of 5 chunks of 2 sequences each, with one CPU unused.
	s = float(len(input)) / cpu_count()
	if int(math.ceil(s)) * (cpu_count() - 1) > len(input):
		return int(s)
	return int(math.ceil(s))


def chunk_maker(n, iterable, fillvalue=None):
	'''
	chunk_maker(3, 'ABCDEFG', 'x') --> ['ABC', 'DEF', 'Gxx']
	where x = fillvalue
	'''
	args = [iter(iterable)] * n
	return [[e for e in t if e is not None] for t in itertools.izip_longest(*args)]


def grouper_nofill(n, iterable):
 	'''
 	list(grouper_nofill(3, 'ABCDEFG')) --> [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
 	'''
 	it = iter(iterable)
 	def take():
		while 1:
			yield list(itertools.islice(it, n))
 	return iter(take().next, [])


def build_cluster_dict(count, vh):
	clusters = {}
	for c in xrange(1, count):
		clusters["lineage_{0}_{1}".format(vh, str(c))] = []
	return clusters


def build_matrix(ichunks, chunksize, size, chunk_count):
	matrix = np.zeros((size, size))
	print 'number of processes:', chunk_count
	print 'matrix:', matrix.shape
	print 'total calculations:', matrix.size
	p = Pool(processes=chunk_count)
	for x, seq in enumerate(grouper_nofill(chunk_count, ichunks)):
		result = p.imap(clonify, seq)
		for i, r in enumerate(result):
			matrix[x * chunksize:x * chunksize + len(r), i * chunksize:i * chunksize + len(r[0])] = r
	p.close()
	p.join()
	return matrix


def squareform(matrix):
	'''
	Transform a squareform distance matrix into a condensed matrix. Input is an array of size N**2
	(N = number of sequences). Output is an array of size (N * (N-1)) / 2.
	'''
	cm = []
	for i, row in enumerate(matrix[:-1]):
		cm.extend(row[i + 1:])
	return cm


def make_clusters(input_seqs, vh):
	chunksize = get_chunksize(input_seqs)
	print 'Chunksize is:', chunksize
	chunks = chunk_maker(chunksize, input_seqs)
	iter_chunks = itertools.product(chunks, repeat=2)
	distMatrix = build_matrix(iter_chunks, chunksize, len(input_seqs), len(chunks))
	print 'condensing the distance matrix...'
	con_distMatrix = squareform(distMatrix)
	print 'clustering...'
	linkageMatrix = fc.linkage(con_distMatrix, method='average', preserve_input=False)
	flatCluster = fcluster(linkageMatrix, args.distance_cutoff, criterion='distance')
	print 'building cluster dict...'
	clusters = build_cluster_dict(max(flatCluster) + 1, vh)
	print 'assigning sequences to clusters...'
	clusters = assign_seqs(flatCluster, clusters, input_seqs, vh)
	return clusters


def assign_seqs(flatCluster, clusters, input_seqs, vh):
	'''
	Partition sequences into lineage clusters
	'''
	for s in xrange(len(flatCluster)):
		s_id = 'lineage_{0}_{1}'.format(vh, str(flatCluster[s]))
		clusters[s_id].append(input_seqs[s])
	return clusters


def analyze_collection(coll):
	'''
	Control function for complete processing of a single MongoDB collection.
	'''
	# query
	startTime = time.time()
	print_query_start(coll)
	seqs = get_seqs(args.db, coll)
	print_query_end(len(seqs))
	queryEnd = time.time()

	# seq splitting
	if len(seqs) < 2:
		return False
	if args.split_by == 'gene':
		split_seqs = split_by_gene(seqs)
	elif args.split_by == 'fam':
		split_seqs = split_by_fam(seqs)
	else:
		split_seqs = {'v0': seqs}

	# clonify
	print_clonify_start()
	clusters = {}
	for vh in sorted(split_seqs.keys()):
		if len(split_seqs[vh]) <= 1:
			continue
		print_vh_info(vh)
		clusters.update(make_clusters(split_seqs[vh], vh))
	print_clonify_end()
	clonifyEnd = time.time()

	# update MongoDB
	if args.update:
		print_update_start()
	else:
		print_no_update()
	cluster_stats = update_db(args.db, coll, clusters)
	print_update_end()
	updateEnd = time.time()

	# write output
	if args.output != '':
		print_write_start()
		write_output(args.output, coll, clusters)
		print_write_end()
	else:
		print_write_null()
	print_run_summary(cluster_stats, startTime, queryEnd, clonifyEnd, updateEnd, len(seqs))




# -----------------------------
#       WRITING/PRINTING
# -----------------------------



def write_output(out_dir, collection, data):
	out_file = os.path.join(out_dir, collection + '_clones.txt')
	rString = ''
	for c in data.keys():
		if len(data[c]) < 2:
			continue
		rString += '#{}\n'.format(c)
		for seq in data[c]:
			rString += '>{0}\n{1}\n'.format(seq.id, seq.junc)
		rString += '\n'
	open(out_file, 'w').write(rString)


def print_query_start(collection):
	print ''
	print ''
	print '========================================'
	print 'processing collection: {0}'.format(collection)
	print '========================================'
	print ''
	print 'Querying MongoDB (db: {0}, collection: {1}) for the input sequences...'.format(args.db, collection)


def print_vh_info(vh):
	print ''
	print '--------'
	print vh
	print '--------'


def print_query_end(seq_count):
	print '...done. Retrieved {} sequences.\n'.format(seq_count)


def print_clonify_start():
	print 'Sorting sequences into clonal families...'


def print_clonify_end():
	print '...done.\n'


def print_update_start():
	print 'Updating MongoDB...'


def print_no_update():
	print 'Calculating cluster statistics...'


def print_update_end():
	print '...done.\n'


def print_write_start():
	print 'Writing clonal families to file...'


def print_write_end():
	print '...done.\n'


def print_write_null():
	print 'No output directory was provided. Lineage assignments are not being written to file.\n'


def print_run_summary(stats, startTime, queryEnd, clonifyEnd, updateEnd, total_count):
	print ''
	print 'Querying MongoDB took {} seconds.'.format(queryEnd - startTime)
	print '{0} sequences were segregated into {1} clonal families in {2} seconds.'.format(total_count, stats[0], clonifyEnd - queryEnd)
	if args.update:
		print 'Updating MongoDB took {} seconds.'.format(updateEnd - clonifyEnd)
	print ''
	print 'The average cluster size was %0.2f.' % (stats[2])
	print 'The largest cluster contains {} seqeunces.'.format(stats[3])
	print '%s sequences were assigned to clonal families (%0.2f%% of all sequences).' % (stats[1], 100.0 * stats[1] / total_count)
	print ''
	print ''



def main():
	for c in get_collections():
		analyze_collection(c)

if __name__ == '__main__':
	main()