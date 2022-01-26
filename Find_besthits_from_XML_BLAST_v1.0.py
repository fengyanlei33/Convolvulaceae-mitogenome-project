#!/usr/bin/python 

################
# Author: Yanlei Feng
# Email: fengyanlei@outlook.com
# Date: 2021-02-14
################

import re
import os
import numpy as np
import time
import glob
from collections import defaultdict
from Bio import SeqIO

def board_struct():
	return defaultdict(user_struct)
def user_struct():
	return dict()

def CombineNumber2( numbers, diff=1 ):
	### CombineNumber2( [1,2,3,5,7,10,11,20,22], 2 ) -> [ [ 1,7 ], [ 10, 11 ], [ 20, 22 ] ]
	## numbers = set( numbers ) # remove duplicates
	## diff min 1
	numbers = list( set( numbers ) )	# Unique
	numbers.sort()	# put numbers in order
	start = numbers[0]
	end = numbers[0]
	out = []
	for i in range( 1, len( numbers ) ):
		if numbers[i] - numbers[i-1] < diff + 1:
			end = numbers[i]
		else:
			out.append( [ start, end ] )
			start = numbers[i]
			end = numbers[i]
	## add the last iteration
	out.append( [ start, end ] )
	return out

def ListRest( List1, List2, cutoff ):
	'''
		ListRest( [ [1, 1500] ], [ [50, 400] ], 100 )
		[401, 1500]
		
		ListRest( [ [1, 1500],[1700,2000] ], [ [1600, 1950] ], 100 )
		[1, 1500]
		
		ListRest( [ [1, 1500],[1700,2000] ], [ [80, 1950] ], 100 )
		[]
	'''
	num1 = []; num2 = []
	for pair in List1:
		num1 += list( range( pair[0], pair[1] + 1 ) )
	for pair in List2:
		num2 += list( range( pair[0], pair[1] + 1 ) )
	num = list( set( num1 ) - set( num2 ) )
	rest = []
	if num:
		out = CombineNumber2( num )
		for o in out:
			if o[1] - o[0] >= cutoff:
				rest.append( o )
	return rest

def ListsIntersect( List1, List2, cutoff ):
	## list1 = [ [1, 5], [7,10] ]; list2 = [3, 9]; ListsIntersect( list1, list2 ) => [ 1,2,3,4,5,6,7,8,9,10,12,13,15 ]
	num1 = []; num2 = []
	for pair in List1:
		num1 += list( range( pair[0], pair[1] + 1 ) )
	for pair in List2:
		num2 += list( range( pair[0], pair[1] + 1 ) )
	num = list( set( num1 ) & set( num2 ) )
	inter = []
	if num:
		out = CombineNumber2( num )
		for o in out:
			if o[1] - o[0] >= cutoff:
				inter.append( o )
	return inter
	
import argparse
if __name__ == '__main__':
	parser = argparse.ArgumentParser( description = 'Parse blast XML format into table' )
	parser.add_argument( '-i', '--input', dest='input', metavar='', help='xml' )
	parser.add_argument( '-o', '--output', dest='output', metavar='', help='output name' )
	
	args = parser.parse_args()
	
	blast = defaultdict( board_struct )
	cutoff = 100
	with open( args.input, 'r' ) as inp, open( args.output, 'w' ) as out:
		inp.readline()  # header
		for i in inp:
			item = i.strip().split( '\t' )
			#		0	  1		 2		  3		4	   5	   6	 7	  8	 9	 10
			# item: Query Length Accession Identity Aln_len Q_start Q_end Evalue Score Genus Title
			
			# filter
			if 'Fagus sylvatica isolate FASYL_29_1 mitochondrion' in item[10] or 'Betula pendula genome assembly, organelle: mitochondrion' in item[10]:
				continue
			
			mito = 1
			if 'mito'.upper() not in item[10].upper():
				mito = 0
			
			for j in [ 1, 3, 4, 5, 6 ]:
				item[j] = int( float( item[j] ) )
			if item[0] in blast:
				if blast[ item[0] ]:
					newRecord = ListsIntersect( blast[ item[0] ], [ [ item[5], item[6] ] ], cutoff )
					if newRecord:
						for i in newRecord:
							out.write( f'{item[0]}\t{i[0]}\t{i[1]}\t{mito}\t{item[2]}\t{item[9]}\t{item[10]}\n' )
							blast[ item[0] ] = ListRest( blast[ item[0] ], [ [ item[5], item[6] ] ], cutoff )
				else:
					continue
			else:
				if item[4] >= cutoff:
					out.write( f'{item[0]}\t{item[5]}\t{item[6]}\t{mito}\t{item[2]}\t{item[9]}\t{item[10]}\n' )
					blast[ item[0] ] = ListRest( [ [ 1, item[1] ] ], [ [ item[5], item[6] ] ], cutoff )
