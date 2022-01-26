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

class GetItems:
	def __init__( self, line ):
		items = line.strip().split( '\t' )
		self.query = items[0]
		self.subject = items[1]
		self.identity = float( items[2] )
		self.length = int( items[3] )
		self.qStart = int( items[6] )
		self.qEnd = int( items[7] )
		self.sStart = int( items[8] )
		self.sEnd= int( items[9] )
		self.evalue = float( items[10] )
		self.score = float( items[11] )
	def qLength( self ):
		qLength = abs( self.qEnd - self.qStart ) + 1
		return qLength
	def sLength( self ):
		sLength = abs( self.sEnd - self.sStart ) + 1
		return sLength

def CombineNumber2( numbers, diff=30 ):
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

def ParseBlast( table ):
	## blast format 6 table => 
	blast = defaultdict( board_struct )
	with open( table, 'r' ) as inp:
		for line in inp:
			hit = GetItems( line )
			###### skip self-blast results
			if hit.query == hit.subject:
				continue
			###### skip same genus
			query_genus = hit.query.split( '_' )[0]
			subject_genus = hit.subject.split( '_' )[0]
			if query_genus == subject_genus:
				continue
			
			if hit.subject in blast[ hit.query ]:
				blast[ hit.query ][ hit.subject ].append( [ hit.qStart, hit.qEnd ] )
			else:
				blast[ hit.query ][ hit.subject ] = []
				blast[ hit.query ][ hit.subject ].append( [ hit.qStart, hit.qEnd ] )
	
	for query in blast:
		for subject in blast[ query ]:
			num = []
			for i in blast[ query ][ subject ]:
				num += list( range( i[0], i[1] + 1 ) )
			blast[ query ][ subject ] = CombineNumber2( num, diff = 1 )
	return blast

def RemoveShort( list1, length ):
	new = []
	for j in list1:
		if length <= j[1] - j[0] + 1:
			new.append( j )
	return new

def CombineLists( *List ):
	numAll = []
	for eachList in List:
		for pair in eachList:
			numAll += list( range( pair[0], pair[1] + 1 ) )
	return numAll

def ReadFastaToDict( file ):
	### fasta -> dict
	#from Bio import SeqIO
	dict = {}
	with open( file, 'r' ) as inp:
		for record in SeqIO.parse( inp, 'fasta' ):
			dict[ record.id ] = record.seq
	return dict

def ExtractFasta( nums, FaName, singleFa, output, cutoff=1000 ):
	pairs = CombineNumber2( nums, 100 )
	with open( output, 'a' ) as out:
		for pair in pairs:
			start, end = pair
			name = '>' + FaName + '_' + str( start ) + '_' + str( end )
			singleFa = str( singleFa )
			seq = singleFa[ ( start - 1 ):end ]
			if len( seq ) > cutoff:
				out.write( name + '\n' + seq + '\n' )

def FastaFileLength( file ):
	length = 0
	with open( file, 'r' ) as inp:
		for record in SeqIO.parse( inp, 'fasta' ):
			length += len( record.seq )
	return length
	

if __name__ == '__main__':
	'''
		!!!(1) merge multi-chromosome into one. Note: the position of the isolated sequences is not the real location if it's multi-chrosomome.
		!!!(2) keep the sequence name and file name same in each species, e.g., file name: Alnus_jgiwogj-MT.fasta, sequence name: >Alnus_jgiwogj-MT
		1) cat file1.fasta file2.fasta file3.fasta .... > all.fasta
		2) blast results: each_single_file blast against all_files -> xxx_blastn.tsv
	'''
	startTime = time.time()
	out_length_cutoff = 300
	
	tables = glob.glob( '*blastn.tsv' )
	
	#### A group of species you want to test
	ingroup = [ 'Lithocarpus_fenestratus-MT', 'Castanea_mollissima-MT', 'Quercuss_robur-MT', 'Quercus_suber-MT', 'Quercus_variabilis-MT' ]

	#### Other species
	outgroup = [ 'Carpinus_cordata-MT', 'Corylus_avellana-MT', 'Ostrya_chinensis-MT', 'Ostryopsis_nobilis-MT', 'Casuarina_equisetifolia-MT', 'Casuarina_glauca-MT', 'Alnus_glutinosa-MT', 'Betula_pendula-MT', 'Betula_platyphylla-MT', 'Cyclocarya_paliurus-MT', 'Fagus_sylvatica-MT', 'Juglans_cathayensis-MT', 'Juglans_hindsii-MT', 'Juglans_microcarpa-MT', 'Juglans_nigra-MT', 'Juglans_regia-MT', 'Juglans_sigillata-MT', 'Platycarya_strobilacea-MT', 'Pterocarya_stenoptera-MT', 'Morella_rubra-MT' ]
	
	### delete old files
	willDel = glob.glob( '*unique.fasta' ) + glob.glob( '*intersect.fasta' )
	for fl in willDel:
		os.remove( fl )
	
	
	for table in tables:
		species = table.replace( '_blastn.tsv', '' )
		fasta = species + '.fasta'
		faseq = ReadFastaToDict( fasta )
		blast = ParseBlast( table )		## blast[ 'G1_Gossypium_bickii-MT' ][ 'Malvastrum_americanum-MT' ] = [ [ 1, 24390 ], [ 30000, 43531 ], ... ]
		
		Ingroup = ingroup.copy()
		Outgroup = outgroup.copy()
		if species in Ingroup:
			Ingroup.remove( species )
		else:
			Outgroup.remove( species )
		
		out_all_common = species + '_all_intersect.fasta'	## homologous length exists in all, differs between species because of the repeats
		out_species_unique = species + '_unique.fasta'	## species-specific sequences
		out_ingroup_unique = species + '_ingroup_unique.fasta'	## only meaningful when the species within $ingroup, specific DNA shared only by the very species and $ingroup
		out_outgroup_unique = species + '_outgroup_unique.fasta'	## only meaningful when the species within $outgroup, specific DNA shared only by the very species and $outgroup
		
		for chr in faseq:
			### common and unique
			union = []
			common = []
			for spec in blast[ chr ]:
				Nums = CombineLists( blast[ chr ][ spec ] )
				if common:
					common = list( set( common ) & set( Nums ) )
				else:
					common = Nums
				union += Nums
			ExtractFasta( common, chr, faseq[ chr ], out_all_common )
			
			unique = list( set( range( 1, len( faseq[ chr ] ) ) ) - set( union ) )
			if len( unique ) > out_length_cutoff:
				ExtractFasta( unique, chr, faseq[ chr ], out_species_unique, out_length_cutoff )
			
			in_num = []
			for cot in Ingroup:
				in_num += CombineLists( blast[ chr ][ cot ] )
			
			out_num = []
			for out in Outgroup:
				out_num += CombineLists( blast[ chr ][ out ] )
				
			onlyIngroup = list( set( in_num ) - set( out_num ) )
			if len( onlyIngroup ) > out_length_cutoff:
				ExtractFasta( onlyIngroup, chr, faseq[ chr ], out_ingroup_unique, out_length_cutoff )
			onlyOut = list( set( out_num ) - set( in_num ) )
			if len( onlyOut ) > out_length_cutoff:
				ExtractFasta( onlyOut, chr, faseq[ chr ], out_outgroup_unique, out_length_cutoff )	
	
	unique_files = glob.glob( '*unique.fasta' )
	for i in unique_files: 
		length = FastaFileLength( i ) 
		print( i + '\t' + str( length ) )
	
	intersect_files = glob.glob( '*intersect.fasta' )
	for i in intersect_files: 
		length = FastaFileLength( i ) 
		print( i + '\t' + str( length ) )

	endTime = time.time()
	runTime = endTime - startTime
	print( '\n\n##########\nrunning time: ' + str( runTime ) )
