##script takes a SAM file and a GTF as input, and outputs a new SAM files with GN and GX tags filled in based on HTSEQ gene overlaps
##script heavily borrows from sciRNA-seq gene assignments: https://github.com/JunyueC/sci-RNA-seq3_pipeline/blob/master/script_folder/sciRNAseq_count.py#L236

import os
import itertools
import collections
import numpy as np
import pandas as pd
import multiprocessing 
import HTSeq
import sys
from functools import partial
import logging

import contextlib
import pysam
import math
from pathlib import Path
import glob

def make_annotation_array(annotation_file):
    
    print('getting annotations')
    
    gtf_file = HTSeq.GFF_Reader(annotation_file,end_included=True)
    
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    genes = HTSeq.GenomicArrayOfSets("auto", stranded=True)

    gene_end = {}

    exon_n = 0
    gene_n = 0
    transcript_n = 0
    gene_count = 0
    
    
    gene_annotate = pd.DataFrame(columns=['gene_id','gene_type','feature','gene_name','number'])

    for feature in gtf_file:
        if feature.type == "exon":
            exon_n += 1
            exons[feature.iv] += feature.attr["gene_id"]
        elif feature.type == "gene":
            gene_n +=1
            genes[feature.iv] += feature.attr["gene_id"]
            
            gene_count += 1

            gene_annotate = gene_annotate.append(dict(zip(gene_annotate.columns,[feature.attr["gene_id"], feature.attr["gene_type"], "exon", feature.attr["gene_name"], str(gene_count)])), ignore_index=True)

            gene_count += 1

            gene_annotate = gene_annotate.append(dict(zip(gene_annotate.columns,[feature.attr["gene_id"]+ "_intron", feature.attr["gene_type"]+ "_intron", "exon", feature.attr["gene_name"]+ "_intron", str(gene_count)])), ignore_index=True)

        elif feature.type == "transcript":
            transcript_n += 1

            if feature.attr["gene_id"] in gene_end.keys():
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)
            else:
                gene_end[ feature.attr["gene_id"] ] = set()
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)



    return gene_annotate, exons, genes, gene_end




def find_nearest_gene(al_end, gene_id_intersect, gene_end):
    gene_id_end = {}
    for gene in gene_id_intersect:
        if gene in gene_end:
            gene_id_end[gene] = (abs(np.array(list(gene_end[gene])) - al_end)).min()
        else:
            continue
    # filter the gene with the least distance. If there are two genes with the least distance, then  "_ambiguous" 
    # would be returned
    gene_end_min = np.min(list(gene_id_end.values()))
    count = 0
    #print ("gene end distance: ", gene_id_end.values())
    for gene in gene_id_end:
        if (gene_id_end[gene] < gene_end_min + 100):
            count += 1
            gene_id = gene
    if count > 1:
        gene_id = "-"
    
    return gene_id




def get_gene_ID(alnmt, gene_annotate, exons, genes, gene_end):

    
    for alnmt in [alnmt]:
        if not alnmt.aligned:
            gene_id = "-"
            continue

        if alnmt.iv.chrom not in genes.chrom_vectors:
            gene_id = "-"
            continue

            
        # First check the intersectin with exons
        gene_id_intersect = set()
        gene_id_combine = set()
        inter_count = 0
        for cigop in alnmt.cigar:
            if cigop.type != "M":
                continue
            for iv,val in exons[cigop.ref_iv].steps():

                gene_id_combine |= val
                if inter_count == 0:
                    gene_id_intersect |= val
                    inter_count += 1

                else:
                    gene_id_intersect &= val

            if len(gene_id_intersect) == 1:
                gene_id = list(gene_id_intersect)[0]

            elif len(gene_id_intersect) > 1:
                gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end)

            else:
                # if there no intersection match, then find the union sets
                if len(gene_id_combine) == 1:
                    gene_id = list(gene_id_combine)[0]

                elif len(gene_id_combine) > 1:
                    gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end)

                else:
                    # if there is no intersection match or union match, then search for genes to find the intronic match
                    gene_id_intersect = set()
                    gene_id_combine = set()
                    inter_count = 0
                    for cigop in alnmt.cigar:
                        if cigop.type != "M":
                            continue
                        for iv,val in genes[cigop.ref_iv].steps():
                            gene_id_combine |= val
                            if inter_count == 0:
                                gene_id_intersect |= val
                                inter_count += 1
                            else:
                                gene_id_intersect &= val

                    if len(gene_id_intersect) == 1:
                        gene_id = list(gene_id_intersect)[0]

                    elif len(gene_id_intersect) > 1:
                        gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_intersect, gene_end)

                    else:
                        # if there no intersection match, then find the union sets
                        if len(gene_id_combine) == 1:
                            gene_id = list(gene_id_combine)[0]


                        elif len(gene_id_combine) > 1:
                            gene_id = find_nearest_gene(alnmt.iv.end_d, gene_id_combine, gene_end)

                        else:
                            gene_id = "-"


        return gene_id
        
        
        
        
#### taken from: https://bioinformatics.stackexchange.com/questions/7052/chunk-alignment-in-a-name-sorted-bam-for-parallel-processing
def splitBAM(input_BAM, chunk_number):  
    print('splitting input BAM into chunks')

    infile = pysam.AlignmentFile(input_BAM,"rb")

    
    BAM_size = int(infile.mapped + infile.unmapped)
    
    chunk_size = math.ceil(BAM_size / int(chunk_number))+1
    
    outfile_pattern = os.environ['TMPDIR']+"/output_segement%d.bam"
    
    chunk = 0
    reads_in_this_chunk = 0
    old_name = None
    outfile = pysam.AlignmentFile(outfile_pattern % chunk, "w", template = infile)

    for read in infile.fetch(until_eof=True):

        if old_name != read.query_name and reads_in_this_chunk > chunk_size:
            reads_in_this_chunk = 0
            chunk += 1
            outfile.close()
            outfile = pysam.AlignmentFile(outfile_pattern % chunk, "w", template = infile)

        outfile.write(read)
        old_name = read.query_name
        reads_in_this_chunk += 1
    
    print('finished splitting BAM')
    outfile.close()
    
    
    
    
def gene_tags_add(input_BAM, gene_annotate, exons, genes, gene_end):
    
    gene_name_lookup = dict(zip(gene_annotate['gene_id'],gene_annotate['gene_name']))
    gene_name_lookup.update({'-':'-'})
    
    in_file = HTSeq.BAM_Reader(input_BAM)
    
    #out_file = Path(input_BAM).stem + '_tag.sam'
    out_file = input_BAM+'_tag.sam'

    with open(out_file, 'a') as new_sam:
        
        print(in_file.get_header_dict(), end='', file=new_sam)
    
        for alnmt in in_file:
        
            gene_id = get_gene_ID(alnmt, gene_annotate, exons, genes, gene_end)

            new_tag_GX = str(gene_id)
            new_tag_GN = str(gene_name_lookup[new_tag_GX])
            
            if alnmt.optional_field('GX') == '-':
                
                line = alnmt.get_sam_line()
                line = line.replace('GX:A:-','GX:Z:'+new_tag_GX)
                line = line.replace('GN:A:-','GN:Z:'+new_tag_GN)
                line = line.replace('GX:Z:-','GX:Z:'+new_tag_GX)
                line = line.replace('GN:Z:-','GN:Z:'+new_tag_GN)

                new_sam.write(line+'\n')
            else:
                line = alnmt.get_sam_line()

                new_sam.write(line+'\n')


        new_sam.close()
        
    return 0
    
        


def cleanup():
    
    print("cleaning up output file")
    
    TMPDIR = os.environ['TMPDIR']
    
    
    with open('temp_SAMs.txt', 'a') as f:
        for line in glob.glob(TMPDIR+"/output_segement*_tag.sam"):
            f.write(line+'\n')
    f.close()
    
   
    pysam.merge("-","-b","temp_SAMs.txt")
    return 0
    

                                    

if __name__ == "__main__":
    gtf_file = sys.argv[1]
    input_SAM = sys.argv[2]
    n_processes = sys.argv[3]

    

    gene_annotate, exons, genes, gene_end = make_annotation_array(gtf_file)
    
    splitBAM(input_SAM, int(n_processes))
    
    pool = multiprocessing.Pool(processes=int(n_processes))
    
    input_BAM_list = glob.glob(os.environ['TMPDIR']+'/output_segement*')
    
    func = partial(gene_tags_add, gene_annotate=gene_annotate, exons=exons, genes=genes, gene_end=gene_end)
    
    print('adding tags to BAM files')
    pool.map(func, input_BAM_list)
    
    pool.close()
    pool.join()
    
    cleanup()
    
    
    
    
