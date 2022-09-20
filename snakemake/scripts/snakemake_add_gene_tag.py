##script takes a SAM file and a GTF as input, and outputs a new SAM files with GN and GX tags filled in based on HTSEQ gene overlaps
##script heavily borrows from sciRNA-seq gene assignments: https://github.com/JunyueC/sci-RNA-seq3_pipeline/blob/master/script_folder/sciRNAseq_count.py#L236


import itertools
import collections
import numpy as np
import pandas as pd
from multiprocessing import Pool
from multiprocessing import *
import HTSeq
import sys
from functools import partial
import logging

import pysam

def make_annotation_array(annotation_file):
    
    
    gtf_file = HTSeq.GFF_Reader(annotation_file,end_included=True)
    
    gene_annotat_file = './gene_name_annotate.txt'
    gene_annotat = open(gene_annotat_file, "w")
    

    exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    genes = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

    gene_end = {}

    exon_n = 0
    gene_n = 0
    transcript_n = 0
    gene_count = 0

    for feature in gtf_file:
        if feature.type == "exon":
            exon_n += 1
            exons[ feature.iv ] += feature.attr["gene_id"]
        elif feature.type == "gene":
            gene_n +=1
            genes[ feature.iv ] += feature.attr["gene_id"]
            gene_count += 1
        # for human and mouse gtf file
            message = (feature.attr["gene_id"] + "," + feature.attr["gene_type"] + "," 
                           + "exon" + "," + feature.attr["gene_name"] + "," + str(gene_count) + "\n")

            gene_annotat.write(message)

            gene_count += 1

         # for human and mouse gtf file
            message = (feature.attr["gene_id"] + "_intron" + "," + feature.attr["gene_type"] + "," 
                           + "intron" + "," + feature.attr["gene_name"] + "_intron" + "," + str(gene_count) + "\n")
            gene_annotat.write(message)

        elif feature.type == "transcript":
            transcript_n += 1
        #print "feature gene name: ", feature.attr["gene_id"]
            if feature.attr["gene_id"] in gene_end.keys():
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)
            else:
                gene_end[ feature.attr["gene_id"] ] = set()
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)



    gene_annotat.close()

    gene_annotat = pd.read_csv(gene_annotat_file, header=None)

    gene_annotat.index =  gene_annotat[0]
    
    return gene_annotat, exons, genes, gene_end




def find_nearest_gene(al_end, gene_id_intersect, gene_end):
    gene_id_end = {}
    for gene in gene_id_intersect:
        if gene in gene_end:
            gene_id_end[gene] = (abs(np.array(list(gene_end[gene])) - al_end)).min()
        else:
            print("****************Found one gene without transcript annotation*****************", "Gene name: ", gene)
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



########################### iterate through the SAM file and look for overlapping genes ####################################

def gene_tags_add(annotation_file, input_BAM, output_BAM):
    
    gene_annotat, exons, genes, gene_end = make_annotation_array(annotation_file)
    
    bam_reader = HTSeq.BAM_Reader(input_BAM)
    
    gene_name_lookup = dict(zip(gene_annotat[0],gene_annotat[3]))
    gene_name_lookup.update({'-':'-'})
    

    
    
    for alnmt in bam_reader:

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
        
        
        new_tag_GX = str(gene_id)
        new_tag_GN = str(gene_name_lookup[new_tag_GX])

        print(gene_id)
        
        with open(output_BAM, 'a') as f:
            
                
            if alnmt.optional_field('GX') == '-':
                line = alnmt.get_sam_line()
                line = line.replace('GX:A:-','GX:Z:'+new_tag_GX)
                line = line.replace('GN:A:-','GN:Z:'+new_tag_GN)
                line = line.replace('GX:Z:-','GX:Z:'+new_tag_GX)
                line = line.replace('GN:Z:-','GN:Z:'+new_tag_GN)

                f.write(line+'\n')
            else:
                line = alnmt.get_sam_line()

                f.write(line+'\n')

                                      
    f.close()
    return 0





                                    

if __name__ == "__main__":
    gtf_file = sys.argv[1]
    input_SAM = sys.argv[2]
    output_SAM = sys.argv[3]
    
    gene_tags_add(gtf_file, input_SAM, output_SAM)
