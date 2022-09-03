##script takes in as input a BAM file, reads the UMI and UMI-tools group unique ID for 3prime marked reads in the BAM file,
##and counts the number of unique UMIs that exist for the group ID. Group IDs with more than 1 associated UMI are considered clashes
##and are marked for removal from the BAM file 


import sys
import HTSeq
import pandas as pd


def generate_table(input_BAM):
##function takes a BAM file as input and returns a dictionary of UMI-tools group IDs and UMI from 3prime reads
	 
	in_file = HTSeq.BAM_Reader(input_BAM)
	
	d = []
	
	for alnmt in in_file:
		
		if alnmt.has_optional_field('UG'):
			if alnmt.optional_field('RG') == '3prime':
				
				d.append(
						{	
							'gene_name': alnmt.optional_field('GN'),
							'gene_ID': alnmt.optional_field('GX'),
							'UMI': alnmt.optional_field('UB'),	
							'cell': alnmt.optional_field('CB'),
							'ID': alnmt.optioanl_field('UG')
							
						}
					)

	collision_table = pd.DataFrame(d)	
				
	perID = collision_table.groupby('ID')['UMI'].nunique()
				
	collisions = perID[perID > 1].index

	collision_table['collision'] = collision_table['ID'].isin(collisions)


	collision_table.to_csv('collision_table.csv')

	return collision_table



def mark_collisions(input_BAM, output_SAM):
	
	in_file = HTSeq.BAM_Reader(input_BAM)
	
	collision_table = generate_table(input_BAM) 
	
	collisionDict = dict(zip(collision_table['ID'],collision_table['collision']))

	out_file = str(output_SAM)
	
	for alnmt in in_file:
		
		with open(out_file, 'a' ) as new_sam:
			
			#SAM header#
			print(in_bam.get_header_dict(), file=new_sam:)

			if alnmt.has_optional_field('UG'):
				
				if alnmt.optional_field('RG') == '3prime':

					ID=alnmt.optional_field('UG')

					collision = collisionDict.get(ID)
					
					line = alnmt.get_sam_line()+'t\co:Z:'+str(collision)

					new_sam.write(line+'\n')

			
		new_sam.close()

	return 0
				
		
		
		

if __name__ == "__main__":
	#input_BAM = sys.argv[1]
	#output_SAM = sys.argv[2]
	input_BAM = snakemake.input[0]
	output_SAM = snakemake.output[0]

	UMI_share(input_BAM, output_SAM)




