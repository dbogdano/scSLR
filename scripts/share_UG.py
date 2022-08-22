
import sys
import HTSeq
from pathlib import Path

from itertools import tee, islice, chain


def next_funct(some_iterable):
	items, nexts = tee(some_iterable, 2)
	nexts = chain(islice(nexts, 1, None), [None])
	return zip(items, nexts)




def share_UG_with_mate(input_SAM, output_SAM):

	in_file = HTSeq.SAM_Reader(input_SAM)

	out_file = str(output_SAM)

	for alnmt, pair in next_funct(in_file):

		with open(out_file, 'a' ) as new_sam:

			if alnmt.has_optional_field('UG'):

				if alnmt.optional_field('RG') == '3prime':

					line = alnmt.get_sam_line()

					new_sam.write(line+'\n')

				else:
					if alnmt.pe_which == 'first':
						UG=alnmt.optional_field('UG')
						BX=alnmt.optional_field('BX')
                    
						line = alnmt.get_sam_line()
                        
						new_sam.write(line+'\n')

						if ((alnmt.get_sam_line().partition('\t')[0] == pair.get_sam_line().partition('\t')[0]) and (pair.pe_which == 'second')):

							line = pair.get_sam_line()+'\tUG:i:'+str(UG)+'\tBX:Z:'+str(BX)
						
							new_sam.write(line+'\n')
		new_sam.close()
	return 0



if __name__ == "__main__":
	#input_SAM = sys.argv[1]
	#output_SAM = sys.argv[2]
	input_SAM = snakemake.input[0]
	output_SAM = snakemake.output[0]
	
	share_UG_with_mate(input_SAM, output_SAM)
	
