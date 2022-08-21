##script takes in as input a BAM file, reads the UMI and UMI-tools group unique ID for 3prime marked reads in the BAM file, creates a dictionary of these values, and shares the UMI values between internal and 5prime reads with the same UMI-tools group unique ID.

import sys
import HTSeq
from pathlib import Path


def make_dicts(input_BAM):
##function takes a BAM file as input and returns a dictionary of UMI-tools group IDs and UMI from 3prime reads
     
    in_file = HTSeq.BAM_Reader(input_BAM)
    
    Barcode_dict = {}
    UMI_dict = {}
    
    for alnmt in in_file:
        
        if alnmt.has_optional_field('UG'):
            if alnmt.optional_field('RG') == '3prime':


                Barcode=alnmt.optional_field('CR')
                UMI=alnmt.optional_field('UR')
                ID=alnmt.optional_field('UG')

                new_entry_Barcode = {ID:Barcode}
                new_entry_UMI = {ID:UMI}

                Barcode_dict.update(new_entry_Barcode)
                UMI_dict.update(new_entry_UMI)
            
            
    
    return Barcode_dict, UMI_dict



def UMI_share(input_BAM):
    
    in_file = HTSeq.BAM_Reader(input_BAM)
    
    Barcode_dict, UMI_dict = make_dicts(input_BAM)
    
    out_file = Path(input_BAM).stem + '_tag.sam'
    
    
    for alnmt in in_file:
        
        with open(out_file, 'a' ) as new_sam:
            
            print(in_file.get_header_dict(), end='', file=new_sam)
            
            if alnmt.has_optional_field('UG'):
                
                if alnmt.optional_field('RG') == '3prime':

                    line = alnmt.get_sam_line()

                    new_sam.write(line+'\n')

                else:
                    ID=alnmt.optional_field('UG')

                    if (ID in Barcode_dict) & (ID in UMI_dict):
                        Barcode=Barcode_dict.get(ID) 
                        UMI=UMI_dict.get(ID)

                        line = alnmt.get_sam_line()+'\tCR:Z:'+Barcode+'\tUR:Z:'+UMI

                    else:

                        line = alnmt.get_sam_line()+'\tCR:Z:-\tUR:Z:-'

                    new_sam.write(line+'\n')
            
        new_sam.close()

    return 0
                
        
        
        

if __name__ == "__main__":
    input_BAM = sys.argv[1]
    
    UMI_share(input_BAM)
