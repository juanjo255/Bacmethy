## This code is to use the same format as the authors of BacMethyl
## They hard coded several extraction from the GFF file by position

import gffpandas.gffpandas as gffpd
import argparse

def transform_gff(args):
    
    df_gff = gffpd.read_gff3('/n/projects/jp2992/MOLNG4331/methylation_analysis/ipdsummary_out/B02/B02_basemods_w_motif.gff')
    bacmethyl_needed_attr = args.order.split(",")

    # Filter and process in a vectorized way
    mask = df_gff.df['type'].isin(['m4C', 'm6A'])
    filtered_df = df_gff.df[mask]

    def reorder_attributes(attr_string):
        attrs = attr_string.split(";")
        # Create a dictionary for quick lookup of positions
        attr_positions = {attr.split('=')[0]: i for i, attr in enumerate(bacmethyl_needed_attr)}
        
        # Sort attributes
        ordered_attrs = sorted(attrs, 
                            key=lambda x: attr_positions.get(x.split('=')[0], len(bacmethyl_needed_attr)))
        return ";".join(ordered_attrs)

    # Apply vectorized operation
    df_gff.df.loc[mask, 'attributes'] = filtered_df['attributes'].apply(reorder_attributes)
    df_gff.to_gff3(args.output)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                    prog='transform_smrt_gff',
                    description='Script to modify the order of the attributes tags to be compatible with Bacmethy',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
                    )
    parser.add_argument("-o",'--output', required=True, help='Output GFF file')
    parser.add_argument('--order', default=['context', 'fracLow', 'motif', 'coverage', 'IPDRatio', 'id', 'fracUp', 'frac', 'identificationQv']
                        ,help='Comma-separated list with the Specific order of atributes')


    
    args = parser.parse_args()
    transform_gff(args)