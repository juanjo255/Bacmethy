#!/usr/bin/env python3
## This code is to use the same format as the authors of BacMethyl
## They hard coded several extraction from the GFF file by position
import gffpandas.gffpandas as gffpd
import argparse


def transform_gff(args):
    
    df_gff = gffpd.read_gff3(args.input)
    bacmethyl_needed_attr = args.order.split(",")

    # Filter and process in a vectorized way
    mask = df_gff.df['type'].isin(['m4C', 'm6A'])
    filtered_df = df_gff.df[mask]

    ## Order established positions and
    ## Add to the original dataset
    df_gff.df.loc[mask, 'attributes'] = filtered_df['attributes'].apply(reorder_attributes, args=(bacmethyl_needed_attr,))

    
    df_gff.to_gff3(args.output)
    
    
def reorder_attributes(attr_string, bacmethyl_needed_attr):
    """Reorder attributes using the position index determined by bacmethyl_needed_attr, every attribute not listed here (extra attribute) is added to the end

    Args:
        attr_string (string): semicolon separated string with format key=value;
        bacmethyl_needed_attr (string): comma separated string with attribute's key.

    Returns:
        string: semicolon separated string with format key=value; following order from bacmethyl_needed_attr
    """
    
    attrs = attr_string.split(";")

    # Create a dictionary for quick lookup of positions
    attr_positions = {attr.split('=')[0]: i for i, attr in enumerate(bacmethyl_needed_attr)}
    
    # Sort attributes
    ordered_attrs = sorted(attrs, 
                        key=lambda x: attr_positions.get(x.split('=')[0], len(bacmethyl_needed_attr)))
    
    return ";".join(ordered_attrs)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                    prog='transform_smrt_gff',
                    description='Script to modify the order of the attributes tags to be compatible with Bacmethy',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
                    )
    parser.add_argument("-i",'--input', required=True, help='Input GFF file')
    parser.add_argument("-o",'--output', required=True, help='Output GFF file')
    parser.add_argument('--order', default='context,fracLow,motif,coverage,IPDRatio,id,fracUp,frac,identificationQv'
                        ,help='Comma-separated list with the Specific order of atributes')

    args = parser.parse_args()
    transform_gff(args)