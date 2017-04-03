import argparse
import csv

# Constants
DICT = 'dict.txt'

def file_to_dict(tsv):
    name_dict = {}

    with open(tsv, 'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t', dialect=csv.excel_tab)
        for row in tsvin:
            name_dict.update(dict([tuple(row[0].split())]))

    return name_dict


def find_replace(vcf, name_dict, vcf_out):
    with open(vcf, 'rb') as vcf:
        vcf_str = vcf.read()
        for key in name_dict:
            vcf_str = vcf_str.replace(key, name_dict[key])

    with open(vcf_out, 'wb') as vcf_out:
        vcf_out.write(vcf_str)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calls HWEfilter with the following arguments')
    parser.add_argument('-DICT', const=DICT, default=DICT, nargs='?', help='Name of the dict file. This always has to be a tab delimited file. First column is barcode/run name ie Xpress_00X. Second column is the sample name ie ARG109_5')
    parser.add_argument('-VCF', type=str, help='Vcf file to filter through')


    args = parser.parse_args()

    vcf_file = str(args.VCF)

    name_dict = file_to_dict(args.DICT)
    rep = find_replace(args.VCF, name_dict, vcf_file)



