import argparse
import re
import csv
import glob, os

# Constants
DICT = 'dict.txt'

def file_to_dict(tsv):
    name_dict = {}

    with open(tsv, 'rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            if row:
                name_dict.update(dict([tuple(row[0].split())]))

    return name_dict

# http://stackoverflow.com/questions/225735/batch-renaming-of-files-in-a-directory
def rename(dir, pattern, new_name):
    for pathAndFilename in glob.iglob(os.path.join(dir, pattern)):
        print 'RENAMED: {} -> {}'.format(pathAndFilename, new_name)
        os.rename(pathAndFilename, os.path.join(dir, new_name))

def rename_in_directory(name_dict, directory):
    for file in name_dict:
        pattern = r'{}*'.format(file)
        new_name = name_dict[file] + '.bam'
        rename(directory, pattern, new_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calls HWEfilter with the following arguments')
    parser.add_argument('-DICT', const=DICT, default=DICT, nargs='?', help='Name of the dict file. This always has to be a tab delimited file. First column is barcode/run name ie Xpress_00X. Second column is the sample name ie ARG109_5')
    parser.add_argument('-DIRECTORY', type=str, help='Explicit path to directory containing files to rename')

    args = parser.parse_args()

    call_directory = os.getcwd()
    directory = call_directory + '/' + str(args.DIRECTORY)

    name_dict = file_to_dict(args.DICT)
    rename_in_directory(name_dict, directory)
