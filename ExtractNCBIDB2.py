#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Extract the annotation in case of blast on imomi database."""

from __future__ import print_function
import os
import sys
import argparse
import csv

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


#===================
# parameters
#===================
def get_arguments():
    """Extract program options
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage="{0} -h [options] [arg]"
                                     .format(sys.argv[0]))
    parser.add_argument('-f', '--BlastResultFile', dest='blast_result_file',
                        type=isfile, required=True,
                        help='Input blast result file, in m8 mode.')
    parser.add_argument('-g', '--gi_taxid_taxonomy', dest='taxonomy_file',
                        type=isfile, required=True,
                        help='gi_taxid_taxonomy file for gi to taxid taxonomy '
                             'correspondancy')
    parser.add_argument('-nb', dest='nbest', type=int, default=0,
                        help='Number of best selected (default:0 - '
                        'based on the number of based aligned)')
    parser.add_argument('-fc', dest='filter_coverage', type=int, default=0,
                        help='Filter the coverage (default >= 0 percent).')
    parser.add_argument('-fi', dest='filter_identity', action='store_false',
                        default=True,
                        help='Filtering on the identity (default False, '
                        'superkingdom >= 0.0, phylum >= 65, class >= 75, '
                        'genus >= 85, specie >=95 percent ).')
    parser.add_argument('-id', dest='identity', type=str, default=None,
                        help="Sample identity")
    parser.add_argument('-o', '--output_file', dest='output_file', type=str,
                        help='Output file')
    parser.add_argument('-r', dest='results', type=isdir,
                        default=os.curdir + os.sep, help='Path to result '
                        'directory.')
    return parser.parse_args()


#===========================================
# parse file
#===========================================
def parse_acc_to_taxid_taxonomy_file(taxonomy_file):
    """
    """
    acc_taxid_taxonomy_dict = {}
    try:
        with open(taxonomy_file, "rt") as taxonomy:
            taxonomy_reader = csv.reader(taxonomy, delimiter='\t')
            taxonomy_reader.next()
            for line in taxonomy_reader:
                #print(line[0])
                acc_taxid_taxonomy_dict[line[0]] = line[2]
    except IOError:
        sys.exit("Error cannot open {0}".format(taxonomy_file))
    return acc_taxid_taxonomy_dict


#===========================================
# Extract blast annotation
#===========================================


def extract_annotation(blast_result_file, acc_taxid_taxonomy_dict):
    """Extract blast annotation
    """
    blast_dict = {}
    try:
        # Remove query dict
        with open(blast_result_file, "rt") as blast_result:
            blast_reader = csv.reader(blast_result, delimiter="\t")
            for line in blast_reader:
                acc = line[1].split('|')[2].split(".")[0]
                #print(acc)
                if acc in acc_taxid_taxonomy_dict:
                    annotation = acc_taxid_taxonomy_dict[acc]
                else:
                    annotation = None
                # id identity coverage
                if line[0] in blast_dict:
                    blast_dict[line[0]] += [[acc, annotation,
                                             float(line[10]), float(line[11]),
                                             float(line[12]), float(line[13])]]
                else:
                    blast_dict[line[0]] = [[acc , annotation,
                                            float(line[10]), float(line[11]),
                                            float(line[12]), float(line[13])]]
            assert(len(blast_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_result_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(blast_result_file))
    return blast_dict

#filter_identity,
def write_annotation(blast_dict, nbest, filter_coverage, filter_identity,
                     output_file, results, identity):
    """Write the result
    """
    idname = ""
    if identity:
        idname = identity + "_"
    split_thres = {0:0, 65:2, 75:3, 85:6, 95:7}
    if not output_file:
        output_file = results + os.sep + 'ncbi_taxonomic_annotation.txt'
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter="\t")
            #output_writer.writerow(["ContigName", "GI", "superkingdom", "kingdom",
            #                        "phylum", "class", "order", "family",
            #                        "genus","species", "PourcID",
            #                        "Coverage", "evalue"])
            for key in blast_dict:
                #print(key)
                #if filter_coverage > 0:
                data = []
                for element in blast_dict[key]:
                    #print(element)
                    if element[3] >= float(filter_coverage):
                        data.append(element)
                #else:
                    #data = blast_data[key]
                #print(data)
                # Sort depending on the identity + coverage
                data.sort(key=lambda x: x[2] + x[3], reverse=True)
                if nbest > 0:
                    short_set = data[0: nbest]
                else:
                    short_set = data
                #print(short_set)
                for hit in short_set:
                    if hit[1]:
                        if hit[2] >= 95.0 or filter_identity:
                            hit[1] = hit[1].split(";")
                        elif hit[2] >= 85.0:
                            hit[1] = hit[1].split(";")[0:7] + ["NA"]
                        elif hit[2] >= 75.0:
                            hit[1] = hit[1].split(";")[0:4] + ["NA"] * 4
                        elif hit[2] >= 65.0:
                            hit[1] = hit[1].split(";")[0:3] + ["NA"] * 5
                        elif hit[1] > 0.0:
                            hit[1] = hit[1].split(";")[0:2] + ["NA"] * 6
                    else:
                        hit[1]= ["NA"] * 9
                for element in short_set:
                    output_writer.writerow([idname + key, element[0]] + element[2:4] + element[1])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))



#===================
# MAIN
#===================
def main():
    """Main program
    """
    args = get_arguments()
    ## Parse gi to taxid
    #print('Parse gi to taxid file')
    #taxid_gi_dict = parse_gi_to_taxid_file(args.gi_taxid_file)
    #print(taxid_gi_dict)
    ## Parse taxid to taxonomy
    #print('Parse taxid to taxonomy file')
    #gi_taxonomy_dict = parse_taxid_to_taxonomy_file(args.taxid_taxonomy_file,
                                                     #taxid_gi_dict)
    acc_taxid_taxonomy_dict = parse_acc_to_taxid_taxonomy_file(args.taxonomy_file)
    #print(gi_taxid_taxonomy_dict)
    # Parse blast result
    blast_dict = extract_annotation(args.blast_result_file,
                                    acc_taxid_taxonomy_dict)
    #print(blast_dict)
    # Write annotation
    #args.filter_identity,
    write_annotation(blast_dict, args.nbest, args.filter_coverage,
                     args.filter_identity, args.output_file, args.results,
                     args.identity)


if __name__ == "__main__":
    main()
# END
