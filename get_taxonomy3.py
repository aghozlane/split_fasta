#!/usr/bin/env python3
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


"""Get NCBI taxonomy with lineage."""
import argparse
import os
import sys

try:
    from taxadb.accessionid import AccessionID
    from taxadb.taxid import TaxID
except ImportError:
    sys.exit("The program requires for the package taxadb")


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2017, Institut Pasteur"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
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


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h (see also "
                                     "ftp://ftp.ncbi.nih.gov/pub/taxonomy/)"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='blast_output_file', type=isfile,
                        required=True, help="blast_output_file")
    parser.add_argument('-d', dest='taxadb_file', type=isfile,
                        required=True, help="Taxadb file")
    parser.add_argument('-o', dest='taxonomy_file', type=str, required=True,
                        help="Output taxonomy_file")
    return parser.parse_args()


def extract_genbank_id(blast_output_file):
    """Extracting Genbank IDS from BLAST output
    """
    len_accession = 0
    accession = []
    try:
        with open(blast_output_file, "rt") as blast_output :
            for aline in blast_output :
                if aline[0] == '#' or aline.strip() == '' :
                    pass
                elif "|" in aline:
                    #print(aline.split()[1].split('|'))
                    accession.append(aline.split()[1].split('|')[2].split('.')[0])
                    #accession.add(aline.split('|')[1])
            assert(len(accession) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_output_file))
    except AssertionError:
        sys.exit("No accession id read from {0}".format(blast_output_file))
    except IndexError:
        sys.exit("{0}".format(aline))
    return accession


def write_results(chunks, accession_db, tax_db, taxonomy_file):
    """Writing results to file
    """
    try:
        with open(taxonomy_file, "wt") as output:
            #output.write("gi\ttaxid\t{0}\n".format(";".join(ranks)))
            output.write("accession\ttaxid\tAnnotation\n")
            for accession in chunks:
                taxids = accession_db.taxid(accession)
                for tax in taxids:
                    lineage = tax_db.lineage_name(tax[1], reverse=True)
                    if lineage:
                        output.write("{0[0]}\t{0[1]}\t{1}\n".format(
                                        tax, ";".join(lineage)))
    except IOError:
        sys.exit("Error cannot open {0}".format(taxonomy_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Step 1
    print("STEP 1: Extracting Genbank IDS from BLAST output...")
    accession = extract_genbank_id(args.blast_output_file)
    chunks = [accession[i:i+999] for i in range (0, len(accession), 999)]
    print("Found {0} ids !".format(len(accession)))
    #print(accession)
    # Step 3
    print("STEP 2: Writing results to file '{0}'...".format(args.taxonomy_file))
    accession_db = AccessionID(dbtype='sqlite', dbname=args.taxadb_file)
    tax_db = TaxID(dbtype='sqlite', dbname=args.taxadb_file)
    taxids = accession_db.taxid(accession)
    #taxids = [accession_db.taxid(acc) for acc in accession]
    #print(taxids)
    write_results(chunks, accession_db, tax_db, args.taxonomy_file)
    print("DONE !")


if __name__ == '__main__':
    main()
