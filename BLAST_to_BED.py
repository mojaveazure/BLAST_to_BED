#!/usr/bin/env python3
"""Run a BLAST search and create a BED file from the resulting hits"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")


import os
import re
import argparse
from itertools import repeat
from copy import deepcopy

try:
    from bs4 import BeautifulSoup
    from bs4 import element
    from bs4 import FeatureNotFound
    from Bio.Blast.Applications import NcbiblastnCommandline
except ImportError as error:
    sys.exit("Please install " + error.name)


DEFAULT_OUTPUT = os.getcwd() + '/output.bed'

VALS = [
    'Hsp_evalue',
    'Hsp_hit-from',
    'Hsp_hit-to',
    'Hsp_hit-frame'
    ]

NO_HIT_MESSAGE = 'No hits found'

HSP_SORT = lambda hsp: (hsp.get_chrom(), hsp.get_start(), hsp.get_end())

#   An error for no reference provided
class NoReferenceError(Exception):
    """A reference was not provided"""


#   An error I probably overuse...
class NoHitError(Exception):
    """A SNP has not been found"""


#   A custom error
class BLASTFailedError(Exception):
    """BLAST seems to have failed..."""


#   A class definition for a BLAST Hsp
class Hsp(object):
    """This is a class for a BLAST Hsp
    It continas the following information:
        Chromosome name
        Query name
        Hsp e-value
        Hsp start relative to subject
        Hsp end relative to subject
        Subject strand (forward or reverse)
        """
    def __init__(self, chrom, name, evalue, hstart, hend, hstrand):
        try:
            assert isinstance(chrom, str)
            assert isinstance(name, str)
            assert isinstance(evalue, float)
            assert isinstance(hstart, int)
            assert isinstance(hend, int)
            assert isinstance(hstrand, int)
        except AssertionError:
            raise TypeError
        try:
            assert hstrand == 1 or hstrand == -1
        except AssertionError:
            raise ValueError
        self._chrom = chrom
        self._name = name
        self._evalue = evalue
        self._start = hstart
        self._end = hend
        self._hstrand = hstrand

    def __repr__(self):
        return self._name + ":" + str(self._evalue)

    def __eq__(self, other):
        if isinstance(other, Hsp):
            return self._name == other._name and self._evalue == other._evalue
        elif isinstance(other, str):
            return self._name == other
        elif isinstance(other, float):
            return self._evalue == other
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, Hsp):
            return self._evalue < other._evalue
        elif isinstance(other, float):
            return self._evalue < other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, Hsp):
            return self._evalue <= other._evalue
        elif isinstance(other, float):
            return self._evalue <= other
        else:
            return NotImplemented

    def __hash__(self):
        return hash(self._name)

    def get_chrom(self):
        """Get the chromosome of our hsp"""
        return self._chrom

    def get_rc(self):
        """Did the query align to the forward (False) or reverse (True) strand"""
        return self._hstrand == -1

    def get_start(self):
        """Get the start of our hsp relative to forward"""
        if self.get_rc():
            return self._end
        else:
            return self._start

    def get_end(self):
        """Get the end of our hsp relative to forward"""
        if self.get_rc():
            return self._start
        else:
            return self._end

    def get_name(self):
        """Get the query name of our hsp"""
        return self._name

    def format_bed(self):
        """Format a BED file"""
        bed = [
            self._chrom
        ]
        if self.get_rc():
            bed += [str(self._end - 1), str(self._start)]
        else:
            bed += [str(self._start - 1), str(self._end)]
        bed.append(self._name)
        return '\t'.join(bed)


#   Ensure we have our BLAST database
def validate_db(db_path):
    """Find a BLAST nucleotide database"""
    try:
        assert isinstance(db_path, str)
    except AssertionError:
        raise TypeError
    db_name = os.path.basename(db_path) # Get the basename of the database
    #   Four regexes to find the four components of a nucleotide BLAST database
    nhr = re.compile(r'(%s\.[0-9\.]*nhr)' % db_name).search
    nin = re.compile(r'(%s\.[0-9\.]*nin)' % db_name).search
    nsq = re.compile(r'(%s\.[0-9\.]*nsq)' % db_name).search
    nal = re.compile(r'(%s\.*nal)' % db_name).search
    db_directory = os.path.abspath(os.path.realpath(os.path.dirname(db_path))) # Get the directory of the database
    if not db_directory: # This is necessary if the BLAST database is in the current directory
        db_directory = os.getcwd()
    print('Searching for proper database files for', db_name, 'in', db_directory, file=sys.stderr)
    db_contents = '\n'.join(os.listdir(db_directory)) # Collect the contents of this directory
    #   If any part of the database doesn't exist...
    if not nhr(db_contents) and not nin(db_contents) and not nsq(db_contents) and not nal(db_contents):
        raise FileNotFoundError("Failed to find BLAST database") # Raise an error


#   A funtion to get the value from a tag
def get_value(tag, value):
    """Get the text from a specific tag from a element.Tag object"""
    try:
        assert isinstance(tag, element.Tag)
        return tag.findChild(value).text
    except AssertionError:
        raise TypeError
    except:
        raise


#   A function to parse the HSP section of a BLAST XML file
def parse_hsp(hsp):
    """Parse the HSP section of a BLAST XML file"""
    try:
        assert isinstance(hsp, element.Tag)
    except AssertionError:
        raise TypeError
    hsp_vals = map(get_value, repeat(hsp, len(VALS)), VALS)
    return tuple(hsp_vals)


#   A function to parse the Hit section of a BLAST XML file
def parse_hit(snpid, hit):
    """Parse the Hit section of a BLAST XML file"""
    try:
        assert isinstance(snpid, str)
        assert isinstance(hit, element.Tag)
    except AssertionError:
        raise TypeError
    chrom = get_value(tag=hit, value='Hit_def')
    #   Holding lists
    vals = list()
    hsps = list()
    for hsp in hit.findAll('Hsp'):
        hsp_vals = parse_hsp(hsp=hsp)
        vals.append(hsp_vals)
    for val in vals:
        (evalue, hsp_start, hsp_end, strand) = val
        hsp = Hsp(
            chrom=chrom,
            name=snpid,
            evalue=float(evalue),
            hstart=int(hsp_start),
            hend=int(hsp_end),
            hstrand=int(strand),
        )
        hsps.append(hsp)
    if hsps:
        return hsps
    else:
        return None


#   A function to run BLASTn
def run_blastn(query, database, evalue, max_seqs, max_hsps):
    """Run BLASTn"""
    try:
        assert isinstance(query, str)
        assert isinstance(database, str)
        assert isinstance(evalue, float)
        assert isinstance(max_seqs, int)
        assert isinstance(max_hsps, int)
    except AssertionError:
        raise TypeError
    #   Create an output name
    query_base = os.path.basename(os.path.splitext(query)[0])
    database_base = os.path.basename(os.path.splitext(database)[0])
    blast_out = os.getcwd() + '/' + query_base + '_' + database_base + '_BLAST.xml'
    try:
        validate_db(database)
    except FileNotFoundError as error:
        sys.exit(error)
    #   Setup BLASTn
    blastn = NcbiblastnCommandline(
        query=query,
        db=database,
        evalue=evalue,
        outfmt=5,
        max_target_seqs=max_seqs,
        max_hsps=max_hsps,
        out=blast_out
    )
    #   Run BLASTn
    print(blastn, file=sys.stderr)
    blastn()
    if not os.path.exists(blast_out):
        raise BLASTFailedError
    return blast_out


#   Make an argument parser
def make_argument_parser():
    """Make an argument parser"""
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    inputs = parser.add_argument_group(
        title='Input options',
        description="Specify the type of input we're working with, '-f | --fasta' requires BLASTn from NCBI to be installed"
    )
    in_mut = inputs.add_mutually_exclusive_group(required=True)
    in_mut.add_argument( # Take FASTA as input
        '-f',
        '--fasta',
        dest='fasta',
        type=str,
        default=None,
        metavar='FASTA FILE',
        help="Input FASTA file to run BLAST on, incompatible with '-x | --xml'"
    )
    in_mut.add_argument( # Take XML as input
        '-x',
        '--xml',
        dest='xml',
        type=str,
        default=None,
        metavar='XML FILE',
        help="Input BLAST XML file to turn into BED file, incompatible with '-f | --fasta'"
    )
    outputs = parser.add_argument_group(
        title='Output options',
        description="Specify how we're writing our output files, create new file or append to preexisting BED file"
    )
    out_mut = outputs.add_mutually_exclusive_group(required=False)
    out_mut.add_argument( # Append to preexisting BED file
        '-b',
        '--bed',
        dest='bed',
        type=str,
        default=None,
        metavar='BED FILE',
        help="BED file to append results to, incompatible with '-o | --outfile'"
    )
    out_mut.add_argument( # Write to new BED file
        '-o',
        '--outfile',
        dest='outfile',
        type=str,
        default=DEFAULT_OUTPUT,
        metavar='OUTPUT FILE',
        help="Name of output file to write to, defaults to '" + DEFAULT_OUTPUT + "', incompatible with '-b | --bed'"
    )
    blast_opts = parser.add_argument_group(
        title='BLAST options',
        description="Options only used when BLASTING, no not need to provide with '-x | --xml'"
    )
    blast_opts.add_argument( # BLAST database to BLAST against
        '-d',
        '--database',
        dest='database',
        type=str,
        required=False,
        default=None,
        metavar='REFERENCE BLAST DATABASE',
        help="Reference BLAST database in nucleotide format, used only when running BLAST"
    )
    blast_opts.add_argument( # E-value threshold
        '-e',
        '--evalue',
        dest='evalue',
        type=float,
        required=False,
        default=1e-1,
        metavar='E-VALUE THRESHOLD',
        help="Evalue threshold for BLAST, defaults to '1e-1', used only when running BLAST"
    )
    blast_opts.add_argument( # Maximum number of hits
        '-s',
        '--max-hits',
        dest='max_hits',
        required=False,
        type=int,
        default=1,
        metavar='MAX HITS',
        help="Maximum hits per query, defaults to '1', used only when running BLAST"
    )
    blast_opts.add_argument( # Maximum number of HSPs
        '-m',
        '--max-hsps',
        dest='max_hsps',
        required=False,
        type=int,
        default=1,
        metavar='MAX HSPS',
        help="Maximum HSPs per hit, defaults to '1', used only when running BLAST"
    )
    blast_opts.add_argument( # Do we keep the XML file from BLAST?
        '-k',
        '--keep-xml',
        dest='keep_xml',
        required=False,
        action='store_const',
        const=True,
        default=False,
        metavar='KEEP XML',
        help="Do we keep the XML results? pas '-k | --keep-xml' to say 'yes', used only when running BLAST"
    )
    filters = parser.add_argument_group(
        title='Filtering options',
        description="Options to filter the resulting BED file, can specify as many or as few chromosome/contig names as you want. You may choose to either exclude or keep, not both. Note: this only affects the immediate output of this script, we will not filter a preexisting BED file"
    )
    filt_opts = filters.add_mutually_exclusive_group(required=False)
    filt_opts.add_argument(
        '-X',
        '--exclude-chrom',
        dest='exclude',
        required=False,
        type=str,
        default=None,
        nargs='*',
        metavar='EXCLUDE CHROMOSOMES',
        help="Chromosomes/Contigs to exclude from BED file, everything else is kept, incompatible with '-K | --keep-chrom'"
    )
    filt_opts.add_argument(
        '-K',
        '--keep-chrom',
        dest='keep',
        required=False,
        type=str,
        default=None,
        nargs='*',
        metavar='KEEP ONLY CHROMOSOMES',
        help="Chromosomes/Contigs to keep from BED file, everything else is excluded, incompatible with '-X | --exclude-chrom'"
    )
    return parser


#   Run the program
def main():
    """Run BLAST_to_BED.py"""
    parser = make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    # args = vars(parser.parse_args())
    args = {key : value for key, value in vars(parser.parse_args()).items() if value is not None}
    try:
        #   Are we running a BLAST search given a FASTA?
        if 'fasta' in args.keys():
            if not args['database']:
                raise NoReferenceError
            print("BLASTing", args['fasta'], "against", args['database'], file=sys.stderr)
            blast_xml = run_blastn(
                query=args['fasta'],
                database=args['database'],
                evalue=args['evalue'],
                max_seqs=args['max_hits'],
                max_hsps=args['max_hsps']
            )
        #   Or are we given a BLAST XML file?
        elif 'xml' in args.keys():
            print("Using", args['xml'], "as input", file=sys.stderr)
            blast_xml = args['xml']
        else:
            raise NotImplementedError("Whatever you're trying to do, we don't do yet...")
        #   Read in the XML as soup
        blast_soup = BeautifulSoup(open(blast_xml, 'r'), 'xml')
    except FeatureNotFound: # No lxml for XML parsing
        sys.exit("Pleast install 'lxml' to properly parse the BLAST results")
    except NoReferenceError: # No reference database passed
        sys.exit(NoReferenceError)
    except FileNotFoundError as error: # Can't open a file
        sys.exit("Cannot find " + error.filename)
    except BLASTFailedError as error:
        sys.exit(error)
    #   Collections for results
    no_hit = set()
    raw_hsps = list()
    #   Parse the XML file
    for query in blast_soup.findAll('Iteration'):
        snpid = get_value(tag=query, value='Iteration_query-def')
        try: # Ask to see if there were no hits
            if get_value(tag=query, value='Iteration_message'):
                print("No hit for", snpid, file=sys.stderr)
                no_hit.add(snpid)
                continue
        except AttributeError: # No message
            pass
        #   For every hit in the iteration
        for hit in query.findAll('Hit'):
            hit_num = get_value(tag=hit, value='Hit_num')
            this_hsps = parse_hit(snpid=snpid, hit=hit)
            try:
                raw_hsps += this_hsps
            except TypeError:
                print('No HSPs for', snpid, 'hit number:', hit_num, file=sys.stderr)
                no_hit.add(snpid)
                continue
    #   Filter our Hsps
    if 'exclude' in args.keys():
        print("Removing SNPs that BLASTed to one of:", args['exclude'], file=sys.stderr)
        hsps = list(filter(lambda hsp: hsp.get_chrom() not in args['exclude'], raw_hsps))
    elif 'keep' in args.keys():
        print("Removing SNPs that did not BLAST to one of:", args['keep'], file=sys.stderr)
        hsps = list(filter(lambda hsp: hsp.get_chrom() in args['keep'], raw_hsps))
    else:
        hsps = deepcopy(raw_hsps)
    filtered = [hsp for hsp in raw_hsps if hsp not in hsps]
    #   Sort our Hsps
    print("Sorting BLAST results by chromosome/contig, start position, and end position")
    hsps.sort(key=HSP_SORT)
    filtered.sort(key=HSP_SORT)
    #   Set our output options
    try:
        if 'bed' in args.keys():
            print("Appending", len(hsps), "searches to", args['bed'], file=sys.stderr)
            # no_hit_name = os.path.basename(args['bed']) + '_failed.log'
            out_base = os.path.splitext(args['bed'])[0]
            outhandle = open(args['bed'], 'a')
        else:
            print("Writing", len(hsps), "searches to", args['outfile'], file=sys.stderr)
            # no_hit_name = os.path.basename(args['outfile']) + '_failed.log'
            out_base = os.path.splitext(args['outfile'])[0]
            outhandle = open(args['outfile'], 'w')
    except FileNotFoundError as error:
        sys.exit("Cannot find " + error.filename)
    no_hit_name = out_base + '_failed.log'
    filtered_name = out_base + '_filtered.bed'
    #   Write the outputs
    for hsp in hsps:
        outhandle.write(hsp.format_bed())
        outhandle.write('\n')
    outhandle.close()
    #   Write any filtered BED lines
    if filtered:
        print("Writing", len(filtered), "filtered hits to", filtered_name, file=sys.stderr)
        with open(filtered_name, 'w') as f:
            for filt in filtered:
                f.write(filt.format_bed())
                f.write('\n')
    #   Write the failed log
    no_hit -= {hsp.get_name() for hsp in hsps}
    if len(no_hit) > 0:
        print("Writing", len(no_hit), "failed searches to", no_hit_name, file=sys.stderr)
        with open(no_hit_name, 'w') as n:
            for fail in sorted(list(no_hit)):
                n.write(fail)
                n.write('\n')
    if 'fasta' in args.keys() and not args['keep_xml']:
        print("Removing intermediate files", file=sys.stderr)
        os.remove(blast_xml) # Remove


if __name__ == '__main__':
    main()

