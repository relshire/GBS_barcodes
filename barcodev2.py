from itertools import product
from operator import itemgetter
import sys
from optparse import OptionParser
import sys,random,Bio
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from optparse import OptionParser
import csv
import re

barcodes = []
nts = [ 'A', 'C', 'G', 'T']
#barcodes = [ 'ACGTATATAT']
unambiguate = { 'R':(['A', 'G']), 'Y':(['C', 'T']), 'S':(['G', 'C']), 'W':(['A', 'T']), \
'K':(['G', 'T']), 'M':(['A', 'C']), 'B':(['C', 'G', 'T']), 'D':(['A', 'G', 'T']),\
'H':(['A', 'C', 'T']), 'V':(['A', 'C', 'G']), 'N':(['A', 'C', 'G', 'T'])}

def main():
    #parse command line options
    try:
        opts, args = parse_options(sys.argv[1:])
    except SyntaxError:
        pass
    #Check if RE is blunt and no restriction in
    check, opts  = check_combination(opts)
    if  check:
        raise Exception(check)
    barcodes, dist_simple = gen_barcodes(opts)
    csv_out(opts, barcodes, dist_simple)
#    xlsx_out(barcode, dist_simple)
#    stat_out = ''
#    for pos, dic in dist.items():
#        stat_out += str(pos) + '\t'
#        for nt, val in dic.items():
#            stat_out += '%.2f\t'%(val)
#        stat_out += '\n'
 #   print stat_out
 
def extend_barcode_RS(opts):
    """Determines: 
    1. Which nt's are added to Barcode to determine stats
    2. How the 5' and 3' strand of the Barcode adapter look like"""
    enzyme = opts.enzyme
    res_site = enzyme.elucidate()
    if res_site.index('^') < res_site.index('_'):
        RS_stat_added = res_site.split('^')[1].replace('_','')
        RS_min_added = res_site[res_site.index('^')+1:res_site.index('_')]
        RS_plus_added = ''
    else:
        RS_stat_added = res_site.split('_')[1].replace('^','')
        RS_plus_added = res_site[res_site.index('_')+1:res_site.index('^')]
        RS_min_added = ''
    return RS_stat_added, RS_plus_added, RS_min_added

def csv_out(opts, barcodes, dist_simple):
    """Generates csv style output"""
    with open(opts.output, 'wb') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerow(["sep=,"])
        csvwriter.writerow(["Barcode script Made by Thomas P. van Gurp see www.deenabio.com\gbs-adapters"])
        csvwriter.writerow([])
        line = '\t\tNucleotide'
        for pos in dist_simple.keys():
            line += '\t%s'%(int(pos)+1)
        csvwriter.writerow(line.split('\t'))
        for nt in dist_simple[0].keys():
            line = '\t\t%%%s\t'%nt
            GC_line = '\t\t%%%s\t'%'GC'
            for pos in sorted(dist_simple.keys()):
                line += '%.2f\t'%(dist_simple[pos][nt]\
                / float(sum(dist_simple[pos].values())))
                GC_line += '%.2f\t'%((int(dist_simple[pos]['G'])+\
                                  int(dist_simple[pos]['C']))/\
                                  float(sum(dist_simple[pos].values())))
            csvwriter.writerow(line.split('\t'))
        csvwriter.writerow(GC_line.split('\t'))
        csvwriter.writerow(['', '', '* for the calculation of these barcode statistics the enzyme recognition site is incorporated.'])
        csvwriter.writerow([])
        csvwriter.writerow(["Number","barcode", "plus_strand_primer 5'--> 3'","negative_strand_primer 5'--> 3'"])
        RS_stat_added, RS_plus_added, RS_min_added = extend_barcode_RS(opts)
        csvwriter.writerow(["Common",'-', RS_min_added + opts.cm_ad_seq,\
                           '%s%s'%(Seq(opts.cm_ad_seq).\
                                   reverse_complement().tostring(),\
                                   RS_plus_added)])
        csvwriter.writerow(["PCR-primer A",'-',"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"])
        csvwriter.writerow(["PCR-primer B",'-',"CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"])
        csvwriter.writerow([])
#    
#        
#        csvwriter.writerow(["Common adapter"])
#        csvwriter.writerow(["plus_strand_primer 5'--> 3'","negative_strand_primer 5'--> 3'"])
#        RS_stat_added, RS_plus_added, RS_min_added = extend_barcode_RS(opts)
#        csvwriter.writerow([RS_min_added + opts.cm_ad_seq,\
#                           '%s%s'%(Seq(opts.cm_ad_seq).\
#                                   reverse_complement().tostring(),\
#                                   RS_plus_added)])
#        csvwriter.writerow([])
#        csvwriter.writerow(["Nummer", "barcode", "plus_strand_primer 5'--> 3'","negative_strand_primer 5'--> 3'"])
#        csvwriter.writerow(["PCR-primer A",'-',"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"])
#        csvwriter.writerow(["PCR-primer B",'-',"CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"])
#        csvwriter.writerow([])
        for n, barcode in enumerate(barcodes):
                    rev_co_bc =  str(Seq(barcode, generic_dna).reverse_complement())
                    barcode_primer_min = str(Seq(opts.ad_seq, generic_dna).reverse_complement())
                    print "plus-strand 5-\t",opts.ad_seq + barcode + RS_plus_added
                    print "neg-strand 5-\t", RS_min_added + rev_co_bc + barcode_primer_min
                    csvwriter.writerow([n+1,barcode, opts.ad_seq + barcode +\
                                         RS_plus_added, RS_min_added + rev_co_bc +\
                                         barcode_primer_min ])
    print '%s barcodes generated'%(n+1)
    f.close()


def check_combination(opts):
    """Checks if the combination of RE and adapters is valid"""
    #Check if RE has overhangs.
    out = ''
    if opts.enzyme.is_blunt():
        out+= "Enzyme %s produces no overhang, it is adviced not to use it" %(opts.enzym)
    #Check if enzyme is nested within the Adapter sequences
    if  opts.enzyme.search(Seq(opts.ad_seq, IUPACAmbiguousDNA())):
        out += """Enzyme %s will cut in Barcode adapter sequence at position %s \n""" %(opts.enzyme, \
                                                                             opts.enzyme.search(Seq(opts.ad_seq, IUPACAmbiguousDNA())))
    if  opts.enzyme.search(Seq(opts.cm_ad_seq, IUPACAmbiguousDNA())):
        out += """Enzyme %s will cut in common adapter sequence at position %s \n""" %(opts.enzyme, \
                                                                             opts.enzyme.search(Seq(opts.cm_ad_seq, IUPACAmbiguousDNA())))
    #Routine for checking if Restriction site will no be recreated by common adapter ligation.
    #insert of C/G or G/C in case of recreating restriction site with common adapter ApeKI
    if opts.enzyme.overhang() == "5' overhang":
        if opts.enzyme.elucidate().split('^')[0] == opts.ad_seq[:len(opts.enzyme.elucidate().split('^')[0])]:
            for nt in nts:
                if nt != opts.enzyme.elucidate().split('^')[0]:
                    opts.ad_seq == nt + opts.ad_seq
                    break
    elif opts.enzyme.overhang() == "3' overhang":
        if opts.enzyme.elucidate().split('_')[0] == \
        str(Seq((opts.ad_seq[:len(opts.enzyme.elucidate().split('_')[0])]), IUPACAmbiguousDNA()).reverse_complement()):
            for nt in nts:
                if nt != opts.enzyme.elucidate().split('_')[0]:
                    opts.ad_seq == opts.ad_seq + nt
                    break
    return out, opts

def parse_options(opts):
    """Parses options provided on command-line"""
    parser = OptionParser()
    parser.add_option("-e", "--enzyme",  metavar = "enzyme",  
                      action = "store", type="string", default= "FseI",
                      dest = "enzyme", help="Choose restriction enzyme")
    parser.add_option("-n", "--number",  metavar = "number", default= 144,  
                      type = "int", action = "store", 
                      help = "Number of output Barcodes")
    parser.add_option("--min",  metavar = "min",  action = "store", 
                      default = 4, dest = "min", type = "int",
                      help = "Minimum Barcode length")
    parser.add_option("--max",  metavar = "max", default = 15,
                       type = "int", action = "store", dest = "max", 
                       help = "Maximum Barcode length")
    parser.add_option("-m", "--mono_nt",  metavar = "mono_nt", default = 2,
                       type = "int", action = "store", dest = "mono_nt", 
                       help = "Maximum number of mononucleotides in Barcode")
    parser.add_option("-a","--barcode_adapter_sequence",  metavar = "ad_seq",
                      action = "store", dest = "ad_seq",
                      default = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",                      
                      help = "Choose alternative 5'-->3'sequence for barcode adapter")
    parser.add_option("--common_adapter_sequence",  metavar = "cm_ad_seq",
                      action = "store", dest = "cm_ad_seq",
                      default = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",                      
                      help = "Choose alternative 5'-->3'sequence for common adapter")
    parser.add_option("-o", "--output",  metavar = "output",  
                      action = "store", type="string",default = 'output.csv',  
                      dest = "output", help = "Output file")
    opts, args = parser.parse_args()
    for enzyme in Restriction.AllEnzymes:
        if "%s"%(enzyme) == opts.enzyme:
            opts.enzyme = enzyme
            break
    return opts, args

def is_dissim(barcode,existing):
    """Indicates number of dissimilarities between 2 sequences"""
    for bc in existing:
        if len(bc)<>len(barcode):
            continue
        dis = 0
        for i in range(len(barcode)):
            if barcode[i] <> bc[i]:
                dis += 1
        if dis < 3:
            return 1
    return 0

def has_trint(barcode):
    """Has trinucleotide"""
    for nt in nts:
        if 3*nt in barcode:
            return 1
    return 0

def is_nested(barcode, existing):
    """Checks if no smaller existing barcodes are nested within this barcode"""
    if existing == []:
        return 0
    bc = existing[0]
    if len(bc) < len(barcode):
        for bc in existing:
            if len(bc) > barcode:
                break
            if bc in barcode:
                return 1
    return 0

def is_low_complex(barcode):
    """Checks if the barcode contains all nucleotides if len>5"""
    count = 0
    for nt in nts:
        if nt in barcode:
            count += 1
    if len(barcode) == 4 and count < 3:
        return 1
    elif len(barcode) > 5 and count < 4:
        return 1
    else:
        return 0

def has_res_site(barcode, opts):
    """Checks if the barcode recreates the restriction site in combination with the bc-adapter"""
    if opts.enzyme.overhang() == "3' overhang":
        #In case of for example PstI bc = Adapter + bc + TGCA.
        seq_pos = Seq(opts.ad_seq + barcode + opts.enzyme.ovhgseq,IUPACAmbiguousDNA())
    else:
        seq_pos = Seq(opts.ad_seq + barcode,IUPACAmbiguousDNA())
    if opts.enzyme.search(seq_pos):
        return 1
    else:
        return 0

def del_dups(items):
    found = set([])
    keep = []

    for item in items:
        if item not in found:
            found.add(item)
            keep.append(item)
    return keep

def exclude_recr(to_test, opts):
    """Exclude nt's which can not be present due 
    to recreating the restriction site from being tested by modifying to_test[-1]"""
    if opts.enzyme.overhang() == "5' overhang":
        excluded_nts = opts.enzyme.elucidate().split('^')[0][::-1]
    elif opts.enzyme.overhang() == "3' overhang":
        excluded_nts = opts.enzyme.elucidate().split('_')[0][::-1]
    for number, nt in enumerate(excluded_nts):
        if nt in to_test[-number-1]:
            #Remove nt's which would recreate a restriction site from barcode sequence.
            to_test[-number-1].pop(to_test[-number-1].index(nt))
    return to_test

def validate_bc(*args):
    """ define the requirement here"""
    comb, dist, existing, opts = args
    to_test = []
    for pos,num in enumerate(comb):
        value = find_keys(dist[pos], num)
        if value:
            to_test += [value]
        else:
            raise Exception('Invalid comb, %s not found in dictionary'%num)
    to_test = exclude_recr(to_test, opts)
    bcs = []
    for bc in product(*to_test):
        bc = ''.join(bc)
        #Start with least expensive functions: Does it have the res_site or is it low_complexity?
        if has_res_site(bc, opts) :
            continue
        if is_low_complex(bc):
           continue
        if has_trint(bc):
            continue
        if is_nested(bc, existing):
            continue
        if is_dissim(bc, existing):
            continue
        else:
            return bc
    return 0

def find_keys(dic, val):
    """return the key of dictionary dic given the value"""
    return [k for k, v in dic.iteritems() if v == val]

def get_nt_dist(input):
    """Yields a distribution of nt's per position for barcodes
    and a key_list indicating the pct of occurence of some nt's"""
    dist = {}
    for barcode in input:
        for pos, nt in enumerate(barcode):
            if not pos in dist:
                dist[pos] = {}
            try:
                dist[pos][nt] += 1
            except KeyError:
                dist[pos][nt] = 1
    for key, nt_dict in dist.items():
        for nt in nts:
            try:
                nt_dict[nt] = nt_dict[nt] / float(len(input))
            except KeyError:
                nt_dict[nt] = 0
    pct_list = []
    for key, value in dist.items():
        keys = []
        for pct in value.values():
            if pct not in keys:
                keys+= [pct]
        pct_list += [keys]
    return dist, pct_list

def check_compositon(nt_list, bc_length):
    """removes non-allowed nt's from beginning and end of preferred nt-list according to the 
    Restriction enzyme recognition site which should not be recreated"""
    #TODO: make function
    return nt_list

def get_stats(dist_simple, bc, opts):
    """Returns updated statistics nt-composition
        based on barcode composition and 'compulsory' nucleotides
        of restriction-site"""
    for pos, nt in enumerate(bc):
        try:
            dist_simple[pos][nt]+=1
        except KeyError:
            dist_simple[pos][nt]=1
    #import enzyme-imposed nt's
    comp_nts = ''.join(re.split(r'[\^,_]',opts.enzyme.elucidate())[1:])
    #TODO: add nucleotides with which sequence starts.
    pos+=1
    for pos_1, nt in enumerate(comp_nts):
        if not pos+pos_1 in dist_simple.keys():
            dist_simple[pos+pos_1] = {}
        try:
            dist_simple[pos+pos_1][nt]+=1
        except KeyError:
            #function for adding fractions to stats of different nt's addressed by ambiguity
            if nt in unambiguate:
                amb_nt = nt
                amb_level = len(unambiguate[amb_nt])
                for nt in unambiguate[amb_nt]:
                    try:
                        dist_simple[pos+pos_1][nt] += 1.0/amb_level
                    except KeyError:
                        dist_simple[pos+pos_1][nt] = 1.0/amb_level
            else:
                dist_simple[pos+pos_1][nt] = 1
    for i in range(len(bc)+len(comp_nts), opts.max + len(comp_nts)):
        if not i in dist_simple.keys():
            dist_simple[i] = {}
        for nt in ['A', 'C', 'G', 'T']:
            try:
                dist_simple[i][nt] += 0.25
            except KeyError:
                dist_simple[i][nt] = 0.25
    dist_list = []
    for dist in dist_simple.values():
        dist_item = []
        for nt, number in dist.items():
            if number in dist_item:
                continue
            else:
                dist_item += [number]
        dist_list += [dist_item]
    return dist_simple, dist_list
    

def gen_barcodes(opts):
    """Generates new barcodes based on criteria"""
    dist_simple = {}
#    barcodes = [ 'ACTA']
    barcodes = []
    extended_ovh = ''.join(re.split(r'[\^,_]',opts.enzyme.elucidate())[1:])
    for pos in range(opts.max + len(extended_ovh)):
        dist_simple[pos] = {'A':0, 'C':0, 'G':0, 'T':0}
    dist_list = []
    for dist in dist_simple.values():
        dist_item = []
        for nt, number in dist.items():
            if number in dist_item:
                continue
            else:
                dist_item += [number]
        dist_list += [dist_item]
#    dist_simple, dist_list = get_stats(dist_simple,barcodes[0], opts)
    #1. Determine how many bc's we want per size-class
    #2. Keep a running counter how the bc's which can not be found get selected from other size-classes
    # bc_size_class = {4:[15],5:[20],6[35]}
    #Let's state we want 150 barcodes. 
    required = opts.number
    batch_size = required / (opts.max+1- opts.min) 
#    print '%s\t'*len(list(barcodes[0]))%(tuple(barcodes[0]))
    for size in range(opts.min, opts.max +1):
        bc_count = 0
        if len(barcodes) == required:
            break
        while True:
            old_len = len(barcodes)
            # nice option, not used.. bc_list = [sorted(v.keys(),key=v.__getitem__) for k,v in dist_simple.iteritems()]
            for index,  comb in enumerate(sorted(product(*dist_list[:size]), key=sum)):
                bc = validate_bc(comb, dist_simple, barcodes, opts)
                if not bc:
                    continue
                else:
                    bc_count += 1
                    dist_simple, dist_list = get_stats(dist_simple,bc, opts)
                    lst = [(v.values()) for k,v in dist_simple.iteritems()]
#                    print '%s\t'*len(list(bc+extended_ovh))%(tuple(bc+extended_ovh))
                    barcodes += [''.join(bc)]
                    break #!!! Otherwise next comb becomes a mess!
            if size == opts.max +1 and len(barcodes) < required:
                continue
            elif size == opts.max  and len(barcodes)  == required:
                bc_count = batch_size
            if size == opts.max and len(barcodes) < required:
                continue
            if bc_count == batch_size:
                #enough barcodes for size class are retained break main loop
                break
            if old_len == len(barcodes):
                #Batch size not fulfilled, going to break. Calculate new required batch size.
                batch_size = (required - len(barcodes)) / (opts.max +1 - (size +1))
                break
    return barcodes, dist_simple

if __name__ == "__main__":
    main()
