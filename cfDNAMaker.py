#!/usr/bin/env python
'''
CFDNAMarker
Author: XiaolongZhang
modify: 2017-03-01
Inputs:
    A position-sorted paired-end BAM file containing reads with a duplex tag in the header.
Outputs:
    1. output bam file(After filter the duplicate reads)
    2. SSCS.tagstate file

    Special arguments description:
    [-p] [-read_type] : A string specifying which types of read to consider.
                        Read types as follows [default: 'dpm']: 
                        n: Neither read 1 or read 2 mapped. 
                        m: Either read 1 or read 2 mapped, but not both. 
                        p: Both read 1 and read 2 mapped, not a propper pair. 
                        d: Both read 1 and read 2 mapped, propper pair. 
                        s: Single ended reads.
    [-r] [-rep_filt] :  Remove tags with homomeric runs of nucleotides of length x. [default: 6]
'''
import sys
import re
import os.path
import pysam
from collections import defaultdict
from optparse import OptionParser

#global regular to match the tag (eg: "|CCAATGCCGCTGATGCGCAAACTA/") in seqname
TagMatch = re.compile('\|[ATCG]+/')


class MatchRead (object):
    ''' Record the base and quality marked by 'M' in cigar value!
        totalqual: calculate the whole quality of the matched base quality!
    '''
    def __init__(self, seq, qual):
        self.seq = seq      # match seq 
        self.qual = qual    # match qual

    def totalqual(self):
        Qual = 0
        for q in self.qual:
            Qual += ord(q)
        return Qual


def TagGet (read, repcount):
    try:
        tag = TagMatch.search(read.qname).group()[1:-1]
        matchnum = 0
        for cig in read.cigar:
            if cig[0] == 0:
                matchnum += cig[1]
        return (tag, matchnum)
    except:
        return ('A' * repcount, repcount)


def GetMatchBase (readlist):
    outlist = []
    for read in readlist:
        seq, qual, index = read.seq, read.qual, 0
        trimseq, trimqual = '', ''

        for cig in read.cigar:
            if cig[0] == 0: # M
                trimseq += seq[index:index+cig[1]]
                trimqual += qual[index:index+cig[1]]
                index += cig[1]
            elif cig[0] == 4 or cig[0] == 1: # S, I
                index += cig[1]
            
        outlist.append(MatchRead(trimseq, trimqual))
    return outlist


def GoodFlagDefine (flagtype):
    goodFlag = []
    if 'd' in flagtype:
        goodFlag.extend((99, 83, 163, 147))
    if 'm' in flagtype:
        goodFlag.extend((181, 117, 137, 133, 73, 89, 69, 153))
    if 'p' in flagtype:
        goodFlag.extend((97, 81, 161, 145, 129, 65, 177, 113))
    if 'n' in flagtype:
        goodFlag.extend((141, 77, 4))
    if 's' in flagtype:
        goodFlag.extend((0, 16))
    return goodFlag


def TagFlagLenCheck (read, flagtype, repfilt, tag, mlength):
    ''' tag = ('ACTTCG', 103)
        'ACTTCG': barcode tag
        103: the matched length of alignment
    '''
    checkstatus = True
    badTag = ('A'*repfilt, 'T'*repfilt, 'G'*repfilt, 'C'*repfilt)
    #pdb.set_trace()
    goodFlag = GoodFlagDefine(flagtype)

    if tag[1] < mlength:
        checkstatus = False

    if (read.flag not in goodFlag) or (tag[0] in badTag):
        checkstatus = False

    return checkstatus


def CreatSSCSRead (readlist, seqlen):
    readnum, matchlist = len(readlist), GetMatchBase(readlist)
    if readnum == 2:
        seq1, seq2 = matchlist[0].seq, matchlist[1].seq
        for i in xrange(seqlen):
            if seq1[i] != seq2[i]:
                if matchlist[0].totalqual() > matchlist[1].totalqual():
                    return readlist[0]
                else:
                    return readlist[1]
        return readlist[0]
    else:
        consensus, outread = '', readlist[0]
        for i in xrange(seqlen):
            nucIdentity = {'A':0,'T':0,'G':0,'C':0,'N':0}

            for j in xrange(readnum):
                nucIdentity[matchlist[j].seq[i]] += 1                
            maxBase = max(nucIdentity, key=nucIdentity.get)

            if nucIdentity[maxBase] / float(readnum) > 0.6:
                consensus += maxBase
            else:
                nucQuality = {'A':0,'T':0,'G':0,'C':0,'N':0}
                for j in xrange(readnum):
                    nucQuality[matchlist[j].seq[i]] += ord(matchlist[j].qual[i])
                consensus += max(nucQuality, key=nucQuality.get)

        outread.seq, outread.mapq = consensus, 255
        outread.cigar = [(0, seqlen)]
        outread.qual = matchlist[0].qual
        return outread

def OverSeqDeal (overlist, overDict, repcount):
    for read in overlist:
        tag = TagGet(read, repcount)
        if tag in overDict:
            overDict[tag].append(read)
        else:
            overDict[tag] = [read]

        
def ConsensusMaker (consensuslist, readDict, readcount):
    for tagkey in readDict.keys():
        if len(readDict[tagkey]) == 1:
            consensuslist.append(readDict[tagkey][0])
            readcount['singlefamilynum'] += 1

        elif len(readDict[tagkey]) >= 2:
            sscsread = CreatSSCSRead(readDict[tagkey], tagkey[1])
            consensuslist.append(sscsread)
            readcount['sscsnum'] += 1


def TagInfoWrite (familycount, tagfile):
    totalreads = 0

    with open (tagfile, 'w') as tagstatfile:
        for size in familycount:
            familycount[size] *= int(size)
            totalreads += int(familycount[size])

        for size in sorted(familycount.keys()):
            tagstatfile.write("%s\t%s\n" %(size, float(familycount[size])/float(totalreads)))


def main():
    usage = 'Usage: %prog -i <input.bam> -o <out.bam> -m <minimum matchlen>'
    parser = OptionParser(usage=usage)
    parser.add_option('-i', "--infile", action="store", dest="infile", 
            help="[required] input BAM file")
    parser.add_option('-o', "--outfile",  action="store", dest="outfile", 
            help="[required] output BAM file")
    parser.add_option('-m', "--matchlen",  action="store", type=int, dest="mlength",
            help="[required] minimum match length for alignment")
    parser.add_option('-t', "--tagfile",  action="store",  dest="tagfile", 
            help="output tagstat file",  default='SSCS.tagstat')
    parser.add_option('-r', "--rep_filt", action="store",  type=int, dest='rep_filt', 
            help="Remove tags with sam nucleotides of length x. [6]", default=6 )
    parser.add_option('-p', '--read_type', type=str, action="store", dest='read_type', default="dpm", 
            help="Read types: n, m, p, d, s. default = ['dpm']")

    (options, args) = parser.parse_args()
    if len(sys.argv) < 7:
        parser.print_help()
        sys.exit(0)

    InBam = pysam.Samfile(options.infile, "rb")
    OutBam = pysam.Samfile(options.outfile,"wb", template = InBam)
    BamEntry = InBam.fetch(until_eof = True)

    previousread, currentread = BamEntry.next(), ''
    readcount = {'allreadnum':1, 'singlefamilynum':0, 'sscsnum':0}
    readDict, overDict, consensuslist, overlist, seqname = {}, {}, [], [], []
    tagDict, familycount = defaultdict( lambda: 0 ), defaultdict( lambda: 0 )
    ''' DataStruct (examples as fellows)
        readDict = {(tag1,len1):[read1,read2,...], (tag2,len2):[read1,read2,...], ...}
        overDict = {(tag1,len1):[read1,read2,...], (tag2,len2):[read1,read2,...], ...}
        consensuslist = [read1,read2,read3 ...]
        tagDict = {tag1:15, tag2: 187, tag3:18, ...}
        familycount = {1: 23454, 2: 4543, 3: 4443 ...}
    '''
    readDict[TagGet(previousread, options.rep_filt)] = [previousread]
    seqname.append(previousread.qname)
    filedone = False

    while (filedone == False):
        try:
            currentread = BamEntry.next()
            readcount['allreadnum'] += 1
        except:
            filedone = True
            ConsensusMaker(consensuslist, readDict, readcount)
            OverSeqDeal(overlist, overDict, options.rep_filt)
            ConsensusMaker(consensuslist, overDict, readcount)

            for sscsread in consensuslist:
                OutBam.write(sscsread)

            for tagvalue in tagDict.values():
                familycount[tagvalue] += 1

            TagInfoWrite (familycount, options.tagfile)
            continue

        if (currentread.pos == previousread.pos):
            if readcount['allreadnum'] % 100000 == 0:
                sys.stdout.write('\r[*]  Reads processed: %d' %(readcount['allreadnum']))
                sys.stdout.flush()

            tag = TagGet(currentread, options.rep_filt)
            if TagFlagLenCheck(currentread, options.read_type, \
                                    options.rep_filt, tag, options.mlength):
                if currentread.qname in seqname:
                    overlist.append(currentread)
                else:
                    seqname.append(currentread.qname)
                    tagDict[tag[0]] += 1

                    if tag not in readDict:
                        readDict[tag] = [currentread]
                    else:
                        readDict[tag].append(currentread)

            previousread = currentread
        else:
            ConsensusMaker(consensuslist, readDict, readcount)
            OverSeqDeal(overlist, overDict, options.rep_filt)
            ConsensusMaker(consensuslist, overDict, readcount)

            for sscsread in consensuslist:
                OutBam.write(sscsread)
            for tagvalue in tagDict.values():
                familycount[tagvalue] += 1

            readDict, overDict, tagDict = {}, {}, defaultdict(lambda:0)
            consensuslist, overlist, seqname = [], [], [currentread.qname]
            
            tag = TagGet(currentread, options.rep_filt)
            tagDict[tag[0]] += 1

            if TagFlagLenCheck(currentread, options.read_type, \
                                    options.rep_filt, tag, options.mlength):
                readDict[tag] = [currentread]
            previousread = currentread


    # Close the files opened
    InBam.close()
    OutBam.close()

    # Write summary statistics
    sys.stdout.write("\n")
    sys.stdout.write("Summary Statistics: \n")
    sys.stdout.write("Total reads processed : %d\n" %readcount['allreadnum'])
    sys.stdout.write("Reads after filtered : %d\n" %(readcount['singlefamilynum'] + readcount['sscsnum']))
    sys.stdout.write("\t(1)Reads of siglefamily : %d\n" %readcount['singlefamilynum'])
    sys.stdout.write("\t(2)Reads os SSCS : %d\n" %readcount['sscsnum'])


if __name__ == '__main__':
    main()
