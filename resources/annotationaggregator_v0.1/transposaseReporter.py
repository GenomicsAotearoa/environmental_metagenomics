'''
python transposaseReporter.py -k Hydrogenimonas/Hydrogenimonas.GCF_900115615_1.TPase.kamoun.txt -o test Hydrogenimonas/Hydrogenimonas.GCF_900115615_1.prod.fna Hydrogenimonas/Hydrogenimonas.GCF_900115615_1.prod.faa
'''
import sys, os
from optparse import OptionParser
from collections import namedtuple
from bioFunctions import VerifyOption, PrepareFastaContents
from bioStructures import ProteinAnnotation
from annotationParsers import ParseHmmTable
from outputFactory import OutputFactory

def main():

    parser = OptionParser()
    usage = "usage: %prog [options] (Prodigal fna file) (Prodigal faa file)"

    # Annotation options
    parser.add_option('-p', '--pfam', help='Pfam annotation file', dest='pfam')
    parser.add_option('-t', '--tigrfam', help='TIGRfam annotation file', dest='tigrfam')
    parser.add_option('-k', '--kamoun', help='Tranposase annotation file base on the Kamoun et al. (2013) set', dest='kamoun')

    parser.add_option('-e', '--evalue', help='Maximum e-value to accept (Default: 0.001)', dest='evalue')
    parser.add_option('-o', '--output', help='Output file name (Default: Input + .transposases', dest='output')

    # Validate the input file and optional flags
    options, arguments = parser.parse_args()
    
    prodigalFileFna, prodigalFileFaa = arguments
    if not (prodigalFileFna and prodigalFileFaa):
        print 'Please provide valud input files...'
        sys.exit()

    prodigalPredictionsFna = PrepareFastaContents(prodigalFileFna)
    prodigalPredictionsFaa = PrepareFastaContents(prodigalFileFaa)

    eValue = 1e-3
    if options.evalue:
        try:
            eValue = float(options.evalue)
        except:
            print 'Unable to parse e-value, using default...'

    # Declare the namedtuple structure and filter lists
    hmmAnnot = namedtuple('hmmAnnot', 'target accession evalue')
    pfamList = ['PF17490', 'PF17489', 'PF14706', 'PF12017', 'PF12596', 'PF01359', 'PF02371', 'PF02992', 'PF02994', 'PF03017', 'PF03004',
                'PF04195', 'PF04754', 'PF00872', 'PF13963', 'PF04236']
    tigrfamList = ['TIGR01765', 'TIGR01766']

    # For each input, read in the full results
    global annotationResults
    annotationResults = {}

    if options.pfam:
        if VerifyOption(options.pfam):
            FilterHits(options.pfam, eValue, hmmAnnot, set(pfamList))

    if options.tigrfam:
        if VerifyOption(options.tigrfam):
            FilterHits(options.tigrfam, eValue, hmmAnnot, set(tigrfamList))

    if options.kamoun: # Make sure to update namedtuple structure (blank acc in output)
        if VerifyOption(options.kamoun):
            FilterHits(options.kamoun, eValue, hmmAnnot)

    ''' Write out the data '''
    outWriter = options.output if options.output else prodigalFileFaa + '.transposases'
    outWriter = open(outWriter, 'w')
    outWriter.write('Query Gene\tGC%\tContig name\tStart position\tStop position\tOrientation\tQuery sequence\tTarget gene\tE-value\tDescription\n')

    for queryGene, annot in annotationResults.iteritems():
        # Parse into useful data based on the ProteinAnnotation structure
        gc, start, stop, aaSeq, contigId, orientation = ParseToProteinAnnotation(queryGene, prodigalPredictionsFaa[queryGene], prodigalPredictionsFna[queryGene])
        outWriter.write( '{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(queryGene, gc, contigId, start, stop, orientation, aaSeq) )
        outWriter.write( '{}\t{}\t{}\n'.format(annot.accession, annot.evalue, annot.target) )

    outWriter.close()

###############################################################################
def FilterHits(resultsFile, eThreshold, hmmAnnot, filterList = None):

    for (qseqid, name, acc, e) in ParseHmmTable(resultsFile): # query, geneName, accession, evalue
        if filterList:
            if float(e) <= eThreshold and name in filterList:
                tHit = hmmAnnot(target=name, accession=acc, evalue=e)
                UpdateAnnotResults(qseqid, tHit)
            
        else: # We're in the Kamoun data and need to declare hits differently
            if float(e) <= eThreshold:
                tHit = hmmAnnot(target='Kamoun transposase', accession=name, evalue=e)
                UpdateAnnotResults(qseqid, tHit)

def UpdateAnnotResults(qseqid, tHit):
    if not qseqid in annotationResults:
        annotationResults[qseqid] = tHit
    else:
        if float(tHit.evalue) < float(annotationResults[qseqid].evalue): annotationResults[qseqid] = tHit

def ParseToProteinAnnotation(queryGene, aaSeq, ntSeq):
    p = ProteinAnnotation(ntSeq, aaSeq)
    return ( p.gc, p.startPos, p.stopPos, p.proteinSequence, p.contigID, p.orientation )
###############################################################################
if __name__ == '__main__':
    main()