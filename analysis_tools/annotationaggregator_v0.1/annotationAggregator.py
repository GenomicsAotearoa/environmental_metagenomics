#!/usr/bin/python
'''
    A script to automatically round up the outputs of a number of annotation programs, and present a single tab-delimited output file.
    Written in python2.7 to avoid python3 module loads for inexperienced members.

    David Waite 14/06/2018
'''

import sys
from collections import namedtuple
from optparse import OptionParser
from outputFactory import NucleotideOutput, ProteinOutput
from bioStructures import ProteinAnnotation
from bioFunctions import VerifyArguments, VerifyBlastParameters, PrepareFastaContents, VerifyOption, SpawnLookupStructures

# Nucleotide/amino acid parsers
from annotationParsers import ProcessMetaxa, ProcessAragorn, ProcessNcRnaAragorn, ProcessCrisprdigger
from annotationParsers import ProcessBlastOutput, ProcessHFamOutput, ProcessSignalpOutput

def main():
    
    # Parse options
    usageString = "usage: %prog [options] [bin or genome file] [prodigal nucleotide predictions] [prodigal protein predictions] [output file]"
    parser = OptionParser(usage=usageString)

    # Protein options, require Prodigal file
    parser.add_option('-r', '--uniref', help='UniRef100 annotation file', dest='uniref')
    parser.add_option('-u', '--uniprot', help='UniProt annotation file', dest='uniprot')
    parser.add_option('-k', '--kegg', help='KEGG annotation file', dest='kegg')
    #parser.add_option('-i', '--interproscan', help='InterProScan annotation file', dest='interproscan')
    parser.add_option('-p', '--pfam', help='Pfam tblout annotation file', dest='pfam')
    parser.add_option('-i', '--tigrfam', help='TIGRfam tblout annotation file', dest='tigrfam')
    parser.add_option('-s', '--signalp', help='SignalP prediction file', dest='signalp')

    # DNA options, sequences exported by software
    parser.add_option('-t', '--aragorn', help='Aragorn tRNA prediction file', dest='aragorn')
    parser.add_option('-m', '--metaxa', help='MeTaxa2 16S/23S prediction file', dest='metaxa')
    parser.add_option('-n', '--ncrna', help='Infernal ncRNA prediction file', dest='ncrna')
    parser.add_option('-c', '--crispr', help='CRISPRdigger prediction file (*.dr file)', dest='crispr')

    # Formating options and flags
    parser.add_option('--restrict_metaxa', help='Use only the MeTaxa2 output file provided. Default is to infer all possible output files and test',
                      dest='restrict_metaxa', action='store_true')
    parser.add_option('--toponly', help='Keep only the top hit for each protein (per method) (default: False)',
                      dest='toponly', action='store_true')
    parser.add_option('--ident', help='Minimum sequence identity to accept a BLAST hit (default 30%)',
                      dest='ident', default=30.0)
    parser.add_option('--coverage', help='Minimum query coverage to accept a BLAST hit (default 50%)',
                      dest='coverage', default=50.0)

    # Validate the arguments and optional flags
    options, arguments = parser.parse_args()
    genomeFile, prodigalFileFna, prodigalFileFaa, outputFile = VerifyArguments(arguments)
    options.ident, options.coverage = VerifyBlastParameters(options.ident, options.coverage)
    
    # Default paths for database lookup tables. Instantiate a namedtuple for the approriate values and check they exist.
    uniprotStruct, unirefStruct, keggStruct = SpawnLookupStructures()

    if not VerifyOption(uniprotStruct.path) or not VerifyOption(uniprotStruct.tax):
        print 'Please contact David Waite about the missing UniProt lookup file.'
        options.uniprot = None

    if not VerifyOption(unirefStruct.path) or not VerifyOption(unirefStruct.tax):
        print 'Please contact David Waite about the missing UniRef100 lookup file.'
        options.uniref = None

    if not VerifyOption(keggStruct.path) or not VerifyOption(keggStruct.tax):
        print 'Please contact David Waite about the missing KEGG lookup file.'
        options.kegg = None

    '''
        Nucleotide results consist of a list of Annotation objects.
        Protein results map to a list that incorporates all options.
        
        This is done because nucleotide annotations generate their own sequences as output, so they are unique between
        programs. Since the protein annotation is deliberately a reannotation of the same sequences over and over, there
        will be competing annotations between the outputs.
    '''
    global proteinResults, nucleotideResults
    prodigalPredictionsFna = PrepareFastaContents(prodigalFileFna)
    prodigalPredictionsFaa = PrepareFastaContents(prodigalFileFaa)

    # This used to be a dict comprehension, but was too messy
    proteinResults = {}
    for geneName, aaSeq in prodigalPredictionsFaa.iteritems():
        prodFnaRecord = prodigalPredictionsFna[geneName]
        proteinResults[geneName] = ProteinAnnotation(prodFnaRecord, aaSeq)

    nucleotideResults = []

    ''' Parse nucleotide results '''
    if options.aragorn:
        if VerifyOption(options.aragorn):
            ProcessAragorn(options.aragorn, genomeFile, nucleotideResults)

    if options.ncrna:
        if VerifyOption(options.ncrna):
            ProcessNcRnaAragorn(options.ncrna, genomeFile, nucleotideResults)

    if options.metaxa:
        if VerifyOption(options.metaxa):
            ProcessMetaxa(options.metaxa, options.restrict_metaxa, nucleotideResults)

    if options.crispr:
        if VerifyOption(options.crispr):
            ProcessCrisprdigger(options.crispr, nucleotideResults)
    
    ''' Parse protein results '''
    if options.uniref:
        if VerifyOption(options.uniref):
            ProcessBlastOutput(options.uniref, proteinResults, options.ident, options.coverage, unirefStruct)

    if options.uniprot:
        if VerifyOption(options.uniprot):
            ProcessBlastOutput(options.uniprot, proteinResults, options.ident, options.coverage, uniprotStruct)

    if options.kegg:
        if VerifyOption(options.kegg):
            ProcessBlastOutput(options.kegg, proteinResults, options.ident, options.coverage, keggStruct)

    ''' Removed, replaced with just flat Pfam/TIGRfam annotations
    if options.interproscan:
        if VerifyOption(options.interproscan):
            ProcessInterproscanOutput(options.interproscan, proteinResults)
    '''

    if options.pfam:
        if VerifyOption(options.pfam):
            ProcessHFamOutput(options.pfam, proteinResults, 'Pfam')

    if options.tigrfam:
        if VerifyOption(options.tigrfam):
            ProcessHFamOutput(options.tigrfam, proteinResults, 'TIGRfam')

    if options.signalp:
        if VerifyOption(options.signalp):
            ProcessSignalpOutput(options.signalp, proteinResults)

    ''' Write out the final file '''
    outputWriterNt = NucleotideOutput(outputFile + '.nt')
    outputWriterAA = ProteinOutput(outputFile + '.aa')

    for ntResult in nucleotideResults:
        outputWriterNt.WriteAnnotation(ntResult)

    for geneID, protAnnotation in proteinResults.iteritems():
        outputWriterAA.WriteAnnotation(geneID, protAnnotation, options.toponly)
    
    outputWriterNt.CloseStream()
    outputWriterAA.CloseStream()
###############################################################################
if __name__ == '__main__':
    main()
