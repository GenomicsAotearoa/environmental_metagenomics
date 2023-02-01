# pylint: disable=print-statement
from bioStructures import NucleotideAnnotation, Sequence, AnnotationResult
from bioFunctions import PrepareFastaContents
from collections import namedtuple, Counter
from itertools import compress
import sys, os

###############################################################################
# Functions to handle MeTaxa2 and CRISPRdigger output
#
def ProcessMetaxa(metaxa, singleFile, nucleotideResults):
    ''' MeTaxa2 does not include coordinate information, so leave it blank for now. Might be required to do additional BLAST step if
        this is required. '''
    commentLookup = {'archaea': 'Archaeal 16S rRNA sequence',
                     'bacteria': 'Bacterial 16S rRNA sequence',
                     'chloroplast': 'Chloroplast-like 16S rRNA sequence',
                     'eukaryota': 'Eukaryotic 16S rRNA sequence',
                     'mitochondria': 'Mitochondria-like 16S rRNA sequence'}

    if singleFile:
            fileType = metaxa.split('.')[-2]
            ParseMetaxaFile(metaxa, commentLookup[fileType], nucleotideResults)
    else:
        stub, ext = os.path.splitext(metaxa)
        stub = '.'.join(stub.split('.')[0:-1])
        for fType in ['archaea', 'bacteria', 'chloroplast', 'eukaryota', 'mitochondria']:
            fName = stub + '.' + fType + ext
            if os.path.isfile(fName):
                ParseMetaxaFile(fName, commentLookup[fType], nucleotideResults)

def ParseMetaxaFile(mtFile, comment, nucleotideResults):
    results = PrepareFastaContents(mtFile)
    if len(results) > 0:
        for k in results.keys():
            seq = results[k]
            s = NucleotideAnnotation.CreateRibosomeAnnotation(comment, seq.sequence, seq.sequenceName)
            nucleotideResults.append(s)

def ProcessCrisprdigger(crisprFile, nucleotideResults):
    content = open(crisprFile, 'r').readlines()[1:]

    # Dynamically build out a per-contig comment, parse to get the counts of all contigs with detected CRISPR spacers
    annotatedContigs = [ x.split('\t')[0] for x in content ]
    annotatedContigsTot = Counter(annotatedContigs)
    annotatedContigsObs = { aC: 1 for aC in annotatedContigsTot.keys() }

    for line in content:
        line = line.strip().split('\t')
        contig, start, stop = line[0:3]
        spacer = line[-1]
        targetNote = 'CRISPR spacer sequence ({} of {})'.format(annotatedContigsObs[contig], annotatedContigsTot[contig])
        annotatedContigsObs[contig] += 1
        nucleotideResults.append( NucleotideAnnotation.CreateCrisprAnnotation(targetNote, spacer, contig, start, stop) )

###############################################################################
# Functions to handle Aragorn and Infernal output
#
def ProcessAragorn(aragorn, genomeFile, nucleotideResults):
    ''' Aragorn only links to the reference sequences by input order, so need to parse that out in a new function '''
    results = PrepareFastaContents(aragorn)
    refGenomeContigs = PullGenomeContigOrder(genomeFile)
    if len(results) > 0:
        for k in results.keys():

            seq = results[k]
            contigName = ReturnContigName(refGenomeContigs, k)

            # Split apart the Aragorn metadata
            tRNA, coordStruct = ParseAragornMetadata(seq.comment)
            nucleotideResults.append( NucleotideAnnotation.CreateTRnaAnnotation(tRNA, seq.sequence, contigName, coordStruct) )

def ProcessNcRnaAragorn(rfam, genomeFile, nucleotideResults):
    ''' Parse the file and extract the coordinates '''
    rfamContent = open(rfam, 'r').readlines()[2:-10]
    refGenomeContigs = PrepareFastaContents(genomeFile)

    for line in rfamContent:
        line = [l for l in line.split(' ') if l ]
        target = '%s (%s)' % (line[0], line[1])
        contigName = line[2]
        start, stop = map(int, line[7:9])
        sequenceSlice = refGenomeContigs[contigName].ReturnSlice(start, stop)
        nucleotideResults.append( NucleotideAnnotation.CreateRfamAnnotation(target, sequenceSlice, contigName, start, stop) )

def PullGenomeContigOrder(genomeFile):
    genomeContents = open(genomeFile, 'r').readlines()
    contigs = [c[1:].strip() for c in genomeContents if '>' in c]
    return contigs

def ReturnContigName(refGenomeContigs, contig):
    cI = int(contig.split('-')[0])
    return refGenomeContigs[ cI-1 ]

def ParseAragornMetadata(seqComment):
    if '(Permuted)' in seqComment:
        tRNA, coordStruct = seqComment.split(')')
        tRNA += ')'
    else:
         tRNA, coordStruct = seqComment.split(' ')
    return tRNA, coordStruct
###############################################################################
# Functions to handle SignalP output
#
def ProcessSignalpOutput(signalpFile, proteinResults):

    signalpContent = open(signalpFile, 'r').readlines()[2:]
    for line in signalpContent:
        signalpAnnot = [l for l in line.strip().split(' ') if l]
        gene, score, threshold, call = list( compress(signalpAnnot, [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]) )

        if float(score) >= float(threshold):
            proteinResults[gene].signal = call
            proteinResults[gene].signalConf = score

###############################################################################
# Functions to handle BLAST output (UniRef100, UniProt, etc.)
#
def ProcessBlastOutput(blastFile, proteinResults, reqIdent, reqCoverage, dbData):

    tempAnnotations = []
    obsHits = set()
    
    ''' Read in the initial annotation data and store as a 7-tuple with the gene ID and siphon off a set of unique protein observations.
        At this stage, discard the bitscore because it doesn't make it to the output anyway. '''
    for line in open(blastFile, 'r'):
        qseqid, sseqid, sstart, send, pident, bitscore, evalue = line.strip().split('\t')
        tempAnnotations.append( (qseqid, sseqid, sstart, send, pident, evalue) )
        obsHits.add(sseqid)

    ''' Import the metadata of all observed hits, index the important information '''
    metadataLookup = ImportAnnotationMetadata(obsHits, dbData.path)
    #taxonomyLookup = ImportTaxonomyData_kegg(obsHits, dbData.tax) if dbData.name == 'KEGG' else ImportTaxonomyData(obsHits, dbData.tax)
    taxonomyLookup = ImportTaxonomyData(obsHits, dbData.tax)

    for (qseqid, sseqid, sstart, send, pident, evalue) in tempAnnotations:

        tLength, tDesc = metadataLookup[sseqid]
        qseqid = ClipQueryName(qseqid)
        coverage = AnnotationResult.ComputeCoverage(sstart, send, tLength)
        pident = float(pident)
        if AnnotationResult.IsAccepted(coverage, pident, reqCoverage, reqIdent):
            if sseqid in taxonomyLookup:
                annot = AnnotationResult(sseqid, pident, coverage, evalue, tDesc, dbData.name, taxonomyLookup[sseqid])
            else:
                annot = AnnotationResult(sseqid, pident, coverage, evalue, tDesc, dbData.name)
            proteinResults[qseqid].AddAnnotationRecord(annot)

def ImportAnnotationMetadata(obsHits, dbPath):
    lookup = {}
    for line in open(dbPath, 'r'):
        tSeq, tLength, tDesc = line.strip().split('\t')
        if tSeq in obsHits:
            lookup[tSeq] = (tLength, tDesc)

    return lookup

# KEGG uses a slightly different mapping in the taxonomy file, so much be handled differently
def ImportTaxonomyData_kegg(obsHits, dbPath):

    obs = {}
    for o in obsHits:
        taxId = o.split(':')[0]
        if not taxId in obs: obs[taxId] = []
        obs[taxId].append(o)

    lookup = {}
    for line in open(dbPath, 'r'):
        name, taxonomy = line.strip().split('\t')

        if name in obs:
            for hitName in obs[name]: lookup[hitName] = taxonomy

    return lookup

''' Need to resolve \n\n issue in uniprot and uniref '''
def ImportTaxonomyData(obsHits, dbPath):
    lookup = {}
    for line in open(dbPath, 'r'):
        try:
            name, taxonomy = line.strip().split('\t')

            if name in obsHits:
                lookup[name] = taxonomy
        except:
            print('Unable to parse line {}'.format( line.strip() ))

    return lookup

def ClipQueryName(qName): # Remove the trailing comment info from prodigal annotation if needed.
    return qName.split(' ')[0] if ' ' in qName else qName
###############################################################################
# Functions to handle InterProScan output (Pfam, TIGRfam)
# Updated with the HMMer parsers
def ProcessHFamOutput(pfamFile, proteinResults, hmmType):
    for (qseqid, description, accession, evalue) in ParseHmmTable(pfamFile):
        annot = AnnotationResult(accession, -1, -1, evalue, description, hmmType)
        proteinResults[qseqid].AddAnnotationRecord(annot)

def ParseHmmTable(hmmerTable):

    content = open(hmmerTable, 'r').readlines()
    results = []

    for line in content[3:-10]:

        ''' The deafult behaviour of hmmer should push the prodigal metadata to a trailing column,
            so don't need to worry about it here. '''
        values = [ l for l in line.split(' ') if l]
        query, geneName, accession, evalue = compress(values, [1, 0, 1, 1, 1])
        results.append( (query, geneName, accession, evalue) )

    return results

def ProcessInterproscanOutput(ipsFile, proteinResults):
    
    for line in open(ipsFile, 'r'):
        
        lineContent = ChompIpsLine(line)

        # Parse the results
        qseqid = lineContent['acc']

        ''' Store dud results for identity and coverage, since they don't have meaning here '''
        annot = AnnotationResult(lineContent['resultAccession'], -1, -1, lineContent['eValue'], lineContent['resultDescription'], lineContent['resultClass'])
        proteinResults[qseqid].AddAnnotationRecord(annot)

def ChompIpsLine(line):
    # Ugly output delimiter...
    ipsData = line.split('\t')

    # Create an object-representational dict to hold the line content.
    ipsObj = {
        'acc': ipsData[0],
        'resultClass': ipsData[3],
        'resultAccession': ipsData[4],
        'resultDescription': ipsData[5],
        'eValue': ipsData[8].split(' ')[0],
    }

    # Check for any KEGG matches
    k = SummariseKeggMatches(ipsData[10:])
    if k:
        ipsObj['resultDescription'] += ' (%s)' % k

    return ipsObj

def SummariseKeggMatches(annotText):
    keggAnnot = []
    for k in ExtractKeggTerms(annotText):
        keggAnnot.append(k)

    return ','.join(keggAnnot) if (len(keggAnnot) > 0) else None

def ExtractKeggTerms(content):
    for elem in content:
        if 'KEGG' in elem:
            for e in elem.split('KEGG: ')[1:]:
                yield 'K' + e.split('+')[0]
###############################################################################
