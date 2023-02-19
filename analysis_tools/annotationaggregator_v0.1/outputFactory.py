# pylint: disable=print-statement

from collections import namedtuple
from itertools import compress
from bioStructures import NucleotideAnnotation, ProteinAnnotation, AnnotationResult

###############################################################################
# Base class
#
class OutputFactory:
    ''' A class to abstract the details of writing out results '''
    #fileName = ''
    #writerHandle = None
    
    def __init__(self, fileName):
        self.fileName = fileName
        self.writerHandle = open(self.fileName, 'w')

    def CloseStream(self):
        self.writerHandle.close()
###############################################################################
# Nucleotide class
#
class NucleotideOutput(OutputFactory):

    def __init__(self, fileName):
        OutputFactory.__init__(self, fileName)
        self.writerHandle.write('Query Gene\tGC%\tContig name\tStart position\tStop position\tOrientation\tTarget\tAnnotation Method\tQuery sequence\n')

    def WriteAnnotation(self, annotation):
        self.writerHandle.write(annotation.contigID + '\t')
        self.writerHandle.write('%f\t' % annotation.gc)
        self.writerHandle.write(annotation.contigID + '\t')

        ''' Alternate flow depending on what information is available '''
        if annotation.target == 'Predicted rRNA subunit':
            self.writerHandle.write('-\t-\t-\t')
        else:
            self.writerHandle.write('%i\t' % annotation.startPos)
            self.writerHandle.write('%i\t' % annotation.stopPos)
            self.writerHandle.write(annotation.orientation + '\t')

        self.writerHandle.write(annotation.target + '\t')
        self.writerHandle.write(annotation.annotationType + '\t')
        self.writerHandle.write(annotation.sequence + '\n')
###############################################################################
# Protein class
#
class ProteinOutput(OutputFactory):

    def __init__(self, fileName):
        OutputFactory.__init__(self, fileName)

        ''' Write out the header information '''
        self.writerHandle.write('Query Gene\tGC%\tContig name\tStart position\tStop position\tOrientation\tQuery sequence')
        
        ''' Write SignalP column '''
        self.writerHandle.write('\tSignalling\tSignal confidence')

        ''' Write the BLAST columns '''
        for c in ['UniProt', 'UniRef100', 'KEGG']:
            self.writerHandle.write('\tTarget gene (%s)\tIdentity\tCoverage\tE-value\tDescription\tTaxonomy' % c)

        ''' Write the HMMer colmns '''
        for c in ['Pfam', 'TIGRfam']:
            self.writerHandle.write('\tTarget gene (%s)\tE-value\tDescription' % c)

        self.writerHandle.write('\n')

    def WriteAnnotation(self, geneID, protAnnotation, topOnly):

        # Write out the start, since this is guaranteed
        self.writerHandle.write(geneID + '\t')
        self.writerHandle.write('%f\t' % protAnnotation.gc)
        self.writerHandle.write(protAnnotation.contigID + '\t')
        self.writerHandle.write('%i\t' % protAnnotation.startPos)
        self.writerHandle.write('%i\t' % protAnnotation.stopPos)
        self.writerHandle.write(protAnnotation.orientation + '\t')
        self.writerHandle.write(protAnnotation.proteinSequence)

        # If present, record the SignalP data
        if protAnnotation.signal:
            self.writerHandle.write('\t' + protAnnotation.signal + '\t')
            self.writerHandle.write(protAnnotation.signalConf)
        else:
            self.writerHandle.write('\t-\t-')

        # Write the annotation information, or fill with N/A values as appropriate
        if protAnnotation.HasAnnotations():
            for c in ['UniProt', 'UniRef100', 'KEGG']:
                annotTuple = ProteinAnnotation.ParseAnnotationString(protAnnotation, c, topOnly)
                if annotTuple:
                    self.writerHandle.write('\t%s\t%s\t%s\t%s\t%s\t%s' % annotTuple)
                else:
                    self.WriteNullBlast()

            for c in ['Pfam', 'TIGRfam']:
                annotTuple = ProteinAnnotation.ParseAnnotationString(protAnnotation, c, topOnly)
                if annotTuple:
                    self.writerHandle.write('\t%s\t%s\t%s' % (annotTuple.target, annotTuple.evalue, annotTuple.description) )
                else:
                    self.WriteNullHmmer()

            self.writerHandle.write('\n')

        else: # Fill in the columns blank...
            self.WriteNullBlast()
            self.WriteNullBlast()
            self.WriteNullBlast()
            self.WriteNullHmmer()
            self.WriteNullHmmer()
            self.writerHandle.write('\n')

    def WriteNullBlast(self): # Annotation method, Identity, Coverage, E-value, Description, Taxonomy
        self.writerHandle.write('\tNot annotated\t-\t-\t-\t-\t-')

    def WriteNullHmmer(self):
        self.writerHandle.write('\tNot annotated\t-\t-') # Annotation method, E-value, Description
