# pylint: disable=print-statement
from collections import Counter, namedtuple
###############################################################################
# Basic class to parse sequence data, and optionally hold comment data, which is useful for some
# filtering criteria in the pipeline
#
class Sequence:
    
    def __init__(self, sequenceName, comment, sequence):
        self.sequenceName = sequenceName
        self.comment = comment
        self.sequence = sequence

    def printSelf(self):
        print self.sequenceName + ' [' + self.comment + ']'
        print self.sequence

    def ReturnSlice(self, start, stop):
        return self.sequence[start:stop] if stop > start else self.sequence[stop:start]

    @staticmethod
    def CreateSequence(sequenceName, sequence, comment = None):
        if comment:
            return Sequence(sequenceName, comment, sequence)
        else:
            return Sequence(sequenceName, '', sequence)

    @staticmethod
    def CalculateGC(geneSequence):
        nt = Counter(geneSequence)
        gc = (float(nt['G']) + nt['C']) / sum(nt.values()) * 100
        return round(gc, 2)

###############################################################################
# The atomic unit of nucleotide annotation. Comes with a few 'constructors' to auto-fill for various annotations
#
class NucleotideAnnotation:

    def __init__(self, target, sequence, annotationType, contigID, startPos, stopPos, orientation = None):
        self.target = target
        self.sequence = sequence
        self.gc = Sequence.CalculateGC(sequence)
        self.annotationType = annotationType
        self.contigID = contigID
        self.startPos = int(startPos)
        self.stopPos = int(stopPos)
        if orientation:
            self.orientation = orientation
        else:
            self.orientation = 'Forward' if (self.startPos < self.stopPos) else 'Reverse'

    def printSelf(self):
        print '---'
        print self.target
        print 'GC: %f' % self.gc
        print 'Contig: %s, start = %i, stop = %i, orientation = %s' % (self.contigID, self.startPos, self.stopPos, self.orientation)
        print self.sequence
        print self.annotationType
        print self.description

    @staticmethod
    def CreateRibosomeAnnotation(target, sequence, contigName):
        return NucleotideAnnotation(target, sequence, 'Predicted rRNA subunit', contigName, -1, -1, 'N/A')

    @staticmethod
    def CreateTRnaAnnotation(target, sequence, contigName, coordStruct):
        if 'c' in coordStruct: # Just account for complementary orientation within my own data structure
            coordStruct = coordStruct.replace('c[', '').replace(']', '')
            stop, start = coordStruct.split(',')[::-1]
        else:
            coordStruct = coordStruct.replace('[', '').replace(']', '')
            start, stop = coordStruct.split(',')
        return NucleotideAnnotation(target, sequence.upper(), 'tRNA', contigName, start, stop)

    @staticmethod
    def CreateRfamAnnotation(target, sequence, contigName, startPos, stopPos):
        return NucleotideAnnotation(target, sequence, 'Non-coding RNA (Rfam)', contigName, startPos, stopPos)

    @staticmethod
    def CreateCrisprAnnotation(target, sequence, contigName, startPos, stopPos):
        return NucleotideAnnotation(target, sequence, 'CRISPR spacer', contigName, startPos, stopPos)

###############################################################################
# The units of protein annotation. Comes with a few static methods for various validations.
#
class ProteinAnnotation:

    def __init__(self, fnaSequence, proteinSequence):
        self.geneSequence = fnaSequence.sequence
        self.gc = Sequence.CalculateGC(fnaSequence.sequence)
        self.proteinSequence = proteinSequence.sequence
        self.signal = None
        self.signalConf = None

        # Deconstruct the prodigal header for this information in nucleotide space
        #start, stop = fnaSequence.comment.split('#')[1:3]
        start = 0
        stop = 0
        contigID = fnaSequence.sequenceName.split('_')[0:-1]

        self.contigID = '_'.join(contigID)
        self.startPos = int(start)
        self.stopPos = int(stop)
        self.orientation = 'Forward' if (self.startPos < self.stopPos) else 'Reverse'
        self.annotations = []

    def AddAnnotationRecord(self, annot):
        self.annotations.append(annot)

    def HasAnnotations(self):
        return len(self.annotations) > 0

    def ReturnAnnotationsOfType(self, annotType):
        return [ annot for annot in self.annotations if annot.annotationType == annotType ]

    def printSelf(self):
        print '---'
        print self.geneSequence
        print self.proteinSequence
        print 'GC: %f' % self.gc
        print 'Contig: %s, start = %i, stop = %i, orientation = %s' % (self.contigID, self.startPos, self.stopPos, self.orientation)
        print 'Has a signalling domain annotation of %s (%s)' % ( self.signal,  self.signalConf)
        print 'Currently holding %i annotation results.' % len(self.annotations)

    @staticmethod
    def ParseAnnotationString(pAnnotation, annotType, topOnly):

        annotSummary = namedtuple('annotSummary', 'target identity coverage evalue description taxonomy')
        annots = pAnnotation.ReturnAnnotationsOfType(annotType)

        if annots:

            if topOnly:
                return annotSummary(target=annots[0].target,
                                    identity=str(annots[0].identity),
                                    coverage=str(annots[0].coverage),
                                    evalue=str(annots[0].evalue),
                                    description=annots[0].description,
                                    taxonomy=annots[0].taxonomy
                                    )
            else:
                # target identity coverage evalue description taxonomy
                vals = ['', '', '', '', '', '']

                ''' Instantiate on the first value, then roll through the remainder '''
                for annot in annots:
                    vals[0] += '; ' + annot.target
                    vals[1] += '; %f' % annot.identity
                    vals[2] += '; %f' % annot.coverage
                    vals[3] += '; %f' % annot.evalue
                    vals[4] += '; ' + annot.description
                    vals[5] += '; ' + annot.taxonomy

                # Strip off the leading '; '
                vals = [v[2:] for v in vals]
                return annotSummary( *vals )
        else:
            return None

class AnnotationResult:

    def __init__(self, target, identity, coverage, evalue, description, annotationType, taxonomy = None):
        self.target = target
        self.identity = float(identity)
        self.coverage = float(coverage)
        self.evalue = float(evalue)
        self.description = description
        self.annotationType = annotationType
        self.taxonomy = taxonomy if taxonomy else '-'

    def printSelf(self):
        print '---'
        print '%s (Ident = %f, coverage = %f)' % (self.target, self.identity, self.coverage)
        print self.annotationType
        print self.description
        print self.taxonomy

    @staticmethod
    def IsAccepted(cov, ident, covLimit, identLimit):
        return (cov >= covLimit) and (ident >= identLimit)

    @staticmethod
    def ComputeCoverage(sstart, send, tlength):
        c = float(send) - float(sstart)
        return c / float(tlength) * 100

###############################################################################
