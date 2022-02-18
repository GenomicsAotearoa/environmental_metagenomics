# pylint: disable=print-statement
import sys, os
from collections import namedtuple
from bioStructures import Sequence

def VerifyOption(optFile):
    if os.path.isfile(optFile):
        return True
    else:
        print 'Warning: Cannot find file ' + optFile + ', skipping...'
        return False

def VerifyArguments(args):

    # Check the right number of inputs
    if len(args) != 4:
        print 'Incorrect number of arguments. Aborting...'
        sys.exit()

    # Do the files exist?
    genome, prodFna, prodFaa, out = args
    if not os.path.isfile(genome):
        print 'Genome file does not exist. Aborting...'
        sys.exit()
    if not os.path.isfile(prodFna):
        print 'Prodigal input file does not exist. Aborting...'
        sys.exit()
    if not os.path.isfile(prodFaa):
        print 'Prodigal input file does not exist. Aborting...'
        sys.exit()

    # Return the values, unpacked
    return genome, prodFna, prodFaa, out

def PrepareFastaContents(fileName):
    # Tokenise the file into per-entry chunks
    content = open(fileName, 'r').read()
    content = content.split('>')[1:]

    # Process the entries into an indexed dict of Sequence objects
    seqs = {}
    for entry in content:
        s = ParseEntryToSequence(entry)
        seqs[s.sequenceName] = s

    return seqs

def ParseEntryToSequence(entry):
    # Prepare the header and sequence variables
    entry = entry.split('\n')
    header = entry[0]
    sequence = ''.join(entry[1:])

    # Parse into a Sequence object, accounting for comments
    if ' ' in header:
        header = header.split(' ')
        return Sequence.CreateSequence(header[0], sequence, ' '.join(header[1:]))
    else:
        return Sequence.CreateSequence(header, sequence)

def VerifyBlastParameters(ident, coverage):
    # Parse the identity
    try:
        ident = float(ident)
    except:
        print 'Cannot convert specified identity threshold to a percentage, using default...'
        ident = 30.0

    # Parse the coverage
    try:
        coverage = float(coverage)
    except:
        print 'Cannot convert specified coverage threshold to a percentage, using default...'
        ident = 50.0
    return ident, coverage

def SpawnLookupStructures():

    # Change here as data changes
    uniprotPath = '/nesi/project/uoa02469/Databases/UniProt/'
    unirefPath = '/nesi/project/uoa02469/Databases/UniRef100/'
    keggPath = '/nesi/project/uoa02469/Databases/KEGG/'

    uniprotDate = '20181026'
    unirefDate = '20181026'
    keggDate = '20180523'

    # Create the paths
    uniprotTax = os.path.join(uniprotPath, 'uniprot.{}.tax'.format(uniprotDate) )
    uniprotMeta = os.path.join(uniprotPath, 'uniprot.{}.txt'.format(uniprotDate) )

    unirefTax = os.path.join(unirefPath, 'uniref100.{}.tax'.format(unirefDate) )
    unirefMeta = os.path.join(unirefPath, 'uniref100.{}.txt'.format(unirefDate) )

    keggTax = os.path.join(keggPath, 'kegg.{}.tax'.format(keggDate) )
    keggMeta = os.path.join(keggPath, 'kegg.{}.txt'.format(keggDate) )

    databaseData = namedtuple('databaseData', 'name path tax')
    uniprot = databaseData(name='UniProt', path=uniprotMeta, tax=uniprotTax)
    uniref = databaseData(name='UniRef100', path=unirefMeta, tax=unirefTax)
    kegg = databaseData(name='KEGG', path=keggMeta, tax=keggTax)
    return uniprot, uniref, kegg
