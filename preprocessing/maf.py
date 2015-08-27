from sys import stdout
import sys
import csv
import fileinput
import time
from subprocess import call
# import vcf
import os


""" This file contains a
"""

# Write the given attributes of the mutation to the file opened in opener.
def writerecord(record, opener, attributes=None):
    if (attributes != None):
        for attribute in attributes:
            if record[attribute] == None:
                opener.write('.')
            else:
                opener.write(str(record[attribute]))
            opener.write("\t")

        opener.write("\n")
    else:
        for attribute in record:
            if record[attribute] == None:
                opener.write('.')
            else:
                opener.write(str(record[attribute]))
            opener.write("\t")

        opener.write("\n")


class GeneFilter:
    def __init__(self, f):
        self.f = f
        self.argtype = 'Hugo_Symbol'
        # self.memo = {}
    def __call__(self, *args):
        return self.f(*args)
        # if not args in self.memo:
        #     self.memo[args] = self.f(*args)
        # return self.memo[args]

@GeneFilter
def OnlyUnknownGene(gene):
    # FilterUnknownGene.argtype = 'Hugo_Symbol'
    return str.lower(gene) == 'unknown'

@GeneFilter
def FilterUnknownGene(gene):
    # FilterUnknownGene.argtype = 'Hugo_Symbol'
    # return not OnlyUnknownGene(gene)
    return str.lower(gene) != 'unknown'

class VariantFilter:
    def __init__(self, f):
        self.f = f
        self.argtype = 'Variant_Classification'
        self.memo = {}
    def __call__(self, *args):
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]

@VariantFilter
def FilterSilentMutation(classification):
    return not str.lower(classification) != 'silent'



# Format for adding an attribute:
# name, argtypes, and memo

class VariantAlleleAttribute:
    def __init__(self, f):
        self.f = f
        self.name = 'Variant_Allele'
        self.argtypes = ['Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
        self.memo = {}
    def __call__(self, *args):
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]

@VariantAlleleAttribute
def VariantAlleleCalculator(ref, allele1, allele2):
    if allele1 != ref: return allele1
    elif allele2 != ref: return allele2
    else:
        print "Error: no variant allele: allele 1 and allele 2 both equal reference allele"
        raise ValueError

class GenotypeAttribute:
    def __init__(self, f):
        self.f = f
        self.name = 'Genotype'
        self.argtypes = ['Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
        self.memo = {}
    def __call__(self, *args):
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]

@GenotypeAttribute
def GenotypeCalculator(ref, allele1, allele2):
    if allele1 != ref and allele2 != ref: return 'BB'
    elif allele1 == ref and allele2 == ref:
        print "Tumor homozygous genotype found"
        return 'AA'
    else: return 'AB'



class MafReader:
    """This object stores all the sample information from the maf file, as a list of mutations.
        Each mutation is a dictionary, of key value pairs.
        All attributes are stored in upper case.
        Each record is represented by a dictionary from attributes to values.
        Records are accessed by entering the correct record ID into the dictionary of records, and getting
        the corresponding dictionary out.
        Attributes are stored as a list.
        The attributevalues structure stores a list of values for each attribute. Each value points to a list
        of the IDS with that value. """

    def __init__(self, filename=None, recordlist=None):
        self.size = 0

        # A dictionary of record ID to the record. Each record is itself a dictionary
        # from the attributes to its corresponding value.
        self.records = {}

        # A list that keeps track of the order in which the attributes were added.
        self.attributes = []

        # A dictionary, with keys as attributes, and values as dictionaries between
        # the value of the attribute, and the set of IDs that have that value.
        self.attributevalues = {}

        if (filename != None):
            self.init_parsemaf(filename)

        elif (recordlist != None):
            self.init_uselist(recordlist)

    def init_parsemaf(self, filename):
        column_header = []
        with open(filename) as tsv:
            for line in csv.reader(tsv, delimiter='\t'):
                if (line != []):
                    if (line[0] != ''):
                        if (line[0][0] == '#'):
                            pass
                        elif (line[0] == 'Hugo_Symbol'):

                            # All attributes are capitalized.
                            column_header = line

                            for attribute in column_header:
                                # There's a problem here.

                                if attribute not in self.attributes:
                                    self.attributes.append(attribute)

                                if attribute not in self.attributevalues:
                                    self.attributevalues[attribute] = {None: set()}

                                    # self.attributes[attribute] = None
                        else:
                            newrecord = {}

                            for i in range(min(len(line), len(column_header))):
                                attribute = column_header[i]
                                value = line[i]

                                if value == '.': value = None

                                newrecord[attribute] = value
                            if newrecord['Hugo_Symbol']:
                                recordID = newrecord['Tumor_Sample_Barcode'] + ':' + \
                                           newrecord['Hugo_Symbol'] + ':' + \
                                           (newrecord['Start_position'] if 'Start_position'in newrecord else newrecord['Start_Position'])

                                for attribute in newrecord:
                                    value = newrecord[attribute]

                                    # Add the recordID to attribute values to enable lookup of the record ID.
                                    if value not in self.attributevalues[attribute]:
                                        self.attributevalues[attribute][value] = set()

                                    self.attributevalues[attribute][value].add(recordID)

                                # Add the record to the dictionary of records.
                                self.records[recordID] = newrecord


                                self.size += 1
        tsv.close()

        # Ensure every record has all the attributes.
        for recordID in self.records:
            record = self.records[recordID]
            for attr in self.attributes:
                if (attr not in record):
                    record[attr] = None
                    self.attributevalues[attr][None].add(recordID)

    #
    # def init_uselist(self, recordlist):
    #     self.dict = recordlist
    #     self.size = len(recordlist)
    #     for key in self.dict[0]:
    #         self.attributes[key] = None

    # def getnewmaf(self, attribute, value):
    #     """Return a newMafReader containing the set of the records of a given value. """
    #
    #     newMafReader = MafReader()
    #     newMafReader.__init__(recordlist=self.getrecords(attribute, value, type='wholerecord'))
    #     return newMafReader


    def getrecords(self, attribute, value, type='ID'):
        """Return a set of the record IDs of those records whose attribute is a given value. """

        if (attribute not in self.attributes):
            print "Sorry, the entered attribute ", attribute, " is not in the dict."
            return

        if (type == 'ID'):
            return self.attributevalues[attribute][value]
        elif (type == 'record'):
            return [self.records[recordID] for recordID in self.attributevalues[attribute][value]]

    def getrecords_multi(self, attributelist, valuelist, type='ID'):
        """Return a set of records whose attributes match the values in valuelist. """

        if (len(attributelist) != len(valuelist)):
            print "ERROR: length of attribute list", attributelist, " not same as length of value list", valuelist
            return

        for attribute in attributelist:
            if (attribute not in self.attributes):
                print "Sorry, the entered attribute ", attribute, " is not in the attribute list."
                return

        # A list of sets of working IDs
        setlist = []
        for i in range(len(attributelist)):
            attribute = attributelist[i]
            value = valuelist[i]
            setlist.append(self.attributevalues[attribute][value])

        # All the sets that match.
        recordIDs = set.intersection(*setlist)

        if type == 'ID':
            return recordIDs
        elif type == 'record':
            return [self.records[recordID] for recordID in recordIDs]


    def getrecords_filter(self, recordfilter, type='ID'):
        # Check if the file has the attributes of the filters.
        attribute = recordfilter.argtype
        if (attribute not in self.attributes):
             print "Sorry, the entered attribute ", attribute, " is not in the attribute list."
             return

        recordIDs = set()

        # Go through the values of that attribute in the file. If filter is true, add to the recordID set.

        print "Number of values to check ", len(self.attributevalues[attribute])

        t = time.time()

        for value in self.attributevalues[attribute]:
            if value != None:
                if recordfilter(value):
                    recordIDs = recordIDs.union(self.attributevalues[attribute][value])

        print "recordfilter done in ", time.time() - t, "for record filter", recordfilter

        if type == 'ID':
            return recordIDs
        elif type == 'record':
            return [self.records[recordID] for recordID in recordIDs]



    def getrecords_filters(self, recordfilters, type='ID'):

        # Check if the file has the attributes of the filters.
        for attribute in [recordfilter.argtype for recordfilter in recordfilters]:
             if (attribute not in self.attributes):
                print "Sorry, the entered attribute ", attribute, " is not in the attribute list."
                return

         # A list of sets of working IDs
        setlist = []
        for recordfilter in recordfilters:
            recordIDs = self.getrecords_filter(recordfilter, type='ID')
            setlist.append(recordIDs)

        #The records that pass all the filters.

        allrecordIDs = set.intersection(*setlist)

        if type == 'ID':
            return allrecordIDs
        elif type == 'record':
            return [self.records[recordID] for recordID in allrecordIDs]


    def writevaluelist(self, recordIDs, filename, attribute='Hugo_Symbol'):
        opener = open(filename, 'w')

        values = set([self.records[recordID][attribute] for recordID in recordIDs])

        for value in values:
            opener.write(str(value) + '\n')



    def getattributes(self):
        return self.attributes

    def getallrecords(self):
        return self.records.values()


    def calculateattribute(self, calculator):
        """Calculate an attribute from the record, and add the new attribute."""

        # Get the attributename and inputs to the calculator.

        attribute = calculator.name
        argtypes = calculator.argtypes

        # Add the attribute to the list of own attributes and to attributevalues.

        if attribute not in self.attributes:
            self.attributes.append(attribute)
        else:
            print "Attribute ", attribute, " already added!"
            return

        if attribute not in self.attributevalues:
            self.attributevalues[attribute] = {None: set()}
        else:
            print "Attribute ", attribute, " already added!"
            return

        # Go through all the records and calculate the new value of the attributes.
        # Add the new value of the attribute to the records.
        # Add the recordID to attributevalues.

        for recordID in self.records:
            args = [self.records[recordID][argtype] for argtype in argtypes]
            value = calculator(*args)
            self.records[recordID][attribute] = value

            # Add the recordID to attribute values to enable lookup of the record ID.
            if value not in self.attributevalues[attribute]:
                self.attributevalues[attribute][value] = set()

            self.attributevalues[attribute][value].add(recordID)




    def addattribute(self, newattributes, newrecords, refattrs):
        # newrecords is a list of records. Each record is a dictionary of key/value pairs.

        # Check if the reference attribute(s) is in the attributes.
        for refattr in refattrs:
            if (refattr not in self.attributes):
                print "Sorry, the reference attribute ", refattr, "is not in the dict."
                return

        # Add the new attributes into the file's list of attributes,
        # And remove the ones that aren't new.

        newattributes = [attr for attr in newattributes if attr not in self.attributes]
        for attribute in newattributes:
            self.attributes.append(attribute)

            # Every single record has no value here.
            allrecordIDs = [recordID for recordID in self.records]
            self.attributevalues[attribute] = {None: set(allrecordIDs)}

            for ID in allrecordIDs:
                self.records[ID][attribute] = None


        print "New attributes are: ", newattributes

        for newrecord in newrecords:
            # Find the IDs  of the records with the given reference value.
            refvalues = [newrecord[refattr] for refattr in refattrs]

            IDs = self.getrecords_multi(refattrs, refvalues, type='ID')

            # Check if the IDs have been found.
            if len(IDs) == 0:
                print "Error: no record of ", refattrs, " equal to ", refvalues, " found"
            else:
                if len(IDs) > 1:
                    print "Note: multiple records of ", refattrs, " equal to ", refvalues, " found"

                # Add the new attributes to the records.
                for ID in IDs:
                    record = self.records[ID]
                    for attr in newattributes:
                        record[attr] = newrecord[attr]

            print "finished record"

    def getgenenum(self, gene):
        # Get the number of records of a given gene.
        return self.getrecordnum('Hugo_Symbol', gene)

    def getrecordnum(self, attribute, value):
        # Get the number of records with attribute equal to value.


        if (attribute not in self.attributes):
            print "Sorry, the entered attribute is not in the dict."
            return

        return (len(self.getrecords(attribute, value)))

    def writetomaf(self, filename=None, attributes=None,
                   useheader=True):

        opener = open(filename, 'w')

        if (attributes == None):
            attributes = self.attributes

        if (useheader):
            for attr in attributes:
                opener.write(attr)
                opener.write("\t")
            opener.write("\n")

        for record in self.records.values():
            writerecord(record, opener, attributes)
        opener.close()
        #
        #
        # def writerecords(self, filename = None, attributes = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"],
        #                    useheader = True):
        #     assert(isinstance(attributes, (list, tuple)))
        #
        #     if (filename != None): opener = open(filename, 'w')
        #     else: opener = stdout
        #
        #     if (useheader):
        #         for attr in attributes:
        #             opener.write(attr)
        #             opener.write("\t")
        #         opener.write("\n")
        #
        #     for record in self.dict:
        #         writerecord(record, opener,attributes)
        #     opener.close()
        #
        #
        # def sortby(self, attribute):
        #     self.dict = sorted(self.dict, key = lambda record: record[attribute])
        #
        #

    # def filter(self, *filters, returnTrue):
    #
    #     for filter in filters:
    #         argtype = filter.argtype
    #         for value in self.attributevalues[argtype].keys():
    #             if not filter(value):
    #                 del self.attributevalues[value]



    def writepatientfiles(self, delimiter='\t', fieldnames=['Chromosome', 'Start_Position', 'Reference_Allele',
                                                            'Variant_Allele', 'Genotype'],
                          basetype='one', suffix='.tsv', recordIDs=None, IDname='Tumor_Sample_Barcode',
                          folderprefix='patientfiles/'):
        """
        For each patient ID, make a new file in the working directory
        and write all the variants specified in records to the file, using the
        corresponding fieldnames.
        :param delimiter: Default: tab.
        :param fieldnames: THe header of the file to write
        :param basetype: Either 'one' or 'zero'. Type of indexing to use.
        :param suffix: Suffix of file name.
        :return: Nothing.
        """

        if recordIDs:
            recordIDs = set(recordIDs)

        for ID in self.attributevalues[IDname]:
            if ID:
                filename = folderprefix + ID + suffix
                with open(filename, 'w') as patient_file:
                    writer = csv.DictWriter(patient_file, fieldnames=fieldnames, delimiter=delimiter, extrasaction='ignore')
                    writer.writeheader()
                    if recordIDs:
                        patient_recordIDs = recordIDs.intersection(self.attributevalues[IDname][ID])
                        recordIDs = recordIDs.difference(patient_recordIDs)
                        for recordID in patient_recordIDs:
                             writer.writerow(self.records[recordID])

                    else:
                        for recordID in self.attributevalues[IDname][ID]:
                            writer.writerow(self.records[recordID])



    def writemutationmatrix(self, filename = None, level = 'patient', justgenes = True, nosilent = True,
                            nounknown = True, varlist = [None], gainloss=False):

        #NO silent and nounknown need to be fixed. -jlu 6/12/15


        if (level == 'patient'):
            endline = 11
        elif (level == 'sample'):
            endline = 14
        elif (level == 'vial'):
            endline = 16
        elif (level == 'portion'):
            endline = 18
        elif (level == 'analyte'):
            endline = 19
        elif (level == 'plate'):
            endline = 25
        elif (level == 'complete'):
            endline = None
        elif (isinstance(level, int)):
            endline = level
            endline -= 1
        else:
            endline = -1

        # The list of variants to be included in the writing.
        varDict = {'missense': 'm', 'splice': 'sp', 'nonsense': 'n', 'frame_shift': 'f',
                   'frameshift': 'f', 'nonstop': 'ns'}

        if (not nosilent): varDict['silent'] = 's'

        if gainloss:
            suffix = 'loss'
        else:
            suffix = ''


        # A dictionary from the sample barcodes of given precision to the corresponding record IDs
        barcodeDict = {}



        # Temporarily store all the Nones.
        NonesampleID = self.attributevalues['Tumor_Sample_Barcode'].pop(None)

        # Temporarily store all the Unknown genes:
        # if nounknown:


        # print "endline is ", endline
        for sampleID in self.attributevalues['Tumor_Sample_Barcode']:
            if endline:

                barcode = sampleID[0:endline+1]
            else:
                barcode = sampleID
            recordIDs = self.attributevalues['Tumor_Sample_Barcode'][sampleID]
            if barcode not in barcodeDict:
                barcodeDict[barcode] = recordIDs
            else:
                barcodeDict[barcode] = barcodeDict[barcode].union(recordIDs)

        with open(filename, 'w') as opener:

            for barcode in barcodeDict:
                geneset = set()

                for recordID in barcodeDict[barcode]:


                    newgene = self.records[recordID]["Hugo_Symbol"]
                    if (not nounknown) or (FilterUnknownGene(newgene)) :
                        if not justgenes:
                            varname = str.lower(self.records[recordID]["Variant_Classification"])
                            for var in varDict:
                                if var in varname:
                                    newgene += varDict[varname]
                                    break



                        if newgene not in geneset:
                            geneset.add(newgene + suffix)

                tab = '\t'
                opener.write(barcode + '\t' + tab.join(geneset) + '\n')

        opener.close()
        #
        # #Check for redundancies
        #
        # for line in fileinput.input(files = [filename], inplace= True):
        #     if (line.split() != []):
        #         entries = line.split()
        #         numremoved = 0
        #         for entry in entries:
        #             if (entries.count(entry) > 1):
        #                 entries.remove(entry)
        #                 #print "entry removed", entry
        #                 numremoved += 1
        #         #print "number of entries removed", numremoved
        #
        #         newline = ''
        #         for i in range(len(entries) - 1):
        #              newline = newline + entries[i] + '\t'
        #
        #         # print entries
        #         # print len(entries) -1
        #         # print entries[len(entries)-1]
        #         newline = newline + entries[len(entries)-1]
        #         print newline
        #     # else:
        #     #     print line




def main():
    # testMafReader = MafReader()
    # testMafReader.__init__(filename = '/Users/jlu96/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf')
    # SilentMut = testMafReader.getnewmaf("Variant_Classification", "Missense_Mutation")
    # # SilentMut.writerecord("Variant_Classification")
    # # SilentMut.printrecord(attributes = ["Entrez_Gene_Id", "Variant_Classification"])
    # testMafReader.writerecords(filename = './testrecord', attributes = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"])
    # testMafReader.sortby("Tumor_Sample_Barcode")
    # testMafReader.writerecords(filename = './testrecordsorted1', attributes = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"])
    # testMafReader.sortby("Variant_Classification")
    # testMafReader.writerecords(filename = './testrecordsorted2', attributes = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification"])
    # testMafReader.writemutationmatrix(filename = './testmutationmatrix.m2')

    # Testing code.
    t = time.time()
    COAD = MafReader()
    COAD.__init__(filename='COAD.maf')
    print "First COAD finished reading!"


    COAD.writemutationmatrix('COAD_test.m2')
    t1 = time.time()

    print "First COAD finished writing to mutation matrix in time: ", t1 - t
    # IDs = COAD.getrecords('Hugo_Symbol', 'TP53', type='ID')
    # for ID in IDs:
    #     print COAD.records[ID]['Hugo_Symbol'],  COAD.records[ID]['Tumor_Sample_Barcode'], ID
    #
    # records = COAD.getrecords('Hugo_Symbol', 'TP53', type='record')
    # for record in records:
    #     print record['Hugo_Symbol'],  record['Tumor_Sample_Barcode']
    #
    # print "number of TP53", COAD.getgenenum('TP53')
    # print "number of TP53 according to IDs:", len(IDs)
    # print "number of TP53 according to records:", len(records)


    COAD2 = MafReader()
    COAD2.__init__(filename='COAD_modifiedheader.maf')
    print "Second COAD finished reading!"

    t2 = time.time()
    print "Time used to read: ", t2 - t1

    newattributes = COAD2.getattributes()
    newrecords = COAD2.getallrecords()[0:10]
    COAD.addattribute(newattributes, newrecords, ['Tumor_Sample_Barcode'])

    print "Finished adding new attributes!"

    t3 = time.time()
    print 'Time used: ', t3 - t2, 'to add new attributes to ', len(newrecords), 'new records'


    COAD.writetomaf('COAD_testnewattr.maf')
    t4 = time.time()
    print 'Time used: ', t4 - t3, ' to write to maf file'


    #
    # idattr = 'Tumor_Sample_Barcode'
    # newattributes = ['TestHa', 'Yolo']
    # newrecords = COAD2.getlist()[0:4]
    # COAD.addattribute(idattr, newattributes, newrecords)


if __name__ == '__main__':
    main()
