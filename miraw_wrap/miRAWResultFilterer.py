import argparse
import sys
import os
import logging
import csv
#import datetime
from datetime import datetime
import hashlib

DEBUG = 1
TESTRUN = 0
PROFILE = 0

MY_NEWLINE = "\n"
if os.name== "Windows":
    MY_NEWLINE ="\r\n"


__all__ = []
__version__ = 0.1
__date__ = '2020-12-19'
__updated__ = '2020-12-19'


EXTRACT_POSITIVE_PREDICTIONS    = "EXPOS"
EXTRACT_NEGATIVE_PREDICTIONS    = "EXNEG"
FILTER_POSITIVE_PREDICTIONS     = "FTPOS"
FILTER_NEGATIVE_PREDICTIONS     = "FTNEG"
FILTER_BY_MIRNA                 = "FTMIR"

TARGET_FILE_EXTENSION           = "csv"
FILTERED_TARGET_FILE_EXTENSION  = "filtered.csv"
SITES_FILE_EXTENSION            = "csv"
FILTERED_SITES_FILE_EXTENSION   = "filtered.csv"
UNIFIED_FILE_EXTENSION          = "csv"
FILTERED_UNIFIED_FILE_EXTENSION = "filtered.csv"   
FILTERED_DETAILS_FILE_EXTENSION = "detailed.csv" 
    
UTR_QUERY_DELIMITER             = ":::"    
bindingEnergyCutoff = -1000000.0
positiveProbabilityCutoff = 0.0
negativeProbabilityCutoff = 0.0
removePositivePredictions = False
extractNegativePredictions = False
extractConflictPredictions = False
filterPredictionsByList = False
removePositivePredictions = False
keepPositivePredictions = False


miRAWtargetPredictionFile = ""
miRAWpositiveSitesFile = ""
miRAWunifiedFile = ""
mirnaFilterFile = ""
filteredTargetPredictionFile = ""
filteredPositiveSitesFile = ""
filteredUnifiedFile = ""
filteredPositiveSitesDetailedFile = ""


filterDropCount = 0
conflictDropCount = 0
keepDropCount = 0
posDropCount = 0
negDropCount = 0
posProbDropCount = 0
negProbDropCount = 0
emptyCount = 0
totalDrops = 0

headerString = "#" + MY_NEWLINE


class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def initLogger(md5string):
    ''' setup log file based on project name'''

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, logFileName + "__" + dt_string + "__" + md5string + ".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)

    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)  # set the new handler
    log.addHandler(handler)
    logging.info("+" + "*" * 78 + "+")
    logging.info("project log file is <" + logfileName + ">")
    logging.info("+" + "*" * 78 + "+")
    logging.debug("debug mode is on")


#GeneName    GeneId    miRNA    Prediction    HighestPredVal    LowestPredVal    PosSites    NegSites    RemovedSites        miRNA    miRNA

GENE_NAME = 0
MIRNA = 1
SITE_START = 2
SITE_END = 3
PREDICTION = 4
PAIR_START_IN_SITE = 5
SEED_START = 6
SEED_END = 7
PAIRS = 8
WC = 9
WOB = 10
MFE = 11
CANONICAL = 12
SITE_TRANSCRIPT_5TO3 = 13
MATURE_MIRNA_TRANSCRIPT = 14
BRACKET_NOTATION = 15
FILTERED = 16
FILTERING  = 17
REASON = 18
ADDITIONAL_PROPERTIES = 19


GENE_NAME_COL          = 0
GENE_ID_COL            = 1
MIRNA_NAME_COL         = 2
PREDICTION_COL         = 3
HIGHEST_PREDICTION_COL = 4
LOWEST_PREDICTION_COL  = 5
POSITIVE_SITES_COL     = 6
NEGATIVE_SITES_COL     = 7
REMOVED_SITES_COL      = 8

GENE_NAME_ID           = "GeneName"
GENE_ID_ID             = "GeneID"
MIRNA_NAME_ID          = "miRNA"
PREDICTION_ID          = "Prediction"
HIGHEST_PREDICTION_ID  = "HighestPredVal"
LOWEST_PREDICTION_ID   = "LowestPredVal"
POSITIVE_SITES_ID      = "PosSites"
NEGATIVE_SITES_ID      = "NegSites"
REMOVED_SITES_ID       = "RemovedSites"

UNIFIED_UTRID_COL      = 2
UNIFIED_UTRSEQ_COL     = 5



HEADER_LINE = GENE_NAME_ID + "\t"\
                + GENE_ID_ID + "\t"\
                + MIRNA_NAME_ID + "\t"\
                + PREDICTION_ID + "\t"\
                + HIGHEST_PREDICTION_ID + "\t"\
                + LOWEST_PREDICTION_ID + "\t"\
                + POSITIVE_SITES_ID + "\t"\
                + NEGATIVE_SITES_ID + "\t"\
                + REMOVED_SITES_ID
                


geneNames = []
geneIDs = []
miRNANames = []
predictions = []
highestPredVal = []
lowestPredVal = []
positiveSites = []
negativeSites = []
removedSites = []
dropPredictions = []          # use this to keep track of which predictions to drop
positiveSiteLines = []
positiveSiteLocations = {}

filteredMiRNAList = []
utrSequences = {}
miRNASequences = {}

logging.getLogger().setLevel(logging.INFO)

def parseArgs(argv):
    '''parse out Command line options.'''

    global parser

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = "filterMiRAWResults"

    program_license = '''%s
    i
      Created by Simon Rayner on %s.
      Copyright 2020 Oslo University Hospital. All rights reserved.

      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0

      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.

    USAGE
    ''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser

        parser = argparse.ArgumentParser(description='parse miRAW prediction results')

        parser.add_argument("-H", "--HelpMe", action="store_true",
                            help="print detailed help")
        parser.add_argument("-t", "--target_pred_file", dest='targetPredFile',
                            help="miRAW summary of results file (.targetPredictionOutput.csv)")
        parser.add_argument("-s", "--target_site_file", dest='targetSiteFile',
                            help="miRAW detailed results file (.positiveTargetSites.csv)")
        parser.add_argument("-u", "--unified_input_file", dest='unifiedInputFile',
                            help="miRAW input file used for target prediction (.unifiedFile.csv)")
        parser.add_argument("-e", "--energy-cutoff", dest='energyCutOff',
                            help="remove predictions ABOVE this energy (i.e. weaker bindings)")
        parser.add_argument("-p", "--pos-prob-cutoff", dest='posProbCutOff',
                            help="remove positive predictions BELOW this prediction probability")
        parser.add_argument("-n ", "--neg-prob-cutoff", dest='negProbCutOff',
                            help="remove negative predictions ABOVE this prediction probability")
        parser.add_argument("-L ", "--extract_from_list", dest='extractListFile',
                            help="remove negative predictions ABOVE this prediction probability")
        parser.add_argument("-R", "--remove_conflicts", action="store_true",
                            help="remove predictions with both positive and negative predictions")
        parser.add_argument("-K", "--keep_conflicts", action="store_true",
                            help="keep predictions with both positive and negative predictions")
        parser.add_argument("-P", "--extract_pos_preds", action="store_true",
                            help="keep positive predictions, remove negative prediction")
        parser.add_argument("-N", "--extract_neg_preds", action="store_true",
                            help="keep negative predictions, remove positive predictions")
        global args
        args = parser.parse_args()
        if args.HelpMe:
            printLongHelpAndExit()
            exit()



    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        print(e)
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2




def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+  miRAWResultFilterer                                                         +")
    logging.info("+    wrapper code for generating files and scripts to run miRAW                +")
    logging.info("+                                                                              +")
    logging.info("+    you need to specify:                                                      +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW summary result file as input        (-t/--target_pred_file)     +")
    logging.info("+        (this is the file which ends in 'targetPredictionOutput.csv')         +")
    logging.info("+                                                                              +")
    logging.info("+      a miRAW detailed results file as input      (-s/--target_site_file)     +")
    logging.info("+        (this is the file which ends in 'positiveTargetSites.csv')            +")
    logging.info("+                                                                              +")
    logging.info("+      the input file used for prediction          (-u/--unified_input_file)   +")
    logging.info("+        (this is the file which ends in 'unifiedFile.csv')                    +")
    logging.info("+                                                                              +")
    logging.info("+      remove conflicts                            (-R/--remove_conflicts)     +")
    logging.info("+         (remove predictions that contain both                                +")
    logging.info("+          positive and negative predictions                                   +")
    logging.info("+          This isn't necessarily a problem, it's                              +")
    logging.info("+          just another way to filter the data)                                +")
    logging.info("+                                                                              +")
    logging.info("+      keep conflicts                              (-K/--keep_conflicts)       +")
    logging.info("+         (keep predictions that contain both                                  +")
    logging.info("+          positive and negative predictions)                                  +")
    logging.info("+                                                                              +")
    logging.info("+      binding energy cutoff:                      (-e/--energy-cutoff)        +")
    logging.info("+        (positive predictions with a probability below                        +")
    logging.info("+         this value will be removed                                           +")
    logging.info("+                                                                              +")
    logging.info("+      positive probability cutoff:                (-p/--pos_prob_cutoff)      +")
    logging.info("+        (positive predictions with a probability below                        +")
    logging.info("+         this value will be removed                                           +")
    logging.info("+         values can be between 0 --> 1                                        +")
    logging.info("+         higher values are better)                                            +")
    logging.info("+                                                                              +")
    logging.info("+      negative probability cutoff:                (-n/--neg_prob_cutoff)      +")
    logging.info("+        (negative predictions with a probability above                        +")
    logging.info("+         this value will be removed                                           +")
    logging.info("+         values can be between 0 --> -1                                       +")
    logging.info("+         lower values are better)                                             +")
    logging.info("+                                                                              +")
    logging.info("+      keep all positive predictions               (-P/--extract_positives)    +")
    logging.info("+      remove all negative predictions                                         +")
    logging.info("+                                                                              +")
    logging.info("+      keep all negative predictions               (-N/--extract_negatives)    +")
    logging.info("+      remove all positive predictions                                         +")
    logging.info("+                                                                              +")
    logging.info("+      extract miRNAs in the supplied list.        (-L/--extract_mirnas)       +")
    logging.info("+          remove miRNAs that are listed in the                                +")
    logging.info("+          supplied file                                                       +")
    logging.info("+          (use this to extract miRNAs you don't want                          +")
    logging.info("+          for example, high numbered/low confidence                           +")
    logging.info("+          or miRNAs in which you are particularly interested)                 +")
    logging.info("+                                                                              +")
    logging.info("+   filtering is performed in the order                                        +")
    logging.info("+      listFilter->Conflicts->PosCutOff->NegCutOff->ExtractPos/Neg Preds       +")
    logging.info("+                ->BindingEnergy                                               +")
    logging.info("+                                                                              +")
    logging.info("+" + "-" * 78 + "+")


def printHelpAndExit():
    parser.print_help()
    logging.info("stopping") 
    sys.exit()

def checkExtractConflicts(): 
    logging.info("extractConflicts")       
    global extractConflictPredictions, keepConflictPredictions
    if args.remove_conflicts & args.keep_conflicts: 
        logging.error("you can't set both '--keep_conflicts' and '--remove_conflicts'")
        logging.info("stopping") 
        sys.exit()
    if args.remove_conflicts:
        if args.remove_conflicts==True:
            extractConflictPredictions=True
            keepConflictPredictions = False
            logging.info("--remove_conflicts is True")
        else:
            logging.info("--remove_conflicts is False")
    if args.keep_conflicts:
        if args.keep_conflicts==True:
            extractConflictPredictions=False
            keepConflictPredictions=True
            logging.info("--keep_conflicts is True")
        else:
            logging.info("--keep_conflicts is False")
    logging.info("--OK")
    
    
def checkBindingEnergyCutOff():
    logging.info("checkBindingEnergyCutOff") 
    
    global bindingEnergyCutoff  
    if args.energyCutOff:
        logging.info(args.energyCutOff)
        if float(args.energyCutOff) >= 0.0:
            bindingEnergyCutoff=float(args.energyCutOff)
        else:
            logging.info(" binding energy cutoff (-e) needs to be less than 0" + str(args.energyCutOff) + ">")
            return()
     

def checkPositiveCutOff():
    logging.info("checkPositiveCutOff") 
    
    global positiveProbabilityCutoff  
    if args.posProbCutOff:
        logging.info(args.posProbCutOff)
        if( float(args.posProbCutOff) >= 0.0) & (float(args.posProbCutOff) <= 1.0):
            positiveProbabilityCutoff=float(args.posProbCutOff)
        else:
            logging.info(" positive probability cutoff (-p) needs to be between 0 and 1" + str(args.posProbCutOff) + ">")
            return()
     

def checkNegativeCutOff():
    logging.info("checknegativeCutOff")   
    
    global negativeProbabilityCutoff
    if args.negProbCutOff:
        logging.info(args.negProbCutOff)
        if (float(args.negProbCutOff) <= 0.0) & (float(args.negProbCutOff) >= -1.0):
            negativeProbabilityCutoff=float(args.negProbCutOff)
        else:
            logging.info(" negative probability cutoff (-n) needs to be between -1 and 0" + str(args.negProbCutOff) + ">")
            return()
    


def checkExtractPositives():
    logging.info("checkExtractPositives")
    global removePositivePredictions
    if args.extract_pos_preds:
        if args.extract_pos_preds==True:
            removePositivePredictions=True
            logging.info("--extract_pos_preds is True")
        else:
            logging.info("--extract_pos_preds is False")
    
    logging.info("--OK")


def checkExtractNegatives():
    logging.info("checkExtractNegatives")
    global extractNegativePredictions
    if args.extract_neg_preds:
        if args.extract_neg_preds==True:
            extractNegativePredictions=True
            logging.info("--extract_neg_preds is True")
        else:
            logging.info("--extract_neg_preds is False")
    logging.info("--OK")
    
    
def checkRemovemiRNAs():
    logging.info("check removeMiRNAsInList File")
    global mirnaFilterFile, filterPredictionsByList
    if args.extractListFile:
        filterPredictionsByList = True
        mirnaFilterFile=args.extractListFile
        logging.info("--extractListFile is " + mirnaFilterFile)
        if not os.path.isfile(mirnaFilterFile):
            logging.error("--can't find miRNA Extract List File at <" + mirnaFilterFile + ">")
            exit()
    logging.info("--OK")
        

def checkTargetPredictionFile():
    global miRAWtargetPredictionFile, filteredTargetPredictionFile, logFileName
    logging.info("check TargetPredictionFile")
    if args.targetPredFile:
        miRAWtargetPredictionFile=args.targetPredFile
        filteredTargetPredictionFile=miRAWtargetPredictionFile.replace(TARGET_FILE_EXTENSION, \
                                                                       FILTERED_TARGET_FILE_EXTENSION)
        logging.info("--Target Prediction File is " + miRAWtargetPredictionFile)
        logging.info("--Filtered Target Prediction File is " + filteredTargetPredictionFile)
        if not os.path.isfile(miRAWtargetPredictionFile):
            logging.error("--can't find Target Prediction File at <" + miRAWtargetPredictionFile + ">")
            exit()
        logFileName =miRAWpositiveSitesFile.split(".")[0]
    logging.info("--OK")


def checkPositiveTargetsFile():
    global miRAWpositiveSitesFile, filteredPositiveSitesFile
    global filteredPositiveSitesDetailedFile
    logging.info("check PositiveTargetsFile")
    if args.targetSiteFile:
        miRAWpositiveSitesFile=args.targetSiteFile
        filteredPositiveSitesFile=miRAWpositiveSitesFile.replace(SITES_FILE_EXTENSION, \
                                                                       FILTERED_SITES_FILE_EXTENSION)
        filteredPositiveSitesDetailedFile=miRAWpositiveSitesFile.replace(SITES_FILE_EXTENSION, \
                                                                       FILTERED_DETAILS_FILE_EXTENSION)
        
        logging.info("--Positive Sites File is " + miRAWpositiveSitesFile)
        logging.info("--Filtered Positive Sites File is " + filteredPositiveSitesFile)
        logging.info("--Detailed Positive Sites File is " + filteredPositiveSitesDetailedFile)
        if not os.path.isfile(miRAWpositiveSitesFile):
            logging.error("--can't find Positive Sites File at <" + miRAWpositiveSitesFile + ">")
            exit()
    logging.info("--OK")


def checkUnifiedInputFile():
    global miRAWunifiedFile, filteredUnifiedFile
    logging.info("check unifiedFile")
    if args.unifiedInputFile:
        miRAWunifiedFile=args.unifiedInputFile
        filteredUnifiedFile=miRAWunifiedFile.replace(UNIFIED_FILE_EXTENSION, \
                                                                       FILTERED_UNIFIED_FILE_EXTENSION)
        logging.info("--Unified Input File is " + miRAWunifiedFile)
        logging.info("--Unified Input File is " + filteredUnifiedFile)
        if not os.path.isfile(miRAWunifiedFile):
            logging.error("--can't find Unified Input File at <" + miRAWunifiedFile + ">")
            exit()
    logging.info("--OK")


def checkForBadFlagCombinations():
    if(removePositivePredictions & extractNegativePredictions):
        logging.warn("you have both --extract_pos_preds and --extract_neg_preds flags set")
        logging.warn("nothing will be filtered with this combination")
        


def readPositiveSitesFile():
    logging.info("read PositiveSiteFile")
    global miRAWpositiveSitesFile, positiveSiteLines
    utrNames = []
    global positiveTargetsHeaderLine
    global gGeneName
    global gMiRNA
    global gSiteStart
    global gSiteEnd
    global gPrediction
    global gPairStartInSite
    global gSeedStart
    global gSeedEnd
    global gPairs
    global gWC
    global gWob
    global gMFE
    global gCanonical
    global gSiteTranscript5to3
    global gMatureMiRNATranscript
    global gBracketNotation
    global gFiltered
    global gFiltering
    global gReason
    global gAdditionalProperties

    gGeneName = []
    gMiRNA = []
    gSiteStart = []
    gSiteEnd = []
    gPrediction = []
    gPairStartInSite = []
    gSeedStart = []
    gSeedEnd = []
    gPairs = []
    gWC = []
    gWob = []
    gMFE = []
    gCanonical = []
    gSiteTranscript5to3 = []
    gMatureMiRNATranscript = []
    gBracketNotation = []
    gFiltered = []
    gFiltering = []
    gReason = []
    gAdditionalProperties = []

    with open(miRAWpositiveSitesFile, 'r') as fP:
        positiveSiteLines = fP.readlines()
        p=0
        for line in positiveSiteLines:
            if p==0:
                positiveTargetsHeaderLine = line
                p+=1
                continue
            if not line.split("\t")[0].isdigit():
                gGeneName.append(line.split("\t")[GENE_NAME].strip())
                gMiRNA.append(line.split("\t")[MIRNA].strip())
                gSiteStart.append(line.split("\t")[SITE_START].strip())
                gSiteEnd.append(line.split("\t")[SITE_END].strip())
                gPrediction.append(line.split("\t")[PREDICTION].strip())
                gPairStartInSite.append(line.split("\t")[PAIR_START_IN_SITE].strip())
                gSeedStart.append(line.split("\t")[SEED_START].strip())
                gSeedEnd.append(line.split("\t")[SITE_END].strip())
                gPairs.append(line.split("\t")[PAIRS].strip())
                gWC.append(line.split("\t")[WC].strip())
                gWob.append(line.split("\t")[WOB].strip())
                gMFE.append(line.split("\t")[MFE].strip())
                gCanonical.append(line.split("\t")[CANONICAL].strip())
                gSiteTranscript5to3.append(line.split("\t")[SITE_TRANSCRIPT_5TO3].strip())
                gMatureMiRNATranscript.append(line.split("\t")[MATURE_MIRNA_TRANSCRIPT].strip())
                gBracketNotation.append(line.split("\t")[BRACKET_NOTATION].strip())
                gFiltered.append(line.split("\t")[FILTERED].strip())
                gFiltering.append(line.split("\t")[FILTERING].strip())
                gReason.append(line.split("\t")[REASON].strip())
                #gAdditionalProperties.append(line.split("\t")[ADDITIONAL_PROPERTIES].strip())
                #positiveSiteLocations[line.split("\t")[0].strip() + UTR_QUERY_DELIMITER + line.split("\t")[1].strip()]=p
                #miRNASequences[line.split("\t")[1].strip()] = line.split("\t")[2].strip()


                dropPredictions.append(0)
            p+=1
                
    logging.info("--read <" + str(len(positiveSiteLines)) + "> lines")            
    logging.info("--read <" + str(len(dropPredictions)) + "> target pairs")
    logging.info("--done")
        
 
    
    
    
def readUnifiedFile():
    logging.info("read UnifiedFile")
    global utrSequences
    
    with open(miRAWunifiedFile, 'r') as fF:
        utrLines = fF.readlines()
        u=0
        for line in utrLines:
            if u!=0:
                utrID = line.split("\t")[UNIFIED_UTRID_COL]
                if not utrID in utrSequences:
                    utrSequences[utrID] = line.split("\t")[UNIFIED_UTRSEQ_COL]
            u+=1       
    logging.info("--read " + str(len(utrSequences)) + " 3'UTR sequence(s)")
    logging.info("--done")     

    


    
    
def readTargetPredictionFile():
    logging.info("read TargetPredictionFile")
    global miRAWtargetPredictionFile      

    with open(miRAWtargetPredictionFile,'r') as fT:
 
        targetLines = fT.readlines()
        l=0
        for line in targetLines:
            #logging.info(line)
            if l!=0:
                geneNames.append(line.split("\t")[GENE_NAME_COL])
                #reader=csv.reader(fT, delimiter='\t')
                #for geneName, geneID, mirName, pred, highPred, lowPred, posSite, negSite, remSite in reader:
                geneIDs.append(line.split("\t")[GENE_ID_COL])
                miRNANames.append(line.split("\t")[MIRNA_NAME_COL])
                predictions.append(line.split("\t")[PREDICTION_COL])
                highestPredVal.append(line.split("\t")[HIGHEST_PREDICTION_COL])
                lowestPredVal.append(line.split("\t")[LOWEST_PREDICTION_COL])    
                positiveSites.append(line.split("\t")[POSITIVE_SITES_COL])
                negativeSites.append(line.split("\t")[NEGATIVE_SITES_COL])
                removedSites.append(line.split("\t")[REMOVED_SITES_COL])     
                #dropPredictions.append(0)
            l+=1
    logging.info("--read " + str(len(geneNames)) + "3'UTR/miRNA prediction pairs")
    logging.info("--done")     
            
            
            
def readFilterFile():
    logging.info("readFilterFile")
    global mirnaFilterFile, filteredMiRNAList
    with open(mirnaFilterFile, 'r') as fM:        
        filteredMiRNAList = fM.readlines()
    logging.info("--read " + str(len(filteredMiRNAList)) + " entries")
    
    
                
def processTargetResults():
    logging.info("process target results")



def filterByList():
    logging.info("filter by list")
    global filterDropCount, filteredMiRNAList, totalDrops
    filterDropCount = 0
    for filterMiRName in filteredMiRNAList:
        filterMiRName = filterMiRName.strip()
        r = 0
        for miRNAName in miRNANames:
            logging.debug(filterMiRName + "\t" + miRNAName)
            if(filterMiRName == miRNAName):
                dropPredictions[r] = 1
                filterDropCount += 1
            r += 1

    logging.info("--Dropped " + str(filterDropCount) + " entries")
    totalDrops += filterDropCount
    
    

def filterByProbability():
    logging.info("filter by prediction probability")
    global posProbDropCount, negProbDropCount, totalDrops
    global positiveProbabilityCutoff, negativeProbabilityCutoff
    
    r = 0
    posProbDropCount = 0
    negProbDropCount = 0
    while r<len(dropPredictions):        
        #logging.info(str(positiveSites[r]) + "|"  \
        #             + str(negativeSites[r]) + "\t" + str(positiveProbabilityCutoff) \
        #             + "|"+ str(negativeProbabilityCutoff) + "\t" + str(predictions[r]) )

        if float(gPrediction[r]) > positiveProbabilityCutoff:
            dropPredictions[r]=1
            posProbDropCount += 1

        r += 1
        
    logging.info("--Dropped " + str(posProbDropCount) + " positive entries")
    totalDrops = posProbDropCount


def filterByMFEAndProbability():
    logging.info("filter by MFE")
    global keepCount

    r = 0
    keepCount = 0
    while r < len(dropPredictions):
        if float(gMFE[r]) < -bindingEnergyCutoff and float(gPrediction[r]) > positiveProbabilityCutoff:
            dropPredictions[r]=1
            keepCount += 1
        r += 1

    logging.info("--kept " + str(keepCount) + " entries")



def filterByPositive():
    logging.info("filter by positives")
    global posDropCount, totalDrops
    
    r = 0
    posDropCount = 0
    while r<len(dropPredictions):
        if int(positiveSites[r]) == 0 :
            dropPredictions[r]=1
            posDropCount += 1
        r += 1
        
    logging.info("--Dropped " + str(posDropCount) + " entries")
    totalDrops += posDropCount


                    
def filterByNegative():
    logging.info("filter by negatives")
    global negDropCount, totalDrops
    
    r = 0
    negDropCount = 0
    while r<len(dropPredictions):
        if int(negativeSites[r]) == 0 :
            dropPredictions[r]=1
            negDropCount += 1
        r += 1
        
    logging.info("--Dropped " + str(negDropCount) + " entries")
    totalDrops += negDropCount


def removeEmptyPredictions():
    global emptyCount, totalDrops
    
    r=1
    while r<len(dropPredictions):        

        if int(positiveSites[r]) == 0 & int(positiveSites[r]) == 0:
            dropPredictions[r]=1
            emptyCount += 1
        r += 1
        
    logging.info("--Dropped " + str(emptyCount) + " empty entries")
    totalDrops += emptyCount
    

    
def keepByConflicts():
    global keepDropCount, totalDrops
    
    r = 0
    keepDropCount = 0
    while r<len(dropPredictions):
        if not (int(positiveSites[r]) != 0) & (int(negativeSites[r]) != 0):
            dropPredictions[r]=1
            keepDropCount += 1
        r += 1
        
    logging.info("--Kept " + str(keepDropCount) + " entries")
    totalDrops += keepDropCount


def filterByConflicts():
    global loggingconflictDropCount, totalDrops
    
    r = 0
    conflictDropCount = 0
    while r<len(dropPredictions):
        if (int(positiveSites[r]) != 0) & (int(negativeSites[r]) != 0):
            dropPredictions[r]=1
            conflictDropCount += 1
        r += 1
        
    logging.info("--Dropped " + str(conflictDropCount) + " entries")
    totalDrops += conflictDropCount


def summarizeFiltering():
    global totalDrops
    totalDrops = 0
    for drops in dropPredictions:
        totalDrops += int(drops)
        
    logging.info("--Dropped a total of " + str(totalDrops) + " entries from " \
                 + str(len(dropPredictions)) + " predictions")





def grabTargetPairLines(queryLine): 
    global positiveSiteLocations, positiveSiteLines
    logging.info(queryLine)
    logging.info(positiveSiteLocations[queryLine])
    
    lineNo = positiveSiteLocations[queryLine]
    queryLines = []
    l = lineNo + 1
    try:
        while(positiveSiteLines[l].split("\t")[1].isdigit()):
            queryLines.append(positiveSiteLines[l])
            l+=1
        return queryLines
    
    except IndexError: # in case of empty line at end of file
        return queryLines
    
    
     



def writeFilteredDetailedData():
    global filteredPositiveSitesDetailedFile
    logging.info(" write filtered detailed Data")
    
    r = 0
    with open(filteredPositiveSitesDetailedFile,'w') as fD:
        while r<len(dropPredictions):
            if dropPredictions[r]==0:
                targetQuery=geneIDs[r]+ UTR_QUERY_DELIMITER + miRNANames[r]
                
                predictionOutputString = geneNames[r] + "\t" + \
                          geneIDs[r] + "\t" + \
                          miRNANames[r] + "\t" + \
                          predictions[r] + "\t" + \
                          highestPredVal[r] + "\t" + \
                          lowestPredVal[r] + "\t" + \
                          positiveSites[r] + "\t" + \
                          negativeSites[r] + "\t" + \
                          removedSites[r] + "\t" 
                # Need to collect all the other data related to this entry
                # from the positive sites file
                # and the sequence from the unifiedFile
                
                targetLines = grabTargetPairLines(targetQuery)
                for targetLine in targetLines:
                    utrStart = int(targetLine.split("\t")[0])-1
                    utrStop = int(targetLine.split("\t")[1])-1
                    logging.info(str(utrStart) + ":" + str(utrStop))
                    fD.write(predictionOutputString + "\t" + \
                             targetLine.strip() + "\t" + \
                             miRNASequences[miRNANames[r]] + "\t" + \
                             utrSequences[geneIDs[r]][utrStart:utrStop] + \
                             MY_NEWLINE)

            r += 1
    logging.info("--finished")
    
    


def writeFilteredPositiveData():
    logging.info(" write filtered Positive Site Data")

    global filteredTargetPredictionFile, headerString, totalDrops
    global conflictDropCount, keepDropCount, emptyCount, filterCount
    global negDropCount, posDropCount, negProbDropCount, posProbDropCount
    global negativeProbabilityCutoff, positiveProbabilityCutoff
    headerString = headerString + "#             " + MY_NEWLINE
    headerString = headerString + "#   filtered positiveTargetSites.detailed.csv file" + MY_NEWLINE
    headerString = headerString + "#   generated by miRAWResultFilterer" + MY_NEWLINE
    headerString = headerString + "#   on  " + str(datetime.now()) + MY_NEWLINE
    headerString = headerString + "#             " + MY_NEWLINE
    if keepDropCount > 0:
        headerString = headerString + "#   kept " + str(keepCount) + " predictions #" + MY_NEWLINE
    else:
        headerString = headerString + \
                       "#   removed " + str(len(dropPredictions) - keepCount) + " entries" + MY_NEWLINE + \
                       "#       with a MFE cut off of <" + str(bindingEnergyCutoff) + ">" + MY_NEWLINE + \
                       "#       and a prediction probability below < " + str(positiveProbabilityCutoff) + ">"+ MY_NEWLINE


    headerString = headerString + "#             " + MY_NEWLINE

    filterLine = "__e" + str(bindingEnergyCutoff) + "_p_" + str(positiveProbabilityCutoff)
    filteredTargetPredictionFile = miRAWpositiveSitesFile.replace(".csv", filterLine + ".csv")
    r = 0
    with open(filteredTargetPredictionFile, 'w') as fT:
        fT.write(headerString + MY_NEWLINE)
        #fT.write(HEADER_LINE + MY_NEWLINE)
        fT.write(positiveTargetsHeaderLine )
        while r < len(dropPredictions):
            if dropPredictions[r] == 1:
                fT.write(gGeneName[r] + "\t" + \
                         gMiRNA[r] + "\t" + \
                         gSiteStart[r] + "\t" + \
                         gSiteEnd[r] + "\t" + \
                         gPrediction[r] + "\t" + \
                         gPairStartInSite[r] + "\t" + \
                         gSeedStart[r] + "\t" + \
                         gSeedEnd[r] + "\t" + \
                         gPairs[r] + "\t" + \
                         gWC[r] + "\t" + \
                         gWob[r] + "\t" + \
                         gMFE[r] + "\t" + \
                         gCanonical[r] + "\t" + \
                         gSiteTranscript5to3[r] + "\t" + \
                         gMatureMiRNATranscript[r] + "\t" + \
                         gBracketNotation[r] + "\t" + \
                         gFiltered[r] + "\t" + \
                         gFiltering[r] + "\t" + \
                         gReason[r] + "\t" + MY_NEWLINE)

            r += 1
    logging.info("--finished")


def checkArgs():
    
    if args.HelpMe :
        printLongHelpAndExit() 
        exit()
        
    checkTargetPredictionFile()
    checkPositiveTargetsFile()
    checkUnifiedInputFile()
    checkRemovemiRNAs()
    checkExtractConflicts()
    checkBindingEnergyCutOff()
    checkPositiveCutOff()
    checkNegativeCutOff()
    checkExtractPositives()
    checkExtractNegatives()
    checkForBadFlagCombinations()


def filterData():
    global extractConflictPredictions, keepConflictPredictions
    logging.info("filterData")
    #removeEmptyPredictions()
    if filterPredictionsByList:
        readFilterFile()
        filterByList()
    # we have removed conflict processing for now.
    # You can have a miRNA targeting different regions of the
    # same 3'UTR with positive and negative predictions,
    # but the model will never predict both
    # for the same miRNA/mRNA fragment.
    #if extractConflictPredictions:
    #    filterByConflicts()
    #if keepConflictPredictions:
    #    keepByConflicts()
    #if (positiveProbabilityCutoff > 0) or (negativeProbabilityCutoff < 0):
    #filterByProbability()
    filterByMFEAndProbability()
    # if removePositivePredictions:
    #    filterByPositive()
    #if extractNegativePredictions:
    #    filterByNegative()
    summarizeFiltering()


def main(argv=None):  # IGNORE:C0111

    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()

    parseArgs(argv)
    checkArgs()
    initLogger(md5String)
    readUnifiedFile()
    readPositiveSitesFile()
    readTargetPredictionFile()

    filterData()
    writeFilteredPositiveData()
    #writeFilteredDetailedData()

if __name__ == "__main__":
    if DEBUG:
        pass
        # sys.argv.append("-h")
        # sys.argv.append("-v")

    if TESTRUN:
        import doctest

        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats

        profile_filename = 'fairpype.virusPipe_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())



