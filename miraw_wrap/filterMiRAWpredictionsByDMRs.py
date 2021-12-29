#!/usr/local/bin/python3.7
# encoding: utf-8
'''


@author:     Simon Rayner

@copyright:  2021 Oslo University Hospital. All rights reserved.

@license:    license

@contact:    simon.rayner@medisin.uio.no
@deffield    updated: Updated
'''

import sys
import csv
import os
from pathlib import Path

from datetime import datetime
import hashlib
import logging

import pandas as pd
import numpy as np

from plotnine import *
from plotnine.data import *


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from Bio.Data.CodonTable import list_possible_proteins


__all__ = []
__version__ = 0.1
__date__ = '2020-12-19'
__updated__ = '2020-12-19'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

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
    projectBaseName = os.path.splitext(os.path.basename(groupedpredsFile))[0]
    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)   
    logfileName = os.path.join(logFolder, projectBaseName + "__" + dt_string + "__" + md5string +".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)
    
    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)
    
    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler 
    log.addHandler(handler) 
    logging.info("+" + "*"*78 + "+")   
    logging.info("project log file is <" + logfileName + ">")         
    logging.info("+" + "*"*78 + "+")   
    logging.debug("debug mode is on")
    

def parseArgs(argv):
    
    '''parse out Command line options.'''

    global parser
    
    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
    i
      Created by Simon Rayner on %s.
      Copyright 2021 Oslo University Hospital. All rights reserved.
    
      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0
    
      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.
    
    USAGE
    ''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-g", "--groupedpredsfile", dest="groupedpredsfile", action="store", 
                            help="set of grouped and filtered predictions generated by filterAndPoolMiRAWpredictions.py [default: %(default)s]")
        parser.add_argument("-m", "--mirbasegff3file", dest="miRBasegff3file", action="store", 
                            help="list of features to select from the predictions [default: %(default)s]")
        parser.add_argument("-u", "--upregulated", dest="upregulated", action="store", 
                            help="list of upregulated proteins [default: %(default)s]")
        parser.add_argument("-d", "--dmrposfile", dest="dmrposfile", action="store", 
                            help="list of upregulated proteins [default: %(default)s]")
        parser.add_argument("-D", "--featuredistance", dest="featuredistance", action="store", 
                            help="minimum distance between miRNA and DMR in nucleotides [default: %(default)s]")

        parser.add_argument("-H", "--HelpMe", action="store_true", 
                            help="print detailed help")

        # Process arguments
        args = parser.parse_args()
        
        if args.HelpMe :
            printLongHelpAndExit() 
            exit()        
        
        
        global groupedpredsFile 
        global miRBaseGFF3File
        global upregulatedProtFile
        global dmrPosFile
        global featureDistance

        
        
        if args.groupedpredsfile:
            groupedpredsFile = args.groupedpredsfile
            print("grouped miRAW prediction file is <" + groupedpredsFile + ">")
        else:
            print("----you need to specify a grouped miRAW prediction file using the -g/--groupedpredsfile parameter") 
            printHelpAndExit()        
        
        
        if args.miRBasegff3file:
            miRBaseGFF3File = args.miRBasegff3file
            print("miRBase GFF3 file is <" + miRBaseGFF3File + ">")
        else:
            print("----you need to specify a miRBase GFF3 file using the -m/--mirbasegff3file parameter") 
            printHelpAndExit()        


        if args.dmrposfile:
            dmrPosFile = args.dmrposfile
            print("DMR position file is <" + dmrPosFile + ">")
        else:
            print("----you need to specify a DMR position file using the -D/--dmrposfile parameter") 
            printHelpAndExit()        

            
        if args.featuredistance:
            featureDistance = args.featuredistance
            print("minimum feature spacing set to <" + featureDistance + ">")
        else:
            print("----you need to specify a minimum feature spacing using the -D/--featuredistance parameter") 
            printHelpAndExit()        
            
            

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



def loadDMRFeatureList():

    global dfDMRFeatureList   
        
    logging.info("loading DMR feature list file <" + dmrPosFile + ">")
    dfDMRFeatureList = pd.read_csv(dmrPosFile, sep="\t", skiprows=3)
    
    logging.info("read <" + str(len(dfDMRFeatureList)) + "> features")



 
def loadMiRBaseFeatureList():

    global miRBaseFeatureList

        
    logging.info("loading DMR feature list file <" + miRBaseGFF3File + ">")
    miRBaseFeatureList = pd.read_csv(miRBaseGFF3File, sep="\t", skiprows=14)
    miRBaseFeatureList.columns = ["chr", "col2", "featureType", "featureStart", 
                                  "featureStop", "col6", "strand", "col8", "Attributes"]
    miRBaseFeatureList=miRBaseFeatureList[miRBaseFeatureList['featureType'] == "miRNA"]
    miRBaseFeatureList["MIMATID"]=miRBaseFeatureList["Attributes"].str.split(";").str[0].str.split("=").str[1]
    
    # we only keep the features that encode miRNAs
    

    logging.info("read <" + str(len(miRBaseFeatureList)) + "> miRNA features")

    
    
def loadGroupPredictionData():
    
    # Format
    #      shortGeneName             shortmiRName                    count
    # 0    ENSG00000005339|CREBBP    hsa-let-7a-2-3p|MIMAT0010195    8
    # 1    ENSG00000005339|CREBBP    hsa-let-7a-3p|MIMAT0004481      3
    
    
    global grpPredData
    # isos/karlsen_isos.MIMAT0000087_hsa-miR-30a-5p_GTAAACA_iso/allTargetSites.csv
    logging.info("loading group prediction data")
    
    grpPredData=pd.read_csv(groupedpredsFile, sep='\t')
    grpPredData['MIMATID']=grpPredData['shortmiRName'].str.split("|").str[1]
        

    logging.info("found <" + str(len(grpPredData))  + "> predictions")
        

def mergeMiRNAData():
    #result = pd.merge(left, right, on=["key1", "key2"])
    global dfFullMiRInfo
    dfFullMiRInfo = pd.merge(grpPredData, miRBaseFeatureList, on=["MIMATID", "MIMATID"])


    
def processSamples():
    # read target list
    # loop through target list
    global targetFiles
    global allPreds
    
    allPreds = pd.DataFrame()    
    
    logging.info("processing target file")
             
             
    # format the filtered predictions for network analysis    
    
    # 1. number of target events / miRNA / 3'UTR
    # Gene name has five parts: ENSG00000106331|ENST00000338516|PAX4|127610811|127610845
    #   Keep ENSG ID and Common Name
    # miRNA has four parts:     hsa-miR-548t-3p_MIMAT0022730_homo_sapiens_miR-548t-3p
    #   keep Common Name and MIMAT id
    allPreds['shortGeneName']=allPreds['GeneName'].str.split("|").str[0] + "|" + allPreds['GeneName'].str.split("|").str[2]
    allPreds['shortmiRName']=allPreds['miRNA'].str.split("_").str[0]+ "|" + allPreds['miRNA'].str.split("_").str[1]    
    groupedAllPreds = allPreds.groupby(['shortGeneName','shortmiRName']).size().reset_index().rename(columns={0:'count'})
    outputFileGrpAllPreds = os.path.splitext(groupedpredsFile)[0] + "_allfilteredGrouped.tsv" 
    groupedAllPreds.to_csv(outputFileGrpAllPreds, sep='\t')
    
    
    # 2. number of targeted 3'UTRs / miRNA
    countsByMiRNAs = groupedAllPreds['shortmiRName'].value_counts()
    dfCounts = pd.DataFrame(countsByMiRNAs)
    dfCounts.columns=['counts']
    outputFileCountsByMiRNAs = os.path.splitext(groupedpredsFile)[0] + "_countsByMiRNAs.tsv" 
    
    count, division = np.histogram(countsByMiRNAs, bins=list(range(0, max(countsByMiRNAs)+2)))
    
    plotHistogramFileCountsByMiRNAs = os.path.splitext(groupedpredsFile)[0] + "_countsByMiRNAs.png"
    countsByMiRNAs.to_csv(outputFileCountsByMiRNAs, sep='\t')

    dataHistogramFileCountsByMiRNAs = os.path.splitext(groupedpredsFile)[0] + "_histByMiRNAs.tsv"
    dfHistmiR = pd.DataFrame([division, count]).T
    dfHistmiR.columns=['bin', 'count']
    dfHistmiR.replace(np.nan, 0)
    dfHistmiR.to_csv(dataHistogramFileCountsByMiRNAs, sep='\t')

    p = ggplot(dfCounts, aes(x='counts')) + geom_histogram(binwidth=1, color="black", fill="white")    
    p.save(filename = plotHistogramFileCountsByMiRNAs, height=5, width=5, units = 'in', dpi=1000)    
    
    # 3. number of miRNAs / 3'UTR
    countsBy3pUTRs = groupedAllPreds['shortGeneName'].value_counts()
    outputFileCountsBy3pUTRs = os.path.splitext(groupedpredsFile)[0] + "_countsBy3pUTRs.tsv" 
    countsBy3pUTRs.to_csv(outputFileCountsBy3pUTRs, sep='\t')
    
    
    
    
    
    

             
def processPredictionFile(predFile, miRName):
    global allPreds
    # grab miRNA name from parent folder name

    # read target predictions and filter by probability and energy

    
    logging.info("done")
          
    
    
def printLongHelpAndExit():
    logging.info("+" + "-" * 78 + "+")
    logging.info("+                                                                              +")
    logging.info("+  filter and pool miRAW prediction results:                                   +")
    logging.info("+                                                                              +")
    logging.info("+    This is a work in progress that is evolving as new filtering needs        +")
    logging.info("+    are identified.  Currently, the user can do the following                 +")
    logging.info("+                                                                              +")
    logging.info("+      Specify a list of miRAW predictions, with each file corresponding       +")    
    logging.info("+      to a miRNA, and a second list of genes to be selected from the          +")
    logging.info("+      prediction set                                                          +")
    logging.info("+                                                                              +")
    logging.info("+      The complete set of filtered predictions are written as a single        +")
    logging.info("+      TSV file that can be used to generate network graphs                    +")            
    logging.info("+                                                                              +")
    logging.info("+    you need to specify:                                                      +")
    logging.info("+                                                                              +")
    logging.info("+      a result file: (-r/--resultsfile)                                       +")
    logging.info("+         This file contains a list of miRAW prediction files                  +")
    logging.info("+         (these are the output files ending with 'allTargetSites.csv')        +")
    logging.info("+                                                                              +")
    logging.info("+                                                                              +")
    logging.info("+    you also need to specify at least one of the following                    +")
    logging.info("+    If you don't, the output will be the same unfiltered input data           +")
    logging.info("+                                                                              +")
    logging.info("+           --probability                                                      +")
    logging.info("+                  cut-off for min probability                                 +")
    logging.info("+                                                                              +")
    logging.info("+           --energy                                                           +")
    logging.info("+                 cut-off for min binding energy                               +")
    logging.info("+                                                                              +")
    logging.info("+           --upregulated                                                      +")
    logging.info("+                a list of up-regulated proteins                               +")
    logging.info("+                  from matching experimental data                             +")
    logging.info("+                                                                              +")
    logging.info("+           --downregulated                                                    +")
    logging.info("+                a list of down-regulated proteins                             +")
    logging.info("+                  from matching experimental data                             +")
    logging.info("+                                                                              +")
    logging.info("+                                                                              +")
    logging.info("+      The output files will be                                                +")
    logging.info("+                                                                              +")    
    logging.info("+                                                                              +")
    logging.info("+" + "-" * 78 + "+")        
    
    
def printHelpAndExit():
    parser.print_help()
    logging.info("stopping") 
    sys.exit()
    
    

def main(argv=None): # IGNORE:C0111
    
    
    #if argv is None:
    #    argv = sys.argv

    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    parseArgs(argv)
    initLogger(md5String)
    loadDMRFeatureList()
    loadGroupPredictionData()
    loadMiRBaseFeatureList()
    mergeMiRNAData()
    processSamples()
    

    logging.info("finished")

if __name__ == "__main__":
    if DEBUG:
        pass
        #sys.argv.append("-h")
        #sys.argv.append("-v")

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
    

    
    
    
