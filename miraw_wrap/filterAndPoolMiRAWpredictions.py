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
    projectBaseName = os.path.splitext(os.path.basename(resultfiles))[0]
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
      Copyright 2020 Oslo University Hospital. All rights reserved.
    
      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0
    
      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.
    
    USAGE
    ''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-r", "--resultsfile", dest="resultfiles", action="store", help="list of result files to process [default: %(default)s]")
        parser.add_argument("-f", "--featurelist", dest="featurelist", action="store", help="list of features to select from the predictions [default: %(default)s]")
        parser.add_argument("-u", "--upregulated", dest="upregulated", action="store", help="list of upregulated proteins [default: %(default)s]")
        parser.add_argument("-d", "--downregulated", dest="downregulated", action="store", help="list of upregulated proteins [default: %(default)s]")
        parser.add_argument("-p", "--probability", dest="probability", action="store", help="cut-off for min probability [default: %(default)s]")
        parser.add_argument("-e", "--energy", dest="energy", action="store", help="cut-off for min binding energy [default: %(default)s]")
        parser.add_argument("-H", "--HelpMe", action="store_true", help="print detailed help")

        # Process arguments
        args = parser.parse_args()
        
        if args.HelpMe :
            printLongHelpAndExit() 
            exit()        
        
        
        global resultfiles 
        global featurelistFile
        global upregulatedProtFile
        global downregulatedProtFile
        global probability
        global energy 
        
        
        if args.resultfiles:
            resultfiles = args.resultfiles
            print("miRAW results file list is <" + resultfiles + ">")
        else:
            print("----you need to specify an miRAW results filelist using the -r/----resultfiles parameter") 
            printHelpAndExit()        
        
        
        if args.featurelist:
            featurelistFile = args.featurelist
            print("feature list file is <" + featurelistFile + ">")
        else:
            featurelistFile = "NONE"
            print("no feature list file specified ")


        if args.downregulated:
            downregulatedProtFile = args.downregulated
            print("up-regulated Protein file is <" + downregulatedProtFile + ">")
        else:
            downregulatedProtFile = "NONE"
            print("no down-regulated Protein file specified ")

        if args.upregulated:
            upregulatedProtFile = args.upregulated
            print("up-regulated Protein file is <" + upregulatedProtFile + ">")
        else:
            upregulatedProtFile = "NONE"
            print("no up-regulated Protein file specified ")

            
        if args.probability:
            probability = args.probability
            print("minimum probability threshold set to <" + probability + ">")
        else:
            probability = 0.0
            print("no minimum probability threshold specified")
            
            
        if args.energy:
            energy = args.energy
            print("minimum energy threshold set to <" + energy + ">")
        else:
            probability = 0.0
            print("no minimum energy threshold specified")
            

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



def loadFeatureList():

    global dfFeatureList
    global featureList
    
    featureList = []
    if featurelistFile == "NONE":
        return
        
    logging.info("loading feature list file <" + featurelistFile + ">")
    with open(featurelistFile) as fFL:
        featureList = fFL.read().splitlines()
    
    logging.info("read <" + str(len(featureList)) + "> features")

    
def loadProteinData():

    global dfMiRsVsProteins
    global proteins
    upProteins = []
    downProteins = []
    proteins = []
    
    if not upregulatedProtFile == "NONE":
        
        logging.info("loading UP protein file <" + upregulatedProtFile + ">")
        with open(upregulatedProtFile) as fUP:
            upProteins = fUP.read().splitlines()
        logging.info("done")
        
    if not downregulatedProtFile == "NONE":
        
        logging.info("loading DOWN protein file <" + upregulatedProtFile + ">")
        with open(downregulatedProtFile) as fDOWN:
            downProteins = fDOWN.read().splitlines()
        logging.info("done")
        
    dfMiRsVsProteins = pd.DataFrame(proteins)
           
    if upregulatedProtFile == "NONE" and downregulatedProtFile == "NONE":     
        return 
                   
        # build data frame, one row / protein
        proteins = upProteins + downProteins
        dfMiRsVsProteins = pd.DataFrame(proteins)
        dfMiRsVsProteins = dfMiRsVsProteins.set_index(dfMiRsVsProteins[0])
        miRNames=[]
        miRNames = miRNAs['miRName'].tolist()
        dfMiRsVsProteins=pd.concat([dfMiRsVsProteins,pd.DataFrame(columns=miRNames)])
        dfMiRsVsProteins = dfMiRsVsProteins.drop([0], axis=1)
        logging.info("done")
    
def loadMiRAWSampleList():
    
    # INDEX    CHANGE    FILE
    # 49    D    /Users/simonray/Dropbox/DResearch/jan_brinchmann/miRAW_targets/full_targetset_76/karlsen_isos.MIMAT0000510_hsa-miR-320a-3p_AAAGCTG_iso/allTargetSites.csv


    
    #
    #  File has 'allTargetSites.csv' ending 
    #  and the following columns:
    # 
    #    FIELD                    SAMPLE VALUE
    #    GeneName                 ENSG00000003137|ENST00000001146|CYP26B1|72129238|72132226    
    #    miRNA                    MIMAT0000087_hsa-miR-30a-5p_GTAAACA_iso    
    #    SiteStart                2879    
    #    SiteEnd                  2919    
    #    Prediction               -0.9872640371322632    
    #    Filtered                 0    
    #    PostfilterPrediction     -0.9872640371322632    
    #    PairsInSeed              7    
    #    FreeEnergy               -11.1    
    #    SiteTranscript           ATAATATTCCATGTTGCATATTAAAAACATGAATGTTGTG    
    #    MatureMiRNATranscript    GTAAACATCCTCGACTGGAAG    
    #    Filtering Reason        
    #    Canonical                0    
    #    Additional properties        
    #
    
    global targetFiles
    global miRNAs
    # isos/karlsen_isos.MIMAT0000087_hsa-miR-30a-5p_GTAAACA_iso/allTargetSites.csv
    logging.info("processing target files")
    targetFiles=[]
    
    miRNAs=pd.read_csv(resultfiles, sep='\t')
    miRNAs['miRName']='unknown'
    for index, row in miRNAs.iterrows():
        miRNAs.at[index, 'miRName'] = row['CHANGE'] + "__" + str(row['INDEX']) + "__" + os.path.basename(Path(row['FILE']).parents[0]).split('.')[1]

        
    
    #with open(resultfiles, 'r') as fd:
    #    csvrTargetFiles = csv.reader(fd)
    #    for targetFile in csvrTargetFiles:
    #        targetFiles.append(targetFile[0])
    #        miRNames.append(os.path.basename(Path(targetFile[0]).parents[0]).split('.')[1])
    logging.info("found <" + str(len(miRNAs))  + "> samples")
        

def processSamples():
    # read target list
    # loop through target list
    global targetFiles
    global allPreds
    
    allPreds = pd.DataFrame()    
    
    logging.info("processing target files")
    for index, row in miRNAs.iterrows():
        logging.info("-- processing file <" + miRNAs.at[index, 'FILE'] + ">")
        processPredictionFile(miRNAs.at[index, 'FILE'], miRNAs.at[index, 'miRName'])
             
    
    proteinsOutputFile = os.path.splitext(resultfiles)[0] + "_proteins.tsv"       
    dfMiRsVsProteins.to_csv(proteinsOutputFile, sep='\t')
    
    outputFileAllPreds = os.path.splitext(resultfiles)[0] + "_allfiltered.tsv" 
    allPreds.to_csv(outputFileAllPreds, sep='\t')
             
             
    # format the filtered predictions for network analysis    
    
    # 1. number of target events / miRNA / 3'UTR
    # Gene name has five parts: ENSG00000106331|ENST00000338516|PAX4|127610811|127610845
    #   Keep ENSG ID and Common Name
    # miRNA has four parts:     hsa-miR-548t-3p_MIMAT0022730_homo_sapiens_miR-548t-3p
    #   keep Common Name and MIMAT id
    allPreds['shortGeneName']=allPreds['GeneName'].str.split("|").str[0] + "|" + allPreds['GeneName'].str.split("|").str[2]
    allPreds['shortmiRName']=allPreds['miRNA'].str.split("_").str[0]+ "|" + allPreds['miRNA'].str.split("_").str[1]    
    groupedAllPreds = allPreds.groupby(['shortGeneName','shortmiRName']).size().reset_index().rename(columns={0:'count'})
    outputFileGrpAllPreds = os.path.splitext(resultfiles)[0] + "_allfilteredGrouped.tsv" 
    groupedAllPreds.to_csv(outputFileGrpAllPreds, sep='\t')
    
    
    # 2. number of targeted 3'UTRs / miRNA
    countsByMiRNAs = groupedAllPreds['shortmiRName'].value_counts()
    dfCounts = pd.DataFrame(countsByMiRNAs)
    dfCounts.columns=['counts']
    outputFileCountsByMiRNAs = os.path.splitext(resultfiles)[0] + "_countsByMiRNAs.tsv" 
    
    plotHistogramFileCountsByMiRNAs = os.path.splitext(resultfiles)[0] + "_countsByMiRNAs.png"
    countsByMiRNAs.to_csv(outputFileCountsByMiRNAs, sep='\t')
    p = ggplot(dfCounts, aes(x='counts')) + geom_histogram(binwidth=1, color="black", fill="white")    
    p.save(filename = plotHistogramFileCountsByMiRNAs, height=5, width=5, units = 'in', dpi=1000)    
    
    # 3. number of miRNAs / 3'UTR
    countsBy3pUTRs = groupedAllPreds['shortGeneName'].value_counts()
    outputFileCountsBy3pUTRs = os.path.splitext(resultfiles)[0] + "_countsBy3pUTRs.tsv" 
    countsBy3pUTRs.to_csv(outputFileCountsBy3pUTRs, sep='\t')
    
    
    
    
    
    

             
def processPredictionFile(predFile, miRName):
    global allPreds
    # grab miRNA name from parent folder name

    # read target predictions and filter by probability and energy
    logging.info("miR <" + miRName + "> " )
    preds = pd.read_csv(predFile, sep='\t')
    logging.info("--read <" + str(len(preds)) + "> lines" )
    predsFilter = preds[(preds['FreeEnergy']<float(energy)) & (preds['Prediction']>float(probability))]
    logging.info("--after processing,  <" + str(len(predsFilter)) + "> lines remain" )
    
    allPreds = pd.concat([allPreds, predsFilter], axis=0)
    allPreds = allPreds.reset_index(drop=True)
    logging.info("--retained preds  <" + str(len(allPreds)) + "> lines remain" )
    
    
    # add up- and down-regulated proteins as columns to the dataframe
    if upregulatedProtFile == "NONE" and downregulatedProtFile == "NONE":     
        return
    
    for protein in proteins:
        #logging.info("----protein:" + protein)
        if len(predsFilter[predsFilter['GeneName'].str.contains(protein)==True].index) > 0:
            dfMiRsVsProteins[miRName][protein]=1
        else:
            dfMiRsVsProteins[miRName][protein]=0
    
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
    loadMiRAWSampleList()
    loadFeatureList()
    loadProteinData()
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
    

    
    
    
