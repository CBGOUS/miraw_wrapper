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
import os

from datetime import datetime
import hashlib
import logging

import pandas as pd
import numpy as np
from matplotlib import pyplot
from matplotlib import cm


from plotnine import *
from plotnine.data import *


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter



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
        parser.add_argument("-p", "--groupedpredsfile", dest="groupedpredsfile", action="store", 
                            help="set of grouped and filtered predictions generated by filterAndPoolMiRAWpredictions.py [default: %(default)s]")
        parser.add_argument("-g", "--mirbasegff3file", dest="miRBasegff3file", action="store", 
                            help="list of features to select from the predictions [default: %(default)s]")
        parser.add_argument("-3", "--3utrsnvfile", dest="threeutrsnvfile", action="store", 
                            help="list of identified 3'UTR SNVs [default: %(default)s]")
        parser.add_argument("-m", "--mirsnvfile", dest="mirsnvfile", action="store", 
                            help="list of identified miRNA SNVs [default: %(default)s]")
        parser.add_argument("-d", "--dmrposfile", dest="dmrposfile", action="store", 
                            help="list of DMR regions [default: %(default)s]")
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
        global threeUTRSNVFile
        global mirSNVFile

        
        
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
            print("----you need to specify a miRBase GFF3 file using the -g/--mirbasegff3file parameter") 
            printHelpAndExit()        


        if args.dmrposfile:
            dmrPosFile = args.dmrposfile
            print("DMR position file is <" + dmrPosFile + ">")
        else:
            print("----you need to specify a DMR position file using the -D/--dmrposfile parameter") 
            printHelpAndExit()        

        if args.threeutrsnvfile:
            threeUTRSNVFile = args.threeutrsnvfile
            print("3'UTR SNV file is <" + threeUTRSNVFile + ">")
        else:
            print("----you need to specify a 3'UTR SNV file using the -3/----3utrsnvfile parameter") 
            printHelpAndExit()        

        if args.mirsnvfile:
            mirSNVFile = args.mirsnvfile
            print("3'UTR SNV file is <" + mirSNVFile + ">")
        else:
            print("----you need to specify a miR SNV file using the -m/----mirsnvfile parameter") 
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

    global miRBaseMiRFeatureList
    global miRBasePreMiRFeatureList

        
    logging.info("loading DMR feature list file <" + miRBaseGFF3File + ">")
    miRBaseMiRFeatureList = pd.read_csv(miRBaseGFF3File, sep="\t", skiprows=14)
    miRBaseMiRFeatureList.columns = ["chr", "col2", "featureType", "featureStart", 
                                  "featureStop", "col6", "strand", "col8", "Attributes"]
    miRBasePreMiRFeatureList=miRBaseMiRFeatureList[miRBaseMiRFeatureList['featureType'] == "miRNA_primary_transcript"]
    miRBasePreMiRFeatureList["MIID"]=miRBaseMiRFeatureList["Attributes"].str.split(";").str[0].str.split("=").str[1]

    # we need the MIID to match the SNV information in 'mergeMiRNAData'
    miRBaseMiRFeatureList=miRBaseMiRFeatureList[miRBaseMiRFeatureList['featureType'] == "miRNA"]
    miRBaseMiRFeatureList["MIMATID"]=miRBaseMiRFeatureList["Attributes"].str.split(";").str[0].str.split("=").str[1]
    miRBaseMiRFeatureList["MIID"]=miRBaseMiRFeatureList["Attributes"].str.split(";").str[3].str.split("=").str[1]
    

    logging.info("read <" + str(len(miRBaseMiRFeatureList)) + "> miRNA features")

    
    
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
        
        
        
def load3pUTRSNVData():
    
    logging.info("loading 3'UTR SNV data from <" + threeUTRSNVFile + ">")    
    
    global df3pUTRSNVData
    df3pUTRSNVData = pd.read_csv(threeUTRSNVFile, sep="\t")
    df3pUTRSNVData["featureStart"]=pd.to_numeric(df3pUTRSNVData["LOC"].str.split(":").str[1].str.split("-").str[0])
    logging.info("found <" + str(len(df3pUTRSNVData)) + "> features")    
    
    logging.info("done")
    
    
    
def loadmiRNASNVData():
    
    logging.info("loading miRNA SNV data from <" + mirSNVFile + ">")    
    
    global dfmiRNASNVData
    dfmiRNASNVData = pd.read_csv(mirSNVFile, sep="\t")
    dfmiRNASNVData["featureStart"]=pd.to_numeric(dfmiRNASNVData["LOC"].str.split(":").str[1].str.split("-").str[0])
    logging.info("found <" + str(len(dfmiRNASNVData)) + "> features")    
    
    logging.info("done")
    
    
    
def mergeMiRNAData():

    global dfMiRsPredsPlusPos
    global dfFullPlusSNVs
    global dfFullMiRInfo

    # 1. merge the miRBase GFF file with the list of miRNAs present in the prediction list so we have
    #    have data frame with a list of all miRNAs present in the prediction list and their genome coordinates   
    uniqueMiRList = grpPredData["MIMATID"].unique()
    dfUniqMiRList = pd.DataFrame(uniqueMiRList)
    dfUniqMiRList.columns=["MIMATID"]
    dfMiRsPredsPlusPos = pd.merge(dfUniqMiRList, miRBaseMiRFeatureList, how="inner", on=["MIMATID", "MIMATID"])
    logging.info("found <" + str(len(dfMiRsPredsPlusPos)) + "> unique miRNAs in the prediction list")
    logging.info("done")
    
    # 2. merge the pre-miRNA SNV information and the miRBase genome position information.
    #    Because ENSEMBL doesn't include the miRBase MIMAT/MI reference IDs, it seems the only way to join
    #    is by first matching the start positions in both dataframes to get the position / SNV information at the
    #    pre-miRNA level
    
    dfPreMiRPosPlusSNVs = pd.merge(miRBasePreMiRFeatureList, dfmiRNASNVData, how="inner", on=["featureStart", "featureStart"])
    dfFullMiRInfo = pd.merge(dfPreMiRPosPlusSNVs, dfMiRsPredsPlusPos, how="inner", on=["MIID", "MIID"])
    
    # 3. merge the SNV information to the miRNAs by matching the MIID columns 
    
    # from this dataframe, i can use the 
    uniqueFullMiRList = dfFullMiRInfo["MIMATID"].unique()
    logging.info("after merging, we have  <" + str(len(uniqueFullMiRList)) + "> unique miRNAs in the prediction list")
    logging.info("done")
    
    # this information still needs to be merged into the grpPredData dataframe
    
    

    
def merge3pUTRData():
# merge the 3'UTR population SNV file with the list of 3'UTRs present in the prediction list so we have
# have data frame with a list of all 3'UTRs and their genome coordinates. 
    global dfFull3pUTRInfo
    unique3pUTRList = grpPredData["shortGeneName"].unique()
    dfunique3pUTRList = pd.DataFrame(unique3pUTRList)
    dfunique3pUTRList.columns=["shortGeneName"]
    dfunique3pUTRList["ENSGID"]= dfunique3pUTRList["shortGeneName"].str.split("|").str[0]

    dfFull3pUTRInfo = pd.merge(dfunique3pUTRList, df3pUTRSNVData, how="inner", on=["ENSGID", "ENSGID"])
    logging.info("found <" + str(len(dfFull3pUTRInfo)) + "> unique 3'UTRs in the prediction list")
    logging.info("done")
    
    
    
def processSamplesByDMR():


    #global targetFiles
    #global dmrHits
    allDMRHits=[]
    # loop through miRNAs in grpPredData list and see if they have an upstream DMR within "featuredistance" nts
    # to do this, iterate through every miRNA present in the predicted target list
    # (we have the coordinates through the mergeMiRNAData method)
    # 
    #uniqueMiRList = dfFullMiRInfo["MIMATID"].unique()
    for index, miR in dfFullMiRInfo.iterrows():        
        dmrHits = dfDMRFeatureList.loc[(dfDMRFeatureList['Chr'] == miR['chr_x']) 
                                       & (dfDMRFeatureList['End position']-miR['featureStart_x'] <= int(featureDistance))
                                       & (dfDMRFeatureList['End position']-miR['featureStart_x'] >= 0)]
        #dmrHits['featureDist']=dmrHits['End position']-miR['featureStart'] 
        if len(dmrHits) > 0:
            
            print(miR['MIMATID'] + "|" + str(len(dmrHits)))
            allDMRHits.append({"MIMATID":miR['MIMATID'], 
                               "number_of_hits":len(dmrHits), 
                               "distance":min(dmrHits['End position'])-miR['featureStart_x'],
                               "EUR":miR['EUR'],
                               "EAS":miR['EAS'],
                               "AMR":miR['AMR'],
                               "SAS":miR['SAS'],
                               "AFR":miR['AFR'],
                               "SNVsPerNT":miR['SNVsPerNT'],
                               "SubPopsPerNT":miR['SubPopsPerNT'],
                               "SupPopsPerNT":miR['SupPopsPerNT']}                            
)
    

    dfDMRHits = pd.DataFrame(allDMRHits)
    logging.info("found <" + str(len(dfDMRHits["MIMATID"].unique())) + "> miRNAs")
    logging.info("and <" + str(len(dfDMRHits["MIMATID"])) + "> DMR miRNA events")

    outputFiledfDMRHitss = os.path.splitext(groupedpredsFile)[0] + "_dmrHits.tsv" 
    dfDMRHits.to_csv(outputFiledfDMRHitss, sep='\t')

    # project the DMR perturbations on to the input prediction set
    # and generate a new list with the miRNAs in the dfDMRHits removed
    logging.info("summarizing DMR data")
    
    grpPredData['DMRcount']=0
    grpPredData['DMRdist']=0    
    grpPredData['EUR']=0    
    grpPredData['EAS']=0    
    grpPredData['AMR']=0    
    grpPredData['SAS']=0    
    grpPredData['AFR']=0    
    grpPredData['SNVsPerNT']=0.0
    grpPredData['SubPopsPerNT']=0.0  
    grpPredData['SupPopsPerNT']=0.0    
    
    for index, dmrMiR in dfDMRHits.iterrows():
        #print(dmrMiR)
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'DMRcount']=dmrMiR['number_of_hits']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'DMRdist']=dmrMiR['distance']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'EUR']=dmrMiR['EUR']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'EAS']=dmrMiR['EAS']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'AMR']=dmrMiR['AMR']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'SAS']=dmrMiR['SAS']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'AFR']=dmrMiR['AFR']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'SNVsPerNT']=dmrMiR['SNVsPerNT']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'SubPopsPerNT']=dmrMiR['SubPopsPerNT']
        grpPredData.loc[(grpPredData['shortmiRName'].str.contains(dmrMiR['MIMATID'])), 'SupPopsPerNT']=dmrMiR['SupPopsPerNT']
                
    outputFileallfilteredGroupedDMRmod = os.path.splitext(groupedpredsFile)[0] + "_allfilteredGroupedDMRmod.tsv" 
    grpPredData.to_csv(outputFileallfilteredGroupedDMRmod, sep='\t')
   
    
    # 1. number of targeted 3'UTRs / miRNA
    countsByMiRNAsNoDMR = grpPredData['shortmiRName'].value_counts()  # <-- this is DMRs / miRNA
    
    dfCountsNoDMR = pd.DataFrame(countsByMiRNAsNoDMR)
    dfCountsNoDMR.columns=['counts']
    outputFileCountsByMiRNAs = os.path.splitext(groupedpredsFile)[0] + "_countsByMiRNAsNoDMR.tsv"     
    countNoDMR, division = np.histogram(countsByMiRNAsNoDMR, bins=list(range(0, max(countsByMiRNAsNoDMR)+2)))    
    countsByMiRNAsNoDMR.to_csv(outputFileCountsByMiRNAs, sep='\t')

    countsByMiRNAsDMR = grpPredData.loc[grpPredData['DMRcount']==0]['shortmiRName'].value_counts()
    
    dfCountsDMR = pd.DataFrame(countsByMiRNAsDMR)
    dfCountsDMR.columns=['counts']
    outputFileCountsByMiRNAs = os.path.splitext(groupedpredsFile)[0] + "_countsByMiRNAsDMR.tsv"     
    countDMR, division = np.histogram(countsByMiRNAsDMR, bins=list(range(0, max(countsByMiRNAsNoDMR)+2)))    
    countsByMiRNAsDMR.to_csv(outputFileCountsByMiRNAs, sep='\t')

    countsByMiRNAsAFRSNVs = grpPredData.loc[grpPredData['AFR']==0]['shortmiRName'].value_counts()
    dfCountsAFRSNVs = pd.DataFrame(countsByMiRNAsAFRSNVs)
    dfCountsAFRSNVs.columns=['counts']
    outputFileAFRSNVsByMiRNAs = os.path.splitext(groupedpredsFile)[0] + "_countsByMiRNAsAFRSNVs.tsv"     
    countAFRSNVs, division = np.histogram(countsByMiRNAsAFRSNVs, bins=list(range(0, max(countsByMiRNAsAFRSNVs)+2)))    
    countsByMiRNAsAFRSNVs.to_csv(outputFileAFRSNVsByMiRNAs, sep='\t')


    dfHistCounts = pd.DataFrame(division[:-1])
    dfHistCounts["Counts - no DMR"] = countNoDMR.tolist()
    dfHistCounts["Counts - DMR"] = countDMR.tolist()
    dfHistCounts["Counts - AFR"] = countAFRSNVs.tolist()    
    
    dfHistCounts.columns=["no of 3UTR targets", "no DMR", "DMR", "AFRSNVs"]
    dfHistCounts.to_csv(os.path.splitext(groupedpredsFile)[0] + "_HistDMRNoDMR.tsv")
    statsFile = os.path.splitext(groupedpredsFile)[0] + "_statsDMRNoDMR.tsv"
    with open(statsFile, 'w') as fStats:
        fStats.write("no DMR: total targeting events" + "\t" + str(countsByMiRNAsNoDMR.sum()) + "\n")
        fStats.write("      : median" + "\t" + str(countsByMiRNAsNoDMR.median()) + "\n")
        fStats.write("      : average" + "\t" + str(countsByMiRNAsNoDMR.mean()) + "\n")
        fStats.write("   DMR: total targeting events" + "\t" + str(countsByMiRNAsDMR.sum()) + "\n")        
        fStats.write("      : median" + "\t" + str(countsByMiRNAsDMR.median()) + "\n")
        fStats.write("      : average" + "\t" + str(countsByMiRNAsDMR.mean()) + "\n")
        fStats.write("   AFR: total targeting events" + "\t" + str(countsByMiRNAsAFRSNVs.sum()) + "\n")        
        fStats.write("      : median" + "\t" + str(countsByMiRNAsAFRSNVs.median()) + "\n")
        fStats.write("      : average" + "\t" + str(countsByMiRNAsAFRSNVs.mean()) + "\n")
    plotHistogramFileCountsByMiRNAs = os.path.splitext(groupedpredsFile)[0] + "_countsByMiRNAsDMRmod.png"
    bins=list(range(0, max(countsByMiRNAsNoDMR)+2))

    pyplot.hist(countsByMiRNAsNoDMR, bins, alpha=0.3, label='noDMR', color='cornflowerblue', edgecolor="navy")
    pyplot.hist(countsByMiRNAsDMR, bins, alpha=0.3, label='DMR', color='plum', edgecolor="slateblue")
    pyplot.hist(countsByMiRNAsAFRSNVs, bins, alpha=0.3, label='AFR', color='green', edgecolor="mediumseagreen")
    pyplot.xlabel("no of targets")
    pyplot.ylabel("no of miRNAs")
    pyplot.title("connectivity DMR vs no DMR")
    pyplot.legend(loc='upper right')
    pyplot.savefig(plotHistogramFileCountsByMiRNAs)
    
      
    
    # 2. number of miRNAs / 3'UTR
    #countsBy3pUTRsNoDMR = grpPredData['shortGeneName'].value_counts()
    countsBy3pUTRDMR = grpPredData.loc[grpPredData['DMRcount']==0]['shortGeneName'].value_counts()
    outputFileCountsBy3pUTRs = os.path.splitext(groupedpredsFile)[0] + "_countsBy3pUTRsDMRmod.tsv" 
    countsBy3pUTRDMR.to_csv(outputFileCountsBy3pUTRs, sep='\t')


             
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
    load3pUTRSNVData()
    loadmiRNASNVData() 
       
    mergeMiRNAData()
    merge3pUTRData()
        
    processSamplesByDMR()
    

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
    

    
    
    
