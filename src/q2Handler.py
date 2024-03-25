#!/usr/bin/env python3

import qiime2
from qiime2 import plugins
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
import psutil
import urllib
import distutils
import IPython
import tempfile
import sys
import time

sns.set()

for plugin in qiime2.plugins.available_plugins():
    qiime2.plugins.importlib.import_module(plugin)

# This function can render .qzv files within a CraneOnDemand JupyterLab environment.
# It creates a temporary directory in the current working directory, copies the
# .qzv file to that directory, then renders the file with an IPython Iframe.
# The code is kinda sketchy, so it might not work in the future if Qiime2 decides to
# radically alter their module layout. Though the general priniciple of rendering
# .qzv files should be the same. In other words, only minor tweaks required (hopefully)
def qzvRender(qzvFile, dimen=(1200, 600), sleepTime=3):
    with tempfile.TemporaryDirectory(dir=pathlib.Path.cwd()) as tempDir:
        distutils.dir_util.copy_tree(qzvFile._archiver.path, tempDir)
        relDir = pathlib.PurePath(tempDir).relative_to(pathlib.Path.cwd())
        #copiedDir = list(pathlib.Path(relDir).iterdir())[0]
        indexPath = pathlib.Path(relDir.joinpath('data/index.html')).as_posix()
        iframe = IPython.display.IFrame(indexPath, dimen[0], dimen[1])
        IPython.display.display(iframe)
        
        # Time in seconds for the temporary directory to exist. After this period,
        # new pages opened by the IFrame could return with a 404 error. The cell
        # can be prematurely stopped before the time is up and the temporary directory
        # will clean itself. You can think of sleepTime as the amount of time you want
        # the imaginary server rendering the image to stay active.
        time.sleep(sleepTime) 

# This function collapses a taxonomy file into individual levels with the
# Qiime2 collapse function. It returns a dict of 7 Qiime2 FeatureTable[Frequency] artifacts,
# with each level in the dict (e.g., level_1) containing a FeatureTable collapsed at
# that taxonomic level
def collapseTaxa(tableFile, taxonomyFile):
    if isinstance(tableFile, str):
        table = qiime2.Artifact.load(tableFile)
        
    if isinstance(taxonomyFile, str):
        taxonomy = qiime2.Artifact.load(taxonomyFile)
        
    collapsedTaxa = {}
    
    for i in range(0, 7):
        level = i + 1
        key = 'level_' + str(level)
        collapsed_table = qiime2.plugins.taxa.methods.collapse(
            table=table,
            taxonomy=taxonomy,
            level=level,
        )
        collapsedTaxa.update({key: collapsed_table.collapsed_table})        
        
    return collapsedTaxa


# ### Save Handling ###
# # This takes a list of tuples, with each individual tuple in the list consisting of:
# # (0): The Qiime2 object to save; (1): The desired name for the saved object
# # (2): (OPTIONAL) saveDir, a directory where the Qiime2 object will be saved
# # Finally, an optional directory string can be passed for all objects to be saved
# # in the same directory
# def saveHandler (saveTup, singleDir=''):
#     for tup in saveTup:
#         qiimeObj = tup[0]
#         saveName = tup[1]

#         if len(tup) == 3 and not singleDir:
#             dirName = tup[2]
#         elif singleDir:
#             dirName = singleDir
#         else:
#             dirName = currDir

#         saveDir = pathlib.Path(dirName)
#         if not saveDir.exists():
#             pathlib.Path.mkdir(saveDir, parents=True)

#         qiimeObj.save(saveDir.joinpath(saveName).as_posix())


# This function takes four required parameters:
# metadataFile: A path to a metadataFile
# dataDir: A path to a directory with sequencing data
# manifestDir: A path to a directory for the manifest to be placed
# manifestStyle: Either 'bash' or 'python', depending on whether Qiime2
# is called from the CLI or the python API
def manifestGen(metadataFile='', dataDir='', manifestDir='', manifestStyle=''):
    currDir = pathlib.Path.cwd()
    metadata = pathlib.Path(metadataFile)
    dataDir = pathlib.Path(dataDir)
    manifestDir = pathlib.Path(manifestDir)

    # This takes an input metadata file. It creates a .tsv file from the metadata
    # if one doesn't already exist (in the current working directory), 
    # then creates a manifest file from the metadata (also in the cwd)
    if metadata.is_absolute() or metadata.exists():
        if dataDir.is_absolute() or dataDir.exists():
           
            if 'xls' in metadata.suffix:
                df = pd.DataFrame(pd.read_excel(metadata))
                tsvFile = manifestDir.joinpath(manifestDir, (metadata.stem + '.tsv'))
                df.to_csv(tsvFile, sep = '\t', index = False)

            elif 'tsv' in metadata.suffix:
                df = pd.DataFrame(pd.read_csv(metadata.as_posix(), delimiter='\t'))

            elif 'csv' in metadata.suffix:
                df = pd.DataFrame(pd.read_csv(metadata, delimiter=','))
                tsvFile = manifestDir.joinpath(manifestDir, (metadata.stem + '.tsv'))
                df.to_csv(tsvFile, sep = '\t', index = False)

            # This removes a header that was on the file. Can be easily changed
            sampIds = tuple(df['sample-id'].astype(str))

            manifest = {}
            for sid in sampIds:
                reads = []

                # If dataDir doesn't have any sequences, the program tries 
                # searching through any directories in dataDir instead
                # (recursive search)
                if not len([match for match in dataDir.glob('*' + sid + '*')]):
                    for direc in (item for item in dataDir.iterdir() if item.is_dir() == True):
                        seqFiles = [match for match in direc.glob('*' + sid + '*')]
                        if seqFiles:
                            for seqFile in seqFiles:
                                reads.append(seqFile)

                else:
                    seqFiles = [match for match in dataDir.glob('*' + sid + '*')]
                    for seqFile in seqFiles:
                        reads.append(seqFile)

                manifest.update({sid: reads})

            # The manifest file is written. This should (in theory) be able
            # to handle both single and paired-end reads. Additionally,
            # because the qiime2 API and CLI use different manifest formats,
            # both are supported.
            with manifestDir.joinpath(metadata.stem + '_' + manifestStyle + '_manifest').open(mode='w') as f:
                # Apparently the qiime2 API can only read .csv manifests 
                # instead of .tsv manifests for python
                if manifestStyle == 'python':
                    f.write('sample-id,absolute-filepath,direction\n')
                    for key, values in manifest.items():
                        r1 = ''
                        r2 = ''
                        for value in values:
                            if value.match('*R1*.fastq*'):
                                r1 = value.as_posix()
                            else:
                                r2 = value.as_posix()

                        f.write(key + ',' + r1 + ',' + 'forward' + '\n') 

                        if r2:
                            f.write(key + ',' + r2 + ',' + 'reverse' + '\n') 

                elif manifestStyle == 'bash':
                    if len(manifest[sampIds[0]]) == 1:
                        f.write('sample-id\tabsolute-filepath\n')
                        for key, values in manifest.items():
                            f.write(key + '\t' + values[0].as_posix() + '\n')
                    else:
                        f.write('sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n')

                        for key, values in manifest.items():
                            r1 = ''
                            r2 = ''
                            for value in values:
                                if value.match('*R1*.fastq*'):
                                    r1 = value.as_posix()
                                else:
                                    r2 = value.as_posix()

                            f.write(key + '\t' + r1 + '\t' + r2 + '\n')

        else:
            print('Please path the absolute path or a valid relative path to ' +
                'a directory with sequencing data')
    else:
        print('Please pass the absolute path or a valid relative path to the metadata file')
