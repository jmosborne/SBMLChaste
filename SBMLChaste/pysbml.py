#!/usr/bin/env python
# import the necessary libraries
import sys
import os
from os.path import basename
from libsbml import *
import translator
import filewriters

def IsSrnModel(path_to_file):
    """ Function to be used in another pipeline (inherits from translator)."""
    reader = SBMLReader()
    document = reader.readSBMLFromFile(path_to_file)
    model = document.getModel()

    return translator.IsSrnModel(model)

def GetOutputFilenames(path_to_file, filename):
    """ Function to get the names of the translated files that need to be moved. Returns in order
    [header file, source file]. """
    outputfiles = []

    #Get the model from the file
    reader = SBMLReader()
    document = reader.readSBMLFromFile(path_to_file)
    model = document.getModel()

    if( translator.IsSrnModel(model) ):
        header_file = filename + "SrnModel.hpp"
        source_file = filename + "SrnModel.cpp"

        outputfiles.append(header_file)
        outputfiles.append(source_file)
    else:
        header_file = filename + "CellCycleModel.hpp"
        source_file = filename + "CellCycleModel.cpp"

        outputfiles.append(header_file)
        outputfiles.append(source_file)

    return outputfiles

def MoveFilesToDirectory(files, new_path):
    """ Moves output files to specified directory. """
    #Get current path
    current_path = os.getcwd()

    for i in range(len(files)):
        os.rename(current_path + "/" + files[i], new_path + "/" + files[i])

def main(args):
    """ Main function. Arguments are:
    Input, Output Filename, Output Directory. """

    path_to_file = args[1]

    #Get the model from the file
    reader = SBMLReader()
    document = reader.readSBMLFromFile(path_to_file)
    model = document.getModel()

    #Get Name of file
    if(len(args) > 2):
        model_name = args[2]
    else:
        filename = basename(path_to_file)
        split_fname = filename.split('.')
        model_name = split_fname[0]

    # If there are events in the model, we define the class as a cell cycle model.
    # Otherwise, the model is assumed to be a subcellular reaction network model.
    if( translator.IsSrnModel(model) ):
        filewriters.WriteSrnModelToFile(model_name, model)
    else:
        filewriters.WriteCcmModelToFile(model_name, model)

    if(len(args) > 3):
        output_directory = args[3]
    else:
        output_directory = "."

    output_files = GetOutputFilenames(path_to_file, model_name)

    #Move created files if output directory has been specified.
    MoveFilesToDirectory(output_files, output_directory)
		
if __name__ == '__main__':
  main(sys.argv) 