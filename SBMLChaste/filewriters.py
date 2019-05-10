#!/usr/bin/env python

import sys
import os.path
from libsbml import *
import translator

# Script with functions that will write all the necessary SBML info into Chaste .cpp 
# and .hpp files. 

############################################################################################################
#####################            Functions common to both SRN and CCM        ###############################
############################################################################################################


def GetOdeSystemName(filename):
    """ Get the OdeSystem of the file. """

    ode_system_name = filename + "OdeSystem"

    return ode_system_name

def GetModelName(filename, model):
    """ Get the SrnModel of the file. """

    is_srn_model = translator.IsSrnModel(model)
    if(is_srn_model):
      model_name = filename + "SrnModel"
    else:
      model_name = filename + "CellCycleModel"

    return model_name

def IncludeString(name):
    """ String defining inclusion of file. """
    include_string = "#include " + name + "\n"

    return include_string

def GetModelDefinitionString(filename, model, forhpp):
    """ Returns string to define the Model using the SbmlSrnWrapperModel class or the SbmlCcmWrapperModel class, 
    depending on whether or not the model has any events. """
    #Get the model by reading from the filename

    #Size is defined by the number of ODEs
    odes_dict = translator.GetOdesDictionary(model)
    rules_dict = translator.GetRulesDictionary(model)
    num_species = len(odes_dict) + len(rules_dict)

    model_defn_string = ""

    ode_name = GetOdeSystemName(filename)
    model_name = GetModelName(filename, model)

    #Define some common names depending on whether or not the model is a SRN or a CCM
    is_srn_model = translator.IsSrnModel(model)
    if (is_srn_model):
      wrapper = "SbmlSrnWrapperModel"

      model_defn_string += translator.GetBlockCommentDefinition(0, "Define SRN model using Wrappers.", True)
    else:
      wrapper = "SbmlCcmWrapperModel"

      model_defn_string += translator.GetBlockCommentDefinition(0, "Define cell cycle model using Wrappers.", True)

    model_defn_string += (IncludeString("\"" + wrapper + ".hpp\"") + 
                          IncludeString("\"" + wrapper + ".cpp\"") + "\n" + 
                          "typedef " + wrapper + "<" + ode_name + "," + str(num_species) + "> " + model_name + ";\n" + 
                          "\n")

    model_defn_string += translator.GetCommentDefinition(0, "Declare identifiers for the serializer", True)
    if ( forhpp ):
        model_defn_string += IncludeString("\"SerializationExportWrapper.hpp\"")
    else:
        model_defn_string += IncludeString("\"SerializationExportWrapperForCpp.hpp\"")

    model_defn_string += ("CHASTE_CLASS_EXPORT(" + ode_name + ")\n" + 
                          "EXPORT_TEMPLATE_CLASS2(" + wrapper + ", " + ode_name + ", " + str(num_species) + ")\n" + 
                          "\n" + 
                          IncludeString("\"CellCycleModelOdeSolverExportWrapper.hpp\"") + 
                          "EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(" + model_name + ")\n\n" )

    return model_defn_string

## Header file functions ##

def GetHeaderFileDefinitionString(filename, model):
    """ Define the header file using #ifndef and #define. """

    is_srn_model = translator.IsSrnModel(model)

    if(is_srn_model):
      ode_and_model =  "ODESYSTEMANDSRNMODEL_HPP_"
    else:
      ode_and_model = "ODESYSTEMANDCELLCYCLEMODEL_HPP_"

    ifndef_str = "#ifndef " + filename.upper() + ode_and_model + "\n"
    define_str = "#define " + filename.upper() + ode_and_model + "\n"

    header_def_str = ifndef_str + define_str + "\n"

    return header_def_str

def GetHeaderFileClosingString(filename, model):
    """ Define the string to close the header file. """

    is_srn_model = translator.IsSrnModel(model)

    if(is_srn_model):
      ode_and_model = "ODESYSTEMANDSRNMODEL_HPP_"
    else:
      ode_and_model = "ODESYSTEMANDCELLCYCLEMODEL_HPP_"

    model_name_block_comment = translator.GetBlockCommentDefinition(0, filename.upper() + ode_and_model, False)

    header_close_string = "#endif " + model_name_block_comment

    return header_close_string

def GetIncludedFilesForHeaderString():
    """ Return string to include the necessary header files to go into the .hpp files. """

    # Don't really need to automate this as it'll be the same for all of them. 
    include_files_string = (IncludeString("\"ChasteSerialization.hpp\"") + 
                            IncludeString("<boost/serialization/base_object.hpp>") + 
                            IncludeString("<boost/serialization/shared_ptr.hpp>") + 
                            "\n" + 
                            IncludeString("<cmath>") + 
                            IncludeString("<iostream>") + 
                            IncludeString("\"AbstractOdeSystem.hpp\"\n") )

    return include_files_string

def GetFriendAndTemplateClassString():
    """This is the same for all of them, so we don't have to automate this."""

    #This essentially defines all the archiving of the ODE system.
    friend_and_template_string = (translator.AddTabs(1) + "friend class boost::serialization::access;\n" + 
                                  translator.AddTabs(1) + "template<class Archive>\n" + 
                                  translator.AddTabs(1) + "void serialize(Archive & archive, const unsigned int version)\n" + 
                                  translator.AddTabs(1) + "{\n" + 
                                  translator.AddTabs(2) + "archive & boost::serialization::base_object<AbstractOdeSystem>(*this);\n" + 
                                  translator.AddTabs(1) + "}\n" + 
                                  "\n")

    return friend_and_template_string

def GetPrivateVariablesForOdeClassString(model):
  """ Returns string defining private classes, including parameters and serialisation etc. """
  private_class_string = ("private:\n" + 
                          "\n" + 
                          translator.GetInitParameterStringForHeader(1, model) + #Initialise parameter values
                          "\n" + 
                          GetFriendAndTemplateClassString())

  return private_class_string

def GetConstructorAndDestructorForOdeClassString(ode_name):
    """ Return string defining public functions etc. """
    construct_string = (translator.GetBlockCommentDefinition(1, "Default constructor.", True) +  #Default constructor comment 
                        translator.AddTabs(1) + ode_name + "(std::vector<double> stateVariables=std::vector<double>());\n" + #Define ODEs system
                        "\n" + 
                        translator.GetBlockCommentDefinition(1, "Destructor.", True) +  #Destructor comment
                        translator.AddTabs(1) + "~" + ode_name + "();\n" + #Define destructor
                        "\n")
    return construct_string

def GetPublicVariablesForOdeClassString(filename, model):
  """ Returns string defining private classes, including parameters and serialisation etc. """
  ode_name = GetOdeSystemName(filename)

  private_class_string = ("public:\n" + 
                          "\n" + 
                          GetConstructorAndDestructorForOdeClassString(ode_name) + 
                          translator.GetFunctionDefinitionsStringForHeader(1, model) + "\n" +
                          translator.AddTabs(1) + "void Init();\n\n" +
                          translator.AddTabs(1) + "void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);\n\n")

  is_srn_model = translator.IsSrnModel(model)

  if (not is_srn_model):
    private_class_string += (translator.AddTabs(1) + "bool CalculateStoppingEvent(double time, const std::vector<double>& rY);\n\n" + 
                             translator.AddTabs(1) + "void CheckAndUpdateEvents(double time, const std::vector<double>& rY);\n\n" + 
                             translator.AddTabs(1) + "bool AreAllEventsSatisfied(double time, const std::vector<double>& rY);\n\n")

  return private_class_string

def GetOdeClassDefinitionString(filename, model):
    """ Return string defining the class in th header file. """
    ode_name = GetOdeSystemName(filename) #Get the name fo the ODE system

    #initialise string
    ode_class_string = ("class " + ode_name + " : public AbstractOdeSystem\n" + 
                        "{\n" + 
                        GetPrivateVariablesForOdeClassString(model) + 
                        GetPublicVariablesForOdeClassString(filename, model) + 
                        "};\n\n")

    return ode_class_string

def GetSaveConstructDataString(ode_name):
    """ Get the string describing the serialization required to construct the ODE system. """
    save_construct_data_string = (translator.GetBlockCommentDefinition(0, "Serialize information required to construct a " + ode_name + ".", True) +
                                  "template<class Archive>\n" + 
                                  "inline void save_construct_data(\n" + 
                                  translator.AddTabs(1) + "Archive & ar, const " + ode_name + " * t, const unsigned int file_version)\n" + 
                                  "{\n" + 
                                  translator.AddTabs(1) + "const std::vector<double> state_variables = t->rGetConstStateVariables();\n" + 
                                  translator.AddTabs(1) + "ar & state_variables;\n" + 
                                  "}\n" )

    return save_construct_data_string

def GetLoadConstructDataString(ode_name):
    """ Get the string describing the de-serialization of the constructor and initialisation of the ODE system. """
    save_load_data_string = (translator.GetBlockCommentDefinition(0, "De-serialize constructor parameters and intiialise a " + ode_name + ".", True) + 
                             "template<class Archive>\n" + 
                             "inline void load_construct_data(\n" + 
                             translator.AddTabs(1) + "Archive & ar, " + ode_name + " * t, const unsigned int file_version)\n" + 
                             "{\n" + 
                             translator.AddTabs(1) + "std::vector<double> state_variables;\n" + 
                             translator.AddTabs(1) + "ar & state_variables;\n" + 
                             translator.AddTabs(1) + "\n" +
                             translator.GetCommentDefinition(1, "Invoke inplace constructor to initialise instance", True) +
                             translator.AddTabs(1) + "::new(t)" + ode_name + "(state_variables);\n" +
                             "}\n")

    return save_load_data_string

def GetSerializationInformationString(filename):
    """ Define the serialization/namespace information required in the
    header file. """
    ode_name = GetOdeSystemName(filename)

    serialization_string = ("namespace\n" + 
                            "{\n" + 
                            "namespace serialization\n" +
                            "{\n" +
                            GetSaveConstructDataString(ode_name) +
                            GetLoadConstructDataString(ode_name) +
                            "}\n" + 
                            "} " + translator.GetCommentDefinition(0, "namespace ...", True) + "\n")

    return serialization_string

## Source files 

def GetIncludedFilesForSourceString(filename, model):
    """ Return string to include the necessary header files to go into the .cpp files. """

    model_name = GetModelName(filename, model)

    # Don't really need to automate this as it'll be the same for all of them. 
    include_files_string = (IncludeString("\"" + model_name + ".hpp\"") + 
                            IncludeString("\"CellwiseOdeSystemInformation.hpp\"") )

    return include_files_string

def GetClassDefinition(filename, model):
  """ Return string to define ODE class in .cpp file. """
  ode_name = GetOdeSystemName(filename)

  ##### EDITED ###################### - in this case, the rules are also variables (as totals of other variables)
  odes_dict = translator.GetOdesDictionary(model)
  rules_dict = translator.GetRulesDictionary(model)
  num_species = len(odes_dict) + len(rules_dict)

  class_defn_str = (translator.GetBlockCommentDefinition(0, "SBML ODE System", True) + 
                    ode_name + "::" + ode_name + " (std::vector<double> stateVariables)\n" +
                    translator.AddTabs(1) + ": AbstractOdeSystem(" + str(num_species) + ")\n" 
                    "{\n" + 
                    translator.AddTabs(1) + "mpSystemInfo.reset(new CellwiseOdeSystemInformation<" + ode_name + ">);\n"
                    "\n" + 
                    translator.AddTabs(1) + "Init();\n\n" + 
                    translator.GetDefaultInitialConditionString(model) + 
                    "\n" + 
                    translator.AddTabs(1) + "if (stateVariables != std::vector<double>())" +
                    translator.AddTabs(1) + "{\n" + 
                    translator.AddTabs(2) + "SetStateVariables(stateVariables);\n" + 
                    translator.AddTabs(1) + "}\n" + 
                    "}\n\n")

  return class_defn_str

def GetClassConstructor(filename):
  """ Instantiate construction of ODE system. """
  ode_name = GetOdeSystemName(filename)

  class_constr_str = (ode_name + "::~" + ode_name + "()\n"
                      "{\n" 
                      "}\n\n")

  return class_constr_str

def GetFunctionDefinitionsForSource(filename, model):
  """ Define SBML function handles. """
  ode_name = GetOdeSystemName(filename)

  funct_defns_str = translator.GetFunctionDefinitionsStringForSource(0, model, ode_name)

  return funct_defns_str

def GetInitForSource(filename, model):
  """ Initialise parameters. """
  ode_name = GetOdeSystemName(filename)

  init_param_str = ("void " + ode_name + "::Init()" + "\n " + 
                      "{\n" + 
                      translator.GetInitParameterStringForSource(1, model) + #Initialise parameter values
                      "}\n\n")
  return init_param_str 

def GetEvaluateYDerivativesVoidString(filename, model):
  """ Define derivatives functions of ODE system. """
  ode_name = GetOdeSystemName(filename)

  ode_derivs_str = ("void " + ode_name + "::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)\n" + 
                    "{\n" + 
                    translator.GetStateVariableString(model) + #Define state variables
                    "\n" + 
                    translator.GetStateParameterString(model) + #Define state parameters
                    "\n " + 
                    # translator.GetCompartmentVariableString(model) + "\n" + #Define compartment variables, e.g. cell
                    # translator.GetConstantParameterString(model) + "\n" + #Define constant parameters
                    translator.GetRulesString(model) + "\n" + #Define rules to be used in reactions
                    translator.GetReactionString(model) + "\n" + #Define reactions
                    translator.GetOdesString(model) + "\n" + #Define odes in terms of reactions
                    "}\n\n")

  return ode_derivs_str 

def GetInitialiseString(filename, model):
  """ Returns templated function to initialise ODE system. """
  ode_name = GetOdeSystemName(filename)

  initialise_void_str = ("template<>\n" + 
                         "void CellwiseOdeSystemInformation<" + ode_name + ">::Initialise()\n" + 
                         "{\n" + 
                         translator.GetInitialInformationString(model) + "\n" + #Initial set up of name, units and initial concentration
                         translator.AddTabs(1) + "this->mInitialised = true;\n" + 
                         "}\n\n")

  return initialise_void_str 

############################################################################################################
#####################                 SRN-specific functiions                ###############################
############################################################################################################

def WriteHeaderFileForSrnModel(filename, model):
    """ Construct the Chaste header file from the SBML file. """ 

    srn_model_name = GetModelName(filename, model) # Get the name of the file we will write 

    #Open to file to write
    header_file = open(srn_model_name + ".hpp", 'w')

    #Define the header files
    header_file_defn = GetHeaderFileDefinitionString(filename, model)
    header_file.write(header_file_defn)

    #Include the appropriate files
    include_files = GetIncludedFilesForHeaderString()
    header_file.write(include_files)

    #Define the ODE System class
    ode_class = GetOdeClassDefinitionString(filename, model)
    header_file.write(ode_class)

    #Define the serialization
    serialization = GetSerializationInformationString(filename)
    header_file.write(serialization)

    #Define the SRN model
    srn_model_defn = GetModelDefinitionString(filename, model, True)
    header_file.write(srn_model_defn)

    #Close the file
    header_close = GetHeaderFileClosingString(filename, model)
    header_file.write(header_close)

    header_file.close()

    print(srn_model_name + ".hpp written!\n")

def WriteSourceFileForSrnModel(filename, model):
  """ Function write source file for SrnModel. """
  srn_model_name = GetModelName(filename, model) # Get the name of the file we will write 

  #Open to file to write
  source_file = open(srn_model_name + ".cpp", 'w')

  #Include header files
  included_files = GetIncludedFilesForSourceString(filename, model)
  source_file.write(included_files)

  #Initialise class
  class_def = GetClassDefinition(filename, model)
  source_file.write(class_def)

  #Constructor for system
  constructor = GetClassConstructor(filename)
  source_file.write(constructor)

  #Functiond efinitions
  funct_defn_str = GetFunctionDefinitionsForSource(filename, model)
  source_file.write(funct_defn_str)

  #Initialise parameters
  init_fn = GetInitForSource(filename, model)
  source_file.write(init_fn)

  #Get the derivative function
  derivs_fn = GetEvaluateYDerivativesVoidString(filename, model)
  source_file.write(derivs_fn)

  #Initialise function
  initialise_fn = GetInitialiseString(filename, model)
  source_file.write(initialise_fn)

  #Define SRN Model
  srn_model_defn = GetModelDefinitionString(filename, model, False)
  source_file.write(srn_model_defn)

  source_file.close()

  print(srn_model_name + ".cpp written!\n")

def WriteSrnModelToFile(filename, model):
  """ Script that calls upon functions to write the header
   and source files needed for SRN models in Chaste. """

  # Write the .hpp file
  WriteHeaderFileForSrnModel(filename, model)

   # Write the .cpp fil
  WriteSourceFileForSrnModel(filename, model)

############################################################################################################
#####################                 CCM-specific functiions                ###############################
############################################################################################################
def GetStoppingEventBooleanString(filename, model):
  """String to define Stopping Event boolean required in Chaste for CCMs to divide. """
  ode_name = GetOdeSystemName(filename)

  stopping_event_str = ("bool " + ode_name + "::CalculateStoppingEvent(double time, const std::vector<double> & rY)\n" + #Define function
                          "{\n" + 
                          translator.GetBlockCommentDefinition(1, "Return true if all events have been triggered.", True) + 
                          translator.AddTabs(1) + "return AreAllEventsSatisfied(time, rY);\n" + #Define stopping event
                          "}\n\n")

  return stopping_event_str

def GetCheckAndUpdateEventsVoidString(filename, model):
  """ Construct the void that checks if any events have been triggered and 
  updates the system accordingly."""
  ode_name = GetOdeSystemName(filename)

  update_events_str = ("void " + ode_name + "::CheckAndUpdateEvents(double time, const std::vector<double> & rY)\n" + #Define function
                      "{\n" + 
                      translator.AddTabs(1) + "std::vector<double> dy(rY.size());\n" + #Initialise derivatives vector
                      translator.AddTabs(1) + "EvaluateYDerivatives(time, rY, dy);\n\n" + #Evaluate derivatives
                      translator.GetAllEventsString(model) + 
                      "}\n\n")

  return update_events_str

def GetAreAllEventsSatisfiedBooleanString(filename):
  """ Boolean to check whether or not all events are satisifed. """
  ode_name = GetOdeSystemName(filename)

  are_events_satisfied_str = ("bool " + ode_name + "::AreAllEventsSatisfied(double time, const std::vector<double>& rY)\n" + #Define function
                             "{\n" + 
                             translator.AddTabs(1) + "CheckAndUpdateEvents(time, rY);\n\n" + 
                              translator.AddTabs(1) + "bool events_satisfied = true;\n\n" + 
                              translator.AddTabs(1) + "if ( std::find(eventsSatisfied.begin(), eventsSatisfied.end(), false) != eventsSatisfied.end() )\n" + 
                              translator.AddTabs(1) + "{\n" +
                              translator.AddTabs(2) + "events_satisfied = false;\n" + 
                              translator.AddTabs(1) + "}\n\n" + 
                              translator.AddTabs(1) + "if (events_satisfied) " + translator.GetCommentDefinition(0, "Reset events vector if division is triggered", True) + 
                              translator.AddTabs(1) + "{\n" +
                              translator.AddTabs(2) + "std::fill(eventsSatisfied.begin(), eventsSatisfied.end(), false);\n" + 
                              translator.AddTabs(1) + "}\n\n" + 
                              translator.AddTabs(1) + "return events_satisfied;\n" +
                              "}\n\n")

  return are_events_satisfied_str

def WriteHeaderFileForCcmModel(filename, model):
    """ Construct the Chaste header file from the SBML file. """ 

    ccm_model_name = GetModelName(filename, model) # Get the name of the file we will write 

    #Open to file to write
    header_file = open(ccm_model_name + ".hpp", 'w')

    #Define the header files
    header_file_defn = GetHeaderFileDefinitionString(filename, model)
    header_file.write(header_file_defn)

    #Include the appropriate files
    include_files = GetIncludedFilesForHeaderString()
    header_file.write(include_files)

    #Define the ODE System class
    ode_class = GetOdeClassDefinitionString(filename, model)
    header_file.write(ode_class)

    #Define the serialization
    serialization = GetSerializationInformationString(filename)
    header_file.write(serialization)

    #Define the SRN model
    srn_model_defn = GetModelDefinitionString(filename, model, True)
    header_file.write(srn_model_defn)

    #Close the file
    header_close = GetHeaderFileClosingString(filename, model)
    header_file.write(header_close)

    header_file.close()

    print(ccm_model_name + ".hpp written!\n")

def WriteSourceFileForCcmModel(filename, model):
  """ Function write source file for SrnModel. """
  ccm_model_name = GetModelName(filename, model) # Get the name of the file we will write 

  #Open to file to write
  source_file = open(ccm_model_name + ".cpp", 'w')

  #Include header files
  included_files = GetIncludedFilesForSourceString(filename, model)
  source_file.write(included_files)

  #Initialise class
  class_def = GetClassDefinition(filename, model)
  source_file.write(class_def)

  #Constructor for system
  constructor = GetClassConstructor(filename)
  source_file.write(constructor)

  #Function definitions
  funct_defn_str = GetFunctionDefinitionsForSource(filename, model)
  source_file.write(funct_defn_str)

  #Initialise parameters
  init_fn = GetInitForSource(filename, model)
  source_file.write(init_fn)

  #Get the derivative function
  derivs_fn = GetEvaluateYDerivativesVoidString(filename, model)
  source_file.write(derivs_fn)

  #Get the stopping event function
  stopping_event_fn = GetStoppingEventBooleanString(filename, model)
  source_file.write(stopping_event_fn)

  #Get the void to check and update SBML events
  events_fn = GetCheckAndUpdateEventsVoidString(filename, model)
  source_file.write(events_fn)

 #Get the void to check and update SBML events
  events_satisfied_fn = GetAreAllEventsSatisfiedBooleanString(filename)
  source_file.write(events_satisfied_fn)

  #Initialise function
  initialise_fn = GetInitialiseString(filename, model)
  source_file.write(initialise_fn)

  #Define SRN Model
  srn_model_defn = GetModelDefinitionString(filename, model, False)
  source_file.write(srn_model_defn)

  source_file.close()

  print(ccm_model_name + ".cpp written!\n")

def WriteCcmModelToFile(filename, model):
    """ Script that calls upon functions to write the header
    and source files needed for CCM models in Chaste. """
    #Write the .hpp file
    WriteHeaderFileForCcmModel(filename, model)

  #Write the .cpp fil
    WriteSourceFileForCcmModel(filename, model)