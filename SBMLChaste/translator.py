#!/usr/bin/env python

import sys
import os.path
from libsbml import *

# Script with useful functions to extract the relevant information from an 
# SBML file needed to define an ODE Model in Chaste

###################################################################################################
############################            Chaste-ey functions         ###############################
###################################################################################################
def AddTabs(n):
    """ Function add tabs however many times we need to. """ 
    return "    " * n

def GetBlockCommentDefinition(num_tabs, var, newline):
    """ Define string as C++ block comment. """
    block_comment_def = AddTabs(num_tabs) + "/* " + var + " */"

    if newline:
        block_comment_def += "\n"

    return block_comment_def

def GetCommentDefinition(num_tabs, var, newline):
    """ Define string as C++ comment. """
    comment_def = AddTabs(num_tabs) + "// " + var

    if newline:
        comment_def += "\n"

    return comment_def

def GetConstDoubleDefinition(num_tabs, name, val, newline):
    """ Get the string defining a constant double in C++. """
    const_double_def = AddTabs(num_tabs) + "const double " + name + " = " + str(val) + ";"

    if newline:
        const_double_def += "\n"

    return const_double_def

def GetDoubleDefinition(num_tabs, name, val, newline):
    """ Get the string defining a double in C++. """
    double_def = AddTabs(num_tabs) + "double " + name +  " = " + str(val) + ";"

    if newline:
        double_def += "\n"

    return double_def

def GetOdeDefinition(num_tabs, n, rhs, comp, newline):
    """ Get a string defining the ODE in C++. """
    ode_def = AddTabs(num_tabs) + "rDY[" + str(n) + "] = (" + rhs + ") / " + comp + ";"

    if newline:
        ode_def += "\n"

    return ode_def


def GetOdeDefinitionFromAssignment(num_tabs, n, rhs, comp, newline):
    """ Get a string defining the ODE in C++. """
    ode_def = AddTabs(num_tabs) + "rDY[" + str(n) + "] = (" + rhs + ") / " + comp + ";"

    if newline:
        ode_def += "\n"

    return ode_def

def GetStateVariableDefinition(num_tabs, n, sv, newline):
    """ Get a string that defines the state variable in C++. """
    state_variable_def = AddTabs(num_tabs) + "double " + sv + " = rY[" + str(n) + "];"

    if newline:
        state_variable_def += "\n"  

    return state_variable_def

def GetStateParameterDefinition(num_tabs, n, sp, newline):
    """ Get a string that defines the state variable in C++. """
    state_parameter_def = AddTabs(num_tabs) + "double " + sp + " = this->mParameters[" + str(n) + "];"

    if newline:
        state_parameter_def += "\n"

    return state_parameter_def


def GetNodeIndex(node_list, node):
    """ Function to get the corresponding node index from an ASTNodeList, as
    these objects are not enumerable, nor are they `lists' in the python sense. """
    indices = [i for i in range(node_list.getSize()) if (node_list.get(i) == node)] #Have to do it as a list
    return indices[0] # Only ever going to be one node

def GetCppEquivalent(ast_name):
    """ Function to get the equivalent of an AST node type. As the list of these operators used is not too large,
    we can essentially hardcode the dictionary of operators and their equivalents. """
    operators = {"and" : " && ", "or" : " || ", "not" : " ! ", "eq" : " == ", "neq": " != ", "gt": " > ",
                    "geq": " >= ", "lt": " < ", "leq": " <= "}

    if ast_name in operators:
        return operators[ast_name]
    else:
        return ast_name

def ConvertVariableNamesIntoChasteEquivalent(string_with_variables, model):
    """ Function to convert a string with variable names into the Chaste equivalent, replacing
    the variable names with rY[INDEX]. """

    # Get the species ids and their Chaste equivalents
    species_ids_dict = GetSpeciesIndices(model)

    # Create the argument to be passed into the holder function
    argument_list = [string_with_variables]

    for key, value in species_ids_dict.iteritems():
        argument_list.append((key, value))

    string_with_chaste_names = reduce(lambda s, p: s.replace(p[0], p[1]), argument_list)

    return string_with_chaste_names


def GetVariableString(node):
    """ Function to get the string equivalent of the node type, whether that is an operator,
    variable name or value. """
    variable = ""

    if(node.isNumber()): #If the node is a number, we convert the value into a string
        node_value = node.getValue()
        variable += str(node_value) + " )"
    else: #We can just grab the name 
        node_name = node.getName()
        if(node.isLogical()): #If the node represents a logical operator
            node_operator = GetCppEquivalent(node_name)
            variable += node_operator
        else:
            if(node.isName()): # If we have a variable name
                variable += "( " + node_name
            else: #Else we have a relational operator, e.g. <, <= etc
                node_relator = GetCppEquivalent(node_name)
                variable += node_relator

    return variable

 ## Constant parameter functions
def GetConstantParameterString(model):
    """ Get all the constant parameters used in the reactions. These include
    local reaction parameters and compartment volumes (assuming single compartment)
    for now. """
    num_parameters = model.getNumParameters()

    const_param_string = GetBlockCommentDefinition(1, "Define constant parameters and their values", True)

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        if(parameter.getConstant()):
            param_id = parameter.getId()
            param_value = parameter.getValue()
            param_def = GetConstDoubleDefinition(1, param_id, param_value, True)
            const_param_string += param_def

    return const_param_string

## Non-constant parameter functions
def GetNonConstantParameterString(model):
    """ Get all the non-constant parameters used in the reactions. These include
    local parameters used to flag division events (e.g Tyson and Novak (2001))."""
    num_parameters = model.getNumParameters()

    non_const_param_string = GetBlockCommentDefinition(1, "Define non-constant parameters and their values", True)

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        if(not parameter.getConstant()):
            param_id = parameter.getId()
            param_value = parameter.getValue()
            param_def = GetConstDoubleDefinition(1, param_id, param_value, True)
            non_const_param_string += param_def

    return non_const_param_string

def GetInOrderTraversal(list_of_nodes, tree, preorder_list):
    """Function to traverse a tree in-order and produce a list of nodes.
    Assumes the tree is an ASTNode, as defined in SBML. Also assumes the list
    passed is initially empty (I don't know how to get around this). """

    if(tree.getNumChildren() != 0):
        left_subtree = tree.getLeftChild()
        GetInOrderTraversal(list_of_nodes, left_subtree, preorder_list)

        node_index = GetNodeIndex(preorder_list, tree)
        list_of_nodes.append(node_index)

        right_subtree = tree.getRightChild()
        GetInOrderTraversal(list_of_nodes, right_subtree, preorder_list)

    else:
        node_index = GetNodeIndex(preorder_list, tree)
        list_of_nodes.append(node_index)
    
    return list_of_nodes

def IsSrnModel(model):
    """Boolean to determine whether a model's an SRN or a CCM. """
    num_events = model.getNumEvents()
    if (num_events != 0):
        return False
    else:
        return True

###################################################################################################
######################        SRN AND Cell Cycle Model Functions         ##########################
###################################################################################################
#Parameter functions
def GetParameter(n, model):
    """ Get the ID of the parameter """
    parameter = model.getParameter(n)
    parameter_id = parameter.getId()
    return parameter_id

def GetParameterValue(n, model):
    """ Get the value of the parameter """
    parameter = model.getParameter(n)
    parameter_id = parameter.getValue()
    return parameter_id

def GetParameterName(n, model):
    """ Get the ID of the parameter"""
    parameter = model.getParameter(n)
    parameter_name = parameter.getName()
    return parameter_name

def IsParameterConstant(n, model):
    """Boolean to check if parameter is constant """
    parameter = model.getParameter(n)

    if(parameter.getConstant()):
        return True
    else:
        return False

def IsParameterValueSet(n, model):
    """Boolean to check if parameter has a pre-set value """
    parameter = model.getParameter(n)

    if(parameter.isSetValue()):
        return True
    else:
        return False

## Function definitions

def GetArgumentsOfFunctionDefinition(funct_def):
    """Function to obtain a string of the arguments in a given functional definition."""
    args_str = ""
    num_args = funct_def.getNumArguments()

    for i in range(num_args):
        funct_arg = funct_def.getArgument(i)
        arg_name = formulaToString(funct_arg)
        args_str += "double " + arg_name
        if ( (num_args != 1)&(i < num_args - 1) ):
            args_str += ", "

    return args_str

def GetSubfunctionEquivalents(funct_def_body):
    """Function to check if any functions used don't match up with their C++ counterparts."""
    functs = {"abs" : "fabs", "arccos" : "acos", "arccsc" : "acsc", "arccsch" : "acsch", "arcsec" : "asec", 
                "arcsech" : "asech", "arcsinh" : "asinh", "arcsin" : "asin", "arctan" : "atan", 
                "arctanh" : "atanh", "max" : "fmax", "min" : "fmin"}

    subfunct_equivalents = {}

    list_of_nodes = []
    preorder_list = funct_def_body.getListOfNodes()
    inorder_list = GetInOrderTraversal(list_of_nodes, funct_def_body, preorder_list)

    for i in range(len(inorder_list)):
        node_index = inorder_list[i]
        node = preorder_list.get(node_index)
        if (node.isFunction()):
            node_formula = formulaToString(node)

            #If the function is root, we have to deal with it in an annoying way...
            if ("root" in node_formula):
                #Split the root by the comma and rearrange
                split_formula = node_formula.split(",", 1)
                first_part = split_formula[0]
                second_part = split_formula[1]

                #Get the exponent of the root
                index = first_part.find("(")
                exponent = first_part[index + 1:len(first_part)]

                #Get the base
                base = second_part[0:len(second_part) - 1] 
                subfunct_equivalents[node_formula] = "pow(" + base + ", " + "1.0 /" + exponent + ")"

            #The rest can be done by simple string replacement
            else:
                new_formula = node_formula

                for key in functs:
                    new_formula = new_formula.replace(key, functs[key])

                subfunct_equivalents[node_formula] = new_formula 

    return subfunct_equivalents

def ConvertFunctionBody(funct_def_body):
    """Function to the C++ equivalent. """
    subfunct_equivalents = GetSubfunctionEquivalents(funct_def_body)

    #Now replace all instances of those sub-functions with their C++ equivalents
    funct_body_formula = formulaToString(funct_def_body)

    for key in subfunct_equivalents:
        funct_body_formula = funct_body_formula.replace(key, subfunct_equivalents[key])

    return funct_body_formula

def GetBodyOfFunctionDefinition(funct_def):
    """Function to obtain the body of a function definition in C++ form. """
    body_str = ""
    funct_body = funct_def.getBody()
    body_str += ConvertFunctionBody(funct_body)
    return body_str

def GetFunctionDefinitionString(n, funct_def, ode_sys_name):
    """Function to define C++ function definition in the source files."""
    funct_def_str = ""
    funct_id = funct_def.getId()
    funct_args = GetArgumentsOfFunctionDefinition(funct_def)
    funct_body = GetBodyOfFunctionDefinition(funct_def)

    funct_def_str += (AddTabs(n) + "double " + ode_sys_name + "::" + funct_id + "(" + funct_args + ")\n" + 
                     AddTabs(n) + "{\n" + 
                     AddTabs(n + 1) + "return " + funct_body + ";\n" + 
                     AddTabs(n) + "}\n")
    return funct_def_str 

def GetFunctionDefinitionsStringForHeader(n, model):
    """Function to define all the function definitions in Initialise of header file. """
    num_funct_defns = model.getNumFunctionDefinitions()

    funct_defns_str = ""

    for i in range(num_funct_defns):
        funct_defn = model.getFunctionDefinition(i)
        funct_id = funct_defn.getId()
        funct_args = GetArgumentsOfFunctionDefinition(funct_defn)

        funct_defns_str += AddTabs(n) + "double " + funct_id + "(" + funct_args + ");\n"

    return funct_defns_str
    
def GetFunctionDefinitionsStringForSource(n, model, ode_sys_name):
    """Function to define all the function definitions in Initialise of source file. """
    num_funct_defns = model.getNumFunctionDefinitions()

    funct_defns_str = ""

    for i in range(num_funct_defns):
        funct_defn = model.getFunctionDefinition(i)
        funct_defns_str += GetFunctionDefinitionString(n, funct_defn, ode_sys_name) + "\n"

    return funct_defns_str

## Reaction rule functions
def GetRule(n, model):
    """ Get the ID of rule n to use as variable name. """
    rule = model.getRule(n)
    rule_id = rule.getId()
    return rule_id

def GetRulesDictionary(model):
    """Dictionary of species defined by reaction rules. """
    rules_dict = {} 
    num_rules = model.getNumRules()

    for i in range(num_rules):
        rule = model.getRule(i)
        rule_id = rule.getId()
        rule_formula = rule.getFormula()
        rules_dict[rule_id] = rule_formula

    return rules_dict

def IsVariableDefinedAsRule(var_id, model):
    """ Boolean to check whether or not a species is defined as a reaction rule. """

    #If species is defined in the dictionary of reaction rules, return true
    rules_dict = GetRulesDictionary(model)
    if (var_id in rules_dict):
        return True
    else:
        return False

def GetRulesString(model):
    """ Get all the rules required to define reactions. """
    num_rules = model.getNumRules()

    rules_string = GetBlockCommentDefinition(1, "Define algebraic rules.", True)
    rules_dict = GetRulesDictionary(model)

    for i in range(num_rules):
        rule_id = GetRule(i, model)
        rule_def = AddTabs(1) + rule_id + " = " + rules_dict[rule_id] + ";\n"
        rules_string += rule_def
        rules_string = GetReactionStringWithSpeciesNamesAndParameterNames(model,rules_string)
        rules_string = rules_string.replace("max","fmax")
        rules_string = rules_string.replace("min","fmin")

    return rules_string

## Reaction functions

def GetReaction(n, model):
    """ Get the ID of a reaction to use as a variable name. """
    reaction = model.getReaction(n)
    reaction_id = reaction.getId()
    reaction_name = reaction.getName()
    return reaction_id

def GetReactionFormula(n, model):
    """ Get the formula associated with the reaction. """
    reaction = model.getReaction(n)
    kinetic_law = reaction.getKineticLaw() # Get the kinetic law
    formula = kinetic_law.getFormula()
    return formula



################################## NEW ######################################################

def GetReactionStringWithSpeciesNamesAndParameterNames(model,string): 

    num_species = model.getNumSpecies()
    num_parameters = model.getNumParameters()
    num_compartments = model.getNumCompartments()

    for i in range(num_species):
        species = model.getSpecies(i)
        species_id= species.getId()
        species_name = species.getName()
        if (species_name != ''):
            string = string.replace(species_id,species_name)

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        parameter_id= parameter.getId()
        parameter_name = parameter.getName()
        if (parameter_name != ''):
            string = string.replace(parameter_id,parameter_name)

    for i in range(num_compartments):
        compartment = model.getCompartment(i)
        compartment_id= compartment.getId()
        compartment_name = compartment.getName()
        if (compartment_name != ''):
            string = string.replace(compartment_id,compartment_name)

    return string


def GetSpeciesNames(model,string): 

    num_species = model.getNumSpecies()

    for i in range(num_species):
        species = model.getSpecies(i)
        species_id= species.getId()
        species_name = species.getName()
        if (species_name != ''):
            string = string.replace(species_id,species_name)

    return string

def GetParameterNames(model,string): 

    num_parameters = model.getNumParameters()

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        parameter_id= parameter.getId()
        parameter_name = parameter.getName()
        if (parameter_name != ''):
            string = string.replace(parameter_id,parameter_name)

    return string

def GetReactionNames(model,string): 

    num_reactions = model.getNumReactions()

    for i in range(num_reactions):
        reaction = model.getReaction(i)
        reaction_id= reaction.getId()
        reaction_name = reaction.getName()
        if (reaction_name != ''):
            string = string.replace(reaction_id,reaction_name)

    return string


def GetCompartmentNameCorrespondingToId(model,string):
    
    compartment_id = string
    num_compartments = model.getNumCompartments()
    compartment_name = ''

    for i in range(num_compartments):
        compartment_ID = GetCompartment(i, model)

        if compartment_ID == compartment_id:
            compartment_name = GetCompartmentName(i,model)

    if compartment_name != '':
        return compartment_name
    else:
        return compartment_id


#############################################################################################  

def GetReactionParameterString(reaction):
    """ Get the reaction parameters string. """
    kinetic_law = reaction.getKineticLaw()
    num_kinetic_law_parameters = kinetic_law.getNumParameters()

    parameters_string = ""

    for i in range(num_kinetic_law_parameters):
        parameter = kinetic_law.getParameter(i)
        parameter_id = parameter.getId()
        parameter_value = parameter.getValue()
        parameter_id_and_value = GetDoubleDefinition(1, parameter_id, parameter_value, True)
        parameters_string += parameter_id_and_value

    return parameters_string


def GetReactionString(model):
    """ Get the reactions in the model. """

    num_reactions = model.getNumReactions()

    reactions_string = GetBlockCommentDefinition(1, "Define the reactions in this model.", True)

    for i in range(num_reactions):
        reaction = model.getReaction(i)
        reaction_name = reaction.getName()
        reaction_id = GetReaction(i, model)
        reaction_formula = GetReactionFormula(i, model)

        reaction_formula = GetReactionStringWithSpeciesNamesAndParameterNames(model,reaction_formula)

        #Get the relevant strings 
        reaction_params_string = GetReactionParameterString(reaction)
        reaction_def = GetDoubleDefinition(1, reaction_id, reaction_formula, True)

        if(reaction_name):
            reaction_name_string = GetCommentDefinition(1, reaction_name, True)
            reactions_string += reaction_name_string

        reactions_string += reaction_params_string
        reactions_string += reaction_def

        reactions_string += "\n"

    reactions_string = GetReactionNames(model,reactions_string)

    return reactions_string


def AddReactionToDictionary(react_dict, reaction):
    """ Decompose a reaction into the sum of products minus sum of reactants and
    add it to the dictionary of reactions. """
    reaction_id = reaction.getId() #Get the reaction variable

    original_model = reaction.getModel() #Need the original model to check for boundary conditions

        # Get the products and reactants, each ODE will essentially be the sum of
    # the products minus the sum fo the reactants divided by the compartment volume
    reaction_products = reaction.getListOfProducts()
    reaction_reactants = reaction.getListOfReactants()

    if (len(reaction_products) != 0):
        for j in range(len(reaction_products)):
            #Get the species concerning the product
            product_species_ref = reaction_products.get(j)
            product_species_id = product_species_ref.getSpecies()

            product_species =  original_model.getSpecies(product_species_id)

            if( (product_species.isSetBoundaryCondition())&(not product_species.getBoundaryCondition()) ):
                if(product_species_id in react_dict):
                    react_dict[product_species_id] += " + " + reaction_id 
                else:
                    react_dict[product_species_id] = reaction_id
            else:
                if(product_species_id in react_dict):
                    react_dict[product_species_id] += " + " + reaction_id 
                else:
                    react_dict[product_species_id] = reaction_id

    if (len(reaction_reactants) != 0):
        for j in range(len(reaction_reactants)):
            reactant_species_ref = reaction_reactants.get(j)
            reactant_species_id = reactant_species_ref.getSpecies()

            reactant_species =  original_model.getSpecies(reactant_species_id)

            if( (reactant_species.isSetBoundaryCondition())&(not reactant_species.getBoundaryCondition()) ):
                if(reactant_species_id in react_dict):
                    react_dict[reactant_species_id] += " - " + reaction_id 
                else:
                    react_dict[reactant_species_id] = " - " + reaction_id 
            else:
                if(reactant_species_id in react_dict):
                    react_dict[reactant_species_id] += " - " + reaction_id 
                else:
                    react_dict[reactant_species_id] = " - " + reaction_id    


def GetOdesDictionary(model):
    """ Get the odes as a dictionary of equations corresponding to each species. """

    # Initialise dictionary
    odes_dict = {}

    num_reactions = model.getNumReactions()

    for i in range(num_reactions):
        reaction = model.getReaction(i) 
        #Add reaction products and reactants to the dictionary
        AddReactionToDictionary(odes_dict, reaction)

    for key,value in odes_dict.items():
        odes_dict[key] = GetReactionStringWithSpeciesNamesAndParameterNames(model,odes_dict[key])

    return odes_dict

def IsSpeciesDefinedAsOde(species_id, model):
    """ Boolean to check whether or not a species is defined as an ODE. """

    #If species is defined in the dictionary of reaction rules, return true
    odes_dict = GetOdesDictionary(model) 

    #We need to get the original species to see if it's set as a boundary condition
    species = model.getSpecies(species_id)

    if (species_id in odes_dict):
        return True
    else:
        return False

def IsTimescaleDifferent(model):
    """Function to check whether or not there's a difference in timescaless, as SBML uses seconds by default
    and Chaste uses hours."""


    num_unit_defns = model.getNumUnitDefinitions()
    if(num_unit_defns != 0):
        for i in range(num_unit_defns):
            unit_defn = model.getUnitDefinition(i)
            unit_id = unit_defn.getId()
            if(unit_id.lower() == "time"): # Do people ever call this something different?
                timescale = unit_defn.getName()
                if("minute" in timescale.lower()):
                    return True
                    break
                elif("hour" in timescale.lower()):
                    return False
    else:
        return True #Generally will be different

def GetTimescaleMultiplier(model):
    #Function to get the multiplier 
    timescale_multiplier = ""

    num_unit_defns = model.getNumUnitDefinitions()

    if(num_unit_defns != 0):
        for i in range(num_unit_defns):
            unit_defn = model.getUnitDefinition(i)
            unit_id = unit_defn.getId()
            if(unit_id.lower() == "time"): # Do people ever call this something different?
                timescale = unit_defn.getName()
                if("minute" in timescale.lower()):
                    timescale_multiplier = "60.0"
    else:
        timescale_multiplier = "3600.0"

    return timescale_multiplier

def AccountForTimescale(model):
    """ Function to account for any differences in timescales. """
    odes_dict = GetOdesDictionary(model) #Get the dictionary describing the ode reactions 
    num_species = len(odes_dict) #This does not include the species defined in algebraic rules

    timescale_multiplier = GetTimescaleMultiplier(model)

    timescale_str = GetBlockCommentDefinition(1, "Account for the differences in timescales.", True)

    for i in range(num_species):
        timescale_str += AddTabs(1) + "rDY[" + str(i) + "] *= " + timescale_multiplier + ";\n"

    return timescale_str

def GetOdesString(model):
    """ Get the string describing the system of ODEs. """ 

    ode_system_string = "" #Initialise string

    odes_dict = GetOdesDictionary(model) #Get the dictionary describing the ode reactions 

    num_species = model.getNumSpecies()
    num_rules = model.getNumRules()

    ode_count = 0
    for i in range(num_species):
        species = model.getSpecies(i)
        species_id= species.getId()
        species_name= species.getName()
        if (species_name == ''):
            species_name = species_id
        if( (species_id in odes_dict)&(not species.getBoundaryCondition()) ): #species defined by algebraic rules are not in odes_dict
            species_ode = odes_dict[species_id] #Get the corresponding ode

            ######### EDITED #############
            species_ode = GetReactionNames(model,species_ode)

        # species = model.getSpecies(key)
            compartment_id = species.getCompartment() #Get the corresponding compartment

            ############# EDITED ##############
            compartment_name = GetCompartmentNameCorrespondingToId(model,compartment_id)

            #Get the C++ definition for the ODE
            ode_def = GetOdeDefinition(1, ode_count, species_ode, compartment_name, False)
            ode_count += 1

            #Get the comment describing the ODE
            comment_str = "d" + species_name + "/dt"
            ode_comment = GetCommentDefinition(0, comment_str, True)

            #Add the ode definition and comment strings together
            ode_string = ode_def + " " + ode_comment

            ode_system_string += ode_string

        ########################## EDITED #############
        elif( (IsVariableDefinedAsRule(species_id,model)) & (not (species_id in odes_dict)) & (not ((species_name == "drag") or (species_id == "drag")))):

            rules_dict = GetRulesDictionary(model)

            rule_string = rules_dict[species_id] ## I assume these are assignments where variables are added together to represent a total: TODO: include other rules

            rule_string = rule_string.replace("max","fmax")
            rule_string = rule_string.replace("min","fmin")

            species_ode = ''

            rule_string2 = rule_string.split(' ')

            for j in range(num_species):
                species_new = model.getSpecies(j)
                species_id_new = species_new.getId()
                species_name_new = species_new.getName()

                for k in rule_string2:
                    if ((species_id_new == k) & (species_ode == '')):
                        species_ode = species_ode + "rDY[" + str(j) + "]"
                    elif ((species_id_new == k) & (not (species_ode == ''))):
                        species_ode = species_ode + "+ rDY[" + str(j) + "]"

            species_ode = GetReactionStringWithSpeciesNamesAndParameterNames(model,species_ode)

            compartment_id = species.getCompartment() #Get the corresponding compartment
            compartment_name = GetCompartmentNameCorrespondingToId(model,compartment_id)

            ode_def = GetOdeDefinitionFromAssignment(1, ode_count, species_ode, compartment_name, False)
            ode_count += 1

            comment_str = "d" + species_name + "/dt"
            ode_comment = GetCommentDefinition(0, comment_str, True)

                    #Add the ode definition and comment strings together
            ode_string = ode_def + " " + ode_comment

            ode_system_string += ode_string

        elif( (species_name == "drag") or (species_id == "drag")): ########################## EDITED #######################
            species_ode = ''
            species_ode = "drag" + " - rY[" + str(i) + "]"

            compartment_id = species.getCompartment() #Get the corresponding compartment
            compartment_name = GetCompartmentNameCorrespondingToId(model,compartment_id)

            ode_def = GetOdeDefinitionFromAssignment(1, ode_count, species_ode, compartment_name, False)
            ode_count += 1

            comment_str = "d" + species_name + "/dt"
            ode_comment = GetCommentDefinition(0, comment_str, True)

                    #Add the ode definition and comment strings together
            ode_string = ode_def + " " + ode_comment

            ode_system_string += ode_string



    if ( IsTimescaleDifferent(model) ):
        ode_system_string += "\n" + AccountForTimescale(model)

    return ode_system_string

## State variables functions 

def GetStateVariable(n, model):
    """ Get the ID of the species to use as state variable """
    species = model.getSpecies(n)
    species_id = species.getId()
    return species_id

def GetStateVariableName(n, model):
    """ Get the ID of the species to use as state variable """
    species = model.getSpecies(n)
    species_name = species.getName()
    return species_name

def GetStateVariableString(model):
    """ Extract state variables and produce string that can be written
    to C++ code. """
    num_species = model.getNumSpecies()

    state_variables_string = GetBlockCommentDefinition(1, "Define state variables", True)

    #We need to define the state variables in a round-about manner. We only define a species
    # to be a state variable if its evolution is described by an ODE.

    state_variable_count = 0
    for i in range(num_species):
        species_id = GetStateVariable(i, model)
#        if (IsSpeciesDefinedAsOde(species_id, model)): ################### EDITED #############################
        state_variable_def = GetStateVariableDefinition(1, state_variable_count, species_id, False)
        state_variable_count += 1
        state_variables_string += state_variable_def

        state_variable_name = GetStateVariableName(i, model)
        if(state_variable_name):
            state_variables_string += " " + GetCommentDefinition(0, state_variable_name, True)
        else:
            state_variables_string += "\n"

    state_variables_string = GetSpeciesNames(model,state_variables_string)

    return state_variables_string

##State parameter functions
def GetStateParameter(n, model):
    """ Get the ID of the species to use as state variable """
    species = model.getSpecies(n)
    species_id = species.getId()
    return species_id

def GetStateParameterName(n, model):
    """ Get the ID of the species to use as state parameter """
    species = model.getSpecies(n)
    species_name = species.getName()
    return species_name

def IsSpeciesDefinedAsParameter(species_id, model):
    """ Boolean to determine whether or not species is defined as parameter for Chaste. """
    if ( (not IsSpeciesDefinedAsOde(species_id, model))&(not IsVariableDefinedAsRule(species_id, model)) ):
        return True
    else:
        return False

def GetStateParameterString(model):
    """ Extract state parameters and produce string that can be written
    to C++ code. """
    state_parameters_string = GetBlockCommentDefinition(1, "Define state parameters", True)
    state_parameter_count = 0

    #Any species not defined by an ODE or rule is defined to be a parameter
    num_species = model.getNumSpecies()
    for i in range(num_species):
        species_id = GetStateVariable(i, model)
        if ( IsSpeciesDefinedAsParameter(species_id, model) ):
            state_parameter_def = GetStateParameterDefinition(1, state_parameter_count, species_id, False)
            state_parameter_count += 1
            state_parameters_string += state_parameter_def

            state_parameter_name = GetStateParameterName(i, model)
            if(state_parameter_name):
                state_parameters_string += " " + GetCommentDefinition(0, state_parameter_name, True)
            else:
                state_parameters_string += "\n"

    #Parameters with unset values should be set to be state parameters as well
    ########## EDITED ############ also parameters with 'wnt' or 'gamma' in their name should be set to state parameters
    num_parameters = model.getNumParameters()
    for i in range(num_parameters):
        parameter_id = GetParameter(i, model)
        parameter_name = GetParameterName(i,model)
        if ( (not IsParameterValueSet(i, model))&(not IsVariableDefinedAsRule(parameter_id, model)) ):
            state_parameter_def = GetStateParameterDefinition(1, state_parameter_count, parameter_id, False)
            state_parameter_count += 1
            state_parameters_string += state_parameter_def

            parameter_name = GetParameterName(i, model)
            if(parameter_name):
                state_parameters_string += " " + GetCommentDefinition(0, parameter_name, True)
            else:
                state_parameters_string += "\n"
        
        elif(parameter_name):
            if (('wnt' in parameter_name) or ('gamma' in parameter_name) or ('ComplexTransit' in parameter_name)):
                state_parameter_def = GetStateParameterDefinition(1, state_parameter_count, parameter_name, False)
                state_parameter_count += 1
                state_parameters_string += state_parameter_def
                state_parameters_string += " " + GetCommentDefinition(0, parameter_name, True)

    return state_parameters_string

def GetSpeciesIndices(model):
    """ Function to return a dictionary of species IDs and their
    Chaste variable name equivalent. """
    odes_dict = GetOdesDictionary(model)

    species_dict = {} # Initialise
    num_species = model.getNumSpecies()

    species_count = 0
    for i in range(num_species):
        species = model.getSpecies(i)
        species_id = species.getId()
        if( (species_id in odes_dict) ):
            species_dict[" " + species_id+ " "] = " rY[" + str(species_count) + "] " 
            species_count += 1

    return species_dict

## Compartment functions 

def GetCompartment(n, model):
    """ Get the id of the compartment to use as a variable """
    compartment = model.getCompartment(n)
    compartment_id = compartment.getId()
    return compartment_id

def GetCompartmentName(n, model):
    """ Get the id of the compartment to use as a variable """
    compartment = model.getCompartment(n)
    compartment_name = compartment.getName()
    return compartment_name

def GetCompartmentValue(n, model):
    """ Get the id of the compartment to use as a variable """
    compartment = model.getCompartment(n)
    compartment_value = compartment.getSize()
    return compartment_value

def GetCompartmentVariableString(model):
    """ Extract compartments and their sizes to produce string that is written to 
    C++ code. """
    num_compartments = model.getNumCompartments()

    compartments_string = GetBlockCommentDefinition(1, "Define compartments and their sizes", True)

    for i in range(num_compartments):
        compartment_id = GetCompartment(i, model)
        ################# EDITED ##################### use name rather than id
        compartment_name = GetCompartmentName(i,model) 
        compartment_val = GetCompartmentValue(i, model)
        compartment_def = GetConstDoubleDefinition(1, compartment_name, compartment_val, True)
        compartments_string += compartment_def

    return compartments_string

## Functions to return information needed to initialise ODEs
def GetInitialInformationForASpecies(species, model):
    """ Get the string describing the initial information for a species; that is,
    its name, units and initial concentrations. """

    #Define the string for the name
    species_name = species.getName()
    species_id = species.getId()

    initial_information_string = ""

    if(IsSpeciesDefinedAsOde(species_id, model)):
        if(species_name == ""): # If no name, set the name to be the species ID
            name_string = AddTabs(1) + "this->mVariableNames.push_back(\"" + species_id + "\");\n"
        else:
            name_string = AddTabs(1) + "this->mVariableNames.push_back(\"" + species_name + "\");\n"

        #Define the string for the units
        species_compartment = species.getCompartment()
        species_units = species.getSubstanceUnits()

        if (species_compartment): # if there's a compartment, then we would have normalised the ODE, so declare it as non-dimensional
            units_string = AddTabs(1) + "this->mVariableUnits.push_back(\"non-dim\");\n"
        else:
            units_string = AddTabs(1) + "this->mVariableUnits.push_back(\"" + species_units + "\");\n"

        #Define the string for the initial condition
        if(species.isSetInitialAmount()):
            species_init_conc = species.getInitialAmount()
        else:
            species_init_conc = species.getInitialConcentration()    
        init_conc_string = AddTabs(1) + "this->mInitialConditions.push_back(" + str(species_init_conc) + ");\n"

        initial_information_string += name_string + units_string + init_conc_string + "\n"
    elif (IsSpeciesDefinedAsParameter(species_id, model)):
        if(species_name == ""): # If no name, set the name to be the species ID
            name_string = AddTabs(1) + "this->mParameterNames.push_back(\"" + species_id + "\");\n"
        else:
            name_string = AddTabs(1) + "this->mParameterNames.push_back(\"" + species_name + "\");\n"

        #Define the string for the units
        species_compartment = species.getCompartment()
        species_units = species.getSubstanceUnits()

        if (species_compartment): # if there's a compartment, then we would have normalised the ODE, so declare it as non-dimensional
            units_string = AddTabs(1) + "this->mParameterUnits.push_back(\"non-dim\");\n"
        else:
            units_string = AddTabs(1) + "this->mParameterUnits.push_back(\"" + species_units + "\");\n"


    elif (IsVariableDefinedAsRule(species_id, model)):  ############################ EDITED ##################################### if species in rule from simbio 
        if(species_name == ""):
            name_string = AddTabs(1) + "this->mVariableNames.push_back(\"" + species_id + "\");\n"

        else:
            name_string = AddTabs(1) + "this->mVariableNames.push_back(\"" + species_name + "\");\n"

        #Define the string for the units
        species_compartment = species.getCompartment()
        species_units = species.getSubstanceUnits()

        if (species_compartment): # if there's a compartment, then we would have normalised the ODE, so declare it as non-dimensional
            units_string = AddTabs(1) + "this->mVariableUnits.push_back(\"non-dim\");\n"
        else:
            units_string = AddTabs(1) + "this->mVariableUnits.push_back(\"" + species_units + "\");\n"

        #Define the string for the initial condition
        if(species.isSetInitialAmount()):
            species_init_conc = species.getInitialAmount()
        else:
            species_init_conc = species.getInitialConcentration()    
        init_conc_string = AddTabs(1) + "this->mInitialConditions.push_back(" + str(species_init_conc) + ");\n"


        #Define string (should be three lines)      
        initial_information_string += name_string + units_string + init_conc_string + "\n"

    return initial_information_string

def GetInitialInformationForAParameter(parameter, model):
    """String describing initial information for a parameter; that is, its name and units. This function should only be called if
    the parameter does not have a set value (and thus must be externally defined). """
    parameter_name = parameter.getName()
    parameter_id = parameter.getId()

    initial_information_string = ""

    if(parameter_name == ""):
        name_string = AddTabs(1) + "this->mParameterNames.push_back(\"" + parameter_id + "\");\n"
    else:
        name_string = AddTabs(1) + "this->mParameterNames.push_back(\"" + parameter_name + "\");\n"

    #Define the string for the units
    if (parameter.isSetUnits()): 
        parameter_units = species.getUnits()
        units_string = AddTabs(1) + "this->mParameterUnits.push_back(\"" + parameter_units + "\");\n"
    else:
        units_string = AddTabs(1) + "this->mParameterUnits.push_back(\"non-dim\");\n"

    initial_information_string += name_string + units_string + "\n"

    return initial_information_string


def GetInitialInformationString(model):
    """Get all the information required to intiialise the system as a string. """
    initial_information_string = "" #Initialise string

    num_species = model.getNumSpecies()

    for i in range(num_species):
        species = model.getSpecies(i)
        species_string = GetInitialInformationForASpecies(species, model)
        initial_information_string += species_string

    num_parameters = model.getNumParameters()

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        parameter_id = parameter.getId()
        parameter_name = parameter.getName()

        if ( (not IsParameterValueSet(i, model))&(not IsVariableDefinedAsRule(parameter_id, model)) ) :
            parameter_string = GetInitialInformationForAParameter(parameter, model)
            initial_information_string += parameter_string
        ############### EDITED ###################
        elif(parameter_name):
            if (('wnt' in parameter_name) or ('gamma' in parameter_name) or ('ComplexTransit') in parameter_name):
                parameter_string = GetInitialInformationForAParameter(parameter, model)
                initial_information_string += parameter_string

    return initial_information_string

def GetDefaultInitialConditionString(model):
    """ Get the default initial conditions. """
    default_init_concs = ""

    num_species = model.getNumSpecies()
    species_count = 0
    for i in range(num_species):
        species = model.getSpecies(i)
        species_id = species.getId()
#        if(IsSpeciesDefinedAsOde(species_id, model)): ###################### EDITED #################################
        if(species.isSetInitialAmount()):
            init_conc = species.getInitialAmount()
        else:
            init_conc = species.getInitialConcentration()
        default_init_concs += AddTabs(1) + "SetDefaultInitialCondition(" + str(species_count) + ", " + str(init_conc) + ");\n"
        species_count += 1

    default_init_concs += "\n"

    num_parameters = model.getNumParameters()
    for i in range(num_parameters):
        parameter = model.getParameter(i) 
        parameter_id = parameter.getId()
        parameter_name = parameter.getName()
        #If the parameter has no set value, nor is it defined by a reaction rule

        ################ EDITED ######################### if wnt or gamma in parameter name, initialise as state parameter
        
        if(parameter_name):
            if (('gamma' in parameter_name)):
                default_init_concs += AddTabs(1) + "this->mParameters.push_back(1.0);" + "//" + parameter_name + "\n"
            if (('wnt' in parameter_name)):
                default_init_concs += AddTabs(1) + "this->mParameters.push_back(0.0);" + "//" + parameter_name + "\n"
            if (('ComplexTransit' in parameter_name)):
                default_init_concs += AddTabs(1) + "this->mParameters.push_back(1.0);" + "//" + parameter_name + "\n"

        elif ( (not IsParameterValueSet(i, model))&(not IsVariableDefinedAsRule(parameter_id, model)) & (not('gamma' in parameter_name)) & (not ('wnt' in parameter_name)) ):

            #Just define the value as 0, it'll be overwritten soon 
            default_init_concs += AddTabs(1) + "this->mParameters.push_back(0.0);" + "//" + parameter_name + "\n"

    return default_init_concs

def GetInitParameterStringForHeader(n, model):
    """ Method to initialise definitions of parameters in the header filei. """
    num_compartments = model.getNumCompartments()
    num_parameters = model.getNumParameters()
    num_rules = model.getNumRules()
    num_events = model.getNumEvents()

    init_string = ""
    #Initialise the compartment values
    init_string += GetBlockCommentDefinition(n, "Initialise compartments and values.", True)
    for i in range(num_compartments):
        compartment_id = GetCompartmentName(i, model)
        if (compartment_id == ''):
            compartment_id = GetCompartment(i,model)
        init_string += AddTabs(n) +  "double " + compartment_id + ";\n"

    init_string += "\n"

    #Initialise the parameters
    ########################### EDITED ##################################
    init_string += GetBlockCommentDefinition(n, "Initialise model parameters.", True)
    for i in range(num_parameters):
        parameter_id = GetParameter(i, model)
        parameter_name = GetParameterName(i,model)
        if (parameter_name == ''):
            parameter_name = parameter_id
        init_string += AddTabs(n) + "double " + parameter_name + ";\n"

    init_string += "\n"

    #Initialise rules ################# EDITED ############################ commented out
    # init_string += GetBlockCommentDefinition(n, "Initialise reaction rules (so that we may call on them later).", True)
    # for i in range(num_rules):
    #     rule = model.getRule(i)
    #     rule_id = rule.getId()
    #     if (not rule.isParameter()):
    #         init_string += AddTabs(n) + "double " + rule_id + ";\n"


    #Add the vector of booleans to check events
    if(num_events != 0): #Add events vector if model is a CCM
        init_string += "\n"
        init_string += GetBlockCommentDefinition(n, "Vector to check if all events have been triggered.", True)
        init_string += AddTabs(n) + "std::vector<bool> eventsSatisfied;\n"

    return init_string

def GetInitParameterStringForSource(n, model):
    """ Method to initialise parameter values in teh source file. """
    num_compartments = model.getNumCompartments()
    num_parameters = model.getNumParameters()
    num_funct_defns = model.getNumFunctionDefinitions()
    num_events = model.getNumEvents()

    init_string = ""
    init_string += GetBlockCommentDefinition(n, "Initialise the parameters.", True)

    ####################### EDITED ####################################
    for i in range(num_compartments):
        compartment_id = GetCompartment(i, model)
        compartment_name = GetCompartmentName(i,model)
        if (compartment_name == ''):
            compartment_name = compartment_id
        compartment_value = GetCompartmentValue(i, model)

        init_string += AddTabs(n) +  compartment_name + " = " + str(compartment_value) + ";\n"

    for i in range(num_parameters):
        parameter_id = GetParameter(i, model)
        parameter_value = GetParameterValue(i, model)
        parameter_name = GetParameterName(i,model)
        if (parameter_name == ''):
            parameter_name = parameter_id
        init_string += AddTabs(n) +  parameter_name + " = " + str(parameter_value) + ";\n"

    #Add the vector of booleans to check events
    if(model.getNumEvents() != 0):
        init_string += "\n"
        init_string += GetBlockCommentDefinition(n, "Initialise vector to check if events have been triggered.", True)
        init_string += AddTabs(n) + "eventsSatisfied.resize(" + str(num_events) + ", false);\n"


    return init_string

###################################################################################################
#####################            Cell Cycle Model Functions         ###############################
###################################################################################################

def GetEventTriggerString(event):
    """ Function to get the condition required to trigger an Event. Output is a boolean string
    that can go straight into C++. """
    trigger_string = "" # Initialise string

    # Get the AST describing the trigger
    event_trigger = event.getTrigger()
    trigger_ast = event_trigger.getMath() 

    preorder_list = trigger_ast.getListOfNodes() # Get the list of nodes in preorder (DFS)

    # Get the list of nodes in-order, so that we can readily generate the condition
    inorder_list = []
    inorder_list = GetInOrderTraversal(inorder_list, trigger_ast, preorder_list)

    # Generate the string using in-order list
    for i in range(len(inorder_list)):
        ast_node = preorder_list.get(inorder_list[i])
        variable_string = GetVariableString(ast_node)
        trigger_string += variable_string

    #We need to replace the species variable names with their equivalents in Chaste.
    original_model = event.getModel()
    
    trigger_string = ConvertVariableNamesIntoChasteEquivalent(trigger_string, original_model)
    return trigger_string

def GetEventAssignmentString(event_assignment):
    """ Get string equivalent for an event assignment. """
    event_assignment_string = ""
    original_model = event_assignment.getModel()

    event_assignment_variable = " " + event_assignment.getVariable() + " " #Get the variable name 
    event_assignment_variable = ConvertVariableNamesIntoChasteEquivalent(event_assignment_variable, original_model) #Convert into Chaste code
    event_assignment_variable = event_assignment_variable.replace("rY", "this->rGetStateVariables()") #Replace rY with modifiable state
    event_assignment_variable = event_assignment_variable.replace(" ", "") #Remove whitespace

    event_assignment_math = event_assignment.getMath() #Get the AST for the maths
    assignment_formula = " " + formulaToString(event_assignment_math)
    assignment_formula = ConvertVariableNamesIntoChasteEquivalent(assignment_formula, original_model) #Convert into Chaste code

    event_assignment_string = AddTabs(2) + event_assignment_variable + " = double ( " + assignment_formula.strip() + " );\n"

    return event_assignment_string

def GetEventAssignments(event):
    """ Get all the event assignments that occur during an event. """
    assignments_string = ""

    for i in range(event.getNumEventAssignments()):
        event_assignment = event.getEventAssignment(i)
        event_assignment_string = GetEventAssignmentString(event_assignment)
        assignments_string += event_assignment_string

    return assignments_string

def GetTriggerAndAssignmentString(n, event):
    """ Get the trigger and assignment for an event. """

    # Get the condition for the trigger
    trigger_string = GetEventTriggerString(event)

    # Get the assignment upon the trigger
    event_assignment = GetEventAssignments(event)

    condition_string = (AddTabs(1) + "if ( " + trigger_string + " )\n" + 
                        AddTabs(1) + "{\n" + 
                        event_assignment +  
                        AddTabs(2) + "eventsSatisfied[" + str(n) + "] = true;\n" + 
                        AddTabs(1) + "}\n")

    return condition_string

def GetAllEventsString(model):
    """ Generate string to describe all the events, their triggers and assignments. """
    num_events = model.getNumEvents()
    all_events = ""

    for i in range(num_events):
        event = model.getEvent(i)
        trigger_and_assignment = GetTriggerAndAssignmentString(i, event)
        all_events += trigger_and_assignment + "\n"

    return all_events#!