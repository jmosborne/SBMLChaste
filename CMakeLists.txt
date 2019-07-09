# Example CMake build configuration file for a Chaste user project.
# See also https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CmakeBuildGuide for full documentation.

find_package(Chaste COMPONENTS cell_based)

chaste_do_project(SBMLChaste)
