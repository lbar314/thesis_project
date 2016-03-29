#!/bin/bash

#Questo script copia le bacro da thesis_project alla cartella corrente

cp $HOME/thesis_project/Topology.* . || printf "Failed to copy Topology.*\n"
cp $HOME/thesis_project/TopDatabase.* . || printf "Failed to copy TopDatabase.*\n"
cp $HOME/thesis_project/MinimTopology.* . || printf "Failed to copy MinimTopology.*\n"
cp $HOME/thesis_project/MinimDatabase.* . || printf "Failed to copy MinimDatabase.*\n"
cp $HOME/thesis_project/Dictionary.* . || printf "Failed to copy Dictionary.*\n"
cp $HOME/thesis_project/BuildDictionary.* . || printf "Failed to copy BuildDictionary.*\n"
cp $HOME/thesis_project/LookUp.* . || printf "Failed to copy Lookup.*\n"
cp $HOME/thesis_project/FastSimulation.* . || printf "Failed to copy FastSimulation.*\n"
cp $HOME/thesis_project/debug.C . || printf "Failed to copy debug.C\n"
cp $HOME/thesis_project/debugMin.C . || printf "Failed to copy debugMin.C\n"
cp $HOME/thesis_project/compiler.C . || printf "Failed to copy compiler.C\n"
cp $HOME/thesis_project/testDecoding5.C . || printf "Failed to copy testDecoding5.C\n"
cp $HOME/thesis_project/testSimulation.C .|| printf "Failed to copy testSimulation.C\n"
