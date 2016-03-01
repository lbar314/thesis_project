#!/bin/bash

#Questo script copia le bacro da thesis_project alla cartella corrente

cp $HOME/thesis_project/Topology.* . || printf "Failed to copy Topology.*\n"
cp $HOME/thesis_project/TopDatabase.* . || printf "Failed to copy TopDatabase.*\n"
cp $HOME/thesis_project/MinimTopology.* . || printf "Failed to copy MinimTopology.*\n"
cp $HOME/thesis_project/MinimDatabase.* . || printf "Failed to copy MinimDatabase.*\n"
cp $HOME/thesis_project/debug.C . || printf "Failed to copy debug.C\n"
cp $HOME/thesis_project/debugMin.C . || printf "Failed to copy debugMin.C\n"
cp $HOME/thesis_project/compiler.C . || printf "Failed to copy compiler.C\n"
