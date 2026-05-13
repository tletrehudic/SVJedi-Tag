#!/usr/bin/env python3

#*****************************************************************************
#  Name: SVJedi-Tag
#  Description: Genotyping of SVs with linked-reads data
#  Copyright (C) 2025 INRIA
#  Author: Anne Guichard, Mélody Temperville
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************

"""
Module 'predict_genotype.py': Process the barcodes signal to genotype the SVs.
"""

from __future__ import print_function
import argparse
import os
import re
import sys
import pickle
import math
import statistics                   #analysis
from decimal import Decimal

from pgGraphs.graph import Graph
from collections import deque
from collections import defaultdict
from pgGraphs.abstractions import Orientation

#pylint: disable=line-too-long, disable=trailing-whitespace, disable=too-many-function-args

#################
# Main function.
#################

def main(args):
    """
    Main method
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a",
        "--gaf",
        metavar="<alignmentGAFFile>",
        type=str,
        nargs=1,
        required=True)

    parser.add_argument(
        "-g",
        "--gfa",
        metavar="<GrapheGFAFile>",
        type=str,
        nargs=1,
        required=True)

    parser.add_argument(
        "-v",
        "--vcf",
        metavar="<inputVCF>",
        type=str,
        nargs=1,
        required=True)
    
    parser.add_argument(
        "-s",
        "--regionSize",
        metavar="<regionSize>",
        type=int,
        required=True)
    
    parser.add_argument(
        "-rs",
        "--regionStart",
        metavar="<regionSize>",
        type=int,
        required=True)

    parser.add_argument(
        "-o", 
        "--output", 
        metavar="<outputVCF>", 
        type=str,
        nargs=1,
        required=True)

    args = parser.parse_args()

    inputGFA = args.gfa[0]
    inputGAF = args.gaf[0]
    inputVCF = args.vcf[0]
    regionSize = args.regionSize
    regionStart = args.regionStart
    outputVCF = args.output[0]


    # Load the dictionary 'chromDict' from pickle file.
    pickleFile = str(inputGAF).rsplit("_vgGiraffe.gaf", maxsplit=1)[0] + "_chromDict.pickle"
    with open(pickleFile, "rb") as pf:
        chromDict = pickle.load(pf)

    gfa_graph:Graph = Graph(
        gfa_file=inputGFA,
        with_sequence=False,
        with_reverse_edges=True,
        low_memory=False,
    )
    gfa_graph.compute_orientations()

    svsDict = {}
    gfaNode2svRegionsDict = defaultdict(list)
    for chr, chrObject in chromDict.items() :
        for sv in chrObject.svs:

            # Regions of interest are created from nodes directly adjacent to the inversion breakpoint.
            # If the node is smaller than the set region size, then it is created using a deep graph traversal.

            ## adjLeft.
            create_region(sv,  sv.gfaNodes[0], Orientation.REVERSE ,regionSize, "adjLeft", gfaNode2svRegionsDict, gfa_graph, regionStart)              

            ## adjRight.
            create_region(sv,  sv.gfaNodes[-1], Orientation.FORWARD ,regionSize, "adjRight", gfaNode2svRegionsDict, gfa_graph,regionStart)
            
            ## nodeSVbegin.
            create_region(sv,  sv.gfaNodes[1], Orientation.FORWARD ,regionSize, "nodeSVbegin", gfaNode2svRegionsDict,gfa_graph,regionStart)

            ## nodeSVend.
            create_region(sv,  sv.gfaNodes[-2], Orientation.REVERSE ,regionSize, "nodeSVend", gfaNode2svRegionsDict, gfa_graph,regionStart)

            svsDict[sv.id] = sv

    # Create a file containing the analysis results.
    analysis_file = str(inputGAF).rsplit("_vgGiraffe.gaf", maxsplit=1)[0] + "_analysis.txt"
    analysisFile = open(analysis_file, "w", encoding='UTF-8')
    analysisFile.write("\t".join(["SV", "Method", "Genotype", "Alt_Allelic_Freq", "NbBarc_tot", "NbAlns_tot", "NbBarc_allele0", "NbBarc_allele1", "NbBarc_alleleNA", "NbAlns_allele0", "NbAlns_allele1", "NbAlns_alleleNA", "NbBarc_adjLeft", "NbBarc_adjRight", "NbBarc_nodeSVbegin", "NbBarc_nodeSVend", "NbAlns_adjLeft", "NbAlns_adjRight", "NbAlns_nodeSVbegin", "NbAlns_nodeSVend", "MinOccAlns_adjLeft", "MinOccAlns_adjRight", "MinOccAlns_nodeSVbegin", "MinOccAlns_nodeSVend", "MaxOccAlns_adjLeft", "MaxOccAlns_adjRight", "MaxOccAlns_nodeSVbegin", "MaxOccAlns_nodeSVend", "MeanOccAlns_adjLeft", "MeanOccAlns_adjRight", "MeanOccAlns_nodeSVbegin", "MeanOccAlns_nodeSVend", "MedianOccAlns_adjLeft", "MedianOccAlns_adjRight", "MedianOccAlns_nodeSVbegin", "MedianOccAlns_nodeSVend"])+"\n")

    #########################
    #B. Process aln results.
    #########################
    compteur_30 = 0
    compteur_20 = 0
    compteur_10 = 0
    compteur_total = 0
    with open(inputGAF, "r", encoding='UTF-8') as gafFile:
        compteurA = 0
        compteurB = 0 
        compteurC = 0
        compteurG = 0
        compteurR = 0 
        for line in gafFile:
            compteur_total += 1
            readID, readLen, __, __, __, path, __, pos_start, pos_end, __, alnLen, mapq, *__ = line.split("\t")
            #readID, __, __, readLen, __, __, __, path, __, pos_start, pos_end, __, alnLen, mapq, *__ = line.split("\t")
            
            #1. Filters to keep only the valid alignments.
            ##############################################
            if path == "*":             #remove unmapped reads
                continue
            cov = int(alnLen) / int(readLen)
            # if cov < 0.9:
            #     continue
            if int(mapq) < 20:
                compteur_20 += 1
            if int(mapq) < 10:
                compteur_10 += 1
            if int(mapq) < 30:
                compteur_30 += 1
                continue

            #2. Get the barcode ID.
            #######################
            if "BX:Z:" in readID:
                barcodeID = "BX:Z:" + ''.join(readID.split("BX:Z:")[1]).split(" ")[0]
                readID = readID.split("BX:Z:")[0]
            else:
                barcodeID = ""


            #3. Get a set of SV regions where at least one vgNode of the path belongs to.
            #############################################################################
            
            list_way_node = extract_nodes(path)       
            for way, node in list_way_node:
                if way == "forward" :
                    compteurA += 1 
                    if node[0] in gfaNode2svRegionsDict.keys() :
                        list_of_sv_regiontype_coords = gfaNode2svRegionsDict[node[0]]            
                    else :
                        compteurB += 1
                        continue

                    for sv_regiontype_coords in list_of_sv_regiontype_coords :
                        sv, region_type, coords_region, node_length = sv_regiontype_coords
                        if ((int(coords_region[0])< int(pos_start)) and (int(coords_region[1]) > int(pos_end))) :
                            compteurC += 1
                            if region_type == "adjLeft" :
                                compteurG += 1
                                sv.adjLeft.addBarcode(barcodeID)
                            elif region_type == "adjRight":
                                compteurR += 1
                                sv.adjRight.addBarcode(barcodeID)
                            elif region_type == "nodeSVbegin" :
                                sv.nodeSVbegin.addBarcode(barcodeID)
                            elif region_type == "nodeSVend" :
                                sv.nodeSVend.addBarcode(barcodeID)

                elif way == "backward" :
                    if node[0] in gfaNode2svRegionsDict.keys() :
                        list_of_sv_regiontype_coords = gfaNode2svRegionsDict[node[0]]
                    else :
                        continue
                    for sv_regiontype_coords in list_of_sv_regiontype_coords :
                        sv, region_type, coords_region, node_length = sv_regiontype_coords
                        start_aln = int(node_length)-int(pos_start)
                        end_aln = int(node_length)-int(pos_end)
                        #Comme on est en backward, l'annotation et comme si le 0 est la fin du noeud
                        if ((int(coords_region[0])< int(end_aln)) and (int(coords_region[1]) > int(start_aln))) :
                            #Comme on est en backward start et end inverse
                            if region_type == "adjLeft" :
                                sv.adjLeft.addBarcode(barcodeID)
                            elif region_type == "adjRight":
                                sv.adjRight.addBarcode(barcodeID)
                            elif region_type == "nodeSVbegin" :
                                sv.nodeSVbegin.addBarcode(barcodeID)
                            elif region_type == "nodeSVend" :
                                sv.nodeSVend.addBarcode(barcodeID)

                # TODO : take into account the information split-reads
    print(f"reads forward : {compteurA}")
    print(compteurB)
    print(compteurC)
    print(compteurG)
    print(compteurR)
    print(f"compteur_total : {compteur_total}")
    print(f"compteur_10 : {compteur_10}")
    print(f"compteur_20 : {compteur_20}")
    print(f"compteur_30 : {compteur_30}")
    print(f"% de reads perdu avec compteur_10 : {compteur_10/compteur_total*100}")
    print(f"% de reads perdu avec compteur_20 : {compteur_20/compteur_total*100}")
    print(f"% de reads perdu avec compteur_30 : {compteur_30/compteur_total*100}")
    ###########################
    #C. Estimate the genotype.
    ###########################
    # start = time()
    compteurD = 0
    compteurE = 0
    with open(inputVCF, "r", encoding='UTF-8') as inVCF, open(outputVCF, "w", encoding='UTF-8') as outVCF:
        sv_id = 0
        for line in inVCF:
            if line.startswith("##"):
                outVCF.write(line)

            elif line.startswith("#C"):
                outVCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                outVCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">\n')
                outVCF.write('##FORMAT=<ID=AD,Number=3,Type=Integer,Description="Depth of each allele by sample (allele0, allele1, alleleNA)">\n')
                outVCF.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternative allelic frequency">\n')
                #outVCF.write(line.rstrip("\n") + "\t" + "\t".join(["FORMAT", "SAMPLE"]) + "\n")
                outVCF.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n")

            else:
                # For each SV of the input VCF file, get the corresponding SV object.
                sv_id += 1

                sv = svsDict[sv_id]
                # Analyse the barcodes signal and estimate the genotype.
                compteurD += 1
                if sv.type == "INV":
                    compteurE += 1

                    #NB. Corresponds to Method 3.

                    # Clean the set of barcodes to keep only the informative ones.
                    adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict = cleanBarcodes(sv)
                    
                    # Estimate the genotype: on #alns with only alns of barcodes specific to one allele.
                    nbBarc_support_alleles, nbAlns_support_alleles = getSupportBarcodes(adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict)
                    #result_GT, allelic_frequency_allele1 = genotype(nbAlns_support_alleles)

                    #Error probability for likelihood (according to inversion length)
                    high, medium, low = 100000, 50000, 25000
                    if int(sv.length) > high : error_proba = 0.1
                    elif int(sv.length) < high and int(sv.length) > medium  : error_proba = 0.1
                    elif int(sv.length) < medium and int(sv.length) > low  : error_proba = 0.1
                    elif int(sv.length) < low  : error_proba = 0.1

                    result_GT, likelihoods = genotype(nbAlns_support_alleles, error_proba)
                                    
                    # Get statistics.
                    nbBarc_total, nbAlns_total, nbBarc_adjLeft, nbBarc_adjRight, nbBarc_nodeSVbegin, nbBarc_nodeSVend, nbAlns_adjLeft, nbAlns_adjRight, nbAlns_nodeSVbegin, nbAlns_nodeSVend = get_statistics(nbBarc_support_alleles, nbAlns_support_alleles, adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict)

                    # Write the analysis results to the 'analysisFile'.
                    #analysisFile.write("\t".join([str(sv.id), "3", str(result_GT), str(allelic_frequency_allele1), str(nbBarc_total), str(nbAlns_total), str(nbBarc_support_alleles[0]), str(nbBarc_support_alleles[1]), str(nbBarc_support_alleles[2]), str(nbAlns_support_alleles[0]), str(nbAlns_support_alleles[1]), str(nbAlns_support_alleles[2]), str(nbBarc_adjLeft), str(nbBarc_adjRight), str(nbBarc_nodeSVbegin), str(nbBarc_nodeSVend), str(nbAlns_adjLeft), str(nbAlns_adjRight), str(nbAlns_nodeSVbegin), str(nbAlns_nodeSVend)])+"\n")
                    analysisFile.write("\t".join([str(sv.id), "3", str(result_GT), str("TODO"), str(nbBarc_total), str(nbAlns_total), str(nbBarc_support_alleles[0]), str(nbBarc_support_alleles[1]), str(nbBarc_support_alleles[2]), str(nbAlns_support_alleles[0]), str(nbAlns_support_alleles[1]), str(nbAlns_support_alleles[2]), str(nbBarc_adjLeft), str(nbBarc_adjRight), str(nbBarc_nodeSVbegin), str(nbBarc_nodeSVend), str(nbAlns_adjLeft), str(nbAlns_adjRight), str(nbAlns_nodeSVbegin), str(nbAlns_nodeSVend)])+"\n")

                    # Output the genotype in the output VCF file.
                    numbers = ",".join(str(y) for y in nbAlns_support_alleles)
                    if len(line.split("\t")) <= 8:
                        new_line = (
                            line.rstrip("\n")
                            + "\t"
                            + "GT:DP:AD:AF"
                            + "\t"
                            + result_GT
                            + ":"
                            + str(round(sum(nbAlns_support_alleles), 3))
                            + ":"
                            + str(numbers)
                            + ":"
                            + str(round(likelihoods[0],3) if likelihoods[0] is not None else 'NA')
                            + ","
                            + str(round(likelihoods[1],3) if likelihoods[1] is not None else 'NA')
                            + ","
                            + str(round(likelihoods[2],3) if likelihoods[2] is not None else 'NA')
                            #+ str(allelic_frequency_allele1)
                        )
                        outVCF.write(new_line + "\n")
                    else:
                        line_without_genotype = line.split("\t")[0:8]
                        new_line = (
                            "\t".join(line_without_genotype)
                            + "\t"
                            + "GT:DP:AD:AF"
                            + "\t"
                            + result_GT
                            + ":"
                            + str(round(sum(nbAlns_support_alleles), 3))
                            + ":"
                            + str(numbers)
                            + ":"
                            + str(round(likelihoods[0],3) if likelihoods[0] is not None else 'NA')
                            + ","
                            + str(round(likelihoods[1],3) if likelihoods[1] is not None else 'NA')
                            + ","
                            + str(round(likelihoods[2],3) if likelihoods[2] is not None else 'NA')
                            #+ str(allelic_frequency_allele1)
                        )
                        outVCF.write(new_line + "\n")

        analysisFile.close()
        print(f"Done. Output genotypes in file {outputVCF}")
        print(compteurD)
        print(compteurE)


#############
# Functions.
#############

def create_region(sv, node, orientation, region_size, region_type, gfaNode2svRegionsDict, gfa_graph, region_start):
    region_end = int(region_size + region_start)
    # If the node is smaller than the set region size, then create the region using a deep graph traversal
    if length_node(node) < region_end :
        print("dfs")
        if region_type == 'nodeSVbegin' or region_type == 'nodeSVend' :
            if region_size > int((sv.length -(2*region_start)) / 2):
                region_end = int((sv.length -(2*region_start)) / 2)
                dico_dfs_region = createSubRegion(node,orientation, gfa_graph, region_start, region_end)
            else :
                dico_dfs_region = createSubRegion(node,orientation, gfa_graph, region_start, region_end)
        else :
            dico_dfs_region = createSubRegion(node,orientation, gfa_graph, region_start, region_end)
        
        dico_dfs_region = clean_region(dico_dfs_region)
        format_region(sv,dico_dfs_region,region_type,gfaNode2svRegionsDict)

    # Otherwise create a region on the node
    else :
        print("associated")
        # associate_GFANode_To_SVRegion(sv, node,region_type, region_size, gfaNode2svRegionsDict, region_start)


def createSubRegion(node, orientation, gfa_graph, start, end):

    full_zone = createRegion_DFS(node, orientation, end, gfa_graph) #we recove all the nodes under the end region for all the path
    full_zone = clean_region(full_zone) #without test, put full_zone in clean region
    if start == 0:
        return full_zone
    forbidden_zone = createRegion_DFS(node, orientation, start, gfa_graph) #we recove all the nodes under the start region for all the path
    forbidden_zone = clean_region(forbidden_zone) #without test, put forbidden_zone in clean region


    final_region = defaultdict(list) 

    for node_id, node_segments in full_zone.items(): # for each nodes that constitute full region
        if node_id not in forbidden_zone:   # if this node not in forbiden region, we add it
            final_region[node_id].extend(node_segments)
        else:                         # if this node is in, how is it (full ?)
            node_len = length_node(node_id)
            max_f = 0                 # initialization of the prohibition borders
            min_r = node_len
            
            for segment_type_f, segment_coord_f in forbidden_zone[node_id]: #we recover the forbiden part of this node and forbidden_zone[node_id] = node_segment_f(orbidden)
                if segment_type_f == "Full": #f for forbidden, if not f, it's full_zone
                    max_f = node_len # Force total exclusion because we enter in the if "maxf >=maxr"
                    break            
                if segment_type_f == "CutF": 
                    max_f = max(max_f, segment_coord_f[1]) # border update
                if segment_type_f == "CutR": 
                    # We convert the CutR overlap to absolute coordinates
                    # If CutR is 100bp on a node of 1000, the barrier is at 900
                    min_r = min(min_r, segment_coord_f[0])

            if max_f >= min_r: #security but normaly never coming because we clean the region of both dict so if it's the case, the node is full, not with cutR and cufF
                continue

            for segment_type, segment_coord in node_segments: #for all node_segment of this node
                start_segment, end_segment = segment_coord    
                new_start = max(start_segment, max_f)         #We calculate the new coordinates of each segment of nodes so that they are only in the non-forbidden zone = clearing. 
                new_end = min(end_segment, min_r)
                if new_end > new_start:                       #Checks that a portion of the sequence remains valid after applying the exclusion bounds (max_f and min_r).
                    final_region[node_id].append(("cut", [new_start, new_end]))
    return final_region


def createRegion_DFS(node, orientation, region_size, gfa_graph):
    ''' Function that allows a graph to be traversed in depth to create a region based on a fixed size. '''

    visited = {}    # We keep the node and the orientation in which it was traversed as a function of the distance covered to reach this node.
    stack = deque()
    region ={}

    stack.append((node,0,orientation))   # Add to stack: node, distance traveled to reach this node (first node therefore 0), orientation of entry into node
    visited[f"{node}-{Orientation.FORWARD}"] = 0    
    visited[f"{node}-{Orientation.REVERSE}"] = 0  

    while stack:
        node, dist, ori = stack.pop()   # dist = "distance covered to reach this node (i.e. sum of nodes covered to reach it)"
        dist_path = 0                   # dist_path = "distance covered to reach the "end" of the node (i.e. previous nodes + node size)"
        visited[f"{node}-{ori}"]=dist
        dist_path = dist + length_node(node)

        if dist_path < region_size :    # Region size = "fixed region size"
            region[node] = []
            region[node].append(('Full',[0,length_node(node)]))

            for info_succ in gfa_graph.segments[node]["out"][ori] : # gfa_graph.segments[node]["out"][ori] -> returns the successor nodes (with outgoing edges) of the node of interest, as well as their direction of entry into the node.
                succ = info_succ[0]
                ori_succ = info_succ[1]     # Ori_succ =direction of entry into the successor node (forward ou reverse)

                if f"{succ}-{ori_succ}" in visited.keys():  # If the node has already been visited in this orientation, we check whether the current route gives a smaller distance, and if so, we keep this distance (because we want to "go as far as possible").
                    if  visited[f"{succ}-{ori_succ}"] > dist_path : # if dist successor > dist_path
                        stack.append((succ,dist_path,ori_succ))
                        visited[f"{succ}-{ori_succ}"] = dist_path

                else :
                    stack.append((succ,dist_path,ori_succ)) 
                    visited[f"{succ}-{ori_succ}"] = dist_path

        elif dist_path > region_size :  # If we reach a size larger than the set region size, we take the part of the node that allows us to reach the set size, according to its reading direction.
            if node not in region.keys():
                region[node] = []
            if ori ==  Orientation.FORWARD:
                region[node].append(('CutF',[0,region_size-dist]))
            elif ori == Orientation.REVERSE:
                region[node].append(('CutR',[length_node(node)-(region_size-dist),length_node(node)]))
        else : #dist_path = region_size
            region[node] = []
            region[node].append(('Full',[0,region_size-dist]))
    return region


def length_node(node):
    ''' Function to retrieve node size from node identifier '''
    node_start = int(str(node).split(":")[1]) - 1            
    node_end = int(str(node).split(":")[2]) 
    node_length = (node_end-node_start)
    return node_length


def clean_region(region):
    for node, list_piece in region.items():
        if len(list_piece) <= 1:
            continue
            
        if any(cut == "Full" for cut, coord in list_piece):
            region[node] = [('Full', [0, length_node(node)])]
            continue

        size_forward = 0
        size_reverse = length_node(node)

        for cut, coord in list_piece:
            if cut == "CutF":
                size_forward = max(size_forward, coord[1])
            elif cut == "CutR":
                size_reverse = min(size_reverse, coord[0])
            
        if size_forward >= size_reverse :
            region[node] = [('Full', [0, length_node(node)])]
        elif size_forward > 0 or size_reverse < length_node(node):
            if size_reverse == length_node(node):
                region[node] = [('CutF', [0, size_forward])]
            elif size_forward == 0:
                region[node] = [('CutR', [size_reverse,length_node(node)])]
            else:
                region[node] = [('CutF', [0, size_forward]), ('CutR', [size_reverse, length_node(node)])]
    return region
            

def format_region(sv,region_dico, region_type,gfaNode2svRegionsDict):
    ''' Function to switch to the expected format from dico_region post DFS '''
    for node,list_piece in region_dico.items():
            node_lenght = length_node(node)
            for piece in list_piece:
                cut,coords = piece
                if region_type == "adjLeft":
                    sv.adjLeft = sv.getAdjLeft(coords,node)
                elif region_type == "adjRight":
                    sv.adjRight = sv.getAdjRight(coords,node)
                elif region_type == "nodeSVbegin":
                    sv.nodeSVbegin = sv.getNodeSVbegin(coords,node)
                elif region_type == "nodeSVend":
                    sv.nodeSVend = sv.getNodeSVend(coords,node)

                gfaNode2svRegionsDict[node].append((sv, region_type,coords,node_lenght))

     
def associate_GFANode_To_SVRegion(sv_object, gfaNode, region_type, regionSize, gfaNode2svRegionsDict,region_start):
    """Method to associate a GFA node to a SV region."""  
    node_length = length_node(gfaNode)     #'node_start' and 'node_end' are 0-based and incl./excl. resp.
    
    # adjLeft.
    if region_type == "adjLeft":
        coords = [(node_length-(regionSize+region_start)), node_length - region_start]
        sv_object.adjLeft = sv_object.getAdjLeft(coords,gfaNode) 


    # adjRight.
    elif region_type == "adjRight":
        coords = [region_start, regionSize+region_start]
        sv_object.adjRight = sv_object.getAdjRight(coords,gfaNode)
    

    # nodeSVbegin
    elif region_type == "nodeSVbegin":
        coords = [region_start, regionSize+region_start]
        sv_object.nodeSVbegin = sv_object.getNodeSVbegin(coords,gfaNode)


    # nodeSVend
    elif region_type == "nodeSVend":
        coords = [(node_length-(regionSize+region_start)), node_length-region_start]
        sv_object.nodeSVend = sv_object.getNodeSVend(coords,gfaNode)

    gfaNode2svRegionsDict[gfaNode].append((sv_object, region_type, coords,node_length))

    if gfaNode not in gfaNode2svRegionsDict:
        gfaNode2svRegionsDict[gfaNode] = [(sv_object, region_type, coords,node_length)]
    else:
        gfaNode2svRegionsDict[gfaNode].append((sv_object, region_type, coords,node_length))

def extract_nodes(path):       
    """Method to extract the nodes contained in a path from an alignment GAF file."""                                        
    list_way_node = []
 
    if path.startswith('>'):
        way_ALN = "forward"
    elif path.startswith('<'):
        way_ALN = "backward"

    node = [node for node in re.split(r'[<>]', path) if node]
    
    list_way_node.append((way_ALN, node))    
    
    return list_way_node


# def genotype(nbAlnBarc_support_alleles):
#     """Method to return the genotype of the SV using the 'allelic_frequency'"""

# 	# Allelic frequency.
# 	####################
#     if nbAlnBarc_support_alleles[0] == 0 and nbAlnBarc_support_alleles[1] == 0 :
#         result_GT = './.'
#         allelic_frequency_allele1 = "NA"    
#     else :
#         allelic_frequency_allele1 = nbAlnBarc_support_alleles[1] / (nbAlnBarc_support_alleles[0] + nbAlnBarc_support_alleles[1])
#         if allelic_frequency_allele1 >= 0.8:	# allelic_frequency_allele1 close to 1 --> supports allele 1.
#             result_GT = "1/1"
#         elif allelic_frequency_allele1 <= 0.2:	# allelic_frequency_allele1 close to 0 --> supports allele 0.
#             result_GT = "0/0"
#         else:                                   # allelic_frequency_allele1 close to 0.5 --> supports both alleles (heterozygous).
#             result_GT = "0/1"

#     return result_GT, allelic_frequency_allele1

def genotype(nbAlnBarc_support_alleles, e):
    """Method to return the genotype of the SV using the genotype likelihood"""
    c1 = nbAlnBarc_support_alleles[0] #number of alignement supporting reference allele
    c2 = nbAlnBarc_support_alleles[1] #number of alignement supporting alternative allele

    if c1+c2 > 9 : #Filter to have minimum of 10 informatifs alignement

        #Likelihood 
        lik0 = Decimal(c1*math.log10(1-e)) + Decimal(c2*math.log10(e)) # 0/0
        lik1 = Decimal((c1+c2)*math.log10(1/2)) # 0/1
        lik2 = Decimal(c2*math.log10(1-e)) + Decimal(c1*math.log10(e)) #1/1
        L = [lik0, lik1, lik2]
        index_of_L_max = [i for i, x in enumerate(L) if x == max(L)]

        if len(index_of_L_max) == 1:
            geno_not_encoded = str(index_of_L_max[0])
            geno = encode_genotype(geno_not_encoded)
        else : geno = './.'

    else :
        geno = './.'
        lik0, lik1, lik2 = [None, None, None]

    out = [-lik0, -lik1, -lik2] if any(x is not None for x in [lik0, lik1, lik2]) else [lik0, lik1, lik2]

    return geno, out


def encode_genotype(g): 
    ''' Encode genotype from 0, 1, 2 to 0/0, 0/1, 1/1 '''
    if g == '0': genotype = "0/0"
    elif g == '1': genotype = "0/1"
    elif g == '2': genotype = "1/1"
    elif g == './.': genotype = "./."
    return genotype


def cleanBarcodes(sv_object):
    """Method to clean the set of barcodes mapping on the regions of the current SV."""
    # Remove barcodes mapping on both 'adjLeft' and 'adjRight'.

    if hasattr(sv_object.adjLeft, 'barcodesDict') and sv_object.adjLeft.barcodesDict:
        adjLeft_barcodesDict = sv_object.adjLeft.barcodesDict.copy()
    else : 
        adjLeft_barcodesDict = {}

    if hasattr(sv_object.adjRight, 'barcodesDict') and sv_object.adjRight.barcodesDict:
        adjRight_barcodesDict = sv_object.adjRight.barcodesDict.copy()
    else :
        adjRight_barcodesDict = {}


    if hasattr(sv_object.nodeSVbegin, 'barcodesDict') and sv_object.nodeSVbegin.barcodesDict :
        nodeSVbegin_barcodesDict = sv_object.nodeSVbegin.barcodesDict.copy()
    else :
        nodeSVbegin_barcodesDict = {}

    if hasattr(sv_object.nodeSVend, 'barcodesDict') and sv_object.nodeSVend.barcodesDict:
        nodeSVend_barcodesDict = sv_object.nodeSVend.barcodesDict.copy()
    else :
        nodeSVend_barcodesDict = {}

    for barcode in list(adjLeft_barcodesDict.keys()):               #pylint: disable=consider-using-dict-items
        if barcode in adjRight_barcodesDict:
            del adjLeft_barcodesDict[barcode]
            del adjRight_barcodesDict[barcode]

    return adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict


def getSupportBarcodes(adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict):
    """
    Method to get the number of alignments (and unique barcodes) that support the allele 0, and the number of alignments (and unique barcodes) that support the allele 1.
    NB: Allele 0 represented by junctions ('adjLeft' and 'nodeSVbegin') and ('nodeSVend' and 'adjRight') for INV.
    NB: Allele 1 represented by junctions ('adjLeft' and 'nodeSVend') and ('nodeSVbegin' and 'adjRight') for INV.
    """
    nbBarc_support_allele0 = 0
    nbBarc_support_allele1 = 0
    nbBarc_undetermined = 0
    nbAln_support_allele0 = 0
    nbAln_support_allele1 = 0
    nbAln_undetermined = 0

    #NB: Count only the barcodes specific to one allele.
    #NB: Count only the alns of barcodes specific to one allele.
    for barcode in list(adjLeft_barcodesDict.keys()):
        # Allele 0.
        if (barcode in nodeSVbegin_barcodesDict) and (barcode not in nodeSVend_barcodesDict):
            nbBarc_support_allele0 += 1
            nbAln_support_allele0 += adjLeft_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode]
        # Allele 1.
        elif (barcode in nodeSVend_barcodesDict) and (barcode not in nodeSVbegin_barcodesDict):
            nbBarc_support_allele1 += 1
            nbAln_support_allele1 += adjLeft_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
        # Undetermined.
        elif (barcode in nodeSVbegin_barcodesDict) and (barcode in nodeSVend_barcodesDict):
            nbBarc_undetermined += 1
            nbAln_undetermined += adjLeft_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
            
    for barcode in list(adjRight_barcodesDict.keys()):
        # Allele 0.
        if (barcode in nodeSVend_barcodesDict) and (barcode not in nodeSVbegin_barcodesDict):
            nbBarc_support_allele0 += 1
            nbAln_support_allele0 += adjRight_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
        # Allele 1.
        elif (barcode in nodeSVbegin_barcodesDict) and (barcode not in nodeSVend_barcodesDict):
            nbBarc_support_allele1 += 1
            nbAln_support_allele1 += adjRight_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode]
        # Undetermined.
        elif (barcode in nodeSVend_barcodesDict) and (barcode in nodeSVbegin_barcodesDict):
            nbBarc_undetermined += 1
            nbAln_undetermined += adjRight_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]

    for barcode in list(nodeSVbegin_barcodesDict.keys()):
        # Undetermined.
        if (barcode in nodeSVend_barcodesDict) and (barcode not in adjLeft_barcodesDict) and (barcode not in adjRight_barcodesDict):
            nbBarc_undetermined += 1
            nbAln_undetermined += nodeSVbegin_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
    
    return [nbBarc_support_allele0, nbBarc_support_allele1, nbBarc_undetermined], [nbAln_support_allele0, nbAln_support_allele1, nbAln_undetermined]


def get_statistics(nbBarc_support_alleles, nbAlns_support_alleles, adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict):
    """Method to get statistics on the results of the genotype estimation."""
    nbBarc_total = sum(nbBarc_support_alleles)
    nbAlns_total = sum(nbAlns_support_alleles)
    nbBarc_adjLeft = len(list(adjLeft_barcodesDict.keys()))
    nbBarc_adjRight = len(list(adjRight_barcodesDict.keys()))
    nbBarc_nodeSVbegin = len(list(nodeSVbegin_barcodesDict.keys()))
    nbBarc_nodeSVend = len(list(nodeSVend_barcodesDict.keys()))
    nbAlns_adjLeft = sum(list(adjLeft_barcodesDict.values()))
    nbAlns_adjRight = sum(list(adjRight_barcodesDict.values()))
    nbAlns_nodeSVbegin = sum(list(nodeSVbegin_barcodesDict.values()))
    nbAlns_nodeSVend = sum(list(nodeSVend_barcodesDict.values()))

    return nbBarc_total, nbAlns_total, nbBarc_adjLeft, nbBarc_adjRight, nbBarc_nodeSVbegin, nbBarc_nodeSVend, nbAlns_adjLeft, nbAlns_adjRight, nbAlns_nodeSVbegin, nbAlns_nodeSVend #, minOccPerBarcode_adjLeft, minOccPerBarcode_adjRight, minOccPerBarcode_nodeSVbegin, minOccPerBarcode_nodeSVend, maxOccPerBarcode_adjLeft, maxOccPerBarcode_adjRight, maxOccPerBarcode_nodeSVbegin, maxOccPerBarcode_nodeSVend, meanOccPerBarcode_adjLeft, meanOccPerBarcode_adjRight, meanOccPerBarcode_nodeSVbegin, meanOccPerBarcode_nodeSVend, medianOccPerBarcode_adjLeft, medianOccPerBarcode_adjRight, medianOccPerBarcode_nodeSVbegin, medianOccPerBarcode_nodeSVend
 
# def local_pos(path, pos_start, pos_end):
#     length_last_node = length_node(path[-1])
#     length_total = sum(length_node(node)for node in path)
#     local_end = length_last_node - (length_total-(pos_start+pos_end))
#     return local_end

def local_pos(path, pos_start, pos_end):
    node = [node for node in re.split(r'[<>]', path) if node]
    length_aln = pos_end - pos_start
    length_first_node_aln = length_node(node[0])-pos_start
    length_aln_without_first = length_aln - length_first_node_aln
    local_end = length_aln_without_first - sum(length_node(node) for node in node[1:-1])
    return local_end

##############################################
if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])

