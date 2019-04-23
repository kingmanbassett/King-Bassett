import sys
import numpy as np 
import os

def aminoAcidComp(seq):
    #numProteins = ProteinCount(seq)
    pFile = open(seq, "r")
    proteins = {}
    name = ""
    proteinCompositions = {}
    
    for line in pFile:  #make a dictionary for proteins
        line = line.strip()
        line = line.strip()
        if(line[0] == ">"):
            name = line
            proteins[line] = 0
        elif(line != ""):
            proteins[name] = line
    pFile.close()
    pFile = open(seq, "r")
    for key in proteins:    #calculates the amount of each amino acid protein
        numMs = 0; numEs = 0; numTs = 0; numLs = 0; numQs = 0
        numNs = 0; numVs = 0; numGs = 0; numKs = 0; numCs = 0
        numRs = 0; numPs = 0; numIs = 0; numSs = 0; numAs = 0
        numHs = 0; numYs = 0; numDs = 0; numFs = 0; numWs = 0
        sequence = proteins.get(key)
        size = len(sequence)
        for char in sequence:
            if(char == "M"):
                numMs += 1
            elif(char == "E"):
                numEs += 1
            elif(char == "T"):
                numTs += 1
            elif(char == "L"):
                numLs += 1
            elif(char == "Q"):
                numQs += 1
            elif(char == "N"):
                numNs += 1
            elif(char == "V"):
                numVs += 1
            elif(char == "G"):
                numGs += 1
            elif(char == "K"):
                numKs += 1
            elif(char == "C"):
                numCs += 1
            elif(char == "S"):
                numSs += 1
            elif(char == "R"):
                numRs += 1
            elif(char == "P"):
                numPs += 1
            elif(char == "I"):
                numIs += 1
            elif(char == "A"):
                numAs += 1
            elif(char == "H"):
                numHs += 1
            elif(char == "Y"):
                numYs += 1
            elif(char == "D"):
                numDs += 1
            elif(char == "F"):
                numFs += 1
            elif(char == "W"):
                numWs += 1
        proteinCompositions[key] = [size,numAs,numCs,numDs,numEs,numFs,numGs,numHs,\
                                    numIs,numKs,numLs,numMs,numNs,numPs,numQs,\
                                    numRs,numSs,numTs,numVs,numWs,numYs]
    
    for key in proteinCompositions:
        data = proteinCompositions.get(key)
        for i in range(1, 21):
            data[i] = data[i]/data[0]
        proteinCompositions[key] = data

    return proteinCompositions

def dipeptide(seq):
    pFile = open(seq, "r")
    proteins = {}
    dipeptideDic = {}
    data = {}
    for line in pFile:  #make a dictionary for proteins
        line = line.strip()
        line = line.strip()
        if(line[0] == ">"):
            name = line
            proteins[line] = 0
        elif(line != ""):
            proteins[name] = line
    for key in proteins:
        sequence = proteins.get(key)
        length = len(sequence)
        for i in range(length-1):
            dipeptideString = sequence[i] + sequence[i+1]
            if dipeptideString not in dipeptideDic:
                dipeptideDic[dipeptideString] = 1
            else:
                dipeptideDic[dipeptideString] = dipeptideDic.get(dipeptideString) + 1
        data[key] = [length,dipeptideDic]
        dipeptideDic = dict.fromkeys(dipeptideDic, 0)
    
    print(data)






    return 0

def ProteinCount(seq):  #returns the amount of protein sequences in a file 
    count = 0
    for line in seq:
        if(line[0] == ">"):
            count = count + 1
    return count

def EuclideanDistance(query, target):   #calculate ED of query and target for likelihood
    d=0
    for i in range(1,21):
        print(query[i], target[i])
        d = d + ((query[i]*100 - target[i]*100)**2)

    d = d**.5
    print(d)


def main():
    
    #infile = "Nucleolus.txt"
    #infile = "Chromatin.txt"
    #infile = "Nuclear_Lamina.txt"
    #infile = "Nuclear_Speckles.txt"
    #infile = "Nucleoplasm.txt"

    infile = "PML_BODY.txt"
    dipep = dipeptide(infile)


    #compisition = aminoAcidComp(infile)
    # strand = compisition.get(">P49959;")
    # strand2 = compisition.get(">Q01196;")
    # EuclideanDistance(strand, strand2)
    # aminoList = ["ala", "cys", "asp", "glu", "phe",\
    #              "gly", "his", "ile", "lys", "leu",\
    #              "met", "asn", "pro", "gln", "arg",\
    #              "ser", "thr", "val", "trp", "tyr"]
    # num = 0 
    # for i in range(1,21):
    #     num = strand[i]*100 + num
    #     #print(aminoList[i-1], strand[i]*100) 
main()