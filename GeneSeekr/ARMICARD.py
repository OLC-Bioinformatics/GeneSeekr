__author__ = 'mike knowles'

import json
import time
import os


class Card():
    '''
    CARD requires a gene # as an input class has three functions:
    Card(antidict, gene)
    resist will :return and list of antibiotics for a revelant gene
        Includes functionality to trace the depenedencies of a gene complex
        Utilizes recurssion to achieve all possible antibiotics
        .resist(genome)
    anti will :return a list of antibiotics for a given gene
        .anti()
    sens will :return a list of sensitivities for a given gene
        .sens()
    '''
    def __init__(self, antidict, index):  # Initialize class and inherit self
        self.index = index  # Defines a gene number that can be passed to functions within the class
        self.antidict = antidict

    def resist(self, genome=None):  # Begin resist function and import initialized self
        resistlist = []  # Initialize list
        if "resist" in self.antidict[self.index]:  # If the key "resist" in gene
            if "complex" in self.antidict[self.index] and genome is not None:
                '''If the key complex in gene defines their are depenedenies'''
                count = 0  # for each resistance set count at zero
                for complex in self.antidict[self.index]["complex"]:  # Allow for multiple dependencies
                    count += 1
                    if len(self.antidict[self.index]["complex"]) >= count:  # if complexes are satisfied
                        resistlist.extend(self.antidict[self.index]["resist"])  # extend the list
            else:  # if no complex then just return the list
                resistlist.extend(self.antidict[self.index]["resist"])
        if "isa" in self.antidict[self.index]:  # Recursion for parent antibiotic traits
            for depend in self.antidict[self.index]["isa"]:
                for amr in Card(self.antidict, depend).resist(genome):  # Call self to recurse through the same class
                    if amr not in resistlist:
                        resistlist.append(amr)
                # resistlist.extend(self.resist(genome))  # Call self to recurse through the same class
        return resistlist  # return the extended list

    def anti(self):
        if "resist" in self.antidict[self.index]:
            return self.antidict[self.index]["resist"]

    def sens(self):
         if "sensitivity" in self.antidict[self.index]:
            return self.antidict[self.index]["sensitivity"]

    def function(self):
        if "function" in self.antidict[self.index]:
            return self.antidict[self.index]["function"]


class dictbuild():
    '''
    Simple class to build a list or dictionary without repeats
    '''
    def __init__(self, index):
        self.key = index

    def add(self, lst):
        if not self.key:
            self.key = lst
        else:
            for drug in lst:
                if drug not in self.key:
                    self.key.append(drug)
        # self.key = sorted(self.key, key=lambda s: s.lower())
        return self.key


def decipher(plusdict, antidict, outputs):
    outputdict = {}
    for genome in sorted(plusdict):  # iterate through plus dict
        outputdict[genome] = {"resist": [], "sensitivity": [], "genes": []}
        resistance = outputdict[genome]["resist"]
        for gene in plusdict[genome]:
            analysis = Card(antidict, gene)
            # refgene = plusdict[genome][gene]
            if plusdict[genome][gene]:
                resist = analysis.resist(genome)  # check resistances
                sens = analysis.sens()  # check sensitivities
                if resist is not None:
                   outputdict[genome]["resist"] = dictbuild(outputdict[genome]["resist"]).add(resist)
                   outputdict[genome]["genes"].append(int(gene))
                if sens is not None:
                    outputdict[genome]["sensitivity"] = dictbuild(outputdict[genome]["sensitivity"]).add(sens)
        outputdict[genome]["genes"].sort()
        outputdict[genome]["resist"].sort()
        outputdict[genome]["sensitivity"].sort()
    json.dump(outputdict,
              open("%s/ARMI_CARD_results_%s.json" % (outputs, time.strftime("%Y.%m.%d.%H.%M.%S")), 'w'),
              sort_keys=True,
              indent=4,
              separators=(',', ': '))
    antilist = []


    for gene in antidict:  # build hearder list
        resistances = Card(antidict, gene).anti()
        if resistances is not None:
            for resist in resistances:
                if resist not in antilist:
                    antilist.append(resist)
    antihead = "Genome"
    drugcounter = {}
    antilist = sorted(antilist, key=lambda s: s.lower())  # sort header case insensitive
    for anti in antilist:
        antihead += ",\"%s\"" % anti
        drugcounter[anti] = 0

    antistr = ""

    ''' Build csv string '''
    for genome in sorted(outputdict):
        genomename = '\n{}'.format(os.path.split(os.path.splitext(genome)[0])[1].replace('_filteredAssembled', ""))
        antistr += genomename
        genomecount = 0
        for drug in antilist:
            if drug in outputdict[genome]["resist"]:
                antistr += ",+"
                drugcounter[drug] += 1
                genomecount += 1
            else:
                antistr += ",-"
        antistr += ",%i" % genomecount
    antihead += "\nCount"
    for drug in antilist:
        antihead += ",%i" % drugcounter[drug]
    antihead += antistr

    with open("%s/ARMI_CARD_results_%s.csv" % (outputs, time.strftime("%Y.%m.%d.%H.%M.%S")), 'w') as f:
        f.write(antihead)


