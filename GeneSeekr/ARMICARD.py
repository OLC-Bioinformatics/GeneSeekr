import json
import time
import os
from collections import defaultdict
__author__ = 'mike knowles'


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
    def __init__(self, antidict, index, plusdict=None):  # Initialize class and inherit self
        self.index = index  # Defines a gene number that can be passed to functions within the class
        self.antidict = antidict
        self.plusdict = plusdict

    def resist(self, genome=None):  # Begin resist function and import initialized self
        resistlist = []  # Initialize dict
        genedict = self.antidict[self.index]
        if "resist" in genedict:  # If the key "resist" in gene
            if "member" in genedict and genome is not None:  # check if dependencies are satisfied
                mdict = {self.index: [memb for memb in genedict['member'] if memb in self.plusdict]}
                if len(mdict[self.index]) == len(genedict['member']):  # check if the list of requirements is complete
                    mdict[self.index].sort()
                    resistlist.extend([dict((resist, mdict) for resist in genedict['resist'])])  # create list of dicts
            elif "complex" in genedict and genome is not None:
                '''If the key complex in gene defines their are depenedenies'''
                # count = 0  # for each resistance set count at zero
                for complex in genedict["complex"]:  # Allow for multiple dependencies
                    resistlist.extend(Card(self.antidict, complex, self.plusdict).resist(genome))
                    # recurse through the same class if complexes are satisfied extend the list
            else:  # if no complex then just return the list
                resistlist.extend([dict((resist, [self.index]) for resist in genedict['resist'])])
        if "isa" in genedict:  # Recursion for parent antibiotic traits
            for depend in genedict["isa"]:
                for amr in Card(self.antidict, depend).resist(genome):  # Call self to recurse through the same class
                    # if amr not in resistlist:
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
        outputdict[genome] = {"resist": defaultdict(list), "sensitivity": [], "genes": []}
        resistance = outputdict[genome]["resist"]
        for gene in plusdict[genome]:
            analysis = Card(antidict, gene, plusdict[genome])
            if plusdict[genome][gene]:

                sens = analysis.sens()  # check sensitivities
                for resist in analysis.resist(genome):  # check resistances
                    if resist is not None:
                        for aro in resist:
                            if resist[aro] not in resistance[aro]:
                                if type(resist[aro]) is dict:
                                    resistance[aro].append(resist[aro])
                                else:
                                    resistance[aro].extend(resist[aro])
                        # resistance[aro].extend([x for x in resist[aro] if x not in resistance[aro]])

                outputdict[genome]["genes"].append(tuple([gene] + plusdict[genome][gene]))
                if sens is not None:
                    outputdict[genome]["sensitivity"] = dictbuild(outputdict[genome]["sensitivity"]).add(sens)
        outputdict[genome]["genes"].sort(key=lambda tup: tup[0])
        # outputdict[genome]["resist"]
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


