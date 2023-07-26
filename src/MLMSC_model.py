import os
import numpy as np
from .species_tree import *
from .locus_tree import *
from .haplotype_tree import *
from .exception import *
from collections import defaultdict
from io import StringIO
from skbio import read
from skbio import TreeNode

class MLMSC_Model:

    def __init__(self, seed=None):
        if seed == None:
            self.__randomState = np.random.RandomState()
        else:
            self.__randomState = np.random.RandomState(seed)
        self.__speciesTree = None
        # self.__haplotypeTree = None
        # self.__locusTree = None
        self.__parameters = {}
        self.__geneSkbioTree = None
        # self.__gene_tips = 0

    @property
    def speciesTree(self):
        return self.__speciesTree

    @property
    def parameters(self):
        return self.__parameters

    @property
    def randomState(self):
        return self.__randomState

    def run(self, inputFile, seedArgs, coalescentArgs, recombinationArgs, duplicationArgs, transferArgs, 
        lossArgs, unlinkArgs, repeatNumber, hemiplasy, verbose):
        # set parameters
        self.setParameters(
            coalescent=coalescentArgs, 
            recombination=recombinationArgs,
            duplication=duplicationArgs,
            transfer=transferArgs, 
            loss=lossArgs, 
            unlink=unlinkArgs,
            repeat=repeatNumber,
            hemiplasy=hemiplasy,
            verbose=verbose)

        # read a species tree from input file
        self.readSpeciesTree(inputFile)
        d_name = str(float(duplicationArgs))
        l_name = str(float(lossArgs))
        c_name = str(1/coalescentArgs)
        name = 'D' + d_name + 'L' + l_name + 'C' + c_name

        outputDir = './output'
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        files = os.listdir(outputDir)
        if 'random_tree_' + name + '.newick' not in files:
            f = open(outputDir + '/random_tree_' + name + '.newick','w')
            f.write('')
            f.close()

        # if 'single_tree_' + name + '.newick' not in files:
        #     f = open(outputDir + '/single_tree_' + name + '.newick','w')
        #     f.write('')
        #     f.close()
        
        # if 'sub_tree_' + name + '.newick' not in files:
        #     f = open(outputDir + '/sub_tree_' + name + '.newick','w')
        #     f.write('')
        #     f.close()

        if 'gene_tree_' + name + '.newick' not in files:
            f = open(outputDir + '/gene_tree_' + name + '.newick','w')
            f.write('')
            f.close()

        # fileName = 'summary_'+ name + '.txt'
        
        # if fileName not in files:
        #     f = open(outputDir + '/' + fileName,'w')
        #     f.write('n_d,n_genes,n_species\n')
        #     f.close()

        fileName = 'summary_'+ name + '.txt'
        
        if fileName not in files:
            f = open(outputDir + '/' + fileName,'w')
            f.write('n_d,n_genes,n_species\n')
            f.close()

        i = 1 

        while i <= repeatNumber:
            if repeatNumber > 1:
                if i%1000 == 0:
                    print('Tree ' + str(i) + ' of ' + str(repeatNumber) + ':')
            # the original locus tree is the same as the species tree
            originalLocusTree = self.constructOriginalLocusTree()
            # construct the original haplotype tree according to the species tree
            originalHaplotypeTree = self.constructOriginalHaplotypeTree()

            events = originalHaplotypeTree.DTLprocess(
                locusTree=originalLocusTree, 
                haplotypeTree=originalHaplotypeTree, 
                initial=True)

            # add new loci
            geneTree, completeCount, incompleteCount = \
                originalHaplotypeTree.addNewLociShell(events=events, 
                haplotypeTree=originalHaplotypeTree, level=0, completeCount=1, incompleteCount=0)
            geneTree.readFromSkbioTree(geneTree.getSkbioTree(), rename = False)
            self.geneSkbioTree = geneTree.getSkbioTree()
            
            # self.gene_tips = 0
            # for tip in self.geneSkbioTree.tips():
            #     self.gene_tips += 1
            self.Lprocess(self.geneSkbioTree.id, distanceAboveRoot=0)

            # geneTree.readFromSkbioTree(self.geneSkbioTree, rename = False)
            # self.geneSkbioTree = geneTree.getSkbioTree()

            geneSkbioTreeTruncated = self.geneSkbioTree.deepcopy()
            geneSkbioTreeTruncated = self.cutTree(geneSkbioTreeTruncated)

            if 'loss' in geneSkbioTreeTruncated.root().name: 
                continue
            
            find_it = True
            while find_it:
                find_it = False
                for node in geneSkbioTreeTruncated.traverse():
                    if 'loss' in node.name:
                        find_it = True
                        node.parent.name = 'l_' + node.parent.name
                        children = []
                        for child in node.parent.children:  
                            if 'loss' in child.name:
                                continue
                            else:
                                children.append(child)
                        node.parent.children = children
                        break
                    
            num_dup1 = 0  
            num_dup2 = 0   
            num_dup3 = 0          
            for node in geneSkbioTreeTruncated.traverse():
                if 'd_' in node.name:
                    num_dup3 += 1
                    if 'l_' in node.name:
                        remainder = node.name.split('_locus')[1]
                        child = node.children[0]
                        if '_locus' + remainder in child.name:
                            num_dup1 += 1
                    else:
                        num_dup1 += 1
                        num_dup2 += 1
  
            geneSkbioTreeTruncated.prune()

            if not geneSkbioTreeTruncated:
                continue
            else:
                Vec = defaultdict(list)
                for node in self.speciesTree.getSkbioTree().tips():
                    Vec[node.name] = []

                survived_all = []
                survived_only = []
                add_original = False
                for node in geneSkbioTreeTruncated.traverse():
                    if node in geneSkbioTreeTruncated.tips():
                        speciesId = int(node.name.split('*')[0])
                        remainder = (node.name.split('*')[1])
                        if remainder:
                            remainder = remainder.split('_locus')[1:]
                            for e in remainder:
                                if e not in survived_all:
                                    survived_all.append(e)
                            if remainder[0] not in survived_only:
                                survived_only.append(remainder[0])
                        elif not add_original:
                            add_original = True
                            survived_all.append('-1')
                            survived_only.append('-1')
                        speciesNode = self.speciesTree.getNodeById(speciesId)
                        node.name = speciesNode.name + str(len(Vec[speciesNode.name])+1)
                        Vec[speciesNode.name].append(node.name)
                    else:
                        node.name = ''

                randomTree = None
                n_genes = 0
                names = []
                single_names = []

                for key, value in Vec.items(): 
                    if Vec[key]:
                        names.append(self.randomState.choice(Vec[key]))
                        n_genes += len(Vec[key])

                # geneSkbioTreeTruncated = geneSkbioTreeTruncated.shear(all_names)

                randomTree = geneSkbioTreeTruncated.deepcopy()
                randomTree = randomTree.shear(names)

                singleTree = None
                single_subtree = None
                # if n_genes == len(names) and len(names) >= 4:
                #     singleTree = geneSkbioTreeTruncated.deepcopy()
                # if len(single_names) >= 4:
                #     single_subtree = geneSkbioTreeTruncated.deepcopy()
                #     single_subtree = single_subtree.shear(single_names)
                
                tipNumber = 0
                full_tip_num = 0
                single_tip_num = 0
                sub_tip_num = 0

                for node in randomTree.tips():
                    tipNumber += 1
                    node.name = node.name[0]
                for node in geneSkbioTreeTruncated.tips():
                    full_tip_num += 1
                    node.name = node.name[0]
                # if singleTree:
                #     for node in singleTree.tips():
                #         single_tip_num += 1
                #         node.name = node.name[0]
                # if single_subtree:
                #     for node in single_subtree.tips():
                #         sub_tip_num += 1
                #         node.name = node.name[0]

                if len(names) < 4:
                    continue 
                else:
                    i = i + 1

                    f = open(outputDir + '/' + fileName,'a')
                    f.write(str(len(survived_only)-1) + ',' + str(n_genes) + ',' + str(len(names)) + '\n') 

                    f = open(outputDir + '/gene_tree_' + name + '.newick','a')
                    string = str(geneSkbioTreeTruncated)
                    for char in string:
                        if char == "'":
                            continue
                        else:
                            f.write(char)
                    f.close()
                    
                    f = open(outputDir + '/random_tree_' + name + '.newick','a')
                    string = str(randomTree)
                    for char in string:
                        if char == "'":
                            continue
                        else:
                            f.write(char)
                    f.close()

                    if repeatNumber == 1:
                        print('gene tree:')    
                        print(geneSkbioTreeTruncated.ascii_art())


                    # if single_tip_num >= 4:
                    #     f = open(outputDir + '/single_tree_' + name + '.newick','a')
                    #     string =  str(singleTree)
                    #     for char in string:
                    #         if char == "'":
                    #             continue
                    #         else:
                    #             f.write(char)
                    #     f.close()

                    # if sub_tip_num >= 4:
                    #     f = open(outputDir + '/sub_tree_' + name + '.newick','a')
                    #     string =  str(single_subtree)
                    #     for char in string:
                    #         if char == "'":
                    #             continue
                    #         else:
                    #             f.write(char)
                    #     f.close()

        print('finished.')
    def setParameters(self, coalescent, recombination, duplication, transfer, loss, unlink,
        repeat, hemiplasy, verbose):
        if coalescent is None:
            raise MLMSC_Error('missing coalescent parameter')
        self.__parameters['coalescent'] = coalescent

        if recombination is None:
            raise MLMSC_Error('missing recombination parameter')
        self.__parameters['recombination'] = recombination

        if duplication is None:
            raise MLMSC_Error('missing duplication parameter')
        self.__parameters['duplication'] = duplication

        if transfer is None:
            raise MLMSC_Error('missing transfer parameter')
        self.__parameters['transfer'] = transfer

        if loss is None:
            raise MLMSC_Error('missing loss parameter')
        self.__parameters['loss'] = loss

        if unlink is None:
            raise MLMSC_Error('missing unlink parameter')
        self.__parameters['unlink'] = unlink

        if repeat is None:
            raise MLMSC_Error('missing repeat times')
        self.__parameters['repeat'] = repeat

        if hemiplasy is None:
            raise MLMSC_Error('missing hemiplasy option')
        self.__parameters['hemiplasy'] = hemiplasy

        if verbose is None:
            raise MLMSC_Error('missing verbose option')
        self.__parameters['verbose'] = verbose

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree(randomState=self.randomState)

        self.speciesTree.initialize(path=path)

        self.speciesTree.setCoalescentRate(
            coalescentPrmt=self.parameters['coalescent'])
        
        self.speciesTree.setRecombinationRate(
            recombinationPrmt=self.parameters['recombination'])

        if self.__parameters['verbose']:
            print('species tree:')	
            print(self.speciesTree)	
            print(self.speciesTree.getSkbioTree().ascii_art())
        # else:
        #     print('species tree:')	
        #     print(self.speciesTree.getSkbioTree().ascii_art())


    def constructOriginalLocusTree(self):
        locusTree = LocusTree(randomState=self.randomState)
        locusTree.initialize(
            nodes=self.speciesTree.getNodes(), 
            skbioTree=self.speciesTree.getSkbioTree())
        locusTree.coalescentRate = self.speciesTree.coalescentRate
        locusTree.recombinationRate = self.speciesTree.recombinationRate
        return locusTree
            
    def constructOriginalHaplotypeTree(self):
        haplotypeTree = HaplotypeTree(
            randomState=self.randomState, 
            speciesTree=self.speciesTree, 
            locusTree=self.speciesTree)
        haplotypeTree.initialize(locusTree=self.speciesTree)

        if self.__parameters['verbose']:
            print('original haplotype tree:')	
            print(haplotypeTree)	
            print(haplotypeTree.getSkbioTree().ascii_art())	

        haplotypeTree.setParameters(
            duplicationPrmt=self.parameters['duplication'],
            transferPrmt=self.parameters['transfer'],
            lossPrmt=self.parameters['loss'],
            unlinkProb=self.parameters['unlink'],
            hemiplasy=self.parameters['hemiplasy'],
            verbose=self.parameters['verbose'])
        return haplotypeTree

    def cutTree(self, untruncatedGeneTree):
        root = untruncatedGeneTree.root()
        self.cutTreeRecurse(root, untruncatedGeneTree)
        return untruncatedGeneTree
    
    def cutTreeRecurse(self, node, untruncatedGeneTree):
        if node.children:
            if 'loss' in node.children[0].name:
                if 'loss' in node.children[1].name:
                    findIt = 1
                else:
                    findIt = self.cutTreeRecurse(node.children[1], untruncatedGeneTree)
            else:
                if 'loss' in node.children[1].name:
                    findIt = self.cutTreeRecurse(node.children[0], untruncatedGeneTree)
                else:
                    findIt = self.cutTreeRecurse(node.children[0], untruncatedGeneTree)*\
                        self.cutTreeRecurse(node.children[1], untruncatedGeneTree)
            if findIt:
                node.name += '_loss'
            return findIt
        else:
            if 'loss' in node.name:
                return 1
            else:
                return 0

    def Lprocess(self, geneSkbioTreeId, distanceAboveRoot=0):
        for node in self.geneSkbioTree.traverse():
            if node.id == geneSkbioTreeId:
                break
        lossRate = self.parameters['loss']
        if lossRate == 0:
            distanceL = float('inf')
        else:
            distanceL = self.randomState.exponential(
                scale=1.0 / lossRate)
        if (distanceL < distanceAboveRoot):   
            node.name += '_loss'
            # if node.children:
            #     tipNumber = 0
            #     for tip in node.tips():
            #         tipNumber += 1
            #     self.gene_tips -= tipNumber
            # else:
            #     self.gene_tips -= 1
        elif node.children:
            for child in node.children:
                distanceToChild = node.distance(child)
                self.Lprocess(
                    geneSkbioTreeId=child.id, 
                    distanceAboveRoot=distanceToChild)
