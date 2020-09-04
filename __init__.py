about = """///////////////////////////////////////////////////////////////////////////////////////////////
// protein-design-plugin.py
//
//  Version:           0.2.2
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     KER        Kyle Roberts          Duke University           ker17@duke.edu
//     JDJ        Jonathan Jou          Duke University           jj@cs.duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////


   Written by Jonathan Jou (2018)
 
   The goal of this pymol plugin is to streamline the process of
   generating protein design configuration files for use in the 
   Open Source Protein REdesign for You (OSPREY) software suite. 
   OSPREY is available online (http://www.cs.duke.edu/donaldlab/osprey.php).

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA

    Contact Info:
        :Bruce Donald
        Duke University
        Department of Computer Science
        Levine Science Research Center (LSRC)
        Durham
        NC 27708-0129 
        USA
        brd@cs.duke.edu
        http://www.cs.duke.edu/~brd/

    If you use or publish any results derived from the use of this
    program please cite:

    Jonathan D. Jou and Bruce R. Donald (2018). Protein Interaction Viewer. 
    http://www.cs.duke.edu/donaldlab/software/proteinInteractionViewer/

    Copyright (C) 2018 Jonathan D. Jou, and Bruce R. Donald

    <signature of Bruce Donald>, 2 Jul, 2013
    Bruce Donald, Professor of Computer Science
"""

import tkinter.simpledialog
import tkinter.messagebox
import ntpath
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.filedialog import askdirectory
from tkinter.simpledialog import askstring
from pymol import cmd
from pymol import util 
from pymol import stored 
from pymol import menu
import sys, urllib.request, urllib.parse, urllib.error, zlib
import tkinter.tix
import tkinter.ttk
from tkinter import *
import Pmw
import subprocess
import os,math,re
import string
from pymol.cgo import *
import queue
import threading
from . import checkboxtreeview
from . import proteininteractionviewer
import re


#DESIGN PLUGIN VARIABLES
DEFAULT_19_STRING = 'ALA MET LEU ILE TRP PHE GLY ARG LYS ASP ASN GLU GLN SER TYR HIS VAL CYS THR'
approot = None
allowed_muts = {}
FULL_AA_TO_3_MAP = {
    'Alanine': 'ALA',
    'Glycine': 'GLY',
    'Isoleucine': 'ILE',
    'Leucine': 'LEU',
    'Proline': 'PRO',
    'Valine': 'VAL',
    'Phenylalanine': 'PHE',
    'Tryptophan': 'TRP',
    'Tyrosine': 'TYR',
    'Aspartate': 'ASP',
    'Glutamate': 'GLU',
    'Arginine': 'ARG',
    'Histidine': 'HIS',
    'Lysine': 'LYS',
    'Serine': 'SER',
    'Threonine': 'THR',
    'Cysteine': 'CYS',
    'Methionine': 'MET',
    'Asparagine': 'ASN',
    'Glutamine': 'GLN'
}

THREE_TO_FULL_AA_MAP= {
    'ALA': 'Alanine',
    'GLY': 'Glycine',
    'ILE': 'Isoleucine',
    'LEU': 'Leucine',
    'PRO': 'Proline',
    'VAL': 'Valine',
    'PHE': 'Phenylalanine',
    'TRP': 'Tryptophan',
    'TYR': 'Tyrosine',
    'ASP': 'Aspartate',
    'GLU': 'Glutamate',
    'ARG': 'Arginine',
    'HIS': 'Histidine',
    'LYS': 'Lysine',
    'SER': 'Serine',
    'THR': 'Threonine',
    'CYS': 'Cysteine',
    'MET': 'Methionine',
    'ASN': 'Asparagine',
    'GLN': 'Glutamine'
}

THREE_TO_ONE_AA_MAP = {
    'ALA': 'A',
    'GLY': 'G',
    'ILE': 'I',
    'LEU': 'L',
    'PRO': 'P',
    'VAL': 'V',
    'PHE': 'F',
    'TRP': 'W',
    'TYR': 'Y',
    'ASP': 'D',
    'GLU': 'E',
    'ARG': 'R',
    'HIS': 'H',
    'HID': 'H',
    'HIE': 'H',
    'HIP': 'H',
    'LYS': 'K',
    'SER': 'S',
    'THR': 'T',
    'CYS': 'C',
    'MET': 'M',
    'ASN': 'N',
    'GLN': 'Q'
}
ONE_TO_THREE_AA_MAP = {
    'A': 'ALA',
    'G': 'GLY',
    'I': 'ILE',
    'L': 'LEU',
    'P': 'PRO',
    'V': 'VAL',
    'F': 'PHE',
    'W': 'TRP',
    'Y': 'TYR',
    'D': 'ASP',
    'E': 'GLU',
    'R': 'ARG',
    'H': 'HIS',
    'H': 'HID',
    'H': 'HIE',
    'H': 'HIP',
    'K': 'LYS',
    'S': 'SER',
    'T': 'THR',
    'C': 'CYS',
    'M': 'MET',
    'N': 'ASN',
    'Q': 'GLN'
}
pdbfiles = {}

def __init__(self):
    self.menuBar.addmenuitem('File', 'command', 'ProteinDesignPlugin',label = 'Save Design Parameters',
                            command = lambda s=self : saveConfigFile())
    self.menuBar.addmenuitem('File', 'command', 'ProteinDesignPlugin',label = 'Load Design Parameters',
                            command = lambda s=self : loadDesignFromConfigFile())
    self.menuBar.addmenuitem('File', 'command', 'ProteinDesignPlugin',label = 'Generate OSPREY CFS files',
                            command = lambda s=self : generateOSPREYCFSFiles())
    self.menuBar.addmenuitem('Plugin', 'command', 'ProteinInteractionViewer',label = 'ProteinInteractionViewer',
                            command = lambda s=self : proteininteractionviewer.ProteinInteractionViewer(s))
    rewritePymolMenu(self)
    proteininteractionviewer.rewritePymolMenu(self)
    approot = self.root

def aa3to1(AA):
    return THREE_TO_ONE_AA_MAP[AA]
cmd.extend('AA3to1', aa3to1)

def aa1to3(AA):
    return ONE_TO_THREE_AA_MAP[AA]
cmd.extend('AA1to3', aa1to3)




def gavilan_cmd(self_cmd, sele):
    rsele = repr(sele)
    return [[ 2, 'Design:'       ,''                        ],
                  [ 1, 'set/check allowed mutations', 'cmd.keyword[\'setmutable\'][0]('+rsele+')'],
                  [ 1, 'set flexible', 'cmd.keyword[\'setflexible\'][0]('+rsele+')'],
                  [ 1, 'define strand', 'cmd.keyword[\'setStrand\'][0]('+rsele+')'],
                  [ 1, 'set backbone flexibility...', [
                        [ 2, 'BB Flex', ''],
                        [ 1, 'CATS', 'cmd.keyword[\'setBBFlex\'][0](\'cats\','+rsele+')'],
                        [ 1, 'DEEPer', 'cmd.keyword[\'setBBFlex\'][0](\'deeper\','+rsele+')'],
                  ]],
                  [ 1, 'add rotation/translation', 'cmd.keyword[\'setRotTrans\'][0]('+rsele+')'],
                  [ 1, 'export search parameters', [
                        [ 2, 'Design type', ''],
                        [ 1, 'K* (affinity)', 'cmd.keyword[\'saveConfigFile\'][0](\'markstar\')'],
                        [ 1, 'GMEC (stability)', 'cmd.keyword[\'saveConfigFile\'][0](\'gmec\')'],
                  ]],
           ]

def generateOSPREYCFSFiles(mutpatten="mut", shellpattern="shell", saveDirectory=None):
    
    if(saveDirectory is None):
        saveDirectory = askdirectory()
    if saveDirectory is None:
        return
    sets = {}
    selectionString = "(error)"
    groupType = "(error)"
    for selection in cmd.get_names('selections'):
        if selection.startswith('mut'):
            groupType = "mut"
            selectionString = selection.replace("mut", "")
        elif selection.startswith('shell'):
            groupType = "flex"
            selectionString = selection.replace("shell", "")
        else: 
            continue
        selectionString = "set"+selectionString
        if selectionString not in list(sets.keys()):
            sets[selectionString] = {'pdb':"(error)", 'mut_res':[], 'flex_res':[], 'strand_defs':{},
                            'strand_flex':{
                                    'protein':{},
                                    'complexProtein':{},
                                    'ligand':{},
                                    'complexLigand':{},
                                },
                            'protein_strand':["(error)", "(error)"], 'ligand_strand':["(error)", "(error)"]}
        cmd.iterate(selection+' and name CA', "'adding '+resi+' to "+groupType+"'", space=sets)
        cmd.iterate(selection+' and name CA', selectionString+'[\''+groupType+'_res\'].append(chain+resi)', space=sets)
        cmd.iterate('last ('+selection+' and name CA)', selectionString+"['model']=model", space=sets)
        print("Found model "+sets[selectionString]['model']+" for "+selectionString)
        sets[selectionString]['pdb'] = sets[selectionString]['model']+".pdb"

    mut_res_string = "(error)"
    for setID, setDict in sets.items():
        mut_resi_list = []
        for residueID in setDict['mut_res']:
            print(residueID)
            setDict['strand_flex']['protein'][residueID] = ['MUT_LIST'] 
            setDict['strand_flex']['complexProtein'][residueID] = ['MUT_LIST'] 
            chain = residueID[0]
            resi = residueID[1:]
            cmd.iterate('first ('+setDict['model']+' and chain '+chain+' and name CA)', setID+"['protein_strand'][0] =chain+resi", space=sets)
            cmd.iterate('last ('+setDict['model']+' and chain '+chain+' and name CA)', setID+"['protein_strand'][1] =chain+resi", space=sets)
            mut_resi_list.append(resi)
        setDict['strand_defs']['protein'] = setDict['protein_strand']
        setDict['strand_defs']['complexProtein'] = setDict['protein_strand']
        mut_res_string = "+".join(mut_resi_list)
        for residueID in setDict['flex_res']:
            print("flex:"+residueID)
            chain = residueID[0]
            if chain in setDict['protein_strand'][0] or chain in setDict['protein_strand'][1]:
                print("Setting flex residue on protein: "+residueID)
                setDict['strand_flex']['protein'][residueID] = ['WT'] 
                setDict['strand_flex']['complexProtein'][residueID] = ['WT'] 
                continue
            cmd.iterate('first ('+setDict['model']+' and chain '+chain+' and name CA)', setID+"['ligand_strand'][0] =chain+resi", space=sets)
            cmd.iterate('last ('+setDict['model']+' and chain '+chain+' and name CA)', setID+"['ligand_strand'][1] =chain+resi", space=sets)
            setDict['strand_flex']['ligand'][residueID] = ['WT'] 
            setDict['strand_flex']['complexLigand'][residueID] = ['WT'] 
        if setDict['ligand_strand'][1].endswith('A'):
            setDict['ligand_strand'][1] = setDict['ligand_strand'][1][:-1]
        setDict['strand_defs']['ligand'] = setDict['ligand_strand']
        setDict['strand_defs']['complexLigand'] = setDict['ligand_strand']
        setDict['file_name'] = mut_res_string+".cfs"
    print(sets)
    for setID, setDict in sets.items():
        outString = "# This cfs file automatically generated with the PDP in PyMOL.\n" \
            + 'structure_dir = "/usr/project/dlab/Users/jj/COVID-19/runs/structures/"\n' \
            + 'complex = structure_dir+"'+setDict['pdb']+'"\n' \
            + 'protein = complex\n' \
            + 'ligand = complex\n' \
            + 'strand_defs = '+str(setDict['strand_defs']).replace("{","{\n").replace("], '","],\n'")+'\n\n' \
            + 'strand_flex = '+str(setDict['strand_flex']).replace("{","{\n").replace("}, '","},\n'").replace("], '","],\n'")+'\n\n'
        print(outString)
        fileName = saveDirectory+"/"+setDict['file_name']
        with open(fileName, 'w+') as configFile:
            configFile.write(outString)
        
cmd.extend('generateOSPREYCFSFiles', generateOSPREYCFSFiles)

def replaceResidues(selection, sourceObj=None, resurfacedSourceObj=None, thirdObj=None, thirdSel=""):
    if sourceObj is None:
        sourceObj = cmd.get_object_list()[0]
    if resurfacedSourceObj is None:
        resurfacedSourceObj = cmd.get_object_list()[1]
    thirdStr = ""
    if thirdObj is not None:
        thirdStr = " or %s and %s" % (thirdObj, thirdSel)
    cmd.save('%s_resurfaced.pdb' % sourceObj, "%s and chain A or %s and chain B and not %s or %s and %s %s" % (sourceObj, sourceObj, selection, resurfacedSourceObj, selection, thirdStr))
cmd.extend('replaceResidues', replaceResidues)

def prep5vay_1r4l():
    prep5vay("1r4l")
cmd.extend('prep5vay_1r4l', prep5vay_1r4l)
def prep5vay_1r42():
    prep5vay("1r42")
cmd.extend('prep5vay_1r42', prep5vay_1r42)

def prep5vay(sourcePDB):
    asnlist = [90, 103, 499, 546]
    objects = cmd.get_object_list()
    cmd.load('J:\\Research\\Coronavirus\\fromRyan\\4-13-2020\\'+sourcePDB+'_prep.pdb')

    for objname in objects:
        cmd.align(sourcePDB+"_prep", objname)
        cmd.select("missing", sourcePDB+"_prep and resi 90+103+481+499")
        cmd.select("helix", objname+" and chain B")
        cmd.select("full", "missing or "+objname+" and chain A or helix")
        cmd.create("merged", "full")
        cmd.save('%s_fixed.pdb' % objname, "merged")

def resurface_2r05():
    resurfaced_resi_2r05 = [375,378,379,382,383,386,387,389,390]
    selection = "resi "+"+".join([str(x) for x in resurfaced_resi_2r05])
    objects = cmd.get_object_list()
    cmd.load("J:\\Research\\Coronavirus\\designs\\resultStructures\\2r05\\resurface\\2r05_resurface.pdb")
    for sourceObject in objects:
        replaceResidues(selection, sourceObj=sourceObject, resurfacedSourceObj="2r05_resurface")
cmd.extend("resurface_2r05", resurface_2r05)

def resurface_3qf7():
    resurfaced_resi_3qf7 = [167,168,171,172,175,179,182,183,186,189,190]
    selection = "resi "+"+".join([str(x) for x in resurfaced_resi_3qf7])
    objects = cmd.get_object_list()
    cmd.load("J:\\Research\\Coronavirus\\designs\\resultStructures\\3qf7\\resurface\\3qf7_resurface.pdb")
    for sourceObject in objects:
        replaceResidues(selection, sourceObj=sourceObject, resurfacedSourceObj="3qf7_resurface")
cmd.extend("resurface_3qf7", resurface_3qf7)

def resurface_5vay():
    resurfaced_resi_5vay = [124,127,128,131,132,135,139]
    selection = "resi "+"+".join([str(x) for x in resurfaced_resi_5vay])
    objects = cmd.get_object_list()
    cmd.load("J:\\Research\\Coronavirus\\designs\\resultStructures\\5vayhelix\\resurface\\5vay_resurface.pdb")
    for sourceObject in objects:
        replaceResidues(selection, sourceObj=sourceObject, resurfacedSourceObj="5vay_resurface")
cmd.extend("resurface_5vay", resurface_5vay)

def prep2r05_1r4l():
    prep2r05("1r4l")
cmd.extend('prep2r05_1r4l', prep2r05_1r4l)

def prep2r05_1r42():
    prep2r05("1r42")
cmd.extend('prep2r05_1r42', prep2r05_1r42)

def prep2r05(sourcePDB):
    frontHalf = "chain B and resi "+"+".join(str(x) for x in range(373,382))
    backHalf = "chain B and resi "+"+".join(str(x) for x in range(382,392))
    asnlist = [90, 103, 499, 546]
    frontHelixPattern = "seq.373*"
    backHelixPattern = "seq.382*"

    objname = cmd.get_object_list()[0]
    objname2 = cmd.get_object_list()[1]
    cmd.load('J:\\Research\\Coronavirus\\fromRyan\\4-13-2020\\'+sourcePDB+'_prep.pdb')
    cmd.align(sourcePDB+"_prep", objname)
    cmd.select("missing", sourcePDB+"_prep and resi 90+103+481+499")
    cmd.select("helix", objname+" and "+frontHalf+" or "+objname2+" and "+backHalf)
    cmd.select("full", "missing or "+objname+" and chain A or helix")
    cmd.create("merged", "full")

def renumber(selection='all', start=1, startsele=None, quiet=1):
    '''
DESCRIPTION

    Set residue numbering (resi) based on connectivity.

ARGUMENTS

    selection = string: atom selection to renumber {default: all}

    start = integer: counting start {default: 1}

    startsele = string: residue to start counting from {default: first in
    selection}
    '''
    start, quiet = int(start), int(quiet)
    model = cmd.get_model(selection)
    cmd.iterate(selection, 'atom_it.next().model = model',
                space={'atom_it': iter(model.atom)})
    if startsele is not None:
        startidx = cmd.index('first (' + startsele + ')')[0]
        for atom in model.atom:
            if (atom.model, atom.index) == startidx:
                startatom = atom
                break
        else:
            print(' Error: startsele not in selection')
            raise CmdException
    else:
        startatom = model.atom[0]
    for atom in model.atom:
        atom.adjacent = []
        atom.visited = False
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])
    minmax = [start, start]

    def traverse(atom, resi):
        print("Atom: "+atom.chain+atom.resi)
        atom.resi = resi
        atom.visited = True
        for other in atom.adjacent:
            if other.visited:
                continue
            if (atom.name, other.name) in [('C', 'N'), ("O3'", 'P')]:
                minmax[1] = resi + 1
                traverse(other, resi + 1)
            elif (atom.name, other.name) in [('N', 'C'), ('P', "O3'")]:
                minmax[0] = resi - 1
                traverse(other, resi - 1)
            elif (atom.name, other.name) not in [('SG', 'SG')]:
                traverse(other, resi)
    traverse(startatom, start)
    cmd.alter(selection, 'resi = atom_it.next().resi',
              space={'atom_it': iter(model.atom)})
    if not quiet:
        print((' Renumber: range (%d to %d)' % tuple(minmax)))

cmd.extend('renumber', renumber)

def prepMolecules(selection="chain B"):
    count = 1
    for model in cmd.get_object_list():
        renumber(model+' and '+selection)
        cmd.set_name(model, 'conf'+str(count))
        count = count + 1
    for model in cmd.get_object_list():
        cmd.save(model+'.pdb', model)
cmd.extend('prepMolecules', prepMolecules)

def dumpCFSString(mutsel, flexsel):
    print("# mutable")
    cmd.iterate(mutsel+" and name ca", "print(\"'\"+chain+resi+\"':['MUT_LIST'],\")")
    print("# flexible ")
    cmd.iterate(flexsel+" and name ca", "print(\"'\"+chain+resi+\"':['WT'],\")")
cmd.extend('dumpCFSString', dumpCFSString)

def generateDockedEnsemble(origin, object_pattern, aligned_residues, dock_selection, start, finish, affix="", save_location=None):
    if(save_location is None):
        save_location = askdirectory()
    if save_location is None:
        return
    for i in range(int(start), int(finish)):
        object_name = object_pattern.replace("*",str(i))
        cmd.align(object_name+" and "+aligned_residues, origin+" and "+aligned_residues)
        cmd.save(save_location+"/"+object_name+"_peptide_docked_"+affix+".pdb", "%s or (%s and %s)" % (object_name, origin, dock_selection))

cmd.extend('generateDockedEnsemble', generateDockedEnsemble)
    

def specifier(strandInfoMap):
    return '/'+strandInfoMap['model_name']+'//'+strandInfoMap['chain_name']+'//CA'
    
def setBBFlex(flexType, selection):
    print("Adding "+flexType+" flexibility to "+selection)
    setflexible(selection)
    mutstr = ""
    if 'mutable' in cmd.get_names('public_selections'):
        mutstr = "mutable and "
    bbselection = selection + " and name c+o+n+a+ca and not ("+mutstr+"name ca)"
    cmd.show("sticks", bbselection)
    if flexType == 'cats' :
        cmd.color("deeppurple", bbselection)
    if flexType == 'deeper' :
        cmd.color("hotpink", bbselection)
    if flexType in cmd.get_names('public_selections'):
        cmd.select(flexType, selection+" or "+flexType)
    else:
        cmd.select(flexType, selection)
    cmd.group('bb_flex', 'cats deeper')
    cmd.deselect()

cmd.extend("setBBFlex", setBBFlex)

def setRotTrans(selection):
    print("Adding rotation and translation to "+selection)
    cmd.show("spheres", selection)
    if 'rot_trans' in cmd.get_names('public_selections'):
        cmd.select("rot_trans", selection+" or rot_trans")
    else:
        cmd.select("rot_trans", selection)
    cmd.deselect()

cmd.extend("setRotTrans", setRotTrans)

def setStrand(selection):
    print("defining strand "+selection)
    cmd.show("ribbon", selection)
    strandsSoFar = cmd.get_names('public_selections')
    index = 1
    strandName = "strand"+str(index)
    while strandName in strandsSoFar:
        index = index+1
        strandName = "strand"+str(index)
    cmd.select(strandName, selection)
    cmd.deselect()

cmd.extend("setStrand", setStrand)


def setmutable(selection, allowed_AAs=None):
    print("setting mutations for "+selection)
    cmd.show("sticks", selection)
    scselection = selection + " and not name c+o+n+a"
    util.cbao(scselection)
    if 'mutable' in cmd.get_names('public_selections'):
        cmd.select("mutable", selection+" or mutable")
    else:
        cmd.select("mutable", selection)
    cmd.deselect()

    my_dict = { 'mutable' : [], 'newmutable': []}
    cmd.iterate("(name ca and "+selection+")", "newmutable.append(chain+resi)", space=my_dict)
    
    cmd.iterate("(name ca and "+selection+")", "mutable.append([chain+resi, resn])", space=my_dict)
    print(my_dict)

    if allowed_AAs is None:
        titleText = "Set allowed mutations for residues"
        if len(my_dict['mutable']) < 2:
            titleText = "Set allowed mutations for residue "+my_dict['mutable'][0][0]
        d = testDialog(approot, titleText, my_dict['mutable'])
        for residue, allowed_AAs in allowed_muts.items():
            if residue in my_dict['newmutable']:
                d.mutations.pre_check(allowed_AAs)
        print(allowed_muts)
    else:
        for residue_pair in my_dict['newmutable']:
            allowed_muts[residue_pair[0]] = allowed_AAs.split(' ')
        print(allowed_muts)

cmd.extend("setmutable", setmutable)

class testDialog(Toplevel):
    def __init__(self, parent, titleText, res_list):
        
        self.res_list = res_list 
        top = self.top = Toplevel(parent)
        top.title("Mutation")

        Label(top, text=titleText).pack()

        mutations = checkboxtreeview.checkboxtreeview(top, show="tree")
        self.mutations = mutations
        mutations.insert("", 0, "1", text="All")
        mutations.insert("1", "end", "11", text="Polar")
        mutations.insert("11", "end", "111", text="Serine")
        mutations.insert("11", "end", "112", text="Threonine")
        mutations.insert("11", "end", "113", text="Asparagine")
        mutations.insert("11", "end", "114", text="Glutamine")
        mutations.insert("11", "end", "115", text="Cysteine")
        mutations.insert("1", "end", "12", text="Nonpolar")
        mutations.insert("12", "end", "121", text="Alanine")
        mutations.insert("12", "end", "122", text="Isoleucine")
        mutations.insert("12", "end", "123", text="Leucine")
        mutations.insert("12", "end", "124", text="Methionine")
        mutations.insert("12", "end", "125", text="Phenylalanine")
        mutations.insert("12", "end", "126", text="Tryptophan")
        mutations.insert("12", "end", "127", text="Tyrosine")
        mutations.insert("12", "end", "128", text="Valine")
        mutations.insert("12", "end", "129", text="Glycine")
        mutations.insert("1", "end", "13", text="Cationic")
        mutations.insert("13", "end", "131", text="Arginine")
        mutations.insert("13", "end", "132", text="Histidine")
        mutations.insert("13", "end", "133", text="Lysine")
        mutations.insert("1", "end", "14", text="Anionic")
        mutations.insert("14", "end", "141", text="Aspartate")
        mutations.insert("14", "end", "142", text="Glutamate")
        mutations.pack()

        buttons = Pmw.ButtonBox(top)
        buttons.add("Apply", command=self.ok)
        buttons.add("Cancel", command=self.top.destroy)
        buttons.pack(fill = 'x',padx=1,pady=1)

    def ok(self):
        mut_list = []
        self.mutations.get_checked(mut_list)
        mut_list_3_letter = [FULL_AA_TO_3_MAP[AA] for AA in mut_list]
        for residue_pair in self.res_list:
            allowed_muts[residue_pair[0]] = mut_list_3_letter
        print(allowed_muts)

        self.top.destroy()
        return mut_list

def setflexible(selection):
    print("setting flexiblity for "+selection)
    cmd.show("sticks", selection)
    mutaffix = ""
    if 'mutable' in cmd.get_names('public_selections'):
        mutaffix = " and not mutable"
    selection = selection + mutaffix
    scselection = selection + " and not name c+o+n+a"
    util.cbac(scselection + mutaffix)
    if 'flexible' in cmd.get_names('public_selections'):
        cmd.select("flexible", selection+" or flexible")
    else:
        cmd.select("flexible", selection)
    cmd.set_bond('stick_radius', 0.1, 'flexible and not name c+o+n+a')
    cmd.deselect()

cmd.extend("setflexible", setflexible)

class prompt_for:


    def __init__(self, parent, titleText, callback, argname, args, choices):
        print("received parent: "+str(parent))
        self.choices = choices
        self.args = args
        self.argname = argname
        self.parent = Tk()
        self.parent.withdraw()
        self.chosen_value = None
        self.callback = callback
        self.radiobuttons  = []
        top = self.top = Toplevel(self.parent)
        top.title(titleText)
        Label(top, text=titleText, justify=LEFT,
        padx=20).pack()

        v = StringVar(self.parent)
        self.v = v
        self.v.set("Nope")
        for choicetext, output in choices.items():
            r = Radiobutton(top,
                text=choicetext,
                padx=20,
                variable=v,
                value=output)
            r.deselect()
            self.radiobuttons.append(r)
            r.pack()
        self.radiobuttons[0].select()

        Button(top, text="Ok", command=self.ok) \
            .pack()

    def show(self):
        self.top.deiconify()
        print("telling "+str(self.parent)+" to wait for "+str(self.top))
        self.parent.wait_window(self.top)
        print("durrhurrhurrhurr")
        value = self.v.get()
        return value

            
    def ok(self):

        print("Index chosen: "+str(self.v.get()))
        self.args[self.argname] = self.v.get()
        self.top.destroy()
        self.callback(arglist=self.args)


def saveConfigFile(arglist=None, run_type=None, save_location=None):        
    print(pdbfiles)
    if(arglist is None):
        arglist= {}
    if("run_type" not in arglist or arglist['run_type'] is None):
        prompt_for(approot, 'Please specify design type:', saveConfigFile, \
                    "run_type", arglist, choices={'K* (recommended)':'K*', 'GMEC':'GMEC'})
        return
    if("output_format" not in arglist or arglist['output_format'] is None):
        prompt_for(approot, 'Please specify output format:', saveConfigFile, \
                    "output_format", arglist, choices={'OSPREY (YAML)':'OSPREY', 'Sylph':'Sylph'})
        return
    if(save_location is None):
        save_location = asksaveasfilename(initialdir="~/", title="Select Config File Save Location", \
                    filetypes=(("cfg", "*.cfg"),("All files", "*.*")))
    if(save_location == ''):
        return
    print("output format is "+arglist['output_format'])
    if(arglist['output_format'] == "OSPREY"):
        generateConfigFileYAML(arglist['run_type'], save_location)    
    if(arglist['output_format'] == "Sylph"):
        generateConfigFile(arglist['run_type'], save_location)    

cmd.extend("saveConfigFile", saveConfigFile)

def generateConfigFileYAML(runType, save_location):
    my_dict = { 'mutable' : [], 'flexible': [], 'strands':{}, 'bb_cats':[], 'bb_deeper':[], 'rot_trans':[]}
    print("Generating design configuration file for "+runType+" run")
    
    # check for relevant parameters
    names = cmd.get_names('public_selections')
    if "mutable" not in names and "flexible" not in names:
        print("ERROR: at least one residue must be flexible or mutable")
        return
    if "mutable" in names:
        cmd.iterate("(name ca and mutable)", "mutable.append([chain+resi, resn])", space=my_dict)
    if "flexible" in names:
        cmd.iterate("(name ca and flexible)", "flexible.append([chain+resi, resn])", space=my_dict)
    if "cats" in names:
        cmd.iterate("(name ca and cats)", "bb_cats.append(chain+resi)", space=my_dict)
    if "deeper" in names:
        cmd.iterate("(name ca and deeper)", "bb_deeper.append(chain+resi)", space=my_dict)
    if "rot_trans" in names:
        cmd.iterate("rot_trans", "rot_trans.append(chain+resi)", space=my_dict)
    index = 1
    strandName = "strand"+str(index)
    while strandName in names:
        my_dict['strands'][strandName] = []
        cmd.iterate("(name ca and "+strandName+")", strandName+".append(chain+resi)", space=my_dict['strands'])
        index = index+1
        strandName = "strand"+str(index)

    print("Mutable:")
    print(my_dict['mutable'])
    print("Flexible:")
    print(my_dict['flexible'])
    print("Backbone flexibility:")
    print("CATS:")
    print(my_dict['bb_cats'])
    print("DEEPer:")
    print(my_dict['bb_deeper'])
    print("Defined strands:")
    print(my_dict['strands'])

    filename = save_location
    print("Saving config file to ["+filename+"]")

    with open(filename, 'w+') as configFile:
        print("runtype "+runType)
        configFile.write("# config file auto-generated by OSPREY Design Plugin for PyMol\n")
        configFile.write("# this config file is formatted to meet YAML specifications\n")
        configFile.write("runtype: "+runType+"\n")
        configFile.write("structures:\n")
        for pdb, pdbfile in pdbfiles.items():
            configFile.write("    \'"+pdb+"\': \'"+pdbfile+"\'\n")
            print(("    \'"+pdb+"\': \""+pdbfile+"\"\n"))

        configFile.write("# strand definitions\n")
        configFile.write("strands:\n")
        for strandDef in sorted(my_dict['strands']):
            configFile.write("    "+strandDef+": ["+", ".join(my_dict['strands'][strandDef])+"]\n")
            print(("    "+strandDef+": ["+", ".join(my_dict['strands'][strandDef])+"]\n"))

        print("# mutable residues")
        configFile.write("# specify mutable and flexible residues\n")
        configFile.write("designed residues:\n")
        configFile.write("# mutable residues\n")
        for (resi, resn) in my_dict['mutable']:
            print("    "+resi+": ["+", ".join(allowed_muts[resi])+", addWTRotamers, continuousRotamers]")
            configFile.write("    "+resi+": ["+", ".join(allowed_muts[resi])+", addWTRotamers, continuousRotamers]\n")

        print("# flexible")
        configFile.write("# flexible residues\n")
        for (resi, resn) in my_dict['flexible']:
            print("    "+resi+": [WT, addWTRotamers, continuousRotamers]")
            configFile.write("    "+resi+": [WT, addWTRotamers, continuousRotamers]\n")

        print("# rotations and translations")
        configFile.write("# rotations and translations\n")
        print("rotation translation: ["+", ".join(list(dict.fromkeys(my_dict['rot_trans'])))+"]")
        configFile.write("rotation translation: ["+", ".join(list(dict.fromkeys(my_dict['rot_trans'])))+"]\n")

        configFile.write("# backbone flexiblity\n")
        configFile.write("bbflexibility:\n")
        if(len(my_dict['bb_cats']) > 0):
            configFile.write("    cats: ["+", ".join(my_dict['bb_cats'])+"]\n")
            print(("    cats: ["+", ".join(my_dict['bb_cats'])+"]\n"))
        if(len(my_dict['bb_deeper']) > 0):
            configFile.write("    deeper: ["+", ".join(my_dict['bb_deeper'])+"]\n")
            print(("    deeper: ["+", ".join(my_dict['bb_deeper'])+"]\n"))

def generateConfigFile(runType, save_location):
    my_dict = { 'mutable' : [], 'flexible': [], 'strands':{}, 'bb_cats':[], 'bb_deeper':[], 'rot_trans':[]}
    print("Generating design configuration file for "+runType+" run")
    
    # check for relevant parameters
    names = cmd.get_names('public_selections')
    if "mutable" not in names and "flexible" not in names:
        print("ERROR: at least one residue must be flexible or mutable")
        #return
    if "mutable" in names:
        cmd.iterate("(name ca and mutable)", "mutable.append([chain+resi, resn])", space=my_dict)
    if "flexible" in names:
        cmd.iterate("(name ca and flexible)", "flexible.append([chain+resi, resn])", space=my_dict)
    if "cats" in names:
        cmd.iterate("(name ca and cats)", "bb_cats.append(chain+resi)", space=my_dict)
    if "rot_trans" in names:
        cmd.iterate("rot_trans", "rot_trans.append(chain+resi)", space=my_dict)
    if "deeper" in names:
        cmd.iterate("(name ca and deeper)", "bb_deeper.append(chain+resi)", space=my_dict)
    index = 1
    strandName = "strand"+str(index)
    while strandName in names:
        my_dict['strands'][strandName] = []
        cmd.iterate("(name ca and "+strandName+")", strandName+".append(chain+resi)", space=my_dict['strands'])
        if "rot_trans" in names:
            cmd.iterate("(rot_trans and "+strandName+")", strandName+".append(chain+resi)", space=my_dict['strands'])
        my_dict['strands'][strandName] = list(dict.fromkeys(my_dict['strands'][strandName]))
        index = index+1
        strandName = "strand"+str(index)

    print("Mutable:")
    print(my_dict['mutable'])
    print("Flexible:")
    print(my_dict['flexible'])
    print("Backbone flexibility:")
    print("CATS:")
    print(my_dict['bb_cats'])
    print("DEEPer:")
    print(my_dict['bb_deeper'])
    print("Defined strands:")
    print(my_dict['strands'])

    filename = save_location
    print("Saving config file to ["+filename+"]")

    with open(filename, 'w+') as configFile:
        print("runtype "+runType)
        configFile.write("# config file auto-generated by OSPREY Design Plugin for PyMol\n")
        configFile.write("# this config file is compatible with Sylph\n")
        configFile.write("runtype "+runType+"\n")
        for pdb, pdbfile in pdbfiles.items():
            configFile.write("pdb "+pdb+" \""+pdbfile+"\"\n")
            print(("pdb "+pdb+" \""+pdbfile+"\"\n"))

        configFile.write("# strand definitions\n")
        for strandDef in sorted(my_dict['strands']):
            configFile.write(strandDef+" "+" ".join(my_dict['strands'][strandDef])+"\n")
            print((strandDef+" "+" ".join(my_dict['strands'][strandDef])))

        print("# mutable")
        configFile.write("# mutable residues\n")
        for (resi, resn) in my_dict['mutable']:
            print(resi+" "+" ".join(allowed_muts[resi])+" addWTRotamers continuousRotamers")
            configFile.write(resi+" "+" ".join(allowed_muts[resi])+" addWTRotamers, continuousRotamers\n")

        print("# flexible")
        configFile.write("# flexible residues\n")
        for (resi, resn) in my_dict['flexible']:
            print(resi+" addWTRotamers continuousRotamers")
            configFile.write(resi+" addWTRotamers continuousRotamers\n")

        print("# rotations and translations")
        configFile.write("# rotations and translations\n")
        print("rotationTranslation "+" ".join(list(dict.fromkeys(my_dict['rot_trans']))))
        configFile.write("rotationTranslation "+" ".join(list(dict.fromkeys(my_dict['rot_trans'])))+"\n")

        if len(my_dict['bb_cats']) > 0 :
            configFile.write("cats "+" ".join(my_dict['bb_cats'])+"\n")
            print(("cats "+" ".join(my_dict['bb_cats'])+"\n"))
        if len(my_dict['bb_deeper']) > 0 :
            configFile.write("deeper "+" ".join(my_dict['bb_deeper'])+"\n")
            print(("deeper "+" ".join(my_dict['bb_deeper'])+"\n"))


def loadDesignFromConfigFile(filename=None):
    if filename == "" or filename is None:
        filename = askopenfilename(initialdir="~/", title="Select Config File to Load", \
            filetypes=(("cfg", "*.cfg"),("All files", "*.*")))
    if filename == "" or filename is None:
        return
    with open(filename, 'r+') as configFile:
        for count, line in enumerate(configFile):
            if line.startswith('#'):
                if 'YAML' in line:
                    print("YAML FILE FOUND.")
                    loadDesignFromConfigFileYAML(filename)
                    return
                continue
            if '#' in line:
                line = line.split('#', 1)[0]
            m = re.match('pdb (\w+) "([^"]+)"', line) 
            if m is not None:
                [pdbid, pdbfile] = m.groups()
                pdbfiles[pdbid] = pdbfile
                cmd.load(pdbfile)
            if line.startswith('cats'):
                cats_res = line.split()[1:]
                for residue in cats_res:
                    setBBFlex('cats', 'chain '+residue[0]+' and resi '+residue[1:])
            if line.startswith('deeper'):
                deeper_res = line.split()[1:]
                for residue in deeper_res:
                    setBBFlex('deeper', 'chain '+residue[0]+' and resi '+residue[1:])
            m = re.match("(\w?\d+) ?((((\w\w\w)|(WT))\s)+)?(addWTRotamers)? ?(continuousRotamers)? ?(addWTRotamers)?", line) 
            if m is not None:
                [residue_number, allowed_AAs, last_AA, last_WT_or_AA, last_non_WT, WT, \
                    add_wild_type_rotamers_string, continuous_rotamers_string, \
                    add_wild_type_rotamers_backup_string] = m.groups()         
                if allowed_AAs is None or ( len(allowed_AAs) < 2 and 'WT' in allowed_AAs ):
                    setflexible('chain '+residue_number[0]+' and resi '+residue_number[1:])
                else:
                    setmutable('chain '+residue_number[0]+' and resi '+residue_number[1:], allowed_AAs=allowed_AAs.strip())
            if line.startswith('strand'):
                strand_residues = line.strip().split(" ")[1:]
                for index, residue_number in enumerate(strand_residues):
                    strand_residues[index] = "(chain "+residue_number[0]+" and resi "+residue_number[1:]+")"
                setStrand("+".join(strand_residues))
            m = re.match("(\w?\d+) rotationTranslation", line) 
            if line.startswith('rotationTranslation'):
                strand_residues = line.strip().split(" ")[1:]
                for index, residue_number in enumerate(strand_residues):
                    residue_number = "(chain "+residue_number[0]+" and resi "+residue_number[1:]+")"
                    setRotTrans(residue_number)

def loadStructures(line):
    print("loading structure from ["+line+"]")
    m = re.match('\s+\'(\w+)\': \'([^\']+)\'', line) 
    if m is not None:
        [pdbid, pdbfile] = m.groups()
        pdbfiles[pdbid] = pdbfile
        cmd.load(pdbfile)

def loadDesignResidues(line):
    print("loading residue ["+line+"]")
    m = re.match("\s+(\w?\d+): ?\[((((\w\w\w)|(WT)),\s)+)?(addWTRotamers,?)? ?(continuousRotamers,?)? ?(addWTRotamers)?\]", line) 
    if m is not None:
        [residue_number, allowed_AAs, last_AA, last_WT_or_AA, last_non_WT, WT, \
            add_wild_type_rotamers_string, continuous_rotamers_string, \
            add_wild_type_rotamers_backup_string] = m.groups()         
        allowed_AAs = allowed_AAs.replace(',','')
        print("Allowed AAs: ["+allowed_AAs+"]")
        if allowed_AAs is None or ( len(allowed_AAs) < 4 and 'WT' in allowed_AAs ):
            setflexible('chain '+residue_number[0]+' and resi '+residue_number[1:])
        else:
            setmutable('chain '+residue_number[0]+' and resi '+residue_number[1:], allowed_AAs=allowed_AAs.strip())

def loadDEEPERFlex(line):
    deeper_res = line.split()
    for residue in deeper_res:
        setBBFlex('deeper', 'chain '+residue[0]+' and resi '+residue[1:])

def loadBBFlex(line):
    print("loading bbflex from ["+line+"]")
    m = re.match("\s+((cats:)|(deeper:)) ?\[((\w?\d+,? ?)+)\]", line) 
    if m is not None:
        [flextype, cats, deeper, residues, first_res] = m.groups()
        print("flextype "+flextype)
        print("residues "+residues)
        if flextype == 'cats:':
            loadCATSFlex(residues)
        if flextype == 'deeper:':
            loadDEEPERFlex(residues)

    

def loadCATSFlex(line):
    print("loading cats flex from ["+line+"]")
    cats_res = line.split()
    for residue in cats_res:
        setBBFlex('cats', 'chain '+residue[0]+' and resi '+residue[1:])

def loadRotationTranslation(line):
    print("loading rottrans from ["+line+"]")
    residues = line.split(':')[1]
    residues = residues.replace('[', '')
    residues = residues.replace(']', '')
    residues = residues.replace(',', '')
    residues = residues.replace('\n', '')
    print("rottrans residues: "+residues)
    strand_residues = residues.strip().split(" ")
    for index, residue_number in enumerate(strand_residues):
        if len(residue_number) < 1:
            continue
        residue_number = "(chain "+residue_number[0]+" and resi "+residue_number[1:]+")"
        setRotTrans(residue_number)

def loadStrands(line):
    print("loading strands from ["+line+"]")
    strand_residues = line.strip().split(" ")[1:]
    for index, residue_number in enumerate(strand_residues):
        strand_residues[index] = "(chain "+residue_number[0]+" and resi "+residue_number[1:]+")"
    setStrand("+".join(strand_residues))

def runType(line):
    return

def loadDesignFromConfigFileYAML(filename=None):
    state_funcs = {
        'structures:':loadStructures,
        'designed residues:':loadDesignResidues,
        'bbflexibility:':loadBBFlex,
        'cats:':loadCATSFlex,
        'deeper':loadDEEPERFlex,
        'strands:':loadStrands,
        'runtype:':runType,
    }
    if filename == "" or filename is None:
        filename = askopenfilename(initialdir="~/", title="Select Config File to Load", \
            filetypes=(("cfg", "*.cfg"),("All files", "*.*")))
    if filename == "" or filename is None:
        return
    state = "normal"
    indentation = ""
    with open(filename, 'r+') as configFile:
        for count, line in enumerate(configFile):
            print("handling line [" +line+"]")
            if line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#', 1)[0]

            m = re.match('(^\S[^:]+:)', line) 
            if m is not None:
                if 'runtype:' in line:
                    continue
                if 'rotation translation:' in line:
                    loadRotationTranslation(line)
                    continue
                state = m.groups()[0]
                print("setting state to "+state)
                continue

            print("Calling state "+state)
            state_funcs[state](line)



cmd.extend("loadDesignFromConfigFile", loadDesignFromConfigFile)



def rewritePymolMenu(self):
    print("Adding design commands to context menu")
    oldAction = menu.sele_action
    oldAction2 = menu.sele_action2
    old_all_action = menu.all_action
    old_mol_action = menu.mol_action
    def newAction(self_cmd, sele):
        output = oldAction(self_cmd, sele)
        output.append([ 1, 'design', gavilan_cmd(self_cmd, sele)])
        return output
    def newAction2(self_cmd, sele):
        output = oldAction2(self_cmd, sele)
        output.append([ 1, 'design', gavilan_cmd(self_cmd, sele)])
        return output
    def new_all_action(self_cmd, sele):
        output = old_all_action(self_cmd, sele)
        output.append([ 1, 'design', gavilan_cmd(self_cmd, sele)])
        return output
    def new_mol_action(self_cmd, sele):
        output = old_mol_action(self_cmd, sele)
        output.append([ 1, 'design', gavilan_cmd(self_cmd, sele)])
        return output
    menu.sele_action = newAction
    menu.sele_action2 = newAction2
    menu.all_action = new_all_action
    menu.mol_action = new_mol_action
    makeCheckList(self)

def makeCheckList(self):
    if True:
        return
    self.cl = tixCheckList(self.root, browscmd=selectItem)
    self.cl.pack()
    self.cl = tixCheckList(self.root, browsecmd=selectItem)
    self.cl.pack()
    self.cl.hlist.add("CL1", text="checklist1")
    self.cl.hlist.add("CL1.Item1", text="subitem1")
    self.cl.hlist.add("CL2", text="checklist2")
    self.cl.hlist.add("CL2.Item1", text="subitem1")
    self.cl.setstatus("CL2", "on")
    self.cl.setstatus("CL2.Item1", "on")
    self.cl.setstatus("CL1", "off")
    self.cl.setstatus("CL1.Item1", "off")
    self.cl.autosetmode()

def selectItem(self, item):
    print(item, self.cl.getstatus(item))

def oneLetterSequence(selection):
    print("Generating one-letter sequence..")
    global iterDict
    iterDict = {'resCodes':[]}
    cmd.iterate(selection+' and name CA', "resCodes.append(resn)", space=iterDict)
    print("Iterated.")
    oneLetterList = [THREE_TO_ONE_AA_MAP[x] for x in iterDict['resCodes']]
    print("Remapped.")
    print("".join(oneLetterList))
cmd.extend("oneLetterSequence", oneLetterSequence)

def colorMutScene(index):
    mutStr = 'mut'+str(index)
    util.cbam(mutStr)
    cmd.show('sticks', mutStr)
    shellStr = 'shell'+str(index)
    util.cbao(shellStr)
    cmd.show('sticks', shellStr)
    cmd.show('lines', 'byres %s around 4' % shellStr)
cmd.extend("colorMutScene", colorMutScene)

def colorScene(index=-1):
    cmd.hide('sticks', 'all')
    cmd.hide('lines', 'all')
    util.cbc()
    util.cnc()
    if(int(index)>0):
        colorMutScene(index)
cmd.extend("colorScene", colorScene)
    
def autoMut(selection):
    nextMutIndex = 1
    nextMutSel = 'mut%d'%nextMutIndex
    while nextMutSel in cmd.get_names('selections'):
        nextMutIndex = nextMutIndex + 1
        nextMutSel = 'mut%d'%nextMutIndex
    cmd.select(nextMutSel, selection)
    cmd.select('shell%d'%nextMutIndex, 'byres %s around 6' % nextMutSel)
    colorScene(nextMutIndex)
    cmd.orient('shell%d'%nextMutIndex)
    cmd.scene(key=nextMutSel, action='store')
cmd.extend("autoMut", autoMut)
    

def devshortcut(approot=None):
    cmd.fetch('1ubq')
    setmutable('resi 43', 'ALA')
    saveConfigFile()

cmd.extend("devshortcut",devshortcut)
def devshortcut(approot=None):
    cmd.fetch('1ubq')
    setmutable('resi 43', 'ALA')
    saveConfigFile()

cmd.extend("devshortcut",devshortcut)

def devshortcut2(approot=None):
    loadDesignFromConfigFile("C:\\Users\\thepa\\Documents\\Research\\PDP\\test.cfg")
    

cmd.extend("devshortcut2",devshortcut2)

oldload = cmd.load 
def cnbc():
    util.cbc()
    util.cnc()
cmd.extend("cnbc", cnbc)

def oneletterseq(selection):
    space={'seqstring':[]}
    cmd.iterate(selection+' and name CA', 'seqstring.append(resn)', space=space)
    oneletterstring = "".join([aa3to1(x) for x in space['seqstring']])
    print(oneletterstring)
    oneletterstring = ",".join([aa3to1(x) for x in space['seqstring']])
    print(oneletterstring)
cmd.extend("getseq", oneletterseq)

def wrapload(filename, object='', state=0, format='', finish=1,
             discrete=-1, quiet=1, multiplex=None, zoom=-1, partial=0,
              mimic=1, object_props=None, atom_props=None, _self=cmd):
    print("loading file "+filename)
    print("full path: "+os.path.abspath(filename))
    pdbid = object
    if pdbid == '':
        pdbid = path_leaf(os.path.abspath(filename)).split('.')[0]
    pdbfiles[pdbid] = os.path.abspath(filename)
    print(pdbfiles)
    oldload(filename, state=state, format=format, finish=finish,
             discrete=discrete, quiet=quiet, multiplex=multiplex, zoom=zoom,
             partial=partial, mimic=mimic, object_props=object_props,
             atom_props=atom_props, _self=_self)

cmd.load = wrapload

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
