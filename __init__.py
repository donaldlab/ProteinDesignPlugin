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
        Bruce Donald
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

import tkSimpleDialog
import tkMessageBox
import ntpath
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfilename
from pymol import cmd
from pymol import util 
from pymol import stored 
from pymol import menu
import sys, urllib, zlib
import Tix
import ttk
from Tkinter import *
import Pmw
import subprocess
import os,math,re
import string
from pymol.cgo import *
import Queue
import threading
import checkboxtreeview
import proteininteractionviewer


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

pdbfiles = {}

def __init__(self):
    self.menuBar.addmenuitem('File', 'command', 'ProteinDesignPlugin',label = 'Save Design Parameters',
                            command = lambda s=self : saveConfigFile())
    self.menuBar.addmenuitem('File', 'command', 'ProteinDesignPlugin',label = 'Load Design Parameters',
                            command = lambda s=self : loadDesignFromConfigFile())
    self.menuBar.addmenuitem('Plugin', 'command', 'ProteinInteractionViewer',label = 'ProteinInteractionViewer',
                            command = lambda s=self : proteininteractionviewer.ProteinInteractionViewer(s))
    rewritePymolMenu(self)
    proteininteractionviewer.rewritePymolMenu(self)
    approot = self.root


            


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
                  [ 1, 'export search parameters', [
                        [ 2, 'Design type', ''],
                        [ 1, 'K* (affinity)', 'cmd.keyword[\'saveConfigFile\'][0](\'markstar\')'],
                        [ 1, 'GMEC (stability)', 'cmd.keyword[\'saveConfigFile\'][0](\'gmec\')'],
                  ]],
           ]

def setBBFlex(flexType, selection):
    print "Adding "+flexType+" flexibility to "+selection
    setflexible(selection)
    bbselection = selection + " and name c+o+n+a+ca and not (mutable and name ca)"
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

def setStrand(selection):
    print "defining strand "+selection
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
    print "setting mutations for "+selection
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
    print my_dict

    if allowed_AAs is None:
        titleText = "Set allowed mutations for residues"
        if len(my_dict['mutable']) < 2:
            titleText = "Set allowed mutations for residue "+my_dict['mutable'][0][0]
        d = testDialog(approot, titleText, my_dict['mutable'])
        for residue, allowed_AAs in allowed_muts.iteritems():
            if residue in my_dict['newmutable']:
                d.mutations.pre_check(allowed_AAs)
        print allowed_muts
    else:
        for residue_pair in my_dict['mutable']:
            allowed_muts[residue_pair[0]] = allowed_AAs.split(' ')
        print allowed_muts

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
        print allowed_muts

        self.top.destroy()
        return mut_list


def setflexible(selection):
    print "setting flexiblity for "+selection
    cmd.show("sticks", selection)
    selection = selection + " and not mutable"
    scselection = selection + " and not name c+o+n+a"
    util.cbac(scselection+" and not mutable")
    if 'flexible' in cmd.get_names('public_selections'):
        cmd.select("flexible", selection+" or flexible")
    else:
        cmd.select("flexible", selection)
    cmd.set_bond('stick_radius', 0.1, 'flexible and not name c+o+n+a')
    cmd.deselect()

cmd.extend("setflexible", setflexible)

class prompt_for:


    def __init__(self, parent, titleText, choicetexts, callback, choices=None):
        print "received parent: "+str(parent)
        self.choicetexts = choicetexts
        self.choices = choicetexts
        if choices is not None:
            self.choices = choices
        self.parent = Tk()
        self.parent.withdraw()
        self.chosen_value = None
        self.callback = callback
        self.radiobuttons  = []
        top = self.top = Toplevel(self.parent)
        print "durr"
        top.title(titleText)
        Label(top, text=titleText, justify=LEFT,
        padx=20).pack()

        self.v = IntVar()
        self.v.set(0)
        for val, choicetext in enumerate(choicetexts):
            r = Radiobutton(top,
                text=choicetext,
                padx=20,
                variable=self.v,
                value=choices[val])
            r.deselect()
            self.radiobuttons.append(r)
            r.pack()
        self.radiobuttons[0].select()

        Button(top, text="Ok", command=self.ok) \
            .pack()
        print "durrhurr"

    def show(self):
        self.top.deiconify()
        print "telling "+str(self.parent)+" to wait for "+str(self.top)
        self.parent.wait_window(self.top)
        print "durrhurrhurrhurr"
        value = self.v.get()
        return value

            
    def ok(self):
        self.top.destroy()
        self.callback(run_type=self.choices[self.v.get()])


def saveConfigFile(run_type=None, save_location=None):        
    print pdbfiles
    if(run_type is None):
        prompt_for(approot, 'Please specify design type:', ['K* (recommended)', 'GMEC'], saveConfigFile, choices=['K*', 'GMEC'])
        return
    if(save_location is None):
        save_location = asksaveasfilename(initialdir="~/", title="Select Config File Save Location", \
                    filetypes=(("cfg", "*.cfg"),("All files", "*.*")))
    if(save_location == u''):
        return
    generateConfigFile(run_type, save_location)    

cmd.extend("saveConfigFile", saveConfigFile)

def generateConfigFile(runType, save_location):
    my_dict = { 'mutable' : [], 'flexible': [], 'strands':{}, 'bb_cats':[], 'bb_deeper':[]}
    print "Generating design configuration file for "+runType+" run"
    
    # check for relevant parameters
    names = cmd.get_names('public_selections')
    if "mutable" not in names and "flexible" not in names:
        print "ERROR: at least one residue must be flexible or mutable"
        return
    if "mutable" in names:
        cmd.iterate("(name ca and mutable)", "mutable.append([chain+resi, resn])", space=my_dict)
    if "flexible" in names:
        cmd.iterate("(name ca and flexible)", "flexible.append([chain+resi, resn])", space=my_dict)
    if "cats" in names:
        cmd.iterate("(name ca and cats)", "bb_cats.append(chain+resi)", space=my_dict)
    if "deeper" in names:
        cmd.iterate("(name ca and deeper)", "bb_deeper.append(chain+resi)", space=my_dict)
    index = 1
    strandName = "strand"+str(index)
    while strandName in names:
        my_dict['strands'][strandName] = []
        cmd.iterate("(name ca and "+strandName+")", strandName+".append(chain+resi)", space=my_dict['strands'])
        index = index+1
        strandName = "strand"+str(index)

    print "Mutable:"
    print my_dict['mutable']
    print "Flexible:"
    print my_dict['flexible']
    print "Backbone flexibility:"
    print "CATS:"
    print my_dict['bb_cats']
    print "DEEPer:"
    print my_dict['bb_deeper']
    print "Defined strands:"
    print my_dict['strands']

    filename = save_location
    print "Saving config file to ["+filename+"]"

    with open(filename, 'w+') as configFile:
        print "runtype "+runType
        configFile.write("# config file auto-generated by OSPREY Design Plugin for PyMol\n")
        configFile.write("runtype "+runType+"\n")
        for pdb, pdbfile in pdbfiles.iteritems():
            configFile.write("pdb "+pdb+" \""+pdbfile+"\"\n")
            print("pdb "+pdb+" \""+pdbfile+"\"\n")
        print "# mutable"
        configFile.write("# mutable residues\n")
        for (resi, resn) in my_dict['mutable']:
            print resi+" "+" ".join(allowed_muts[resi])+" addWTRotamers continuousRotamers"
            configFile.write(resi+" "+" ".join(allowed_muts[resi])+" addWTRotamers continuousRotamers\n")
        print "# flexible"
        configFile.write("# flexible residues\n")
        for (resi, resn) in my_dict['flexible']:
            print resi+" addWTRotamers continuousRotamers"
            configFile.write(resi+" addWTRotamers continuousRotamers\n")
        configFile.write("# strand definitions\n")
        for strandDef in my_dict['strands']:
            configFile.write(strandDef+" "+" ".join(my_dict['strands'][strandDef])+"\n")
            print(strandDef+" "+" ".join(my_dict['strands'][strandDef]))
        configFile.write("cats "+" ".join(my_dict['bb_cats'])+"\n")
        print("cats "+" ".join(my_dict['bb_cats'])+"\n")
        configFile.write("deeper "+" ".join(my_dict['bb_deeper'])+"\n")
        print("deeper "+" ".join(my_dict['bb_deeper'])+"\n")


def loadDesignFromConfigFile(filename=None):
    if filename == "" or filename is None:
        filename = askopenfilename(initialdir="~/", title="Select Config File to Load", \
            filetypes=(("cfg", "*.cfg"),("All files", "*.*")))
    if filename == "" or filename is None:
        return
    with open(filename, 'r+') as configFile:
        for count, line in enumerate(configFile):
            if line.startswith('#'):
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


cmd.extend("loadDesignFromConfigFile", loadDesignFromConfigFile)



def rewritePymolMenu(self):
    print "Adding design commands to context menu"
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
    print item, self.cl.getstatus(item)

def devshortcut(approot=None):
    cmd.fetch('1ubq')
    setmutable('resi 43', 'ALA')
    saveConfigFile()

cmd.extend("devshortcut",devshortcut)

def devshortcut2(approot=None):
    loadDesignFromConfigFile("D:\\Research\\pymol\\pymoltool\\test.cfg")
    

cmd.extend("devshortcut2",devshortcut2)

oldload = cmd.load 

def wrapload(filename, object='', state=0, format='', finish=1,
             discrete=-1, quiet=1, multiplex=None, zoom=-1, partial=0,
              mimic=1, object_props=None, atom_props=None, _self=cmd):
    print "loading file "+filename
    print "full path: "+os.path.abspath(filename)
    pdbid = object
    if pdbid == '':
        pdbid = path_leaf(os.path.abspath(filename)).split('.')[0]
    pdbfiles[pdbid] = os.path.abspath(filename)
    print pdbfiles
    oldload(filename, state=state, format=format, finish=finish,
             discrete=discrete, quiet=quiet, multiplex=multiplex, zoom=zoom,
             partial=partial, mimic=mimic, object_props=object_props,
             atom_props=atom_props, _self=_self)

cmd.load = wrapload

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

