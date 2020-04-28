import tkinter.simpledialog
import tkinter.messagebox
import ntpath
import sys, urllib.request, urllib.parse, urllib.error, zlib
import tkinter.tix
import tkinter.ttk
from tkinter import *
import Pmw
import subprocess
import os,math,re
import string
import queue
import threading

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

class checkboxtreeview(tkinter.ttk.Treeview):
    """
        Treeview widget with checkboxes left of each item.
        The checkboxes are done via the image attribute of the item, so to keep
        the checkbox, you cannot add an image to the item.
    """


    def __init__(self, master=None, **kw):
        tkinter.ttk.Treeview.__init__(self, master, **kw)
        # checkboxes are implemented with pictures
        print((os.getcwd()))
        print((os.path.dirname(os.path.abspath(__file__))))
        scriptdir = os.path.dirname(os.path.abspath(__file__))+"/"

        print("Looking for images in "+scriptdir)
        self.im_checked = PhotoImage(file=scriptdir+'checked.gif')
        self.im_unchecked = PhotoImage(file=scriptdir+'unchecked.gif')
        self.im_tristate = PhotoImage(file=scriptdir+'tristate.gif')
        self.tag_configure("unchecked", image=self.im_unchecked)
        self.tag_configure("tristate", image=self.im_tristate)
        self.tag_configure("checked", image=self.im_checked)
        # check / uncheck boxes on click
        self.bind("<Button-1>", self.box_click, True)
        self.selectmode = None
        style = tkinter.ttk.Style(master)

        style.layout('nodotbox.Treeview.Item', 
                     [('Treeitem.padding',
                       {'children': [('Treeitem.indicator', {'side': 'left', 'sticky': ''}),
                         ('Treeitem.image', {'side': 'left', 'sticky': ''}),
                         ('Treeitem.text', {'side': 'left', 'sticky': ''})],
                       'sticky': 'nswe'})])

        self.configure(style='nodotbox.Treeview')

    def insert(self, parent, index, iid=None, **kw):
        """ same method as for standard treeview but add the tag 'unchecked'
            automatically if no tag among ('checked', 'unchecked', 'tristate')
            is given """
        if not "tags" in kw:
            kw["tags"] = ("unchecked",)
        elif not ("unchecked" in kw["tags"] or "checked" in kw["tags"]
                  or "tristate" in kw["tags"]):
            kw["tags"] = ("unchecked",)
        tkinter.ttk.Treeview.insert(self, parent, index, iid, **kw)

    def check_descendant(self, item):
        """ check the boxes of item's descendants """
        children = self.get_children(item)
        for iid in children:
            self.item(iid, tags=("checked",))
            self.check_descendant(iid)

    def check_ancestor(self, item):
        """ check the box of item and change the state of the boxes of item's
            ancestors accordingly """
        self.item(item, tags=("checked",))
        parent = self.parent(item)
        if parent:
            children = self.get_children(parent)
            b = ["checked" in self.item(c, "tags") for c in children]
            if False in b:
                # at least one box is not checked and item's box is checked
                self.tristate_parent(parent)
            else:
                # all boxes of the children are checked
                self.check_ancestor(parent)

    def tristate_parent(self, item):
        """ put the box of item in tristate and change the state of the boxes of
            item's ancestors accordingly """
        self.item(item, tags=("tristate",))
        parent = self.parent(item)
        if parent:
            self.tristate_parent(parent)

    def uncheck_descendant(self, item):
        """ uncheck the boxes of item's descendant """
        children = self.get_children(item)
        for iid in children:
            self.item(iid, tags=("unchecked",))
            self.uncheck_descendant(iid)

    def uncheck_ancestor(self, item):
        """ uncheck the box of item and change the state of the boxes of item's
            ancestors accordingly """
        self.item(item, tags=("unchecked",))
        parent = self.parent(item)
        if parent:
            children = self.get_children(parent)
            b = ["unchecked" in self.item(c, "tags") for c in children]
            if False in b:
                # at least one box is checked and item's box is unchecked
                self.tristate_parent(parent)
            else:
                # no box is checked
                self.uncheck_ancestor(parent)

    def box_click(self, event):
        """ check or uncheck box when clicked """
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        if "image" in elem:
            # a box was clicked
            item = self.identify_row(y)
            tags = self.item(item, "tags")
            if ("unchecked" in tags) or ("tristate" in tags):
                self.check_ancestor(item)
                self.check_descendant(item)
            else:
                self.uncheck_descendant(item)
                self.uncheck_ancestor(item)

    def pre_check(self, allowed_AAs, item=""):
        children = self.get_children(item)
        if children is not None and len(children) > 0:
            for iid in children:
                self.pre_check(allowed_AAs, item=iid)
            return
        if FULL_AA_TO_3_MAP[self.item(item, "text")] in allowed_AAs:
            self.check_ancestor(item)
            self.check_descendant(item)
        return

    def get_checked(self, AAList,item=""):
        children = self.get_children(item)
        if children is not None and len(children) > 0:
            for iid in children:
                self.get_checked(AAList, item=iid)
            return
        if "checked" in self.item(item, "tags"):
            text = self.item(item, "text")
            AAList.append(text)
        return
