#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
'''
(C) 2008, 2009 Kerim Mansour
For licensing information please refer to license.txt
'''
from xml.dom import minidom
from xml.dom import Node
from .animate import *
from .filter import *
from .gradient import *
from .linking import *
from .script import *
from .shape import *
from .structure import *
from .style import *
from .text import *

def calculateMethodName(attr):
    name=attr
    name=name.replace(':','_')
    name=name.replace('-','_')
    name='set_'+name
    return name
    
def setAttributes(attrs,obj):
    for attr in list(attrs.keys()):
        if hasattr(obj, calculateMethodName(attr)):
            eval ('obj.'+calculateMethodName(attr))(attrs[attr].value)
        else:
            print(calculateMethodName(attr)+' not found in:'+obj._elementName)
        
def build(node_, object):
    attrs = node_.attributes
    if attrs != None:
        setAttributes(attrs, object)
    for child_ in node_.childNodes:
        nodeName_ = child_.nodeName.split(':')[-1]
        if child_.nodeType == Node.ELEMENT_NODE:
            try:
                objectinstance=eval(nodeName_) ()                
            except:
                print('no class for: '+nodeName_)
                continue
            object.addElement(build(child_,objectinstance))
        elif child_.nodeType == Node.TEXT_NODE:
            #print "TextNode:"+child_.nodeValue
            #if child_.nodeValue.startswith('\n'):
            #    print "TextNode starts with return:"+child_.nodeValue
            #else:
#            print "TextNode is:"+child_.nodeValue
            #object.setTextContent(child_.nodeValue)
            if child_.nodeValue != None:
                object.appendTextContent(child_.nodeValue)
        elif child_.nodeType == Node.CDATA_SECTION_NODE:  
            object.appendTextContent('<![CDATA['+child_.nodeValue+']]>')          
        elif child_.nodeType == Node.COMMENT_NODE:  
            object.appendTextContent('<!-- '+child_.nodeValue+' -->')          
        else:
            print("Some node:"+nodeName_+" value: "+child_.nodeValue)
    return object

#TODO: packageprefix ?
def parse(inFileName):
    doc = minidom.parse(inFileName)
    rootNode = doc.documentElement
    rootObj = svg()
    build(rootNode,rootObj)
    # Enable Python to collect the space used by the DOM.
    doc = None
    #print rootObj.getXML()
    return rootObj


