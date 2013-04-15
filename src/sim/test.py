#!/usr/bin/python

import libxml2
import sys

xmlstr = sys.stdin.read()
doc = libxml2.parseMemory(xmlstr, len(xmlstr))
ctxt = doc.xpathNewContext()
for node in ctxt.xpathEval("/*[local-name() = \"cml\"]"):
	print dir(node)
	print node.name
