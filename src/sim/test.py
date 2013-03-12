#!/usr/bin/python

import libxml2
import sys

xmlstr = sys.stdin.read()
doc = libxml2.parseMemory(xmlstr, len(xmlstr))
ctxt = doc.xpathNewContext()
ctxt.xpathRegisterNs('cml', "http://www.xml-cml.org/schema")
ATOM="O1.2"
xpath = {
	"EFG":
	"//cml:propertyList[@ref = \"%s\"]//cml:property[@dictRef = \"nwchem:efgPrincipalComponents\"]"%(ATOM),
	"shielding":
	"//cml:propertyList[@ref = \"%s\"]//cml:property[@dictRef = \"nwchem:shieldingPrincipalComponents\"]"%(ATOM),
	"EulerAlpha":
	"//cml:propertyList[@ref = \"%s\"]//cml:property[@dictRef = \"nwchem:eulerAngleAlpha\"]"%(ATOM),
	"EulerBeta":
	"//cml:propertyList[@ref = \"%s\"]//cml:property[@dictRef = \"nwchem:eulerAngleBeta\"]"%(ATOM),
	"EulerGamma":
	"//cml:propertyList[@ref = \"%s\"]//cml:property[@dictRef = \"nwchem:eulerAngleGamma\"]"%(ATOM),
}
values = {}
for k in xpath:
	x = xpath[k]
	for node in ctxt.xpathEval(x):
		values[k] = node.getContent().strip()
		print x, k
		print node.get_properties()
		print node.getContent()

print values['shielding']
print values['EulerAlpha'], values['EulerBeta'], values['EulerGamma']
print "-0.02578e-24", [ "%e"%(float(a)*3.24136E15) for a in values['EFG'].split() ]
