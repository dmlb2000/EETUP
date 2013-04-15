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
	# this one needs to be shifts not shielding then subtract the isotropic variant and negate the result.
	"isoshift":
	"//cml:propertyList[@ref = \"%s\" and @dictRef = \"nwchem:nmrChemicalShift\"]//cml:property[@dictRef = \"nwchem:isotropicChemicalShift\"]"%(ATOM),
	"shift":
	"//cml:propertyList[@ref = \"%s\" and @dictRef = \"nwchem:nmrChemicalShift\"]//cml:property[@dictRef = \"nwchem:ChemicalPrincipalComponents\"]"%(ATOM),
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

shifts = [ "%e"%(float(a)-float(values['isoshift'])) for a in values['shift'].split() ]
print "%s %s %s" % (shifts[0], shifts[1], shifts[2])
print values['EulerAlpha'], values['EulerBeta'], values['EulerGamma']
efg = [ "%e"%(float(a)*3.24136E15) for a in values['EFG'].split() ]
print "-0.02578e-24 %s %s %s" % (efg[0], efg[1], efg[2])
