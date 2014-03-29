coupled-quasi-harmonic-bases
============================

Implementation of paper "Coupled Quasi-harmonic Bases"

<b>Usage:</b>

<b>setMu():</b> set the strength of coupling.

<b>setPackage():</b> choose between NLopt and Opt++.

<b>setOfftermType():</b> set the type of off-diagonal type, determine whether 
the output eigenbases are ordered or not.

<b>setInnerProductType():</b> choose the type of inner product between correspondence function and eigenbases.
If the input is point-wise correspondence, set inner product to standard.
If the input is region-based correspondence, set inner product to area-weighted.

<b>setObjectiveType():</b> choose between changing both eigenbases or only changing eigenbase a with respect to eigenbase b.

<b>enableScaleBase()/disableScaleBase():</b> optionally scale eigenbases of two shapes, the output is scaled back.

<b>enabldScaleCorrespondingFunction()/disableScaleCorrepondingFunction():</b>
If the input is point-wise correspondence, call disableScaleBase().
If the input is region-based correspondence, call enableScaleBase().

<b>setCouplingScale():</b> set different coupling strength for different eigenbases. Sometimes we want strong coupling for the first few eigenbases.

<b>enableOutput()/disableOutput():</b> enable/disable output of optimization.

<b>Dependency:</b>

NLopt: http://ab-initio.mit.edu/wiki/index.php/NLopt

Opt++: https://software.sandia.gov/opt++/
