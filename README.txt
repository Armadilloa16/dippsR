Author: Lyron Juan Winderbaum
Date  : 18/8/2015

R code I have written for the analysis
of MALDI imaging peaklist data is in 
localFunctions.R, I want to get around 
to wrapping it properly into a package
at some point.

tute.Rnw provides both:
 - a simple example of the use of knitr
     to produce a nicely formatted and
     reproducible LaTeX document, 
     including bibtex referencing.
 - a step-by-step explanation of how to 
     use most of the code in 
     localFunctions.R 

Note that chunk output is stored in 
./cache -- of particular note is the
intial reading in of the peaklists in 
chunk `read_data', in order to reproduce 
this one would need to unpack the 
`A1.7z' archive containing the original
peaklist data.

TODO:

Add a section on making ROC curves maybe? 
Meh.