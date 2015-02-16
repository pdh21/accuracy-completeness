# accuracy-completeness
Replicate the likelihood ratios, accuracy and completeness in arxiv:1312.0552 

THESE CODES WILL CHANGE INCLUDING THE OUTPUTS/RESULTS AS IT'S A WORK IN PROGRESS

To run:
Change the filenames in LRfiles.py
Run in this order and update the file LRfiles.py with the relevant outputted values:
findingsigmar.py
findingpds.py
findingqds.py
findingLRthresh.py
findingLR.py
--from here the codes haven't been tested on a working dataset so there might be issues.
findingacc.py
findingcompl.py
--codes for other uses 
radialflux2Dhist.py: plots fig 3 in arxiv:1312.0552, showing r and dS are independent.
comparepdsqds.py: plots qds/pds so you can see if there's any weird behaviour there 
kdTreecheck.py: checks that the number of sources per leaf in the kdTree isn't affecting the output of the results (everything's fine but kept it in the upload as might as well)


The output of findingLR.py are the indices of the input catalogue that are matched and the associated dS and r's, acting as diagnostic tools. Not quite sure if everything works correctly yet (getting weird matches and negatives in qds). Haven't written 
