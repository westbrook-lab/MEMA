# MEMA

The MEMA package includes functions to preprocess, QA, normalize and explore Microenvironment Microarray (MEMA) experiments. 
MEMAs are spotted cellular microarrays that are imaged at the population or cell level. Each spot can have different microenvironment proteins. A MEMA typically has 700-4000 spots. The arrays are in wells that expose the spots to conditions such as cell type, media and drug. Wells can be in 384 well plates which typically hold 25 spots
or in 8 well plates that can hold an array of ~700 spots.

The cells are stained with up to four antibodies that typically target DNA and specific proteins. The spots are imaged at the 
population or cell level using quantitative immunofluroescence.

The MEMA package functions merge the experiment metadata with the imaging data, perform a bivariate loess-based QA and support various normalization methods.
