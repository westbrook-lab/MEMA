# MEMA

The MEMA package includes functions to preprocess, QA, normalize and explore Microenvironment Microarray (MEMA) experiments. 
MEMAs are spotted cellular microarrays that are imaged at the population or cell level. Each spot can have different microenvironment proteins. A MEMA typically has 700-4000 spots. The arrays are in wells that expose the spots to conditions such as cell type, media and drug. Wells can be in 96 well plates which typically hold 64 spots
or in 8 well plates that can hold an array of ~700 spots.

The cells are stained with antibodies that typically target DNA and specific proteins. The spots are imaged at the population or cell level using quantitative immunofluroescence.

The MEMA package functions merge the experiment metadata with the imaging data, perform a bivariate loess-based QA and implement several normalization methods.

-   The MEMA package can be downloaded and installed from this repo with the command:

    ``` r
    devtools::install_github("markdane/MEMA")

    ```
