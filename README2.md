Welcome to Corretto! Now that the excitement has subsided- here are a couple of things you'll need to know to understand this project within Adam. 
What you'll need to be able to run Corretto:
1) Maven
2) Spark
3) Adam
4) Corretto

Now that you're ready to use this amazing piece of code- here's the basic input for Corretto:
./bin/corretto-submit [INPUT_FILE] [OUTPUT_FILE]

The input file is a file of reads in any format that is supported by ADAM.
The output file will be a parquet file containing the corrected reads. The parquet file can also be outputted in any format that Adam supports.

Now you might be wondering- my that's aamazing but what does corretto actually do?
Given a set of reads, corretto goes in and detects bases that are faulty, corrects them, and returns the new and improved set of reads which can now be aligned with a higher accuracy than alignment perfromed without this type of error base fixing. 
The process for identifying the bases that are faulty and accurately fixing these bases is as follwos:
The input file of reads is split into q-mers, which is then split into k-mers. The k-mer distribution forms a Poisson Mixture 
distribution with error k-mers forming the left-tail of the distribution and the good k-mers forming the crux of the right distribution. 
