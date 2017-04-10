Welcome to Corretto- a denovo error correction program! Now that the excitement has subsided- here are a couple of things you'll need to know to understand this project within Adam. 
What you'll need to be able to run Corretto:
* Maven
* Spark
* Adam
* Corretto

Now that you're ready to use this amazing piece of code- here's the basic input for Corretto:

```
./bin/corretto-submit [INPUT_FILE] [OUTPUT_FILE]
```

The input file is a file of reads in any format that is supported by ADAM.
The output file will be a parquet file containing the corrected reads. The parquet file can also be outputted in any format that Adam supports. For more information on ADAM input/outputs, please see: https://github.com/meghanaravikumar/adam/blob/master/README.md 

Now you might be wondering- my that's aamazing but what does corretto actually do?
Given a set of reads, corretto goes in and detects bases that are faulty, corrects them, and returns the new and improved set of reads which can now be aligned with a higher accuracy than alignment perfromed without this type of error base fixing. 
The process for identifying the bases that are faulty and accurately fixing these bases is as follwos:
The input file of reads is split into q-mers, which is then split into k-mers. The k-mer distribution forms a Poisson Mixture 
distribution with error k-mers forming the left-tail of the distribution and the good k-mers forming the crux of the right distribution. Essentially, the error k-mers form a Poisson distribution with a smaller expected value for the distribution than the good k-mers. 
A probability threshold is calculated to determine whether a given k-mer falls within the error or good Poisson distribution. 
The probability that a given k-mer is in the good distribution/is valid is:
        K = L(k is good)/L(k is bad) for all k in K
So modeling the correction becomes:
        P(b|K,k,m) = P(b_k|K,k,m)
        if b_k in k:
          P(b_k|K,k,m) = P(i|K) (probability k-mer was good under Poisson model)
        else:
          P(b_k|K,k,m) = m 
where b = k-mers covering a base
      K = k-mer probability vector
      k = set of k-mers observed
      m = probability that missing kmer is corrected
Softmax is used to check the likelihood probabilites for substitution of a base. Bases marked for substitution are determined by the quality score for each base. Coordinate descent is used to move the entire sequence to a state of larger likelihood estimate. 
When there are two or more errors present in the read/k-mer, the first error is corrected to find the next positions in the sequence that are deemed as errors. 
