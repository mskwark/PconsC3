# PconsC3
Faster, more accurate and entirely open source method for predicting contacts in proteins

If you use PconsC3 please cite:
 *  Carlo Baldassi, Marco Zamparo, Christoph Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weigt and Andrea Pagnani, (2014) PLoS ONE 9(3): e92721. doi:10.1371/journal.pone.0092721
 *  Christoph Feinauer, Marcin J. Skwark, Andrea Pagnani, and Erik Aurell. (2014) PLoS Comp Bio: e1003847. doi:10.1371/journal.pcbi.1003847

# Prerequisites

* Julia interpreter (ver. 0.3 and up is supported). Present in most Linux repositories (Ubuntu , otherwise download it from [Julia](http://julialang.org/) website.
* Python interpreter (2.7+)
* CD-HIT. Available in most Linux distributions, otherwise downloadable from [GitHub](https://github.com/weizhongli/cdhit)
* A way to generate multiple sequence alignments (or a FASTA formatted MSA).
* A way to generate PSIPRED-like secondary structure predictions
* A way to generate NetSurfP-like solvent accessibility predictions
* An external source of contact information (e.g. PhyCMAP, CMAPpro...), capable of producing contact estimates in CASP RR format

If Julia, Python and CD-HIT are in your search path, you are set to go. Otherwise, you need to either add them to the path, or modify the necessary scripts.

For Julia: `rungdca.py` and `runplm.py`
For CD-HIT: `alignmentstats.py`

# Installation

1. Check out PconsC3 from GitHub
    ```
    git checkout https://github.com/mskwark/PconsC3.git
    ````
2. Install Julia packages. Start Julia and install:
    * [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl).
    ```
    julia> Pkg.add("NLopt")
    ```
    * [GaussDCA](https://github.com/carlobaldassi/GaussDCA.jl)
    ```
    julia> Pkg.clone("https://github.com/carlobaldassi/GaussDCA.jl")
    julia> Pkg.clone("https://github.com/carlobaldassi/ArgParse.jl")
    ```
    * [PlmDCA](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003847). Install a (very slightly) modified version of [PlmDCA.jl](https://github.com/pagnani/PlmDCA)
    ```
    julia> Pkg.clone("https://github.com/mskwark/PlmDCA")
    ```
3. Make sure the prerequisites are installed.

# Running the software

Before the first run, download and unpack the trained Random Forests in the same directory as the PconsC3 code. You should have six subdirectories named `tforest0, tforest1,...tforest5`. You can get them from [Google Drive](https://drive.google.com/open?id=0BxpeugdrylmAV2pwVXpPcW5JR3c) or a [local mirror](https://share.ics.aalto.fi/project/pconsc2/pconsc3.forests.all.tar.xz) (378MiB). If you just want to give PconsC3 a try or for some other reason need a smaller archive, feel free to download the mini-version either from [Google Drive](https://drive.google.com/open?id=0BxpeugdrylmAS1UzNG1oemh3Y3c) or [local mirror](https://share.ics.aalto.fi/project/pconsc2/minitrees.tar.xz) (38MB), being advised that this version may not perform as well as the fully-fledged one (but will be roughly 10x faster!).

```
> tar -xJf pconsc3-forests.tar.xz
```

You may want to put them on a fast filesystem (on a relatively recent Linux machine `/dev/shm/` is a good choice and by default PconsC3 will look for them there (i.e. it will check if `/dev/shm/tforest0` etc. exist and are sane). As a fallback it will look in the same directory `./predict.py` is located. If you want to change it, you need to modify `forestlocation` variable in the head of `./predict.py`. 

To run PconsC3, you need to have at hand:
 * Your input alignment in FASTA format, retaining only these columns that you want to run the prediction for (most often: all the match states/amino acids in your target sequence). For A3M alignments this can be attained by filtering all the lowercase letters (inserts) from the sequences.
 * Predicted secondary structure in a PSIPRED ss2 format. Download or run [PSIPRED here](http://bioinf.cs.ucl.ac.uk/psipred/)
 * Predicted relative solvent accessibility (RSA) in NetSurfP format. Download or run [NetSurfP here](http://www.cbs.dtu.dk/services/NetSurfP/)
 * A source of external estimates of contact propensities in a [CASP RR format](http://predictioncenter.org/casp8/index.cgi?page=format#RR) to act as a contact prior. We have tested the method with CMAPpro and PhyCMAP, but other methods should work as well.

You can name these files any way you want, but assuming your alignment is named `myprotein.fas`, your contact priors are named `external.RR`, secondary structure prediction file is named `psipred.ss2` and RSA is named `netsurf.rsa`, to run the prediction do the following. 

 1. Infer evolutionary couplings with GaussDCA:
    ```
    ./rungdca.py myprotein.fas
    ```
    It will produce a file named `myprotein.gdca`

 2. Infer evolutionary couplings with plmDCA.jl:
    ```
    ./runplm.py myprotein.fas
    ```
    It will produce a file named `myprotein.0.02.plm20`

 3. Compute alignment statistics:
    ```
    ./alignmentstats.py myprotein.fas
    ```
    It will produce a file named `myprotein.stats`

 4. Run PconsC3:
    ```
    ./predict.py myprotein.gdca myprotein.0.02.plm20 external.RR netsurf.rsa psipred.ss2 myprotein.stats myprotein.fas outputfile
    ```
    This will run for a while, but will provide you with estimates of running time. It will result in a number of intermediate files being generated: `outputfile.l0, outputfile.l1...outputfile.l5` and an `outputfile.RR` containing final predictions in RR format (by default only non-local prediction are output).

# Making PconsC3 run faster

There are a few parameters in `./predict.py` that can be tweaked, notably:
 * `maxtime` -- the maximum time in seconds spent on a single prediction layer (total prediction time will be at most 6x maxtime + time spent on i/o). For the larger proteins setting maxtime too low may result in sub-par performance
 * `treefraction` -- the fraction of trees to be used. While we recommend leaving `treefraction` set to 1., benchmarks have demonstrated satisfactory performance at values as low as `0.3` and for proteins with a lot of sequence information, as low as `0.1`. 

# Help and Support

If you run into any problems with the software or observe it performing poorer than expected, we would appreciate an email to Marcin J. Skwark (firstname@lastname.pl or firstname.middleinitial.lastname@vanderbilt.edu).
