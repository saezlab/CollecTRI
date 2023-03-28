## Environments
We provide the dependencies required to run the following steps:
- CollecTRI GRN generation
- Benchmarking of different GRNs
- Case study using CollecTRI and decoupleR

These dependencies are provided in the form of a `.yaml` environment file which can be deployed to your system automatically. We refer to the [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) for detailed instructions on how to do this.

----

### Binaries for FIMO
The motif enrichment analysis using FIMO, requires a working installation of the [meme-suite](https://meme-suite.org/meme/).
Please find the instructions on how to install it on the [bioconductor page of the memes packages](https://bioconductor.org/packages/release/bioc/vignettes/memes/inst/doc/install_guide.html)

Once the MEME Suite is installed, you can check the installation path in R:
```R
library(memes)

check_meme_install('/path/to/bin') #change to correct path
```

Before running the following script (`04.2_weights_FIMO.R`) , please point to the correct binaries, by changing the following code snippet in the header of the file:
```R
options(meme_bin = "/opt/local/bin/") # change to correct path on your system
```

