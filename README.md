# CAGEcleaner

[![DOI](https://zenodo.org/badge/904110273.svg)](https://doi.org/10.5281/zenodo.14726119)

**>>> `CAGEcleaner` will be integrated into [`cblaster`](https://github.com/gamcil/cblaster)! <<<**

## Outline

`CAGEcleaner` removes genomic redundancy from gene cluster hit sets identified by [`cblaster`](https://github.com/gamcil/cblaster). The redundancy in target databases used by `cblaster` often propagates into the result set, requiring extensive manual curation before downstream analyses and visualisation can be carried out.

Given a session file from a `cblaster` run (or from a [`CAGECAT`](https://cagecat.bioinformatics.nl/) run), `CAGEcleaner` retrieves all hit-associated genome assemblies, groups these into assembly clusters by ANI and identifies a representative assembly for each assembly cluster using `skDER`. In addition, `CAGEcleaner` can reinclude hits that are different at the gene cluster level despite the genomic redundancy, and this by different gene cluster content and/or by outlier `cblaster` scores. Finally, `CAGEcleaner` returns a filtered `cblaster` session file as well as a list of retained gene cluster IDs for easier downstream analysis.

For installation instructions, usage, explanations and more, head over to the [`CAGEcleaner` wiki](https://github.com/LucoDevro/CAGEcleaner/wiki)!

![workflow](workflow.png)

## Citations
If you found `CAGEcleaner` useful, please cite our manuscript:

```
De Vrieze, L., Biltjes, M., Lukashevich, S., Tsurumi, K., Masschelein, J. (2025) CAGEcleaner: reducing genomic redundancy in gene cluster mining. bioRxiv https://doi.org/10.1101/2025.02.19.639057
```

`CAGEcleaner` relies heavily on the `skDER` genome dereplication tool and its main dependency `skani`, so please give these proper credit as well.

```
Salamzade, R., & Kalan, L. R. (2023). skDER: microbial genome dereplication approaches for comparative and metagenomic applications. https://doi.org/10.1101/2023.09.27.559801`
Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nature Methods, 20(11), 1661–1665. https://doi.org/10.1038/s41592-023-02018-3
```

## License

`CAGEcleaner` is freely available under an MIT license.

Use of the third-party software, libraries or code referred to in the References section above may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.
