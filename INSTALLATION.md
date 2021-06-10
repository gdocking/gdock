### gdock

Download the latest relase or clone a stable branch of the repository.

### Python environment

An `environment.yml` file is provided and Anaconda is recomended.

```bash
$ conda env create -f environment.yml
$ conda activate gdock
(gdock) $
```

### Third-party Dependencies

`gdock` uses [DComplex](https://sparks-lab.org/Publications_files/zhou061.pdf) as the scoring function, [FCC](https://github.com/haddocking/FCC) as the clustering engine and [PROFIT](http://www.bioinf.org.uk/software/profit) to calculate CAPRI metrics.


**After cloning the repository and setting up the python environment with anaconda**, you can use the very rough install script provided to install the third-party dependencies and configure `gdock.ini`:

Example:
```bash
$ bash install.sh /home/rodrigo/software/gdock
```

* * *

:octopus: