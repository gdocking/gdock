### `gdock` Install instructions


```bash
$ git clone https://github.com/rvhonorato/gdock.git
$ cd gdock
$ python setup.py develop
$ bash install.sh `pwd`
$ gdock -h
```

To run the example, [follow this instructions](example/README.md)

* * *

### Third-party Dependencies

**These have DIFFERENT LICENSES, handle them accordingly.**

* [DComplex](https://sparks-lab.org/Publications_files/zhou061.pdf) as the scoring function,
* [FCC](https://github.com/haddocking/FCC) as the clustering engine,
* some script from [haddock-tools](https://github.com/haddocking/haddock-tools) to calculate the fitness and
* [ProFit](http://www.bioinf.org.uk/software/profit) to calculate CAPRI metrics.

The paths for these are defined in `/etc/gdock.ini`
```toml
[third_party]
profit_exe        = /Users/rodrigo/repos/gdock/src/ProFit_V3.3/src/profit
fcc_path          = /Users/rodrigo/repos/gdock/src/fcc
pdbtools_path     = /Users/rodrigo/repos/gdock/src/pdb-tools
haddocktools_path = /Users/rodrigo/repos/gdock/src/haddock-tools
; dcomplex this needs to be "specially compiled", change the fort21 path and the charges
dcomplex_exe = /Users/rodrigo/repos/gdock/src/dcomplex_single_file/dcomplex
```

* * *

:octopus:
