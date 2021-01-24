### gdock

Download the latest relase or clone a stable branch of the repository.

### Python Dependencies

An `environment.yml` file is provided and Anaconda is recomended.

```bash
$ cd gdock
$ conda env create -f environment.yml
$ conda activate gdock
(gdock) $
```

### Third-party Dependencies

**gdock** uses [DComplex](https://sparks-lab.org/Publications_files/zhou061.pdf) as the scoring function, [FCC](https://github.com/haddocking/FCC) as the clustering engine and [DockQ](https://github.com/bjornwallner/DockQ) to calculate CAPRI metrics.

-   DComplex
    -   Download from [this link](http://servers.sparks-lab.org/downloads/dcomplex2.tar.gz)
    -   For it to work as intended in **gdock** you need to change line 115 and 148 of `dcomplex.c` and provide the full path to the DComplex installation,   then compile normally.
        ```c++
        // L28
        #define MAXA		54000  //total atom number in one protein.
        // L115
        read_charge_file("/Users/rodrigo/software/dcomplex_single_file/charge_inp.dat");
        // L148
        fp=(FILE *)fopen("/Users/rodrigo/software/dcomplex_single_file/fort.21_alla","r"); //monomic 1.61
        ```

-   FCC
    -   Clone the repository from [this link](https://github.com/haddocking/FCC)
    -   Follow the install instructions
    -   Set the default branch to `python3`
        ```bash
        ~/repos/fcc $ git checkout python3
        ```

-   DockQ
    -   Download it from [this link](https://github.com/bjornwallner/DockQ)
    -   Follow the install instructions

Edit `etc/gdock.ini` with the correct paths