## Example

After installing `gdock` simply:

```
$ cd examples/
$ pyton ../gdock.py run.toml
```

## How-to-use

All run-specific parameters are defined in the `run.toml` file and all the fields provided in the example are needed.

```toml
[main]
identifier = 'dev'
number_of_processors = 4

[restraints]
A = [69, 242, 259, 260, 292, 294, 295, 321, 323, 324]
B = [76, 71, 67, 39, 41, 42, 44, 69]

[molecules]
A = '1a2k_r_u.pdb'
B = '1a2k_l_u.pdb'
native = '1a2k_complex_bound.pdb'
```

* * *

:octopus: