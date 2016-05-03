# fcsp
An encoder for fragmental code of substructures superposition  (FCSS). FCSP is historical name obtained by transilteration

## Building

fcsp is build with [SCons](http://scons.org) tool. There are no extra library dependencies except for [Boost C++](http://www.boost.org).

On Ubuntu therefore it should be enough to install the following:

```sh
sudo apt-get install scons libboost-dev build-essential
```

And then invoke scons in the folder of the project:
```sh
scons
```



## Command-line options

Some commonly useful options are:
`-d|--descriptors` - directory that contains DB of chemical patterns i.e. files descr1.csv, descr2.sdf and replacement.sdf
`--format` - output format, currently supported 'txt' - plain text, 'csv' - pairs of file name + text of FCSS codes, and the most complete 'json' format that also includes location of each decriptor in the molecule.

These options are followed by a list of MOL files to process, the result is outputtted to stdout in the format specified by `--format` flag. Alternatively is no MOL files are given, reads single MOL file from stdin.


