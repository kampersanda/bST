# integer\_sketch\_search

A C++17 implementation of data structures described in the paper: Shunsuke Kanda and Yasuo Tabei, "[b-Bit Sketch Trie: Scalable Similarity Search on Integer Sketches](https://arxiv.org/abs/1910.08278)," IEEE BigData 2019. This library supports to benchmark the single- and multi-index approaches through the b-bit sketch trie or the hash table.

## Build instructions

You can download and compile this library with the following commands:

```
$ git clone https://github.com/kampersanda/integer_sketch_search.git
$ cd integer_sketch_search
$ mkdir build && cd build
$ cmake ..
$ make -j
```

After the commands, the executables will be produced in `build/bin` directory.

### Requirements

- C++17 supported compiler such as g++ >= 7.0 or clang >= 4.0.
- CMake >= 3.0
- [sdsl-lite](https://github.com/simongog/sdsl-lite) (should be installed at your home directory when the default `CMake` setting is used)

## Example to run

You can use executable `bin/search` to benchmark the data structures.

```
$ ./bin/search 
usage: ./bin/search --name=string --index_fn=string --base_fn=string --query_fn=string [options] ... 
options:
  -n, --name          index name (hash | trie) (string)
  -i, --index_fn      input/output file name of index (string)
  -d, --base_fn       input file name of database sketches (string)
  -q, --query_fn      input file name of query sketches (string)
  -m, --dim           dimension (<= 64) (int [=32])
  -b, --bits          #bits of alphabet (<= 8) (int [=2])
  -B, --blocks        #blocks (B=1 means to use single index) (int [=1])
  -e, --errs_range    range of errs (min:max:step) (string [=1:5:1])
  -v, --validation    validation (bool [=0])
  -s, --suf_thr       suf_thr (float [=2])
  -?, --help          print this message
```

You can try it using toy datasets `data/news20.scale_{base|query}.cws.bvecs`. (The input data format will be explained later.)

### 1) Testing single-trie index

The following command performs to construct a single-trie index (indicated by options `-n trie` and `-B 1`) from database `data/news20.scale_base.cws.bvecs` (indicated by `-d`) and search for queries `data/news20.scale_query.cws.bvecs` (indicated by `-q`), where the testing sketch dimension is `16` (indicated by `-m`), the testing number bits for integers is `2` (indicated by `-b`), and the testing error thresholds are `[1,2,3]` (indicated by `-e`). The parameter `\lambda` for sparse layer (defined in the paper) can be indicated by `-s`.

```
$ ./bin/search -n trie -i news20 -d ../data/news20.scale_base.cws.bvecs -q ../data/news20.scale_query.cws.bvecs -m 16 -b 2 -B 1 -e 1:3:1 -v 0 -s 2
### sketch_trie ###
Now loading keys...
--> 15935 keys
--> 0 sec
Now constructing index
--> 0 sec
Now writing news20.16m2b1B.trie
--> 71816 bytes; 0.0684891 MiB
Statistics of sketch_trie
--> perf_height: 5
--> suff_dim: 10
--> suf_thr: 2
--> rep_type: HYBRID
Now loading queries...
--> 100 queries
Now simlarity searching...
--> 1 errs; 0.1 ans; 0 cands; 0 ms
--> 2 errs; 0.13 ans; 0 cands; 0 ms
--> 3 errs; 0.25 ans; 0 cands; 0.03 ms
```

For each threshold, the average number of answers, the average number of answer candidates (for multi-index approaches), and average search time (in ms) are reported. After this, the index file `news20.16m2b1B.trie` will be written whose prefix is indicated by `-i`. When the same parameters are tested again, the index file will be read.

### 2) Verifying the correctness

When option `-v 1` is set, you can verify the correctness of answers by using the middle value of error thresholds.

```
$ ./bin/search -n trie -i news20 -d ../data/news20.scale_base.cws.bvecs -q ../data/news20.scale_query.cws.bvecs -m 16 -b 2 -B 1 -e 1:3:1 -v 1 -s 2
### sketch_trie ###
Now loading keys...
--> 15935 keys
--> 0 sec
Now loading index
--> 71816 bytes; 0.0684891 MiB
Statistics of sketch_trie
--> perf_height: 5
--> suff_dim: 10
--> suf_thr: 2
--> rep_type: HYBRID
Now loading queries...
--> 100 queries
Now validating with 2 errs...
--> No problem!!
```

### 3) Testing multi-trie index

By setting `-B` to a value no less than 2, the multi-trie index can be tested as follows.

```
$ ./bin/search -n trie -i news20 -d ../data/news20.scale_base.cws.bvecs -q ../data/news20.scale_query.cws.bvecs -m 32 -b 4 -B 4 -e 1:5:2 -v 0 -s 2
### multi_index<sketch_trie> ###
Now loading keys...
--> 15935 keys
--> 0 sec
Now constructing index
--> 0 sec
Now writing news20.32m4b4B.trie
--> 541789 bytes; 0.51669 MiB
Statistics of sketch_trie
--> perf_height: 2
--> suff_dim: 5
--> suf_thr: 2
--> rep_type: HYBRID
Statistics of sketch_trie
--> perf_height: 2
--> suff_dim: 5
--> suf_thr: 2
--> rep_type: HYBRID
Statistics of sketch_trie
--> perf_height: 2
--> suff_dim: 5
--> suf_thr: 2
--> rep_type: HYBRID
Statistics of sketch_trie
--> perf_height: 2
--> suff_dim: 5
--> suf_thr: 2
--> rep_type: HYBRID
Now loading queries...
--> 100 queries
Now simlarity searching...
--> 1 errs; 0.04 ans; 0.08 cands; 0 ms
--> 3 errs; 0.05 ans; 0.12 cands; 0 ms
--> 5 errs; 0.11 ans; 0.14 cands; 0 ms
```

### 4) Testing multi-hash index

By setting `-n hash`, the multi-hash index can be tested as follows.

```
$ ./bin/search -n hash -i news20 -d ../data/news20.scale_base.cws.bvecs -q ../data/news20.scale_query.cws.bvecs -m 32 -b 4 -B 4 -e 1:5:2 -v 0 -s 2
### multi_index<hash_table> ###
Now loading keys...
--> 15935 keys
--> 0 sec
Now constructing index
--> 0 sec
Now writing news20.32m4b4B.hash
--> 1717553 bytes; 1.63799 MiB
Now loading queries...
--> 100 queries
Now simlarity searching...
--> 1 errs; 0.04 ans; 0.08 cands; 0 ms
--> 3 errs; 0.05 ans; 0.12 cands; 0 ms
--> 5 errs; 0.11 ans; 0.14 cands; 0.01 ms
```


### Input data format

Integer sketches should be stored in [TEXMEX's bvecs format](http://corpus-texmex.irisa.fr/). That is, the dimension number and feature values for each sketch are interleaved (in little endian), where the number is 4 bytes of size and each feature is 1 byte of size. To convert vectors of real numbers into sketches, you can use [consistent\_weighted\_sampling](https://github.com/kampersanda/consistent_weighted_sampling) to generate such a dataset using the [GCWS](https://doi.org/10.1145/3097983.3098081) algorithm.

## Licensing

This program is available for only academic use, basically. For the academic use, please keep MIT License. For the commercial use, please make a contact to Shunsuke Kanda or Yasuo Tabei.

If you use the library, please cite the following paper:

```latex
@inproceedings{kanda2019bst,
  author = {Kanda, Shunsuke and Tabei, Yasuo},
  title = {b-Bit Sketch Trie: Scalable Similarity Search on Integer Sketches},
  booktitle = {Proceedings of the 2019 IEEE International Conference on Big Data},
  pages={810--819},
  year = {2019}
}
```
