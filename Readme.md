# Two-step graph mapper
## About

## Install
BWA mem and Minimap2 are requirements. Make sure they are installed first.
Then install by cloning:
```bash
git clone 
cd two_step_graph_mapper 
pip3 install .
```

## How to use 
This assumes you already have some reads mapped to a graph. If not, you can use rough graph aligner first.

Predict path:
```bash
two_step_graph_mapper predict_path -d graphdir/ -a some.graphalignments -c 1,2,3 -o my_predicted_path
```
NB: The above command also runs bwa index, which will take an hour or two extra. You can skip that step by adding `-s`.

Map reads to that path:
```bash
two_step_graph_mapper map_to_path -r my_predicted_path.fa -f my_reads.fa -o aligned_reads.sam
```
You now have a `aligned_reads.sam` which is a sam-file relative to your predicted path. 
We can convert the coordinates to linear reference genome coordinates by running:
```bash
two_step_graph_mapper convert_to_reference_positions -s aligned_reads.sam -d graph_dir/ -l my_predicted_path -c 1,2,3 -o converted.sam
```

## Benchmarking and comparison with vg
