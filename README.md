# Light Field Codec

## Clone and build
To build this repository, you should have cmake 3.15+ installed in your system.

First, clone the repo.

```
$ mkdir repo
$ git clone https://github.com/IDSRosler/Light-Field-Codec.git repo
$ cd repo
```

Then, make a directory to build the project binaries

```
$ mkdir build
$ cd build
$ cmake ..
$ make 
```

From the `build` directory, run the binary `./LF_Codec`.


## Encoder Parameters (CLI)
| Parameters                  | Type      | Value | Description |
|-----------------------------|:---------:|:------:|-------------|
| -input                      | String    | Required | Path to folder containing the *.ppm files to be encoded |
| -output                     | String    | Required | Path to folder where the reconstructed *.ppm files will be stored alongside with any execution file|
| -blx                        | Int       | Required | Block size in x dimension
| -bly                        | Int       | Required | Block size in y dimension
| -blu                        | Int       | Required | Block size in u dimension
| -blv                        | Int       | Required | Block size in v dimension
| -lfx                        | Int       | Required | Light field size in x dimension
| -lfy                        | Int       | Required | Light field size in y dimension
| -lfu                        | Int       | Required | Light field size in u dimension
| -lfv                        | Int       | Required | Light field size in v dimension 
| -qp                         | Float     | Default: 1 | Linear quantization parameter for the quantization 4D-volumn
| -qx                         | Int       | Default: 1 | Quantization weight in x dimension
| -qy                         | Int       | Default: 1 | Quantization weight in y dimension
| -qu                         | Int       | Default: 1 | Quantization weight in u dimension
| -qv                         | Int       | Default: 1 | Quantization weight in v dimension
| -prediction                 | String    | Default: `angular` | Prediciton mode. Valid modes are `none`, `angular`, `all`, `DC`, `IBC(Intra Block Copy)`
| -lambda                     | Float     | Default: 100 | Lagrangian multiplicative for the RD-COST 
| -use-transforms             | List      | Default: `DCT_II` | Defines which transforms will be used during the encoding step. A space separated list is expected. Valid values are `DCT_II`, `DST_I`, `DST_VII`
| -transform-min-angular-size | Int       | Default: 4| Minumum size of a block on the dimensions u,v
| -transform-min-spatial-size | Int       | Default: 4| Minumum size of a block on the dimensions x,y
| -partition-tree-max-depth   | Int       | Default: 2| Maximum depth for the partition tree. 

| Flag                  | Description |
|-----------------------|-------------|
| -lytro                | If set, the Light Field will be treated as if captured by a Lytro camera
| -verbose              | If set, shows a verbose output for each step of the encoding process
| -experimental         | If set, enables experimental features
| -show-progress-bar    | If set, displays a progress bar with the estimated remaining time
| -disable-transforms   | If set, disables the transform/quantization step
| -lossless             | If set, no quantization is applyied
| -uniform-quantization | If set, quantize all coefficients with the value passed to `-qp` parameter and ignore others




### Running the binary
```bash
# Please note that both paths must end with a trailing slash!
DATASET_DIR="/path/to/dataset/ppm/"
RESULT_DIR="/path/to/result/"
LFCODEC_BIN="./build/LF_Codec"

# Result folder must exist before calling the encoder!
mkdir -p "${RESULT_DIR}"


$LFCODEC_BIN -input "${DATASET_DIR}"              \
             -output "${RESULT_DIR}"              \
             -lfx 625 -lfy 434 -lfu 13 -lfv 13    \
             -blx 15  -bly 15  -blu 13 -blv 13    \
             -qx 3    -qy 3    -qu 3   -qv 3      \
             -qp 3                                \
             -prediction angular                  \
             -lambda 1                            \
             -use-transforms DCT_II DST_I DST_VII \
             -partition-tree-max-depth 2          \
             -lytro                               \
             -show-progress-bar                   \
             -experimental

```





<!-- Maybe remove this part?  -->
### Scripts 
- Scripts to run the Light Field Codec for different **QPs** and different dataset
##### Run simulation.py:
1. In the simulation.py, change the input_file variable and the output_file variable
    - **input_file:** Directory with .ppm files to be encoded
    - **output_file:** Output directory for the encoded *.ppm
2. Execute
    - taskset -c 0 python3 simulation.py


