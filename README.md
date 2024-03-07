# bamsta

Simple tool to calculate various per-base stats over a bam file, and returns bedgraph files. Currenly only collects stats for depth and mean MAPQ.

## Installation

Either download prebuilt binaries (Linux x86-64 and MacOS ARM64 only) from the release page: https://github.com/bjhall/bamsta/releases

...or compile from main branch. Requires that the rust toolchain is installed.

```
git clone https://github.com/bjhall/bamsta.git
cargo build --release
```

## Usage

```
bamsta input.bam output_prefix
```

## Output

`PREFIX.dp.bedgraph` - Sequencing depth per base

`PREFIX.mapq.bedgraph` - Mean MAPQ per base
