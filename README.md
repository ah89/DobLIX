# DobLIX: A Dual-Objective Learned Index for Log-Structured Merge Trees

This work implements DobLIX, a dual-objective learned index specifically designed for Log-Structured Merge (LSM) tree-based key-value stores. Traditional learned indexes primarily focus on optimizing index lookups, often overlooking the critical role of data access from storage, which can become a significant performance
bottleneck. 

In LSM-tree-based systems, a considerable portion of the index is stored on disk, making lookups highly dependent on the efficient coordination between in-memory structures and disk-resident data. Poorly optimized access patterns can lead to excessive I/O operations, negatively impacting read latency and overall system performance. DobLIX addresses this by incorporating a second objective, data access optimization, into the learned index training process. 

This dual-objective approach ensures that both index lookup efficiency and data access costs are minimized, leading to significant improvements in read performance while maintain- ing write efficiency in real-world LSM-tree systems. Additionally, DobLIX features a reinforcement learning agent that dynamically tunes the system parameters, allowing it to adapt to varying workloads in real-time. Experimental results using real-world datasets demonstrate that DobLIX reduces indexing overhead and improves throughput by 1.19× to 2.21× compared to state-of-the-art methods within RocksDB, a widely used LSM-tree-based storage engine.




## Build Instructions

### Prerequisites

Ensure you have the following dependencies installed:

```bash
sudo apt update
sudo apt install -y build-essential cmake libgflags-dev libsnappy-dev zlib1g-dev libbz2-dev liblz4-dev libzstd-dev
```

### Clone and Build DobLIX

```bash
# Clone the repository
git clone https://github.com/your-username/DobLIX.git
cd DobLIX

# Initialize and update RocksDB submodule (if applicable)
git submodule update --init --recursive

# Create a build directory and navigate to it
mkdir -p build && cd build

# Run CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Compile the project
make -j$(nproc)

# Optionally install
sudo make install
```

### Running DobLIX

After a successful build, you can run DobLIX using:

```bash
./doblix
```

### Cleaning the Build

To remove all compiled files and reset the build directory:

```bash
rm -rf build
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.
