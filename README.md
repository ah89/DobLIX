## DobLIX: A Dual-Objective Learned Index for Log-Structured Merge Trees

This work implements DobLIX, a dual-objective learned index specifically designed for Log-Structured Merge (LSM) tree-based key-value stores. Traditional learned indexes primarily focus on optimizing index lookups, often overlooking the critical role of data access from storage, which can become a significant performance
bottleneck. 

In LSM-tree-based systems, a considerable portion of the index is stored on disk, making lookups highly dependent on the efficient coordination between in-memory structures and disk-resident data. Poorly optimized access patterns can lead to excessive I/O operations, negatively impacting read latency and overall system performance. DobLIX addresses this by incorporating a second objective, data access optimization, into the learned index training process. 

This dual-objective approach ensures that both index lookup efficiency and data access costs are minimized, leading to significant improvements in read performance while maintain- ing write efficiency in real-world LSM-tree systems. Additionally, DobLIX features a reinforcement learning agent that dynamically tunes the system parameters, allowing it to adapt to varying workloads in real-time. Experimental results using real-world datasets demonstrate that DobLIX reduces indexing overhead and improves throughput by 1.19× to 2.21× compared to state-of-the-art methods within RocksDB, a widely used LSM-tree-based storage engine.

### How to run
Run  `cmake .` in main DobLIX path
