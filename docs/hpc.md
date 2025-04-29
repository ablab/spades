# Running SPAdes across multiple nodes of HPC cluser

Running software on high-performance computational clusters usually involves
many options including, for example, the grid engine used, queue, resource allocation
among others. It is notoriously hard to support all these options in a
generic way. Essentially each cluster requires its unique configuration and
settings that cannot be predicted reasonably.

The grid mode of SPAdes handles all these details via a so-called 'Executor'. By
default, we provide executors for SLURM-based clusters (expecting a more or less
fresh SLURM and OpenMPI installation under the hood) and direct run over bare
`mpiexec`. While we have done our best to make them as generic as possible, there
still might be many unhandled corner cases unique to a particular computational
cluster. Users are encouraged to change the default implementation of executors
and/or write their own suitable for a particular cluster and job submission
system. Please contact us if some guidance is required.

There are a few important notes that should be taken into account here:

* Assembly with SPAdes involves the whole computational pipeline with different
  steps that need to be run in a particular order, as they use the results from
  the previous steps.
* Not all steps of the SPAdes pipeline can be made parallel across multiple
  nodes of a cluster.
* Even those that can be made parallel incur some non-negligible overhead
  required to communicate the inputs, synchronize the global state, and combine
  the results.
* The mentioned overheads usually grow with the number of cluster nodes, so it
  might easily happen that running an assembly job across eight nodes would be
  slower than running across four nodes.
* All the nodes periodically need to see the common state. As a result, the
  memory consumption of a single node could be comparable to running SPAdes
  alone on a single node. Essentially, while the total assembly time could be
  reduced by spreading the workload across multiple nodes, the memory
  consumption normally becomes the same.
* The first node is the main one, combining all the data from worker
  nodes. Usually, it is the master node that consumes more memory, etc., and
  could benefit from, for example, `high-memory` nodes that are sometimes present
  on certain clusters. However, binding nodes in a non-uniform way across
  different available hardware resources usually requires some arcane settings
  for the job submission system, and this is fully left to the users.

## Compilation

To enable SPAdes HPC mode run

```
./spades_compile -SPADES_ENABLE_PROJECTS=hpcspades
```

As SPAdes HPC mode relies on MPI to split tasks across jobs sufficiently
sane MPI implementation is required. We only tested using OpenMPI. OpenMPI 3.x series
are known to contain important bugs resulting in deadlocks. OpenMPI 4.x and 
5.x worked fine in our tests.

After the compilation is complete, `spades-hpc` executable will be located in 
the `bin` folder.

## Command line options

The options specific to SPAdes support for computational clusters are:

`--grid-engine <engine>` The grid engine to use. `engine` might be one of:

* `local`: run everything locally, no grid engine / job submission is involved
* `mpi`: same as `local`, but `mpiexec` is used to run the applications that
  support MPI
* `slurm`: job submission via SLURM
* `save_yaml`: dry run, save all the steps of SPAdes pipeline one after another into
  `run_spades.yaml` file.

All executors also save `run_spades.sh` shell script. The execution of this
script is equal to running via `spades.py` (though some options are certainly
not supported here, like `--continue`). The shell script could be edited to tune
to a particular computational cluster and job submission system.

`--grid-queue <name>` Name of queue to submit jobs to.

`--grid-nnodes <nodes>` Number of nodes to use

`--grid-wait` By default SPAdes submits all jobs and exits. This option would
make it wait for completion of the results

The following addtition options are specific to SLURM executor:

`--grid-time <value>` Time limit for the job (corresponds to `--time` option of
`sbatch`)

`--grid-extra <options>` Any extra options to pass to job submission system

## Troubleshooting

`--only-generate-config` option might be useful to troubleshoot the job
submission and option. For SLURM executor it generates
`run_spades_on_cluster.sh` script with all `sbatch` / `srun` invocations used to
set up the execution of different steps of SPAdes pipeline via SLURM.


## References

TBA
