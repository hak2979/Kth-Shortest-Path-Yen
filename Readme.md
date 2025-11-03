# Parallel Kth-Shortest Path Finder using MPI and OpenMP

This project implements a parallelized version of **Yen’s Algorithm** for finding the Kth-shortest paths in a graph.  
It leverages **MPI (Message Passing Interface)** for distributed parallelism and **OpenMP** for shared-memory parallelism, enabling efficient computation across multiple processors and machines.

---

## 🧩 Features

- Implements **Yen’s Algorithm** for Kth-shortest paths.
- Utilizes **MPI** for inter-process communication.
- Supports **OpenMP** for multi-threaded execution.
- Can be executed on a **single machine** or a **cluster of multiple devices**.
- Highly configurable — easily adjust number of processes and machines.

---

## ⚙️ Compilation

To compile the program, use the following command:

```bash
mpic++ w.cpp -o e -fopenmp
```
## ⚙️ Execution

```bash

mpiexec -n 10 ./e
```

If execute on multiple device:
```bash
mpiexec  -n 10 -f machinefile ./e
```


## Explaination:
Explanation of the compilation command:

- **`mpic++`** — The MPI C++ compiler used to compile programs that utilize MPI (Message Passing Interface).
- **`-fopenmp`** — Compiler flag that enables OpenMP support for multi-threaded parallelism.
- **`w.cpp`** — The source file containing the program’s C++ implementation.
- **`-o e`** — Specifies the name of the output executable file (`e`).


