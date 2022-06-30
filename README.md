# RapidFlow
## Introduction
Continuous subgraph matching (CSM) is an important building
block in many real-time graph processing applications. Given a
subgraph query 𝑄 and a data graph stream, a CSM algorithm reports
the occurrences of 𝑄 in the stream. Specifically, when a new
edge 𝑒 arrives in the stream, existing CSM algorithms start from the
inserted 𝑒 in the current data graph 𝐺 to search 𝑄. However, this
rigid matching order of always starting from 𝑒 can lead to a massive
number of partial results that will turn out futile. Also, if 𝑄 contains
automorphisms, there will be a lot of redundant computation in
the matching process. To address these two problems, we propose
RapidFlow, an effective approach to CSM. First, we design a query
reduction technique, which reduces CSM to batch subgraph matching
(BSM) where we enumerate all results in a region of 𝐺 that will
be affected by the update. The well-established BSM techniques can
determine effective matching orders, not necessarily starting from
the newly inserted edge. Second, to eliminate redundant computation
caused by automorphisms in 𝑄, we propose dual matching,
which leverages the duality of 𝑄 and 𝐺 in the matching process.
Extensive experiment results show that RapidFlow outperforms
state-of-the-art algorithms, including TurboFlux and SymBi, by up
to two orders of magnitude on various workloads.

For the details, please refer to our VLDB'2022 paper
"RapidFlow: An Efficient Approach to Continuous Subgraph Matching"
by [Shixuan Sun](https://github.com/shixuansun), [Xibo Sun](https://github.com/xibosun), [Bingsheng He](https://www.comp.nus.edu.sg/~hebs/),
[Qiong Luo](http://www.cse.ust.hk/~luo/).
If you have any further questions, please feel free to contact us.


## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Test
Execute the following commands to test the correctness of the binary file. The script will
execute 100 test cases on insert and delete streams, respectively. The test data is stored
in the insert and delete folders. The test process will take around 15 minutes. As we have
tested the correctness of the binary, execute the command unless you modified the source code.

```zsh
cd test
python test.py ../build/streaming/RapidFlow.out
```

## Execute
After compiling the source code, you can find the binary file 'RapidFlow.out'
under the 'build/streaming' directory.  Execute the binary with the following
command './RapidFlow.out -d data_graphs -q query_graphs -u update_streams
-num number_of_embeddings -time_limit time_in_seconds',
in which '-d' specifies the input of the data graphs, '-q' specifies the
input of the query graphs and '-u' specifies the input of the graph update stream.
The '-num' parameter sets the maximum number of embeddings that you would like to find for each edge update.
If the number of embeddings enumerated reaches the limit or all incremental results for the update have been found,
then the program will process next update. Set '-num' as 'MAX' to find all incremental results for each update.
The '-time_limit' parameter configures the time budget for the query. If the query cannot be completed within the time limit,
then the program will terminate the query and return the number of results found. The default value is 3600 seconds (1 hour).
[This figure](/log-annotation.png) illustrates the output log of the binary.

Example: The time limit is 3600 seconds for each query. The target number of incremental matches is 429496729, which is the limit of uint32_t (Note that the target number is
the limit of uint64_t if set -num to MAX. In our experiments, we set -num to 429496729 to obtain the query time).

```zsh
./RapidFlow.out -d ../../test/insert/data_graph/data.graph -q ../../test/insert/query_graph/Q_0 -u ../../test/insert/data_graph/insertion.graph -num 429496729 -time_limit 3600
```

Example: The time limit is 3600 seconds for each query. The target number of incremental matches is 1 (In our experiments,
we set -num to 1 to obtain the response time).

```zsh
./RapidFlow.out -d ../../test/insert/data_graph/data.graph -q ../../test/insert/query_graph/Q_0 -u ../../test/insert/data_graph/insertion.graph -num 1 -time_limit 3600
```

The query time of processing a stream includes the global indexing time, which is the elapsed time of updating the global index,
the local indexing time, which is the elapsed time of updating the local index, and the enumeration time, which is the elapsed
time of enumerating the incremental results. If you want to measure the three metrics separately, add the macro
'#define MEASURE_INDEXING_COST' to Line 18 in streaming_engine.cpp.

## Experiment Datasets and Baseline Methods in our Experiments

The datasets and baseline methods in our experiments can be found in this [repository](https://github.com/RapidsAtHKUST/ContinuousSubgraphMatching),
which is an in-depth study of existing continuous subgraph matching methods by our research team.
