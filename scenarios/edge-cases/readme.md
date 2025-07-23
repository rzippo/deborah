# Edge-cases

This folder documents examples of what `deborah` does *not* handle well.

We provide examples with the "wrong" version, that gives incorrect results, and the "fixed" version.

In most instances, the "wrong" version can be identified by the fact that LUDB (`--ludb`) produces a lower result than the lower bound (`--lb`).

> Reminder: you can run both studies with the same command, i.e. `deborah config.conf --ludb -lb`.
> If you do so, `deborah` will also compare them and highlight issues / note if the bound is tight.

Note: these edge-cases may be considered bugs to be fixed.
Pull requests are welcome.

## One flow per source-destination pair

`deborah` assumes that the configuration provides one flow per source-destination pair.

Thus, the following is not correct:

```
FLOW	1	4	579	215
FLOW	1	4	553	102
FLOW	1	4	557	269
```

The flows should instead be given as their aggregate:
```
FLOW	1	4	1689	586
```

This assumption is without loss of generality: 
in FIFO tandems, both worst-case delay and backlog of an aggregated flow are the same as those of each of its parts.

## Nodes should be aggregated when possible

This issue arises if there is a sequence of nodes $i, i+1, ... i+k$ such that all flows transiting through $i$ also transit through all the others.
This means that we can replace this sequence of nodes with a single one with $\beta_{aggr} = \beta_{i} \otimes \beta_{i+1} \otimes ... \otimes \beta_{i+k}$.

`deborah` assumes this processing is done before-hand.

For example, the following is not correct, because nodes 2 and 3 can be aggregated:

```
TANDEM 4 3

NODE 1  0       2054
NODE 2  0       10845
NODE 3  0       5171
NODE 4  0       1910

FLOW    1       4       595     400
FLOW    2       4       2235    1329
FLOW    1       3       85      945
```

It should instead be fixed like this:

```
TANDEM 3 3

NODE 1  0       2054
NODE 2  0       5171
NODE 3  0       1910

FLOW    1       3       595     400
FLOW    2       3       2235    1329
FLOW    1       2       85      945
```
