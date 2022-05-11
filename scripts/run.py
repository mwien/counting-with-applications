from dct_policy import dct_policy
from causaldag import DAG
import sys

n, m = map(int, sys.stdin.readline().split())
sys.stdin.readline()
arcs = []
for i in range(m):
    a, b = map(int, sys.stdin.readline().split())
    arcs.append((a,b))


g = DAG(arcs = arcs)
intervened_nodes = dct_policy(g)
print(len(intervened_nodes))
print()
for v in intervened_nodes:
    print(v)
