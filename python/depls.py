#!/usr/bin/env python

def depsort(deps):
  ''' returns a sorted list s.t. dependencies specified in 'deps' are satisfied
  INPUT:  deps[i]:    a list items that i depends on
  OUTPUT: (d, cycl)
          d[0..n-1]:  a permutation of n = len(deps)
          such that d[i] only depends on d[0], ..., d[i-1]
          return None if such an arrangement is impossible
          cycl: a chain of cyclic dependencies '''
  d = sort(deps)
  if d: return d, None
  else: return None, findcycl(deps)


def sort(deps):
  ''' returns a sorted list s.t. dependencies specified in 'deps' are satisfied
  INPUT:  deps[i]:    a list items that i depends on
  OUTPUT: d[0..n-1]:  a permutation of n = len(deps)
          such that d[i] only depends on d[0], ..., d[i-1]
          return None if such an arrangement is impossible 
  NOTE: the algorithm works best if the input array already sorted '''
  n = len(deps)
  d = range(n)
  for i in range(n): # figure out d[i] every time
    for j in range(i, n):
      ls = deps[ d[j] ] # list of d[j]'s dependencies
      if len(ls) == 0: break # an independent item
      pos = map(d.index, ls) # map to dependencies' positions
      if max(pos) < i: break # only depends on d[0], ..., d[i-1]
    else: return None
    if j != i: d[i], d[j] = d[j], d[i] # swap d[i] and d[j]
  return d

def findcycl(deps):
  ''' find cyclic dependencies
  INPUT:  deps[i]:    a list items that i depends on
  OUTPUT: 0 if there's no cyclic dependencies
          otherwise a chain of dependencies '''
  n = len(deps)
  good = [0] * n # good[i] = 1 if i has no dependency problem
  while 1:
    # find a dirty item with at least one dependency
    for i in range(n):
      if len(deps[i]) and not good[i]: break
    else: return 0 # all nodes are indep.

    # grow a dependency tree starting from i
    st = [0] * (n + 1)  # st[i] is the item
    jt = [0] * (n + 1)  # jt[i] is the current index during looping over deps[ st[i] ]
    top = 0
    st[top] = i # push i to the root of tree
    jt[top] = 0 # index of child nodes of st[top]
    while top >= 0:
      i = st[top]
      j = jt[top]
      dl = deps[i]
      # about to examine dl[j]
      if j >= len(dl): # examined all dependencies of i
        good[i] = 1 # mark i as a good item
        top -= 1
        if top >= 0: jt[top] += 1
      elif good[ dl[j] ]:
        jt[top] += 1 # skip it
      else: # search dependencies of dl[j]
        top += 1
        x = dl[j]
        bad = st[:top] # a list of forbidden items
        if x in bad: # return a list of dependencies
          return bad[bad.index(x):] + [x]
        st[top] = x # proceed to the next level
        jt[top] = 0
  return 0

def main():
  ''' test depsort() '''
  deps = [ 
      [],     # 0
      [2, 4], # 1
      [3, 4], # 2
      [0],     # 3
      [0],    # 4
        ]
  d = depsort(deps)
  if d:
    print d
  else:
    print '-->'.join( findcycl(deps) )

if __name__ == "__main__":
  main()
