#!/usr/bin/python3

# Implemented by J.T. Brandenburg for using the script from Brandenburg et al. (2012)


import re
import sys

def parse(newick):
    tokens = re.findall(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(nextid = 0, parentid = -1): # one node
        thisid = nextid;
        children = []

        name, length, delim, ch = tokens.pop(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid+1, thisid)
                children.append(node)
            name, length, delim, ch = tokens.pop(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None,
                "parentid": parentid, "children": children}, delim, nextid

    return recurse()[0]

def print_newick(newick, Cmt=0, length=0):
  if len(newick['children'])==0:
      return newick['name']+":"+str(newick['length']+length)
  if len(newick['children'])==1:
      if newick['length']==None:
         lng=0
      else : 
          lng=newick['length']
      return print_newick(newick['children'][0],Cmt+1, lng)

  strout='('
  cmtchil=0
  if len(newick['children'])==1 and Cmt==1:
      print(newick)
      print(Cmt)
      sys.exit('children is 1')
  for children in newick['children'] :
      if cmtchil==2 :
          strout="("+strout
          strout+='):0,'
          cmtchil=1
      strout+=print_newick(children,Cmt+1)
      if cmtchil==0 :
         strout+=',' 
      cmtchil+=1
  if newick['length']!=None :
     return strout+'):'+str(length+newick['length'])
  else :
     return strout+')'

          
Test=False
if Test :
  AAA=parse('((A:01,C:23,D:5):2,G:89)')
  print(AAA)
  print(print_newick(AAA))
  sys.exit()

FileIn=sys.argv[1]
FileOut=sys.argv[2]
ReadIn=open(FileIn)
Writ=open(FileOut, 'w')
ListVec=ReadIn.readlines()
for lines in ListVec :
    Cmt=0
    treeform=parse(lines.replace('\n', ''))
    treeformnet=print_newick(treeform)+';'+'\n'
    Writ.write(treeformnet)
