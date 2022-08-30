# -*- coding: utf-8 -*-
#Edited by christian bou chahine
#Complementary fonctions

def minmin(L):
    M=[]
    for i in range(len(L)):
        M.append(min(L[i]))
    return min(M)

def maxmax(L):
    M=[]
    for i in range(len(L)):
        M.append(max(L[i]))
    return max(M)
