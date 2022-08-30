# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#MILLARD MODEL
#From "Flexing Computational Muscle: Modeling and Simulation of Musculotendon Dynamics"
#Written by: Mattew Millard, Thomas Uchida, Ajay Seth, Scott L.Delp
#--> The original code:
    #https://github.com/mjhmilla/Millard2012EquilibriumMuscleMatlabPort/blob/master/src/calcIndex.m
#--------------------------------------------------------------------------------------
#Edited by christian bou chahine

def calcIndex(x, ptsM,moduloRange,tol): #options=[moduloRange,tol]
 
    """Given a series of sequential ranges in the 2xRealN matrixReal
    % ptsM, this function will calculate the row indexReal that
    % has an interval ptsM(n,1)-ptsM(n,2) that contains the
    % point xReal. If xReal lies on a border that defines a sub-interval
    % the function will return the indexReal such that xReal is at the 
    % start of hte interval.
    %
    % @param xReal: a double 
    % @param ptsM: a 2 xReal n matrixReal that defines the range of n
    %              sub intervals. For exRealample
    %
    %
    %        ptsM = [0 1 2 3 4; ...
    %                1 2 3 4 5]
    %
    %        defines 5 adjacent sub intervals 
    %        [[0,1],[1,2],[2,3],[3,4],[4,5]]
    %
    % @param tol: how close a value xReal is allowed to be to a
    %             sub interval border before it is declared to be
    %             on the border.
    %
    % @returns idx: the column indexReal of ptsM that contains an
    %               interval that includes xReal."""
    

    xReal = x.real
    idx = 0
    flag_found = 0
    
    rows = ptsM.shape[0]   
    
    #xReal is either monotonically increasing or decreasing
    dxSign = 1
    
    if ptsM[1,0]-ptsM[0,0] < 0:
        dxSign = -1
    
    for i in range(ptsM.shape[1]):
        
        if dxSign == 1:
            if (xReal >= ptsM[0,i] and xReal < ptsM[rows-1,i]) or (abs(xReal-ptsM[0,i]) < tol or abs(xReal-ptsM[rows-1,i]) < tol) :
                idx = i
                flag_found = 1

        else: 
            if (xReal <= ptsM[0,i] and xReal > ptsM[rows-1,i]) or (abs(xReal-ptsM[0,i]) < tol or abs(xReal-ptsM[rows-1,i]) < tol ):
                idx = i
                flag_found = 1

    #Check if the value xReal is identically the last point
    if flag_found == 0 and xReal == ptsM[rows-1,ptsM.shape[1]-1]: 
        idx = ptsM.shape[1]-1
        flag_found = 1
    
    assert flag_found == 1,'Error: A value of xReal was used that is not within the Bezier curve set.'
    
    return idx
