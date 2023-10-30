#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 11:50:30 2022

@author: shannon

Code Written by Noah Bagazinski
=============================================================================================

The Target of this Script is to Define the functions and parameters to define a ship hull. 

This Code is a continuation of Noah Bagazinski's work in ship hull genneration using a 
parametric design scheme.  

=============================================================================================
The hull parameterization is defined in five chunks:
    
1) General Hull Form. 
    
    This includes metrics such as the beam  at the deck, the 
    hull taper, and the stern taper.
    
2) The parallel mid body cross section. 
    
    This includes parameters for the flare, chine radius, deadrise, and keel radius
    
3) The Bow Form. 
    
    This includes functions that define the bow rise (rake), the profile drift angle with respect to depth,
    the profile of rocker at the front of the ship, and the profile of the location where the full breadth 
    is achieved for a given depth.
    
4) The Stern Form
    
    This includes ADD DETAIL HERE WHEN FIGURED OUT
    
5) Bulb Forms
    
    This includes HOW TO FORM A BULBOUS BOW AND A BULBOUS STERN
"""


#import all the goodies:
import numpy as np
from matplotlib import pyplot as plt

'''
#load the design vector:
#Vec = np.load('DesVec.npy')
Vec = [100.0, #LOA
       25.0,  #Lb
       25.0,  #Ls
       15.0,  #Bd
       10.0,  #Dd       
       12.0,  #Bs
       3.0,   #WL
       9.5,   #Bc
       6.0,   #Dc
       20.0,  #Beta
       0.0,   #Rc
       1.0,   #Rk
       0.02,  #Abow
       0.01,  #Bbow
       0.05,  #BK_z
       0.9,   #Kappa
       2.0,   #Adel
       2.0,   #Bdel
       30.0,  #Adrft
       -15.0, #Bdrft
       10.0   #Cdrft
       
       ]

Vec = np.array(Vec)

'''
class Hull_Parameterization:
  
    #Define parameters of targethull
    def __init__(self, inputs):
        '''
        inputs is a numpy vector that represents the parameterization of a hull.
        the instantiation function generates all of the constants and factors used 
        in defining the parameterization of the hull. 
        '''
        
        self.LOA = inputs[0]
        self.Lb = inputs[1]
        self.Ls = inputs[2]
        self.Bd = inputs[3]/2.0 #half breadth
        self.Dd = inputs[4]
        self.Bs = inputs[5]
        self.WL = inputs[6]
        self.Bc = inputs[7]/2.0 #half breadth
        self.Dc = inputs[8]
        self.Beta = inputs[9]
        self.Rc = inputs[10]
        self.Rk = inputs[11]
        self.BOW = np.zeros((3,))
        self.BOW[0:2] = inputs[12:14]
        self.BK = np.zeros((2,))
        self.BK[1] = inputs[14] #BK_z is an input - BK_x is solved for
        self.Kappa = inputs[15]
        self.DELTA_BOW = np.zeros((3,))
        self.DELTA_BOW[0:2] = np.array(inputs[16:18])
        self.DRIFT = np.array(inputs[18:21])
   
        

        
        
        
        
        
    '''
    =======================================================================
                        Section 1: General Hull Form
    =======================================================================
    
    The General hull form is characterized by 5 characteristics:
        
        0) LOA -> length overall of the vessel in [m] or = 1
        1) Lb  -> length of the bow taper in [m] or fraction of LOA
        2) Ls  -> length of the stern taper in [m] or fraction of LOA
        3) Bd  -> Beam at the top deck of the vessel in [m] or fraction of LOA
        4) Dd  -> Depth of the vessel at the deck in [m] or fraction of LOA
        5) Bs  -> Beam at the stern in [m] or fraction of LOA
        6) WL  -> Waterline depts in [m] or fraction of LOA
        
    Constraints / NOTES to ensure realistic sizing/ shape of a hull: 
        0) The length of the parallel mid body is equal to LOA-Lb-Ls = Lm
        1) Lb + Ls <= LOA
        2) Bd is not necessarily the maximum beam of the vessel.  It is only the breadth of the 
            main deck. BOA is calculated in the Section 2: Cross Section 
        3) 0 <= Bs <= BOA
        4) Lb is used to define the limits of the bow taper from the forwardmost point on the 
            bow rake to the point where the parallel mid-body starts. The profile of the ship 
            at different waterlines is dependent of the other parameters defined later in the 
            parameterization.
        5) Ls is used to define the limits of the stern taper from the aftmost point on the 
            stern rake to the point where the parallel mid-body ends. The profile of the ship 
            at different waterlines is dependent of the other parameters defined later in the 
            parameterization. 
        6) WL < Dd
        7) All variables are positive or 0
    '''
    def GenGeneralHullform(self):
        '''
        This funciton computes the other form factors of the general hullform
        that can be calculate from the inputs
        '''
        self.Lm = self.LOA - self.Ls - self.Lb    
    
    def GenralHullformConstraints(self):
        '''
        This function checks that constraints are satisfied for the hullfrom.
        If no constraint violations are found, 
        '''
        C = np.array([self.LOA >= self.Ls+self.Lb, self.WL < self.Dd])
        return C
    '''
    =======================================================================
                        Section 2: Cross Section
    =======================================================================
    
    The Cross Section is defined by the following inputs:
        0) Bd   -> The Beam at the Deck in [m] or fraction of LOA
        1) Dd   -> The Depth of the Deck in [m] or fraction of LOA
        2) Bc   -> The Beam at the Chine  (intersection) in [m] or fraction of LOA
        3) Dc   -> The Depth of the Chine (intersection) in [m] or fraction of LOA
        4) Beta -> The deadrise angle in degrees 
        5) Rc   -> The Chine Radius in [m] or fraction of LOA
        6) Rk   -> The keel Radius in [m] or fraction of LOA
    
    Constraints/ NOTES to ensure realistic sizing/ shape of a hull: 
        0) 0 <= Dc < Dd
        1) 0 <= Beta <= 90
        2) Rc and Rk are agebraically limited to ensure that the radius can exist with the 
            given Bd,Dd,BcdC, and Beta values. 
    
    '''
    ## cross section of pmb

    def GenCrossSection(self):
        '''
        This function calculates the constants and other form factors that will allow future
        analysis of the cross section.

        '''
        # Upper Gunwhale Line: A*z + B*y + C = 0, where UG = [A,B,C]
        A = np.array([[self.Dc, self.Bc, 1.0],
                      [self.Dd, self.Bd, 1.0],
                      [1.0, 1.0, 1.0]])
        
        b = np.array([0.0,0.0,1.0])
        
        self.UG = np.linalg.solve(A,b)
        
        # Lower Gunwhale line: A*z + B*y + C = 0 where LG = [A,B,C]. -B/A = tan(Beta)
        #NOTE: Machine Precision makes this feasible for us since tan(90) = 10^16 here
        A = np.array([[self.Dc, self.Bc, 1.0],
                      [np.sin(self.Beta*np.pi/180.0), np.cos(self.Beta*np.pi/180.0), 0.0],
                      [1.0, 1.0, 1.0]])
        
        b = np.array([0.0,0.0,1.0])
        
        self.LG = np.linalg.solve(A,b)

        
        #make math more readable in the next step
        A1 = self.UG[0]
        B1 = self.UG[1]
        theta = np.arctan2(-B1,A1)
        #print(theta)
        A2 = self.LG[0]
        B2 = self.LG[1]
        
        
        #Find the position of the center of the chine radius and intersections with the chine radius
        
        self.Rc_Center = np.zeros((2,)) #(y,z) pair for center
        self.Rc_UG_int = np.zeros((2,)) #(y,z) pair for intersection of chine radius and UG line
        self.Rc_LG_int = np.zeros((2,)) #(y,z) pair for intersection of chine radius and LG line
        
        # This matrix is used to solve for the center and intersections of chine radius: v = inv(A)*b
        A = np.array([[0.0, 0.0, B1, A1, 0.0, 0.0], 
                      [0.0, 0.0, 0.0, 0.0, B2, A2],
                      [0.0, 1.0, 0.0, 0.0, 0.0, -1.0], 
                      [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                      [0.0, 1.0, 0.0, -1.0, 0.0, 0.0],
                      [-1.0, 0.0, 1.0, 0.0, 0.0, 0.0]])
        #print(np.shape(A))  
        
        b = np.array([-self.UG[2], 
                      -self.LG[2], 
                      self.Rc*np.cos(self.Beta*np.pi/180.0), 
                      self.Rc*np.sin(self.Beta*np.pi/180.0), 
                      self.Rc*np.cos(theta), 
                      self.Rc*np.sin(theta)])
        
        c = np.linalg.solve(A,b)
        
        self.Rc_Center = c[0:2] # center of the chine radius
        self.Rc_UG_int = c[2:4] # intersection btwn chine radius and upper gunwhale
        self.Rc_LG_int = c[4:6] # intersection of chine radius and lower gunwhale
        
        
        #Find the position of the center of the chine radius and intersections with the Keel radius
        self.Rk_Center = np.zeros((2,))
        self.Rk_LG_int = np.zeros((2,))
        self.Z_bottom = 0.0
        
        del A, b, c
        
        A = np.array([[-1.0, 0.0, 1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, -1.0, 0.0],
                      [0.0, 0.0, B2, A2, 0.0],
                      [0.0, 1.0, 0.0, 0.0, -1.0],
                      [1.0, 0.0, 0.0, 0.0, 0.0]])
                     
        #print(np.shape(A))  
          
        b = np.array([self.Rk*np.sin(self.Beta*np.pi/180.0),
                      self.Rk*np.cos(self.Beta*np.pi/180.0),
                      -self.LG[2],
                      self.Rk*(0.5 + 0.5*np.sign(self.Rk)),
                      -self.Rk*(0.5 - 0.5*np.sign(self.Rk))])
            
        c = np.linalg.solve(A,b)
         
        self.Rk_Center = c[0:2]
        self.Rk_LG_int = c[2:4]
        self.Z_bottom = c[4]
        
        self.Translate_Z() 
    
    def Translate_Z(self):
        # This funciton translate the z components of the design vector to correct that z_bottom = 0. 
        # This will also adjust any external design vectors in the future as well.
        # new_Z = org_z - Z_bottom
        Z = self.Z_bottom
        
        self.Dc = self.Dc - Z
        self.Dd = self.Dd - Z
        self.Rc_Center[1] = self.Rc_Center[1] - Z        
        self.Rk_Center[1] = self.Rk_Center[1] - Z  
        self.Rc_UG_int[1] = self.Rc_UG_int[1] - Z
        self.Rc_LG_int[1] = self.Rc_LG_int[1] - Z
        self.Rk_LG_int[1] = self.Rk_LG_int[1] - Z
        self.UG[2] = -(self.UG[0]*self.Dc + self.UG[1]*self.Bc)
        self.LG[2] = -(self.LG[0]*self.Dc+ self.LG[1]*self.Bc)
        self.Z_bottom = self.Z_bottom - Z
        
     
        
    def CrossSectionConstraints(self):
        
        C = [self.Rc_UG_int[1] >= self.Dc, self.Rc_LG_int[0] <= self.Bc, self.Rk_LG_int[0] <= self.Rc_LG_int[0]]
        return C
    
    def halfBeam_MidBody(self, z):
        # This funtion calculates the half beam of the cross section at a given height, z
        # If 0 > z or Dd < z, then the function returns -1 as an error
        if z < 0.0 or z > self.Dd:
            return -1
        elif z >= 0.0 and z < self.Rk_LG_int[1]:
            return np.sign(self.Rk)*np.sqrt((self.Rk**2) - (z-self.Rk_Center[1])**2) + self.Rk_Center[0]
        elif z >= self.Rk_LG_int[1] and z < self.Rc_LG_int[1]:
            return -(self.LG[0] * z + self.LG[2])/self.LG[1]
        elif z >= self.Rc_LG_int[1] and z < self.Rc_UG_int[1]:
            return np.sqrt((self.Rc**2) - (z-self.Rc_Center[1])**2) + self.Rc_Center[0]
        else:
            return -(self.UG[0] * z + self.UG[2])/self.UG[1]


    def plotCrossSection(self):
        
        # Plot intersection points in blue
        # Plot chine pt in green
        # Plot Center of Rc and Rk in red
        # half Beam(z) in black
        
        z = np.linspace(0.0, self.Dd, num = 200)
        y = np.zeros((200,))
        for i in range(0,len(z)):
            y[i] = self.halfBeam_MidBody(z[i])
        
        
        plt.subplot(1,1,1)
        plt.axis('equal')
        #plt.axis([0,10,0,10])
        
        plt.plot([self.Bd, self.Rc_UG_int[0], self.Rc_LG_int[0], self.Rk_LG_int[0], 0.0], 
                 [self.Dd, self.Rc_UG_int[1], self.Rc_LG_int[1], self.Rk_LG_int[1], self.Z_bottom], 'o', color = 'blue')
        plt.plot([self.Rc_Center[0], self.Rk_Center[0]], [self.Rc_Center[1], self.Rk_Center[1]],'o' ,color = 'red')  
        plt.plot([self.Bc], [self.Dc],'o' ,color = 'green')  
        plt.plot(y,z,'-', color = 'black', linewidth = 0.75)
        
    
    
    
    '''
    =======================================================================
                        Section 3: Bow Form
    =======================================================================
    
    The Bow Form is defined by the following inputs:
        0) Dd   -> The Depth of the Deck in [m] or fraction of LOA
        1) Abow  -> The z^2 term for Bow(z) that defines the profile of the bowrise
        2) Bbow  -> The z term for Bow(z) that defines the profile of the bowrise
        3) BK_z -> The Z Point of the intersection of the Bow rise and keel rise as percentage of Dd
        4) Kappa-> The X position where the Keel rise begins. percentage of Lb
        5) Adel -> z^2 term for delta(z), the x position where the max Beam is achieved for a given height,  x is normalized to be between 0 and 1 as a percentage of Lb
        6) Bdel -> z term for delta(z), the x position where the max Beam is achieved for a given height
        7) Adrft-> z^2 term for drift(z), the drift angle along the bowrise and keel rise
        8) Bdrft-> z term for drift(z), the drift angle along the bowrise and keel rise
        9) Cdrft-> const term for drift(z), the drift angle along the bowrise and keel rise
    
    These Parameters solve for 4 functions:
        0) Bow(z)   -> gives the X position of the bow rise in the form Az^2 + Bz + C
        1) Keel(x)  -> gives the z height of the keel rise with respect to X in the form A*(X-Kappa*Lb)^2
        2) Delta(z) -> gives the x position between 0 and Lb where the full breadth is achieved for a given z: A(z/Dd)^2 + B(z/Dd) + C = x/Lb
        3) Drift(z) -> gives the drift angle of the bow for a given z: Az^2 + Bz + C
    
    These four functions define the following curve for each z:
        halfBeam_Bow(x) = Y(x) = A*sin(w*x + a) + B for all z between 0 and Dd
        Since we know two points and the derivatives of  those two points
    
    Constraints/ NOTES to ensure realistic sizing/ shape of a hull: 
        0) Kappa*Lb < delta(z=0)
        1) 0 < drift(z) < 90 for 0 <= z <= Dd (only need to check at z = 0, Dd, and -B/(2*A) if within range of z )
        2) 0 <= BK_x < Kappa*Lb
        3) 0 <= BK_z < Dd
        4) delta(z) > Bow(z) and Keel(z) for 0 <= z <= Dd 
    
    '''
    
    def GenBowform(self):
        '''
        This funciton computes the other form factors of the Bowform
        that can be calculated from the inputs
        '''
        # Update BK_z after Dd is readjusted:
            
        self.BK[1] = self.BK[1]*self.Dd
        
        if self.BOW[0] == 0:
            Zv = -1.0
        else:
           Zv = -self.BOW[1]/(2*self.BOW[0]) #Find Z of vertex of bowrise(z)
        
        C = np.array([self.BOW[0]*self.Dd**2.0 + self.BOW[1]*self.Dd,   #Bow rise protrusion at Deck
                      self.BOW[0]*self.BK[1]**2.0 + self.BOW[1]*self.BK[1], #Bow rise protrusion at Bow-keel  intersection
                      self.BOW[0]*Zv**2.0 + self.BOW[1]*Zv]) #Bowrise protrusio at vertex of bow rise eqn
                      
        if (Zv >= self.BK[1]*self.Dd and Zv <= self.Dd):
            self.BOW[2] = -np.amin(C) # If the vertex is between the BK intersect and the Deck, then it is included in the min search
        else:
            self.BOW[2] = -np.amin(C[0:2])
            
        
        # X Position of BK intersect
        self.BK[0] = self.bowrise(self.BK[1])
        
        
        # Calculate the Keelrise equation: it is of the form X = sqrt(Z/A) + Kappa*Lb or Z = A(X-K*Lb)**2, where self.Keel = A
        self.KEEL = self.BK[1]/((self.BK[0]-self.Kappa*self.Lb)**2.0)
        
        
        #Calculate the C for the Delta equation, where C is the constant such that max(Delta(z)) = 0 between 0 and Dd
        if self.DELTA_BOW[0] == 0:
            Zv = -1.0
        else:
            Zv = -self.DELTA_BOW[1]/(2*self.DELTA_BOW[0]) #Find Z of vertex of Delta(z)
        
        C = np.array([self.DELTA_BOW[0]*self.Dd**2.0 + self.DELTA_BOW[1]*self.Dd,   #BDelta at Deck
                      0.0, #As is, Delta(0) = 0
                      self.DELTA_BOW[0]*Zv**2.0 + self.DELTA_BOW[1]*Zv]) #Bowrise protrusio at vertex of bow rise eqn
                      
        if (Zv >= 0.0 and Zv <= self.Dd):
            self.DELTA_BOW[2] = -np.amax(C) # If the vertex is between z = 0  and the Deck, then it is included in the search
        else:
            self.DELTA_BOW[2] = -np.amax(C[0:2])      
    
    #The following funcitons return the 
        
    def bowrise(self, z):
        #returns the x position of the bowrise for a given z for BK_z <= z <= Dd
        return self.BOW[0]*z**2.0 + self.BOW[1]*z + self.BOW[2]
    
    def keelrise_bow(self, z):
        #returns the x position of the keelrise  at the bow for a given z for 0 <= z <= Bk_z
        return -np.sqrt(z/self.KEEL) + self.Kappa*self.Lb
    
    def delta_bow(self, z):
        #returns the x position where the full cross section width is achieved for a given z for 0 <= z <= Dd
        return self.Lb + self.DELTA_BOW[0]*z**2.0 + self.DELTA_BOW[1]*z + self.DELTA_BOW[2]
    
    
    def drift(self, z):
        return self.DRIFT[0]*z**2.0 + self.DRIFT[1]*z + self.DRIFT[2]
        
        
    def solve_waterline_bow(self,z):
        #this function solves for a cubic function: y(half beam) = Ax^3 + Bx^2 + CX + D for the half beam of the profile between the bow/keel rise and delta for a given z for 0 <= z <= Dd

        X1 = self.bow_profile(z)
        
        X2 = self.delta_bow(z)
        
        Y2 = self.halfBeam_MidBody(z)
        
        A = np.array([[X1**3.0, X1**2.0, X1, 1.0],
                      [3.0*X1**2.0, 2*X1, 1.0, 0.0],
                      [X2**3.0, X2**2.0, X2, 1.0],
                      [3.0*X2**2.0, 2.0*X2, 1.0, 0.0]])
            
        b = np.array([0.0,
                      np.tan(np.pi*self.drift(z)/180.0),
                      Y2,
                      0.0])
        return np.linalg.solve(A,b)
    
    def bow_profile(self, z):
        if z <= self.BK[1]:
            X1 = self.keelrise_bow(z)
        else:
            X1 = self.bowrise(z)
        return X1
    
    def halfBeam_Bow(self, x, PROF):
        #returns the halfbeam along the bow taper between the bow/keel rise and delta(z), PROF is the output of solve)waterline_bow(z)
        return PROF[0]*x**3.0 + PROF[1]*x**2.0 + PROF[2]*x + PROF[3]
    
    
    def gen_waterline_bow(self,x, z):
   
        prof = self.solve_waterline_bow(z)
        y = self.halfBeam_Bow(x, prof)
            
        return y
        
    def BowformConstraints(self):
        #
        C = []
        return C
    
    
    '''
    =======================================================================
                        Section 4: Stern Form
    =======================================================================
    '''
    
    
    
    
    
    
    
    '''
    =======================================================================
                        Section 5: Bulb Forms
    =======================================================================
    '''


