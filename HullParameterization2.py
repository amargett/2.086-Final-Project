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
from scipy.optimize import fsolve
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
        self.Bs = inputs[5]/2.0 #half breadth
        self.WL = inputs[6]
        self.Bc = inputs[7]/2.0 #half breadth
        #self.Dc = inputs[8]
        self.Beta = inputs[8]
        self.Rc = inputs[9]
        self.Rk = inputs[10]
        self.BOW = np.zeros((3,))
        self.BOW[0:2] = inputs[11:13]
        self.BK = np.zeros((2,))
        self.BK[1] = inputs[13] #BK_z is an input - BK_x is solved for
        self.Kappa_BOW = inputs[14]
        self.DELTA_BOW = np.zeros((3,))
        self.DELTA_BOW[0:2] = np.array(inputs[15:17])
        self.DRIFT = np.array(inputs[17:20])
        self.TRANS = np.zeros((2,))
        self.TRANS[0] = inputs[20]
        self.SK = np.zeros((2,))
        self.SK[1] = inputs[21]
        self.Kappa_STERN = inputs[22]
        self.DELTA_STERN = np.zeros((3,))
        self.DELTA_STERN[0:2] = np.array(inputs[23:25])
        self.RY_STERN = np.array(inputs[25:27])
        self.RX_STERN = np.array(inputs[27:29])
        self.Beta_trans = inputs[29]
        self.Bc_trans = inputs[30]/2.0 # half breadth
        self.Rc_trans = inputs[31]
        self.Rk_trans = inputs[32]
        self.CONVERGE = np.array(inputs[33:36])

               
        #Generate and Check the Forms of the Overall Hull
        
        self.GenGeneralHullform()
        C1 = print(self.GenralHullformConstraints())
        
        self.GenCrossSection()
        C2 = print(self.CrossSectionConstraints())
        
        self.GenBowForm()
        C3 = print(self.BowformConstraints())
        
        self.GenSternForm()
        C4 = print(self.sternFormConstraints())
        
        
        
        
        
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
        C = np.array([self.LOA >= self.Ls+self.Lb,
                      self.WL < self.Dd])
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
        1) Rc and Rk are agebraically limited to ensure that the radius can exist with the 
            given Bd,Dd,BcdC, and Beta values. 
    
    '''
    def GenCrossSection(self):
        '''
        This function calculates the constants and other form factors that will allow future
        analysis of the cross section.

        '''
        
        #(y,z) pair for center of keel radius
        self.Rk_Center = np.array([-self.Rk*(0.5 - 0.5*np.sign(self.Rk)), 
                                    self.Rk*(0.5 + 0.5*np.sign(self.Rk))])
        #(y,z) pair for intersection of keel radius and LG line at the transom
        self.Rk_LG_int = np.array([self.Rk_Center[0] + self.Rk*np.sin(np.pi*self.Beta/180.0),
                                        self.Rk_Center[1] - self.Rk*np.cos(np.pi*self.Beta/180.0)])
       
        
        #solve for the lower gunwhale line: A*z + B*y + C = 0
        A = np.array([[1.0, 1.0, 1.0],
                      [self.Rk_LG_int[1], self.Rk_LG_int[0], 1.0],
                      [-(self.Rk_LG_int[0]-self.Rk_Center[0]), (self.Rk_LG_int[1]-self.Rk_Center[1]), 0.0]])
        b = np.array([1.0, 0.0, 0.0])

        self.LG = np.linalg.solve(A,b)
        
        del A, b     
        
        self.Dc = -(self.LG[1]*self.Bc + self.LG[2])/self.LG[0]
     
        # Upper Gunwhale Line: A*z + B*y + C = 0, where UG = [A,B,C]
        A = np.array([[self.Dc, self.Bc, 1.0],
                      [self.Dd, self.Bd, 1.0],
                      [1.0, 1.0, 1.0]])
        
        b = np.array([0.0,0.0,1.0])
        
        self.UG = np.linalg.solve(A,b)
        
        del A, b
        
        # Calculate terms for the half beam of the cross section of the transom:
        self.Rc_Center = np.zeros((2,)) #(y,z) pair for center of chine radius at the transom
        self.Rc_UG_int = np.zeros((2,)) #(y,z) pair for intersection of chine radius and UG line at the transom
        self.Rc_LG_int = np.zeros((2,)) #(y,z) pair for intersection of chine radius and LG line at the transom
        
        #make math more readable to solve the chine
        A1 = self.UG[0]
        B1 = self.UG[1]
        theta = np.arctan2(-B1,A1)
        beta = self.Beta*np.pi/180.0
        A2 = self.LG[0]
        B2 = self.LG[1]
        
        
        A = np.array([[B1, A1, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, B2, A2, 0.0, 0.0],
                      [1.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                      [0.0, -1.0, 0.0, 0.0, 0.0, 1.0],
                      [0.0, 0.0, 1.0, 0.0, -1.0, 0.0],
                      [0.0, 0.0, 0.0, -1.0, 0.0, 1.0]])
      
        b = np.array([-self.UG[2],
                      -self.LG[2],
                      self.Rc*np.sin(theta),
                      self.Rc*np.cos(theta),
                      self.Rc*np.sin(beta),
                      self.Rc*np.cos(beta)])
        
        C = np.linalg.solve(A,b)
        
        self.Rc_UG_int = C[0:2]
        self.Rc_LG_int = C[2:4]
        self.Rc_Center = C[4:6]
   
        
    def CrossSectionConstraints(self):
        
        C = [self.Rc_UG_int[1] >= self.Dc,
             self.Rc >= 0.0,
             self.Bc > 0.0,
             self.Dc >= 0.0,
             self.Rc_LG_int[0] <= self.Bc,
             self.Rk_LG_int[0] <= self.Rc_LG_int[0]]
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


    def plot_MidBody_CrossSection(self):
        
        # Plot intersection points in blue
        # Plot chine pt in green
        # Plot Center of Rc and Rk in red
        # half Beam(z) in black
        
        z = np.linspace(0.0, self.Dd, num = 200)
        y = np.zeros((200,))
        for i in range(0,len(z)):
            y[i] = self.halfBeam_MidBody(z[i])
        
        
        fig2, ax2 = plt.subplots()
        ax2.axis('equal')
        #plt.axis([0,10,0,10])
        
        ax2.plot([self.Bd, self.Rc_UG_int[0], self.Rc_LG_int[0], self.Rk_LG_int[0], 0.0], 
                 [self.Dd, self.Rc_UG_int[1], self.Rc_LG_int[1], self.Rk_LG_int[1], 0.0], 'o', color = 'blue')
        ax2.plot([self.Rc_Center[0], self.Rk_Center[0]], [self.Rc_Center[1], self.Rk_Center[1]],'o' ,color = 'red')  
        ax2.plot([self.Bc], [self.Dc],'o' ,color = 'green')  
        ax2.plot(y,z,'-', color = 'black', linewidth = 0.75)
        
    
    
    
    '''
    =======================================================================
                        Section 3: Bow Form
    =======================================================================
    
    The Bow Form is defined by the following inputs:
        0) Dd   -> The Depth of the Deck in [m] or fraction of LOA
        1) Lb   -> The length of the bow taper in [m] or fraction of LOA
        2) Abow  -> The z^2 term for Bow(z) that defines the profile of the bowrise
        3) Bbow  -> The z term for Bow(z) that defines the profile of the bowrise
        4) BK_z -> The Z Point of the intersection of the Bow rise and keel rise as percentage of Dd
        5) Kappa_bow-> The X position where the Keel rise begins. percentage of Lb
        6) Adel -> z^2 term for delta(z), the x position where the max Beam is achieved for a given height
        7) Bdel -> z term for delta(z), the x position where the max Beam is achieved for a given height
        8) Adrft-> z^2 term for drift(z), the drift angle along the bowrise and keel rise
        9) Bdrft-> z term for drift(z), the drift angle along the bowrise and keel rise
        10) Cdrft-> const term for drift(z), the drift angle along the bowrise and keel rise
    
    These Parameters solve for 4 functions:
        0) Bow(z)   -> gives the X position of the bow rise in the form Az^2 + Bz + C
        1) Keel_BOW(x)  -> gives the z height of the keel rise with respect to X in the form A*(X-Kappa_BOW*Lb)^2
        2) Delta_BOW(z) -> gives the x position between 0 and Lb where the full breadth is achieved for a given z: A(z/Dd)^2 + B(z/Dd) + C = x/Lb
        3) Drift(z) -> gives the drift angle of the bow for a given z: Az^2 + Bz + C
    
    These four functions define the following curve for each z:
        halfBeam_Bow(x) = Y(x) = A*x^3 + Bx^2 + Cx + D for all z between 0 and Dd
        Since we know two points and the derivatives of  those two points
    
    Constraints/ NOTES to ensure realistic sizing/ shape of a hull: 
        0) Kappa_BOW*Lb < delta(z=0)
        1) 0 < drift(z) < 90 for 0 <= z <= Dd (only need to check at z = 0, Dd, and -B/(2*A) if within range of z )
        2) 0 <= BK_x < Kappa_BOW*Lb
        3) 0 <= BK_z < Dd
        4) delta(z) > Bow(z) and Keel(z) for 0 <= z <= Dd 
    
    '''
    
    def GenBowForm(self):
        '''
        This funciton computes the other form factors of the Bowform
        that can be calculated from the inputs
        '''
        # Update BK_z after Dd is readjusted:
            
        self.BK[1] = self.BK[1]
        
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
        
        
        # Calculate the Keelrise equation: it is of the form X = sqrt(Z/A) + Kappa_Bow*Lb or Z = A(X-K*Lb)**2, where self.Keel = A
        self.KEEL_BOW = self.BK[1]/((self.BK[0]-self.Kappa_BOW*self.Lb)**2.0)
        
        
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
        return -np.sqrt(z/self.KEEL_BOW) + self.Kappa_BOW*self.Lb
    
    def delta_bow(self, z):
        #returns the x position where the full cross section width is achieved for a given z for 0 <= z <= Dd
        return self.Lb + self.DELTA_BOW[0]*z**2.0 + self.DELTA_BOW[1]*z + self.DELTA_BOW[2]
    
    
    def drift(self, z):
        #returns the drift angle in radians
        return np.pi*(self.DRIFT[0]*z**2.0 + self.DRIFT[1]*z + self.DRIFT[2])/180.0
        
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
                      np.tan(self.drift(z)),
                      Y2,
                      0.0])
        return np.linalg.solve(A,b)
    
    def bow_profile(self, z):
        # This assumes that z >= 0 and z <= Dd
        
        if z <= self.BK[1]:
            X1 = self.keelrise_bow(z)
        else:
            X1 = self.bowrise(z)
        return X1
    
    def halfBeam_Bow(self, x, PROF):
        #returns the halfbeam along the bow taper between the bow/keel rise and delta(z), PROF is the output of solve)waterline_bow(z)
        return PROF[0]*x**3.0 + PROF[1]*x**2.0 + PROF[2]*x + PROF[3]
    
    
    def gen_waterline_bow(self, z, NUM_POINTS = 100, SPACING = [0,1], BOOL_PTS_OR_SPACE = 1):
        '''
        This fuction generates a set of points [[X1,Y1] .... [X2,Y2]] that detail the curvature of the bow taper for a given z, for 0 <= z <= Dd
        
        it can either be created as with a number of set of points, or an even spacing based on the global x spacing (better for station plotting) 

        BOOL_PTS_OR_SPACE controls whether a set number of points will produce the waterline (1) or the spacing vector will(0)
        
        '''
       
        x1 = self.bow_profile(z)
        
        x2 = self.delta_bow(z)
        
        prof = self.solve_waterline_bow(z)
               
        
        
        if BOOL_PTS_OR_SPACE:
            x = np.linspace(x1,x2,NUM_POINTS)
            
            XY = np.zeros((NUM_POINTS,2))
            for i in range(0,NUM_POINTS):
                XY[i,:] = [x[i], self.halfBeam_Bow(x[i], prof)]
            
            return XY
        
    def BowformConstraints(self):
        #This fuction returns booleans if the bow constraints are satisfied as detailed above:
            
        #Check that the vertex (Zv) of the drift angle equation satisfies the constraint of the drift angle if
        # if it lies within the bounds of 0 and Dd. If drift(z) is a line, or the vertex is outside the bounds,
        # Then True is returned 
        if self.DRIFT[0] == 0.0:
            Zv = -1.0
        else:
            Zv = -self.DRIFT[1]/(2.0*self.DRIFT[0])
        
        
        if Zv >= 0.0 and Zv <= self.Dd:
            vert_drift = [self.drift(-self.DRIFT[1]/(2.0*self.DRIFT[0])) < np.pi/2.0,
                          self.drift(-self.DRIFT[1]/(2.0*self.DRIFT[0])) >= 0.0]
        else:
            vert_drift = [True, True]
            
            
        #Check that Delta_Bow(z) is always greater than the leading edge of the ship (keelrise(z) and bow(z))
        #Check at z = 0, vertex of delta(z), vertex of bow(z), BKz, Dd
        if self.DELTA_BOW[0] == 0.0:
            Zv = -1.0
        else: 
            Zv = -self.DELTA_BOW[1]/ (2.0*self.DELTA_BOW[0])
        
        if Zv >=0.0 and Zv <= self.Dd:
            vert_delta_bow = (self.delta_bow(Zv) > self.bow_profile(Zv))
        else:
            vert_delta_bow = True
        
        
        if self.BOW[0] == 0.0:
            Zv = -1.0
        else: 
            Zv = -self.BOW[1]/ (2.0*self.BOW[0])
        
        if Zv >=0.0 and Zv <= self.Dd:
            vert_bow = (self.delta_bow(Zv) > self.bow_profile(Zv))
        else:
            vert_bow = True
        
        
        
        
        C = [self.Kappa_BOW*self.Lb < self.delta_bow(0.0),
            self.drift(0.0) < np.pi/2.0,
            self.drift(0.0) >= 0.0,
            self.drift(self.Dd) < np.pi/2.0,
            self.drift(self.Dd) >= 0.0,
            vert_drift[0],
            vert_drift[1],
            self.BK[0] >= 0.0,
            self.BK[0] < self.Kappa_BOW*self.Lb,
            self.BK[1] >= 0.0,
            self.BK[1] < self.Dd,
            self.delta_bow(self.Dd) > self.bow_profile(self.Dd),
            self.delta_bow(self.BK[1]) > self.BK[0],
            vert_delta_bow,
            vert_bow]
        return C
    
    
    '''
    =======================================================================
                        Section 4: Stern Form
    =======================================================================
        The Stern Form is defined by the following inputs:
        0) Bs   -> The width of the stern at the deck of the ship in [m] or fraction of LOA
        1) Ls   -> The length of the stern taper in [m] or fraction of LOA
        2) A_trans -> The A term that defines the transom slope X = Az + B
        3) SKz -> The Z Point of the intersection of the Stern rise and transom as percentage of Dd
        4) Kappa_STERN -> The X position where the Stern rise begins aft of the end of the parallel midbody as a fraction of Ls
        5) Adel -> z^2 term for delta_stern(z), the x position where the max Beam is achieved for a given height,
        6) Bdel -> z term for delta_stern(z), the x position where the max Beam is achieved for a given height
        7) A_Ry-> z term for  Ry(z), the y-raduis of the ellipse at the stern of the ship
        8) B_Ry-> const for Ry(z), the y-raduis of the ellipse at the stern of the ship
        9) A_Rx-> z term for  Rx(z), the x-raduis of the ellipse at the stern of the ship
        10) B_Rx-> const for Rx(z), the x-raduis of the ellipse at the stern of the ship
        11) Bc_trans -> The beam of the chine point at the transom in [m] or fraction of LOA
        12) Dc_trans -> The depth of the chine point at the transom in [m] or fraction of LOA
        13) Rc_trans -> The Chine radius of the chine at the transom in [m] or fraction of LOA
        14) Rk_trans -> the keel radius of the chine at the transom in [m] or fraction of LOA
        15) AconvT -> the z^2 term for Converge Angle(z) the tangent angle of the gunwhale at the transom
        16) BconvT -> the z term for Converge Angle(z) the tangent angle of the gunwhale at the transom
        17) CconvT -> the const term for Converge Angle(z) the tangent angle of the gunwhale at the transom
    
    
    These Parameters solve for 7 functions:
        0) Transom(z)   -> gives the X position of the transom in the form  Az + B
        1) Sternrise(x)  -> gives the z height of the stern rise with respect to X in the form A*(X-Kappa*Ls)^2
        2) Delta_Stern(z) -> gives the x position between LOA-Ls and LOA where the full breadth is achieved for a given z: A(z)^2 + B(z) + C = X
        3) Converge(z) -> gives the convergence tangent angle of the gunwhale at the transom for a given z: Az^2 + Bz + C
        4) Ry(z) -> gives the y radius of the stern ellipse in the form Ry = Az + B
        5) Rx(z) -> gives the x radius of the stern ellipse in the form Rx = Az + B
        6) halfBeam_transom(z) -> gives the halfbeam of the transom for z between SKz and Dd
    
    These four functions define the following curve for each z:
        halfBeam_Stern(x) = Y(x) = Parabola + Ellipse for all z between 0 and Dd
    
    Constraints/ NOTES to ensure realistic sizing/ shape of a hull: 
        0) Lb+Lm + Kappa_Stern*Ls > delta_Stern(z=0)
        1) 0 < converge(z) < 90 for 0 <= z <= Dd (only need to check at z = 0, Dd, and -B/(2*A) if within range of z )
        2) 0 <= SK_x > Lb+Lm+ Kappa*Ls
        3) 0 <= SK_z < Dd
        4) delta(z) < Transom(z) and Sternrise(z) for 0 <= z <= Dd 
    
    
    
    '''
    def GenSternForm(self):
        
        # Recalculate SK to be a value instead of a percentage
        self.SK[1] = self.SK[1]*self.Dd
        
        # Solve for the B value such that max(Transom(z)) = LOA
        if self.TRANS[0] >= 0.0:
            self.TRANS[1] = self.LOA - self.TRANS[0]*self.Dd
        else:
            self.TRANS[1] = self.LOA - self.TRANS[0]*self.SK[1]
        
        #calculate the x value for the SK intersect
        self.SK[0] = self.transom(self.SK[1])
        
        # find the constant term in the sternrise equation: z = A(x-Lb+Lm+Ls*Kappa_stern)^2
        self.STERNRISE = self.SK[1]/(self.SK[0] - (self.Lb + self.Lm + self.Ls*self.Kappa_STERN))**2.0
        
        #Calculate the C for the Delta_stern equation, where C is the constant such that min(Delta_stern(z)) = 0 for z between 0 and Dd
        if self.DELTA_STERN[0] == 0:
            Zv = -1.0
        else:
            Zv = -self.DELTA_STERN[1]/(2*self.DELTA_STERN[0]) #Find Z of vertex of Delta(z)
        
        C = np.array([self.DELTA_STERN[0]*self.Dd**2.0 + self.DELTA_STERN[1]*self.Dd,   #Stern Delta at Deck
                      0.0, #As is, Delta_Stern(0) = 0
                      self.DELTA_STERN[0]*Zv**2.0 + self.DELTA_STERN[1]*Zv]) #vertex of Delta_STERN equation
                      
        if (Zv >= 0.0 and Zv <= self.Dd):
            self.DELTA_STERN[2] = np.amin(C) # If the vertex is between z = 0  and the Deck, then it is included in the search
        else:
            self.DELTA_STERN[2] = np.amin(C[0:2])
        
        
        #(y,z) pair for center of keel radius
        self.Rk_Center_trans = np.array([-self.Rk_trans*(0.5 - 0.5*np.sign(self.Rk_trans)), 
                                         self.SK[1] + self.Rk_trans*(0.5 + 0.5*np.sign(self.Rk_trans))])
        #(y,z) pair for intersection of keel radius and LG line at the transom
        self.Rk_LG_int_trans = np.array([self.Rk_Center_trans[0] + self.Rk_trans*np.sin(np.pi*self.Beta_trans/180.0),
                                        self.Rk_Center_trans[1] - self.Rk_trans*np.cos(np.pi*self.Beta_trans/180.0)])
       
        
        #solve for the lower gunwhale line: A*z + B*y + C = 0
        A = np.array([[1.0, 1.0, 1.0],
                      [self.Rk_LG_int_trans[1], self.Rk_LG_int_trans[0], 1.0],
                      [-(self.Rk_LG_int_trans[0]-self.Rk_Center_trans[0]), (self.Rk_LG_int_trans[1]-self.Rk_Center_trans[1]), 0.0]])
        b = np.array([1.0, 0.0, 0.0])

        self.LG_trans = np.linalg.solve(A,b)
        
        del A, b     
        
        self.Dc_trans = -(self.LG_trans[1]*self.Bc_trans + self.LG_trans[2])/self.LG_trans[0]
     
        # Upper Gunwhale Line: A*z + B*y + C = 0, where UG = [A,B,C]
        A = np.array([[self.Dc_trans, self.Bc_trans , 1.0],
                      [self.Dd, self.Bs, 1.0],
                      [1.0, 1.0, 1.0]])
        
        b = np.array([0.0,0.0,1.0])
        
        self.UG_trans = np.linalg.solve(A,b)
        
        del A, b
        
        # Calculate terms for the half beam of the cross section of the transom:
        self.Rc_Center_trans = np.zeros((2,)) #(y,z) pair for center of chine radius at the transom
        self.Rc_UG_int_trans = np.zeros((2,)) #(y,z) pair for intersection of chine radius and UG line at the transom
        self.Rc_LG_int_trans = np.zeros((2,)) #(y,z) pair for intersection of chine radius and LG line at the transom
        
        #make math more readable to solve the chine
        A1 = self.UG_trans[0]
        B1 = self.UG_trans[1]
        theta = np.arctan2(-B1,A1)
        beta = self.Beta_trans*np.pi/180.0
        A2 = self.LG_trans[0]
        B2 = self.LG_trans[1]
        
        
        A = np.array([[B1, A1, 0.0, 0.0, 0.0, 0.0],
                      [0.0, 0.0, B2, A2, 0.0, 0.0],
                      [1.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                      [0.0, -1.0, 0.0, 0.0, 0.0, 1.0],
                      [0.0, 0.0, 1.0, 0.0, -1.0, 0.0],
                      [0.0, 0.0, 0.0, -1.0, 0.0, 1.0]])
      
        b = np.array([-self.UG_trans[2],
                      -self.LG_trans[2],
                      self.Rc_trans*np.sin(theta),
                      self.Rc_trans*np.cos(theta),
                      self.Rc_trans*np.sin(beta),
                      self.Rc_trans*np.cos(beta)])
        
        C = np.linalg.solve(A,b)
        
        self.Rc_UG_int_trans = C[0:2]
        self.Rc_LG_int_trans = C[2:4]
        self.Rc_Center_trans = C[4:6]
        
        
        
    
    def transom(self, z):
        #returns the x position of the transom for a given z fr SK_z <= z <= Dd
        return self.TRANS[0]*z + self.TRANS[1]    
    
    def sternrise(self, z):
        #returns the x position of the sternrise for a given z for 0 <= z <= SK_z
        return np.sqrt(z/self.STERNRISE) + self.Lb + self.Lm + self.Ls*self.Kappa_STERN
    
    def stern_profile(self, z):
        if z <= self.SK[1]:
            return self.sternrise(z)
        else:
            return self.transom(z)
    
    def delta_stern(self, z):
        #returns the starting position of the stern taper at a given heigt
        return self.Lb + self.Lm + self.DELTA_STERN[0]* z**2.0 + self.DELTA_STERN[1]*z + self.DELTA_STERN[2]
   
    def Rx_stern(self, z):
        #Assumes z is between 0 and Dd
        # Returns the x radius of curvature for the stern ellipse
        if z < self.SK[1]:
            return self.RX_STERN[0]*z +self.RX_STERN[1]
        else:
            return self.RX_STERN[0]*self.SK[1] + self.RX_STERN[1]
        
    def Ry_stern(self,z):
        #Assumes z is between 0 and Dd
        # Returns the y radius of curvature for the stern ellipse
        if z < self.SK[1]:
            return self.RY_STERN[0]*z +self.RY_STERN[1]
        else:
            return self.RY_STERN[0]*self.SK[1] + self.RY_STERN[1]        
     
    def converge(self,z):
        #returns the convergence angle  in radians between the gunwhale and the transom in degrees (-90 deg = gunwhale curved to be tangent to transom. 0 deg = gunwhale is perpendicular to transom)
        if z < self.SK[1]:
            conv = -np.pi/2.0
        else:
            conv = np.pi*(self.CONVERGE[0]* z**2.0 + self.CONVERGE[1]*z + self.CONVERGE[2])/180.0
        
        return conv
   
    def halfBeam_Transom(self, z):
        #Returns the x,y pair of the transom at a height z. This assumes that SK_z <= z <= Dd. otherwise y returns -1
        x = self.stern_profile(z)
        
        if z < self.SK[1] and z > 0.0:
            y = 0.0
        elif z < 0.0 or z > self.Dd:
            y = -1.0
        
        elif z >= self.SK[1] and z < self.Rk_LG_int_trans[1]:
            y =  np.sign(self.Rk_trans)*np.sqrt((self.Rk_trans**2) - (z-self.Rk_Center_trans[1])**2) + self.Rk_Center_trans[0]
        elif z >= self.Rk_LG_int_trans[1] and z < self.Rc_LG_int_trans[1]:
            y =  -(self.LG_trans[0] * z + self.LG_trans[2])/self.LG_trans[1]
        elif z >= self.Rc_LG_int_trans[1] and z < self.Rc_UG_int_trans[1]:
            y =  np.sqrt((self.Rc_trans**2) - (z-self.Rc_Center_trans[1])**2) + self.Rc_Center_trans[0]
        else:
            y =  -(self.UG_trans[0] * z + self.UG_trans[2])/self.UG_trans[1]
         
        return [x,y]
    
    def plot_Transom_CrossSection(self):
        
        # Plot intersection points in blue
        # Plot chine pt in green
        # Plot Center of Rc and Rk in red
        # half Beam(z) in black
        
        z = np.linspace(self.SK[1], self.Dd, num = 200)
        y = np.zeros((200,2))
        for i in range(0,len(z)):
            y[i] = self.halfBeam_Transom(z[i])
        
        
        fig1,ax1 = plt.subplots()
        ax1.axis('equal')
        #plt.axis([0,10,0,10])
        
        ax1.plot([self.Bs, self.Rc_UG_int_trans[0], self.Rc_LG_int_trans[0], self.Rk_LG_int_trans[0], 0.0], 
                 [self.Dd, self.Rc_UG_int_trans[1], self.Rc_LG_int_trans[1], self.Rk_LG_int_trans[1], self.SK[1]], 'o', color = 'blue')
        ax1.plot([self.Rc_Center_trans[0], self.Rk_Center_trans[0]], [self.Rc_Center_trans[1], self.Rk_Center_trans[1]],'o' ,color = 'red')  
        ax1.plot([self.Bc_trans], [self.Dc_trans],'o' ,color = 'green')  
        ax1.plot(y[:,1],z,'-', color = 'black', linewidth = 0.75)
    
    def halfBeam_Stern(self, x, PROF):
        #returns the halfbeam along the stern taper between delta(z) and stern_profile(z), PROF is the output of solve_waterline_stern(z)
        # assumes x is in bounds of the stern taper
        if x < PROF[7]:
            return PROF[0]*x**2.0 + PROF[1]*x + PROF[2] #parabola
        else:
            return np.sqrt(PROF[4]**2.0 * (1.0 - ((x - PROF[5])/PROF[3])**2.0)) + PROF[6] #ellipse  
        
    
    
    def solve_waterline_stern(self, z):
        #returns PROF, a parabola [A,B,C], an ellipse [Rx, Ry, Cx, Cy], and an intersection [Xint, Yint] of the two curves such that they are tangent at the intersection
        # Accounts for the convergence slope and the half_breath of the transom, assumes that 0 < z < Dd
        # PROF = [A,B,C, Rx, Ry, Cx, Cy, Xint, Yint]
        
        Rx = self.Rx_stern(z)
        Ry = self.Ry_stern(z)
        
        conv = self.converge(z)
        [x_stern, Bs_z] = self.halfBeam_Transom(z)
        
        
         #Algebraically solve for the center of the ellipse because it could not be solved linearly
        A = np.sqrt((-np.tan(conv)**2.0 * Rx**4.0) / (Ry**2.0 - np.tan(conv)**2.0 *Rx**2.0) )
        
        B = np.sqrt(1 - (A/Rx)**2.0) * Ry
        
        
        Cx = x_stern - A
        Cy = Bs_z - B
    
        del A, B
        
        
        #Solve for the parabola and the intersection points
        x0 = np.array([-1,Cx+1, Cy+1])
        
        X_delta = self.delta_stern(z)
        Bm_z = self.halfBeam_MidBody(z)
        
        
        solve_parabola = lambda b: [b[0]*(b[1] - X_delta)**2.0 + Bm_z - b[2], #parabola is at the intersection point
                                    2.0*b[0]*(b[1] - X_delta) + (b[1] - Cx)*Ry**2.0 / (Rx**2.0 * (b[2] - Cy)), #parabola is tangent to ellipse
                                    ((b[1] - Cx)/Rx)**2.0 + ((b[2] - Cy)/Ry)**2.0 - 1.0] #intesect point is on the ellipse
                                    
        
        b = fsolve(solve_parabola, x0)
        
        
        A = b[0]
        B= -2.0*b[0]*X_delta
        C= Bm_z + b[0]*X_delta**2.0
        Xint = b[1]
        Yint = b[2]
        
        PROF =  [A,B,C,Rx,Ry,Cx,Cy,Xint,Yint]
        #print(PROF)
        return PROF
          
    
    def gen_waterline_stern(self, z, NUM_POINTS = 100, SPACING = [0,1], BOOL_PTS_OR_SPACE = 1):
        '''
        This fuction generates a set of points [[X1,Y1] .... [X2,Y2]] that detail the curvature of the bow taper for a given z, for 0 <= z <= Dd
        
        it can either be created as with a number of set of points, or an even spacing based on the global x spacing (better for station plotting) 

        BOOL_PTS_OR_SPACE controls whether a set number of points will produce the waterline (1) or the spacing vector will(0)
        
        '''
       
        x1 = self.stern_profile(z)
        
        x2 = self.delta_stern(z)
        
        prof = self.solve_waterline_stern(z)
               
        
        
        if BOOL_PTS_OR_SPACE:
            x = np.linspace(x1,x2,NUM_POINTS)
            
            XY = np.zeros((NUM_POINTS,2))
            for i in range(0,NUM_POINTS):
                XY[i,:] = [x[i], self.halfBeam_Stern(x[i], prof)]
            
            return XY
      
    
    
    def sternFormConstraints(self):
        # this is an incomplete list of geometric constrains for the hull form
        
        C = [self.delta_stern(0.0) < (self.Lb + self.Lm + self.Ls*self.Kappa_STERN),
             (self.Lb + self.Lm + self.Ls*self.Kappa_STERN) < self.SK[0],
             self.Bc_trans < self.halfBeam_MidBody(self.Dc_trans)]
        
        
        return C
    
    
    
    
    '''
    =======================================================================
                        Section 5: Bulb Forms
    =======================================================================
    '''









