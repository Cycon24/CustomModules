# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 09:45:15 2023

@author: cycon
"""
import sys
import numpy as np
import EngineErrors as EngineErrors
import EnginePerformanceFunctions as EPF

# Use to define the general states/functions shared by each/most stages
class Stage():
    def __init__(self, **kwargs):
        '''
        A general stage class that serves as a baseline for every stage.
        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.R = 287 # J/kg*K
        self.gam_a = 1.4
        self.gam_g = 4/3
        self.cp_a = 1.005 # kJ/kg*K
        self.cp_g = 1.148 # kJ/kg*K
        # General properties that could be used by all stages 
        # so all components know the atm conditions
        self.Ta  = kwargs.get('Ta') 
        self.Pa  = kwargs.get('Pa')
        self.Vinf = kwargs.get('Vinf')
        self.Minf = kwargs.get('Minf')
        
        self.Toi = kwargs.get('Toi')
        self.Poi = kwargs.get('Poi')
        self.Ti  = kwargs.get('Ti')
        self.Pi  = kwargs.get('Pi')
        self.Mi  = kwargs.get('Mi')
        self.Vi  = kwargs.get('Vi')
        
        self.m_dot = kwargs.get('m_dot') # Stays constant through component
        self.mdot_ratio = 1 # used to track mass flow ratio through sections
        self.ni   = kwargs.get('ni', 1) # Isentropic efficiency
        self.BPR = kwargs.get('BPR',1)
        
        self.Toe = kwargs.get('Toe')
        self.Poe = kwargs.get('Poe')
        self.Te  = kwargs.get('Te')
        self.Pe  = kwargs.get('Pe')
        self.Me  = kwargs.get('Me')
        self.Ve  = kwargs.get('Ve')
        
        self.StageName = ""
        self.Power = None
        
    def forward(self, next_Stage):
        next_Stage.Toi = self.Toe
        next_Stage.Poi = self.Poe
        next_Stage.Ti  = self.Te
        next_Stage.Pi  = self.Pe
        next_Stage.Mi  = self.Me
        next_Stage.Vi  = self.Ve
        next_Stage.m_dot = self.m_dot
        next_Stage.mdot_ratio = self.mdot_ratio
    
    def printOutputs(self):
        form ='{:9.3f}'
        print('Stage: ', self.StageName)
        if self.Toe != None:
            print('\t Toe = {} K'.format(form).format(self.Toe))
        if self.Poe != None:
            print('\t Poe = {} Pa'.format(form).format(self.Poe))
        if self.Te != None:
            print('\t Te  = {} K'.format(form).format(self.Te))
        if self.Pe != None:
            print('\t Pe  = {} Pa'.format(form).format(self.Pe))
        if self.m_dot != None:
            print('\tmdot = {} kg/s'.format(form).format(self.m_dot))
        if self.Me != None:
            print('\t Me  = {}'.format(form).format(self.Me))
        if self.Ve != None:
            print('\t Ve  = {} m/s'.format(form).format(self.Ve))
        if self.Power != None:
            print('\t Pow = {} W'.format(form).format(self.Power))
        if self.SpecPower != None:
            print('\t Specific Pow = {} J/kg'.format(form).format(self.specPower))
        self.extraOutputs()
    
    def extraOutputs(self):
        # Overwrite this and put any extra outputs here within individual stages
        return None
        
class Intake(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Intake"
        # NOTE: Ram efficiency ~= Isentropic Efficiency
        
    def calculate(self):
        # Always assume Pi/Pa and Ti/Ta are given (atmos conditions)
        self.Pi = self.Pa
        self.Ti = self.Ta
        self.Vi = self.Vinf
        self.Mi = self.Minf
        # If no vel or mach num inputted, assume stationary
        if self.Mi == None:
            if self.Vi == None:
                self.Mi = 0
            else:
                self.Mi = self.Vi/np.sqrt(self.gam_a*self.R*self.Ti)
        else:
            if self.Vi == None:
                self.Vi = self.Mi*np.sqrt(self.gam_a*self.R*self.Ti)
        
        # Now we should have mach num no matter what
        # and the static props (atm props)
        self.Toe = self.Ti * (1 + (self.gam_a-1)*(self.Mi**2)/2)
        self.Poe = self.Pi * (1 + self.ni*(self.Mi**2)*(self.gam_a-1)/2)**(self.gam_a/(self.gam_a-1))
        
class Compressor(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Compressor"
        # Adding PR and BPR
        self.r = kwargs.get('rc') # Pressure Ratio of stage
        self.BPR = kwargs.get('BPR', 1) # Bypass Ratio: total mass flow (air)/mass flow through core
        self.np = kwargs.get('np') # Polytropic efficiency
        self.mdot_ratio = 1 # Starts as 1 for fan, will be updated by prior comp if
                            # different from 1 from forward section
    def calculate(self):
        # Should always have input To and Po, need to calculate power
        # and output To and Po. r will always be given, BPR will affect output 
        # to next stage
        if self.r == None:
            raise EngineErrors.MissingValue('R-Press. Ratio','Compressor')
        elif self.np == None:
            self.np = ((self.gam_a-1)/self.gam_a)*np.log(self.r) / \
                        np.log( (self.r**((self.gam_a-1)/self.gam_a) - 1)/self.ni + 1)
        
        n_frac =  (self.gam_a-1)/(self.gam_a*self.np)
        self.Toe = self.Toi + self.Toi*(self.r**n_frac - 1)
        self.Poe = self.r*self.Poi
        
        if self.m_dot == None:
            self.specPower = self.mdot_ratio*self.cp_a*(self.Toe-self.Toi)
        else:
            self.Power = self.m_dot*self.cp_a*(self.Toe-self.Toi)
        # Done
        
    
    def forward(self, next_Stage_hot, next_Stage_cold=None):
        next_Stage_hot.Toi = self.Toe
        next_Stage_hot.Poi = self.Poe
        next_Stage_hot.Ti  = self.Te
        next_Stage_hot.Pi  = self.Pe
        next_Stage_hot.Mi  = self.Me
        next_Stage_hot.Vi  = self.Ve
        next_Stage_hot.mdot_ratio = self.mdot_ratio
        
        if next_Stage_cold == None:
            next_Stage_hot.m_dot = self.m_dot
        else:
            if self.BPR == None:
                raise EngineErrors.MissingValue('BPR','Compressor')
            else:
                if self.m_dot != None:
                    m_dot_h = self.m_dot/(1 + self.BPR)
                    m_dot_c = self.m_dot - m_dot_h
                    
                    next_Stage_hot.m_dot = m_dot_h
                    next_Stage_cold.m_dot = m_dot_c
                else:
                    # No inputted mdot
                    mdot_ratio_h = 1/(1+self.BPR)
                    
                    next_Stage_hot.m_dot = None
                    next_Stage_hot.mdot_ratio = mdot_ratio_h
                    next_Stage_cold.m_dot = None
                    # Dont need to send mdot ratio to cold section
                    
                next_Stage_cold.Toi = self.Toe
                next_Stage_cold.Poi = self.Poe
                next_Stage_cold.Ti  = self.Te
                next_Stage_cold.Pi  = self.Pe
                next_Stage_cold.Mi  = self.Me
                next_Stage_cold.Vi  = self.Ve
               
                
        
        
    def calculate_nc(self, np, gamma=1.4):
        '''

        Parameters
        ----------
        np : Float, 0-1
            Polytropic efficiency.
        gamma : Float, optional
            Gamma. The default is 1.4.

        Returns
        -------
        Isentropic Efficiency.

        '''
        nc = ( self.r**((self.gam_a-1)/self.gam_a) - 1 ) / ( self.r**((self.gam_a-1)/(self.gam_a*np)) - 1 )
        return nc
        
        
class Combustor(Stage):
    def __init__(self, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Combustor"
        self.dTo = kwargs.get('dTb')
        self.dPo = kwargs.get('dPb_dec', 0) # the pressure loss within the compressor as a decimal (0.05 = 5% loss)
        self.f  = kwargs.get('f')
        self.Q  = kwargs.get('Q_fuel')
        self.nb = kwargs.get('nb', 1) # Combustor efficiency
        
    def calculate(self):
        # Assuming we have the Ti and Pi from compressor/prev stage
        # We need to have the exit 
        if self.Toe == None: 
            # No Turbine inlet temp given
            if self.dTo == None: 
                # No combustor increase temp given
                if self.f == None and self.Q == None:
                    # No air-fuel ratio given, cant calculate temps
                    raise EngineErrors.MissingValue('Toe, dTo, or f&Q','Combustor')
                else: 
                    # We have f and Q to calculate exit temp
                    f_ideal = self.f*self.nb # inputted f would be actual
                    self.Toe = (f_ideal*self.Q + self.cp_a*self.Toi)/(self.cpg(1+f_ideal))
            else:
                # We dont have exit temp, but do have temp increase
                self.Toe = self.Toi + self.dTo
         # else: Dont need to use since we have what we need
             # We have turbine inlet temp (Te)
             
        self.Poe = self.Poi*(1-self.dPo)
        self.dTo = self.Toe - self.Toi # will use later for f calcs
        
        if self.f == None:
            if self.Q != None: 
                # Assuming non-ideal, will calculate f and use in mass fuel flow
                self.f = (self.cp_g*self.Toe - self.cp_a*self.Toi) / (self.ni*(self.Q - self.cp_g*self.Toe))
                
        if self.m_dot != None and self.f != None:
            self.m_dot += self.f*self.m_dot
        elif self.m_dot == None and self.f != None:
            self.mdot_ratio = (1+self.f)/(self.BPR + 1)
            
class Turbine(Stage):
    def __init__(self, Comp_to_power, **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Turbine"
        self.np = kwargs.get('np') # Polytropic efficiency
        self.nm = kwargs.get('nm',1)
        self.Compressor = Comp_to_power # Could be list
        # Will have inlet temp, compressor power
        self.r  = kwargs.get('rt') # Add for later, not used now
        # this will be for generators or when turbine pressure ratio is specified
        
    def calculate(self):
        if self.m_dot != None:
            if type(self.Compressor) == list:
                com_power = 0
                for i in range(0,len(self.Compressor)):
                    com_power += self.Compressor[i].Power
                self.Power = com_power/self.nm 
            else:
                self.Power = self.Compressor.Power/self.nm 
            # Calculate exit temp
            self.Toe = self.Toi - self.Power/(self.m_dot*self.cp_g)
        else:
            # No m_dot is given, need to power balance based on 
            # BPR ratios instead
            if type(self.Compressor) == list:
                com_power = 0
                for i in range(0,len(self.Compressor)):
                    com_power += self.Compressor[i].Power_ratio
                self.specPower = com_power/self.nm 
            else:
                self.specPower = self.Compressor.Power_ratio/self.nm 
            # Calculate exit temp
            self.Toe = self.Toi - self.specPower/(self.mdot_ratio*self.cp_g)
        
            
        if self.np == None:
            if self.r != None:
                # Calculate np
                self.np = np.log(1- self.ni*(1 - self.r**((self.gam_g-1)/self.gam_g)))
                self.np /= np.log(self.r)*(self.gam_g-1)/self.gam_g
            else:
                print('Warning: insufficient parameters given to turbine')
                print('Continuing assuming polytropic efficiency = 1')
                self.np = 1
                
        m_frac = self.np*(self.gam_g-1)/self.gam_g
        self.Poe = self.Poi*(1- (self.Toi-self.Toe)/self.Toi )**(1/m_frac)
        
class Nozzle(Stage):
    def __init__(self, air_type='hot', **kwargs):
        Stage.__init__(self, **kwargs)
        self.StageName = "Nozzle"
        if air_type == 'hot':
            self.gam = self.gam_g
        else:
            self.gam = self.gam_a
        
    def calculate(self):
        # Check if choked
        Tc = self.Toi*(2/(self.gam_g+1))
        Pc = self.Poi*(1 - (1/self.ni)*(1-Tc/self.Toi))**(self.gam/(self.gam-1))
        
        P_rat = self.Poi/self.Pa
        P_crit = self.Poi/Pc
        if P_rat > P_crit:
            # Nozzle is choked
            if self.Pe == None:
                self.Pe = Pc
            else:
                # We are given exit pressure. For now assuming this
                # is because Pe = Pa for a fully expanded CD Nozzle
                self.Te = self.Toi*(1-self.ni*(1-(self.Pe/self.Poi)**((self.gam-1)/self.gam)))
        else:
            # Nozzle is not choked
            self.Pe = self.Pa
        
        self.Te = self.Toi*(1-self.ni*(1-(self.Pe/self.Poi)**((self.gam-1)/self.gam)))
        self.Me = np.sqrt((2/(self.gam-1))*(self.Toi/self.Te - 1))
        self.Ve = self.Me*np.sqrt(self.gam*self.R*self.Te)
        
        # Stag props at exit
        self.Toe = self.Toi
        self.Poe = self.Pe * (1 + (self.gam -1)*(self.Me**2)/2)**(self.gam/(self.gam -1))
       
            
class Engine():
    def __init__(self):
        print(None)
        # Will use to define and connect all of the stages so the 
        # outlets of one stage is the inputs for the next stages
        
        
class Turbofan_SingleSpool():
    def __init__(self, **kwargs):
        '''
        A signgle spool turbofan which has one turbine to power the
        fan and compressor. It has a cold-air bypass which is after the fan
        and goes straight to a nozzle. The core-flow goes to a second compressor
        and then to the combustor, followed by a single turbine and lastly the
        core nozzle. 
        Parameters
        ----------
        **kwargs : Dictionary
            Contains all needed and optional parameters with the keys listed below.
            Required:
            'Ta': Atmospheric static temperature
            'Pa': Atmospheric static pressure
            'rfan': Fan Pressure Ratio
            'rc':   Compressor Pressure Ratio
            'BPR':  Bypass Ratio (If not passed, assumed 1 so no bypass occurs)
            'T_turb_in': Turbine inlet temp
            
            Optional
            'Vinf': None, # Or Minf
            'Minf': 0.85, # Or Vinf, if none its assumed stationary
            'mdot_a': # Mass flow rate of air into engine (kg/s)
            'Q_fuel': Heat energy of fuel, pass if real to calculate f and include fuel flow in power/velecity calcs
            'F': Thrust of the engine produced
            Efficiencies (Assumed to be 1 if not passed)
            'ni': Inlet Isentropic Efficiency
            'nj': Nozzle Isentropic Efficiency
            'nf': Fan Isentropic Efficiency
            'nc': Compressor Isentropic Efficiency
            'nt': Turbine Isentropic Efficiency
            'nb': Cobustor Efficincy
            'nm': Mechanical Efficiency
            'npf': Fan Polytropic Efficiency (overrides isentropic)
            'npc': Compressor Polytropic Efficiency (overrides isentropic)
            'npt': Turbine Polytropic Efficiency (overrides isentropic)
            'dP_combustor': Decimal pressure drop in combustor (ex: 0.05 for 5% P loss, 0 for ideal)
            
        Returns
        -------
        None.

        '''
        # Stages
        # Atm moving
        # Inlet
        # Fan  (is a compressor)
        # Bypass Nzzle
        # LP Compressor
        # HP Compressor
        # Combustor
        # HP Turbine
        # LP Turbine
        # Nozzle
        self.inputs = kwargs.copy()
        gen_kwargs = {
            'Ta': kwargs.get('Ta'),
            'Pa': kwargs.get('Pa'),
            'Vinf': kwargs.get('Vinf'),
            'Minf': kwargs.get('Minf'),
            'BPR': kwargs.get('BPR',1), # Need this here on case mdot=None
            'Q_fuel':  kwargs.get('Q_fuel')}# kJ/kg
        # Efficiencies
        ni = kwargs.get('ni',1) # Inlet
        nj = kwargs.get('nj',1) # Nozzle
        nf = kwargs.get('nf',1) # Compressor - Isentropic
        nc = kwargs.get('nc',1) # Compressor - Isentropic
        nt = kwargs.get('nt',1) # Turbine - Isentropic
        nb = kwargs.get('nb',1) # Cobustor
        nm = kwargs.get('nm',1) # Mechanical
        npf = kwargs.get('npf') # Fan - Polytropic
        npc = kwargs.get('npc') # Compressor - Polytropic
        npt = kwargs.get('npt') # Turbine - Polytropic
        # Pressure Ratios/Relations
        dP_b = kwargs.get('dP_combustor') # Decimal pressure drop in combustor
        rfan = kwargs.get('rfan') # Fan PR
        rc   = kwargs.get('rc')   # Compressor PR
        # Turbine Inlet
        To_ti = kwargs.get('T_turb_in') # K - Turbine inlet temp
        # Air Mass flow
        mdot = kwargs.get('mdot_a') # kg/s
        
        self.F = kwargs.get('F')
        
        # Define each stage and pass in parameters
        self.inlet = Intake(**gen_kwargs,ni=ni,m_dot=mdot)
        self.fan = Compressor(**gen_kwargs, rc=rfan, np=npf, ni=nf)
        self.BP_nozzle = Nozzle('cold',**gen_kwargs, ni=nj)
        self.HP_comp = Compressor(**gen_kwargs, rc=rc, np=npc, ni=nc)
        self.combustor = Combustor(**gen_kwargs, Toe=To_ti, dPb_dec=dP_b, ni=nb)
        self.HP_turb = Turbine([self.fan, self.HP_comp], **gen_kwargs, nm=nm, ni=nt, np=npt)
        self.nozzle = Nozzle(**gen_kwargs, ni=nj) # Nozzle/Exhaust?
        
        # Set names for easier readout checks
        self.fan.StageName = 'Fan'
        self.BP_nozzle.StageName = 'Cold Nozzle'
        self.nozzle.StageName = 'Hot Nozzle'
        
        # Define all stages in engine to iterate through
        # Two dimensional since there is a bypass, ie one stage
        # passes params to two different stages
        self.AllStages = [[self.inlet, None ],
                          [self.fan, None], 
                          [self.HP_comp,  self.BP_nozzle],
                          [self.combustor,None],
                          [self.HP_turb, None],
                          [self.nozzle, None]]
        
    def calculate(self, printVals=True):
        '''
        Calculates the properties of the air throughout the engine.

        Parameters
        ----------
        printVals : Bool, optional
            When true, it will print out each components exit conditions. The default is True.

        Returns
        -------
        None.

        '''
        for i in range(0,len(self.AllStages)):
            # Calculate each row and print outputs
            self.AllStages[i][0].calculate()
            if printVals: self.AllStages[i][0].printOutputs()
            # Check if current stage has a parallel (ie, prev stage passes air to 2 stages)
            if self.AllStages[i][1] != None:
                self.AllStages[i][1].calculate()
                if printVals: self.AllStages[i][1].printOutputs()
                
            # Move forward/propogate
            if i != len(self.AllStages)-1: # It is not at the end, so forward
                if self.AllStages[i+1][1] != None: 
                    # Means that this stage delivers to two stages: fan -> HPC & BP Noz
                    self.AllStages[i][0].forward(self.AllStages[i+1][0],self.AllStages[i+1][1])
                else:
                    # Stage delivers to one stage
                    self.AllStages[i][0].forward(self.AllStages[i+1][0])
                    
    def getOutputs(self):
        '''
        Returns a dictionary containing important values fromwithin the engine
        using typical engine notaion. Must be run after calculate in order to contain values.

        Returns
        -------
        outs : Dictionary
            Conatins the following items:
                'mdot_c' - bypass flow [kg/s]
                'mdot_h1' - core flow  [kg/s]
                'mdot_h2' - core flow + fuel [kg/s]
                'mdot' - air flow into the engine [kg/s]
                'Ca'  - Initial vel of air (rel to engine) [m/s]
                'C9'  - Exhast vel of core    [m/s]
                'C19' - Exhuast vel of bypass [m/s]
                'Pa'   - Static atmospheric pressure [Pa]
                'P9'   - Exit pressure of core noz   [Pa]
                'P19'  - Exit pressure of bypass noz [Pa]
                'To3': - Combustor inlet stagnation temp [K]
                'To4': - Combustor outlet stagnation temp [K]
                'nb':  - combustor efficiency
                'dH':  - fuel energy  [kJ/kg]
                'f':   - air to fuel ratio

        '''
        outs = {
            'mdot_c': self.BP_nozzle.m_dot, # bypass flow
            'mdot_h1': self.HP_comp.m_dot, # core flow
            'mdot_h2': self.nozzle.m_dot, # core flow + fuel
            'mdot': self.inlet.m_dot, # air flow in
            'Ca': self.inlet.Vi, # Initial vel of air
            'C9': self.nozzle.Ve, # Exhast vel of core
            'C19': self.BP_nozzle.Ve, # Exhuast vel of bypass
            'Pa': self.inputs.get('Pa'),
            'P9': self.nozzle.Pe,    # Exit pressure of core noz
            'P19': self.BP_nozzle.Pe,# Exit pressure of bypass noz
            'To3': self.combustor.Toi, # combustor inlet temp
            'To4': self.combustor.Toe, # combustor outlet temp
            'nb': self.combustor.ni, # combustor efficiency
            'dH': self.combustor.Q,  # fuel energy
            'f': self.combustor.f    # air to fuel ratio
            }
        return outs
    
    def printInputs(self):
        '''
        Prints out all of the kwargs entered on intialization

        Returns
        -------
        None.

        '''
        print('Inputs')
        for key,val in self.inputs.items():
            print('\t {}  =  {}'.format(key,val))
            
    def calculatePerformanceParams(self):
        outs = self.getOutputs()
        
        # Calculate thrust or mdot if one is given and not the other
        if self.F == None:
            if outs['mdot'] != None:
                # Calcualte thrust
                self.F = EPF.Thrust_1(outs['mdot'], self.inputs['BPR'], outs['C9'], outs['C19'], outs['Ca'])
        else:
            if outs['mdot'] == None:
                # calculate mdot
                mdot = EPF.mdot_2(self.F, self.inputs('BPR'), outs['C9'], outs['C19'], outs['Ca'])