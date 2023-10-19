This folder will house the final Liquid Solid Hybrid Model for HADES engine

Primary Notes:
Liquid Model   -   Zillac-Karabeyoglu Model
Solid Model    -   Fuel Regression Marxman-Gilbert relation

FullModel.py is a attempt at ODE45

FUllModelEuler.py is conerting Rhodes Model from English to Metric to prepare for Liquid Blowdown ZK Model
    
chemicalProperties.py is a class that is used in the ZK model.

FullModelEulerZK is a refined model with ZK implemented onto Rhodes Metric model
     Assumptions:
        Peng-Robinson Equation of State 
        Fuel Regression Marxman-Gilbert Relation

What we need to refine thoughts:
    O/F Ratio in GDL program for each P_chmb to determine chamber Temperature and Combustion Molar Mass
    Fuel Regression Models
        a and n coefficients of fuel regression
    CD of injector and CD of nozzle
    Nozzle equations for better thrust modeling
    Make quality of life code changes