# 
# Keywords:
# 
opti conp prop compare dist angle phonon inten                                  
# 
# Options:
# 
title
alpha quartz                                                                    
end
cell
   4.940911   4.940911   5.448922  90.000000  90.000000 120.000000
fractional    2
Si    core 0.4647780 0.4647780 0.0000000 2.39999999 1.00000 0.00000             
O     core 0.1553033 0.4268321 0.1248693 -1.1999999 1.00000 0.00000             
space
P 32 2 1        
totalenergy          -175.0161148629 eV
species   2
Si     core    2.400000            
O      core   -1.200000            
buck     
O     core O     core  1388.7730     0.362319 175.00      0.000 10.000
buck     
O     core Si    core  18003.757     0.205205 133.54      0.000 10.000
dump sio2.res                                                    
output movie arc quartz_opt.arc                                              
