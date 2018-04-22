# 
# Keywords:
# 
opti conp rfo molmec phonon comp full noden prop  
# 
# Options:
# 
name model_0                                                                    
cell
   8.950759   5.452717   7.153444  90.000000  90.000000  90.000000
fractional  1   
Ba    core 0.1868657 0.2500000 0.1588604 2.00000000 1.00000 0.00000             
S     core 0.4441955 0.7500000 0.1908703 1.36000000 1.00000 0.00000             
O     core 0.5956544 0.7500000 0.1157554 -0.8400000 1.00000 0.00000             
O     core 0.3282445 0.7500000 0.0463340 -0.8400000 1.00000 0.00000             
O     core 0.4242330 0.9660968 0.3099143 -0.8400000 1.00000 0.00000             
space
P N M A         
totalenergy          -181.4880304544 eV
species   3
Ba     core    2.000000                  
S      core    1.360000                  
O      core   -0.840000                  
buck     
O     core Ba    core  4204.79370     0.290700  0.000000      0.00 10.00
buck inter     
O     core O     core  103585.020     0.200000  25.98000      0.00 15.00
morse intra bond
O     core S     core 5.0000000     1.2000      1.52029  0.0000
three bond intra regular  
S     core O     core O     core 7.1524     109.47    
print    1
dump model_0.res                                                 
