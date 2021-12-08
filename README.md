# G4Damage

### Pipeline for cross mapping DNA damage with regard to G-quadruplex structures.

DNA damages from various context are mapped onto G4 stuctures with varied stability. The damage propensity is plotted and discussed in detail in the publication.

### Workflow to get similar plots in the publication

Change the variables as follow:

1. UV damage plots                                                                                                                                                            
   dmg.names = c('CPD', 'PP')                                                                                                                                            
   dmg.type = "UV"                                                                                                                                                       
   dmg.pattern = c("CT", "TC", "TT")                                                                                                                                     
   strand.sensitive = T                                                                                                                                                  
   combine.plot = F                                                                                                                                                      
2. cisplatin damage                                                                                                                                                      
   dmg.names = "cisplatin"                                                                                                                                               
   dmg.type = "cisplatin"                                                                                                                                                
   dmg.pattern = "GG"                                                                                                                                                    
   strand.sensitive = T                                                                                                                                                  
3. 8oxoG damage                                                                                                                                                          
   dmg.names = "oxoG"                                                                                                                                                    
   dmg.type = "oxoG"                                                                                                                                                     
   dmg.pattern = "G"                                                                                                                                                     
   strand.sensitive = T                                                                                                                                                  
4. breakage damage                                                                                                                                                       
   dmg.names = c("sonication", "enzymatic", "ancient")                                                                                                                   
   dmg.type = "breakage"                                                                                                                                                 
   dmg.pattern = "NN"                                                                                                                                                    
   strand.sensitive = F                                                                                                                                                  
   combine.plot = T
