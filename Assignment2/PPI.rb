#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: PPI Object
# author: Lucía Prieto Santamaría
#-----------------------------------------------------



require 'net/http'
require 'json'
require './Protein.rb'

class PPI
  
  @@ppis = nil
  
  @@init = false
  
#-----------------------------------------------------  
  def initialize (params = {})
    
    @@ppis = Array.new # Array containing all interactions  and the corresponding networks
    # For example: [[P1, P2, Network1], [P1, P3, Network1], [P2, P5, Network2] ...]
    
    @@init = true
    
  end
#-----------------------------------------------------

  
#-----------------------------------------------------  
  def self.all_ppis
    # Class method to get the hash with the whole list
    
    return @@ppis
  
  end
#-----------------------------------------------------


#-----------------------------------------------------
  def self.get_ppis (intact_code)
    # Class method that returns the proteins interating with a given protein, applying some filters
    # It is necessary to provide the protein IntAct accession code 
    
    res = uri_fetch ("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{intact_code}?format=tab25")
    
    if res
         
      body = res.body
      lines = body.split("\n") # Each line defines an interaction, we need look for all the proteins the protein interacts with
      
      ppis = Array.new
      
      lines.each do |line|
        fields = line.split("\t") # We split the information of each interaction in the columns defined by PSI-MI TAB 2.5 Format
        
        p1 = fields[0].sub(/uniprotkb:/, "") # First protein in the interaction
        p2 = fields[1].sub(/uniprotkb:/, "") # Second protein in the interaction
        
        intact_miscore = fields[14].sub(/intact-miscore:/, "").to_f # Confidence value, converted to float
        
        
        #Filtering the interactions
        case
          
          when (fields[6]=~/psi-mi:"MI:0018"/ || fields[6]=~/psi-mi:"MI:0397"/)
            # Filter 1: we remove interactions detected by "two hybrib" and "two hybrid array" (too many false positives)
            next
        
        
          when p1 == p2
            # Filter 2: we remove self interactions
            next
          
          when intact_miscore < 0            
            # Filter 3: we remove iteractions that have a score below the value we want
            next
       
          
          when (fields[9]=~/taxid:3702/ && fields[10]=~/taxid:3702/)         
            # Filter 4: we remove those interactions from proteins which are not from Arabidopsis Thaliana
            # Aclaration. In the gene list, there were genes from Arabidopsis lyrata, but there are not interactions of this specie in IntAct
            

            # To add the interaction in the, we must check the order, so that the first protein is before the second (alfabetically)
            # We do this in order to avoid inserting duplicates ([P1, P2] is the same thing as [P2, P1])
           
            if p1 < p2
              
              if not @@ppis.include?([p1, p2])
                @@ppis << [p1, p2]
                ppis << [p1, p2]
              end
            
            else
              
              if not @@ppis.include?([p2, p1])
                @@ppis << [p2, p1]
                ppis << [p2, p1]
              end
              
            end
          
          
          else
            next
          
        end # End filtering
               
                
      end # End loop lines
      
      @@ppis.sort!
      return ppis.sort!
    
    else
      return nil
      
    end
    
  end
#-----------------------------------------------------


#-----------------------------------------------------
  def self.create_ppis (intact_id, level)
    # Class method that call Protein.create_prot to create Protein objects with the interacting proteins
    
    if not @@init
      PPI.new()
    end
     
    ppis = get_ppis(intact_id) # We ask for the interacting proteins
    level += 1 # We deep another level
    
    if ppis
 
      ppis.each do |id1, id2, gene|
  
        if ((id2 == intact_id) && (not Protein.exists(id1)))
          # We check if the new protein is the first one and if it has not been created, we create the object
          # To get the prot_id (from UniProt) from the IntAct ID we call function "intactcode2protid"
          Protein.create_prot(Protein.intactcode2protid(id1), level, gene_id = nil, intact = id1)
          # Here is the recursivity: Protein.create_prot calls PPI.create_ppis and vice versa
    
        elsif ((id1 == intact_id) && (not Protein.exists(id2)))
          # We check if the new protein is the second one and if it has not been created, we create the object
          # To get the prot_id (from UniProt) from the IntAct ID we call function "intactcode2protid"
          Protein.create_prot(Protein.intactcode2protid(id2), level, gene_id = nil, intact = id2)
          # Here is the recursivity: Protein.create_prot calls PPI.create_ppis and vice versa

        end
        
      end

    end   
    
  end
#-----------------------------------------------------


end