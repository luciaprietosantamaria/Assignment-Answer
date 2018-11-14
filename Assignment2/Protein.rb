#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: PROTEIN Object
# author: Lucía Prieto Santamaría
#-----------------------------------------------------


require 'net/http'
require 'json'
require './PPI.rb'

class Protein
  
  attr_accessor :prot_id # UniProt ID
  attr_accessor :intact_id # IntAct ID, if the protein interacts with another protein
  attr_accessor :network # Network ID, if the protein is member of one

  @@total_protein_objects = Hash.new # Class variable that will save in a hash all the instances of Protein created (key: prot_id)
  @@total_protwithintact_objects = Hash.new # Class variable that will save in a hash all the instances of Protein created (key: intact_id)
  
#-----------------------------------------------------  
  def initialize (params = {})
    
    @prot_id = params.fetch(:prot_id, "XXXXXX") 
    @intact_id = params.fetch(:intact_id, nil)
    @network = params.fetch(:network, nil)
    
    @@total_protein_objects[prot_id] = self
    
    if intact_id
      @@total_protwithintact_objects[intact_id] = self
    end  
    
  end
#-----------------------------------------------------  


#----------------------------------------------------- 
  def self.all_prots
    # Class method to get the hash with all the instances of Protein object
    
    return @@total_protein_objects
  
  end
#-----------------------------------------------------


#----------------------------------------------------- 
  def self.all_prots_withintact
    # Class method to get the hash with all the instances of Protein object with key the IntAct ID
    
    return @@total_protwithintact_objects
  
  end
#-----------------------------------------------------  
  

#-----------------------------------------------------
  def self.create_prot (prot_id, level, gene_id = nil, intact = nil)
    # Class method that creates a Protein object, given its protein ID and the current level of interaction depth
    # The gene ID (when we are at level 0) and the IntAct ID (when we are creating Protein object from a PPI) can be given
    

    if not intact # If the IntAct ID is not provided, we look if the protein has it
      intact = self.get_prot_intactcode(gene_id) # Retrieve whether the protein is present in IntAct Database, and which is its accession code
    end

    if intact && (level < $MAX_LEVEL)
      # Once we know the protein has an IntAct ID, and if we are in a lower level than the maximum, we look for the interacting proteins
      PPI.create_ppis(intact, level)
    end
    
    Protein.new(
            :prot_id => prot_id,
            :intact_id => intact,
            :network => nil # We first set all networks relationships to empty
            )
      

    level += 1 # We deep one level
    
  end
#-----------------------------------------------------


#-----------------------------------------------------
  def self.exists(intact_id)
    # Class method that given an IntAct ID returns TRUE if a Protein instance has already been created with it,
    # FALSE if not
    
    if @@total_protwithintact_objects.has_key?(intact_id)
      return true
    
    else
      return false
    
    end
    
  end
#-----------------------------------------------------


#-----------------------------------------------------
  def self.get_prot_intactcode(gene_id)
    # Class method that searchs whether a given protein is present in IntAct database or not
    # If it is, the function will return th IntAct ID
    
    address = URI("http://togows.org/entry/ebi-uniprot/#{gene_id}/dr.json")
    response = Net::HTTP.get_response(address)

    data = JSON.parse(response.body)

    if data[0]['IntAct']
      return data[0]['IntAct'][0][0] # If the protein is present, the result returned will be the protein accession code. If it not, the result will be empty
    else
      return nil
    end


  end
#-----------------------------------------------------


#-----------------------------------------------------
  def self.intactcode2protid(intact_accession_code)
    # Class method to obtain the UniProt ID given a protein IntAct accession code
    
    intact_accession_code.delete!("\n")
    
    if intact_accession_code =~ /[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}/
      # We check the accession UniProt code is correct.
      # Regular expression obtained from "https://www.uniprot.org/help/accession_numbers"
      
      begin
        
        address = URI("http://togows.org/entry/ebi-uniprot/#{intact_accession_code}/entry_id.json")
        response = Net::HTTP.get_response(address)
  
        data = JSON.parse(response.body)

        return data[0]
      
      rescue
        return "UniProt ID not found"      
      
      end
    
    else
      
      return "UniProt ID not found"
    
    end
    
  end  
#-----------------------------------------------------




end