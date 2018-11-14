#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: GENE Object
# author: Lucía Prieto Santamaría
#-----------------------------------------------------


require 'net/http'
require 'json'
require './Protein.rb'

class Gene
  
  attr_accessor :gene_id 
  attr_accessor :prot_id
  attr_accessor :kegg # Hash that contains the KEGG annotations (key: KEGG Pathway ID, value: KEGG Pathway Name)
  attr_accessor :go # Hash that contains the GO annotations (key: GO ID, value: GO Term name)
  
  @@total_gene_objects = Hash.new # Class variable that will save in a hash all the instances of Gene created (key: protein ID)
  
  

#-----------------------------------------------------  
  def initialize (params = {})
    
    genecodif = params.fetch(:gene_id, "AT0G00000")
    # We have to check whether the given gene ID has the right format or not, we use a regular expression:
    if genecodif =~ /A[Tt]\d[Gg]\d\d\d\d\d/ 
      @gene_id = genecodif
    else # If it does not have the right format we stop the program
      abort "Sorry... The gene ID does not have the right format.\nPlease try again with ATXGXXXXX (X stands for a digit)."
    end
    
    @prot_id = params.fetch(:prot_id, "XXXXXX") #UniProt ID
    
    #Annotations, only if the protein codificated by the gene is a member of an interaction network
    @kegg = params.fetch(:kegg, Hash.new) # If the gene is not KEGG annotaded, the hash will be empty
    @go = params.fetch(:go, Hash.new) # If the gene is not KEGG annotaded, the hash will be empty
    
   
    @@total_gene_objects[prot_id] = self # Everytime a Gene object is initialized we add it to the hash that contains all the instances of this object
    
  end
#-----------------------------------------------------  
  
  
#-----------------------------------------------------  
  def self.all_genes
    # Class method to get the hash with all the instances of Gene object
    
    return @@total_gene_objects
  
  end
#-----------------------------------------------------


#-----------------------------------------------------
  def self.get_prot_id(gene_id)
    # Class method that will return de UniProt ID of a protein codified by a gived gene.
    
    address = URI("http://togows.org/entry/ebi-uniprot/#{gene_id}/entry_id.json")
    response = Net::HTTP.get_response(address)

    data = JSON.parse(response.body)
    
    return data[0]

  end
#-----------------------------------------------------
  
    
#-----------------------------------------------------  
  def self.load_from_file(filename)
    # Class method to create instances of Gene object based on an input file
       
    fg = File.open(filename, "r")
    
    fg.each_line do |line|
      
      line.delete!("\n") # We remove end of line \n 
      
      protid = self.get_prot_id(line)
      
      Gene.new(
            :gene_id => line,
            :prot_id => protid
            )
      
      Protein.create_prot(protid, 0, gene_id = line) # We create Protein objects once we have the Gene object 
  
    end
    
    fg.close
  
  end
#-----------------------------------------------------


#-----------------------------------------------------
  def annotate
    # Instance method that obtains (if present) KEGG ID, KEGG Pathway,GO ID and GO Term, and updates the object with them
    
    addressKEGG = URI("http://togows.org/entry/kegg-genes/ath:#{self.gene_id}/pathways.json")
    addressGO = URI("http://togows.org/entry/ebi-uniprot/#{self.gene_id}/dr.json")
    
    responseKEGG = Net::HTTP.get_response(addressKEGG)
    responseGO = Net::HTTP.get_response(addressGO)

    dataKEGG = JSON.parse(responseKEGG.body)
    dataGO = JSON.parse(responseGO.body)
    
    
    # Annotate with KEGG Pathways
    if dataKEGG[0]
      dataKEGG[0].each do |path_id, path_name|
        self.kegg[path_id] = path_name # We insert a new key (KEGG Pathway ID) with its corresponding value (pathway name)
      end
    end
    
    
    # Annotate with GO
    if dataGO[0]["GO"]
      dataGO[0]["GO"].each do |num|
        if num[1] =~ /^P:/ # We must check the GO refers to a biological proccess (it will start with a 'P')
          self.go[num[0]] = num[1].sub(/P:/, "") # We insert a new key (GO ID) with its corresponding value (GO term name) 
        end 
      end  
    end
    
    
  end  
#-----------------------------------------------------
  


end