#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: INTERACTION NETWORK Object
# author: Lucía Prieto Santamaría
#-----------------------------------------------------


require './Gene.rb' #Import of Gene object
require './Protein.rb' #Import of Protein object

class InteractionNetwork
  
  attr_accessor :network_id # The number that identifies the network
  attr_accessor :num_nodes # Number of nodes that the networks has
  attr_accessor :members # Array containing the Gene objects from the given file that belong to the network
  
  @@total_network_objects = Hash.new # Class variable that will save in an array all the instances of InteractionNetwork created
  @@number_of_networks = 0
  @@genes_in_networks = Array.new

#-----------------------------------------------------  
  def initialize (params = {})
    
    @network_id = params.fetch(:network_id, "X")
    @num_nodes = params.fetch(:num_nodes, "0")
    @members = params.fetch(:members, Hash.new)
   
    @@total_network_objects[network_id] = self # Everytime a InteractionNetwork object is initialized we add it to the array that contains all the instances of this object
    @@number_of_networks += 1 # Each time we create a new network
    
    
  end
#-----------------------------------------------------


#-----------------------------------------------------
 def self.all_networks
   # Class method to get the hash with all the instances of InteractionNetwork object
  
   return @@total_network_objects
 
 end
#-----------------------------------------------------


#-----------------------------------------------------
  def self.create_network
    # Class method to create new networks
    
    networkid = @@number_of_networks + 1
    
    InteractionNetwork.new(
                          :network_id => networkid,
                          :num_nodes => 2,
                          :members => Hash.new
                          )

    return networkid
    
  end  
#-----------------------------------------------------


#-----------------------------------------------------
  def self.add_nodes2network(networkid)
    # Class method to add nodes to a network
    
    @@total_network_objects[networkid].num_nodes += 1
    
  end  
#-----------------------------------------------------


#-----------------------------------------------------
  def self.add_gene2network(networkid, gene_object)
    # Class method to add genes to a network. It also annotates the Gene objects
    
    @@total_network_objects[networkid].members[gene_object.gene_id] = gene_object
    
    gene_object.annotate

  end  
#-----------------------------------------------------


#-----------------------------------------------------
  def self.assign(protein_object, network_id)
   # Class method to recursively to go through all the branches of the networks
   # We need the complete list of PPIs
   
   protein_object.network = network_id
   
   # We go all over the PPIs
   $PPIS.each do |ppi|
     
     if (ppi[0] == protein_object.intact_id) && (Protein.all_prots_withintact[ppi[1]].network == nil)
       self.add_nodes2network(network_id)
       self.assign(Protein.all_prots_withintact[ppi[1]], network_id) # Call the routine itself
     end
     
     if (ppi[1] == protein_object.intact_id) && (Protein.all_prots_withintact[ppi[0]].network == nil)
       self.add_nodes2network(network_id)
       self.assign(Protein.all_prots_withintact[ppi[0]], network_id) # Call the routine itself
     end
    
   end
   
   if Gene.all_genes[protein_object.prot_id] # If the protein has an associated gene, add it to the network
     self.add_gene2network(network_id, Gene.all_genes[protein_object.prot_id])
   end
     
  end
#-----------------------------------------------------

  




end