#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment2: MAIN SCRIPT
# author: Lucía Prieto Santamaría
#-----------------------------------------------------------------

# Input: file with a co-expressed gene list

# Function: the program creates interaction networks between the
# genes, annotating them (Gene, Protein and InteractionNetwork
# Objects are used to achieve this)

# Output: report of which members of the gene list interact with
# one another, together with the KEGG/GO functional annotations of
# those interacting members

#-----------------------------------------------------------------

puts ""

#Import the different modules needed
require './Gene.rb'
require './Protein.rb'
require './InteractionNetwork.rb'
require 'net/http'


#-----------------------------------------------------------------

# Check the input file

unless ARGV[0] && ARGV[1] # We check user provides the gene list file
    abort "USAGE: main.rb geneList.txt output.txt"
end

unless File.exists?(ARGV[0]) # We check the given file exists
    abort "Error: File #{ARGV[0]} does not exist"
end

#-----------------------------------------------------------------

# Function that will be called from some objects to get data from an URI 

def uri_fetch(uri_str)
    # Routine that retrives data from an URI doing some basic error-handling.  

    address = URI(uri_str)  
    response = Net::HTTP.get_response(address)
    
    case response
      when Net::HTTPSuccess then  # when response Object is of type Net::HTTPSuccess
        # successful retrieval of web page
        return response  # return that response object to the main code
      else
        abort "Something went wrong... the call to #{uri_str} failed; type #{response.class}"
    end 
end

#-----------------------------------------------------------------
 
def write_record(networks, filename)
  #Class method to record the report of the members of an InteractionNetwork, together with the KEGG/GO functional annotations of those interacting members
  
  if File.exists?(filename) 
    File.delete(filename) # We remove the file in case it exits to update it
  end
  
  fnr = File.open(filename, "a+")
  fnr.puts "--------------------------------"
  fnr.puts "--------------------------------"
  fnr.puts "Assigment2 RECORD"
  fnr.puts "author: Lucía Prieto Santamaría"
  fnr.puts "--------------------------------"
  fnr.puts "--------------------------------\n\n"
  fnr.puts "For every network created we display:"
  fnr.puts "  1) Network ID"
  fnr.puts "  2) Number of nodes in the network"
  fnr.puts "      (Same as number of proteins - will vary depending on the depth of interactions)"
  fnr.puts "      (Depth of interactions selected: #{$MAX_LEVEL})"
  fnr.puts "  3) Genes form the list that are included in the network"
  fnr.puts "\n\n\n"
  
  fnr.puts "\n-----------------------------------------------\n"
  networks.each do |id_net, network_obj|  
    fnr.puts "NETWORK NUMBER: #{id_net}"
    fnr.puts "\tNodes_number: #{network_obj.num_nodes}"
    fnr.puts "\tGenes associated and annotations:"
    network_obj.members.each do |id, gene|
      fnr.puts "\t\t-- #{id}"
      gene.kegg.each do |kegg_id, kegg_name|
        fnr.puts "\t\t\t\t\t KEGG Pathways ID: #{kegg_id};\tKEGG Pathways Name: #{kegg_name}"
      end
      gene.go.each do |go_id, go_term_name|
        fnr.puts "\t\t\t\t\t GO ID: #{go_id};\tGO Term Name: #{go_term_name}"
      end
    end
    fnr.puts "\n-----------------------------------------------\n"
  end
  
end
  
#-----------------------------------------------------------------

$MAX_LEVEL = 2
# We define in a global constant how deep we want to search in the interacting proteins
# In this case we choose 2, what it means that, given the list of genes, the program will get
#       1) The proteins codified by those genes
#       2) The first proteins that directly interact with proteins in 1)
#       3) The proteins that directly interact with proteins in 2)
# We can set the level of interacting proteins anywhere we want, as we could like to explore
# deeper relationships.



Gene.load_from_file(ARGV[0])
# We create the Gene objects from the genes gived in the file.
# This will lead to the creation of Protein objects that will propagate generating Protein
# objects of the interacting proteins in turn.
# A list with all the interactions will be created.


$PPIS = PPI.all_ppis
# We retrieve this list with all the interacting proteins 


Protein.all_prots_withintact.each do |id, prot_object|
  
  if not prot_object.network
    
    new_network = InteractionNetwork.create_network
    InteractionNetwork.assign(prot_object, new_network)
    # We call the recursive routine to explore the posible branches of the network
    
  end
  
end


write_record(InteractionNetwork.all_networks, ARGV[1])
# Output the result in a file.
  