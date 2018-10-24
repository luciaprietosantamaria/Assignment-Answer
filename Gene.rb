#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment1: GENE Object
# author: Lucía Prieto Santamaría
#-----------------------------------------------------


class Gene
  
  attr_accessor :gene_id 
  attr_accessor :name
  attr_accessor :mutant_phenotype
  attr_accessor :linked 
  
  @@total_gene_objects = Hash.new # Class variable that will save in a hash all the instances of Gene created (key: gene ID)
  
  

#-----------------------------------------------------  
  def initialize (params = {})
    
    genecodif = params.fetch(:gene_id, "AT0G00000")
    # We have to check whether the given gene ID has the right format or not, we use a regular expression:
    if genecodif =~ /A[Tt]\d[Gg]\d\d\d\d\d/ 
      @gene_id = genecodif
    else # If it does not have the right format we stop the program
      abort "Sorry... The gene ID does not have the right format.\nPlease try again with ATXGXXXXX (X stands for a digit)."
    end
    
    @name = params.fetch(:name, "xxx")
    @mutant_phenotype = params.fetch(:mutant_phenotype, "Phenotype not available")
    @linked = params.fetch(:linked, false) # In case the gene is linked to other gene, @linked will contain the Gene object corresponding to the linked gene
   
    @@total_gene_objects[gene_id] = self # Everytime a Gene object is initialized we add it to the hash that contains all the instances of this object
    
  end
#-----------------------------------------------------  
  
  
#-----------------------------------------------------  
  def self.all_genes
    # Class method to get the hash with all the instances of Gene object
    
    return @@total_gene_objects
  
  end
#----------------------------------------------------- 
  

#----------------------------------------------------- 
  def self.get_gene (id)
    # Class method that returns the instance of the Gene object identified by a given ID
    
    if @@total_gene_objects.has_key? (id) # We first check that the instance corresponding to the given ID has been intialized
      return @@total_gene_objects[id]
    else
      abort "Error: You have not created a Gene object with #{id} yet"
    end
    
  end  
#-----------------------------------------------------  
  
  
#-----------------------------------------------------  
  def self.load_from_file(filename)
    # Class method to create instances of Gene object based on an input file
    
    unless File.exists? (filename) # We check that the filename stands for an existing file
      abort "Error: The file #{filename} does not exist" # If not, we stop the program
    end
    
    fg = File.open(filename, "r")
    fg.readline # We discard the first line (it has the header)
    
    fg.each_line do |line|
      id, name, phenotype = line.split("\t") # We split the line in each of the columns

      # We create the instance with the given properties in each column
      Gene.new(
            :gene_id => id,
            :name => name,
            :mutant_phenotype => phenotype
            )
  
    end
    
    fg.close
  
  end
#-----------------------------------------------------  
  


end