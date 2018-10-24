#-----------------------------------------------------
# Bioinformatic programming challenges
# Assignment1: HYBRID CROSS Object
# author: Lucía Prieto Santamaría
#-----------------------------------------------------


require './SeedStock.rb' #Import of SeedStock object



class HybridCross
  
  attr_accessor :parent1
  attr_accessor :parent2
  attr_accessor :f2_wild
  attr_accessor :f2_p1
  attr_accessor :f2_p2
  attr_accessor :f2_p1p2
  
  @@total_hcross_objects = Array.new # Class variable that will save in an array all the instances of HybridCross created 
  
  
  
  
#-----------------------------------------------------  
  def initialize (params = {})

    @parent1 = params.fetch(:parent1, "Instance of SeedStock object not found") # This instance variable will save the SeedStock object which ID corresponds to 'Parent1'in the cross file
    @parent2 = params.fetch(:parent2, "Instance of SeedStock object not found") # This instance variable will save the SeedStock object which ID corresponds to 'Parent2'in the cross file
    @f2_wild = params.fetch(:f2_wild, "0")
    @f2_p1 = params.fetch(:f2_p1, "0")
    @f2_p2 = params.fetch(:f2_p2, "0")
    @f2_p1p2 = params.fetch(:f2_p1p2, "0")
    
    @@total_hcross_objects << self # Everytime an HybridCross object is initialized we add it to the array that contains all the instances of this object
  
  end
#-----------------------------------------------------  
  
  
  
#-----------------------------------------------------  
  def self.all_hybridcross
    # Class method to get the array with all the instances of HybridCross object
    
    return @@total_hcross_objects
  
  end
#-----------------------------------------------------

  
  
#-----------------------------------------------------  
  def self.load_from_file(filename)
    # Class method to create instances of HybridCross object based on an input file
    
    unless File.exists? (filename) # We check that the filename stands for an existing file
      abort "Error: The file #{filename} does not exist" # If not, we stop the program
    end
    
    fhc = File.open(filename, "r")
    fhc.readline # We discard the first line (it has the header)
    
    fhc.each_line do |line|
      p1_id, p2_id, f2_w, f2_p1, f2_p2, f2_p1p2 = line.split("\t") # We split the line in each of the columns

      # We create the instance with the given properties in each column
            HybridCross.new(
              :parent1 => SeedStock.get_seed_stock(p1_id), # We call the SeedStock method that retrieves the instance by giving the seed stock ID
              :parent2 => SeedStock.get_seed_stock(p2_id), # We call the SeedStock method that retrieves the instance by giving the seed stock ID
              :f2_wild => f2_w.to_i,
              :f2_p1 => f2_p1.to_i,
              :f2_p2 => f2_p2.to_i,
              :f2_p1p2 => f2_p1p2.to_i
              )
  
    end
    
    fhc.close
  
  end
#-----------------------------------------------------



#-----------------------------------------------------
  def self.linked (cross_obj)
    # Class method that determines by a Chi Square function whether two genes are linked or not
    
    total = cross_obj.f2_wild + cross_obj.f2_p1 + cross_obj.f2_p2 + cross_obj.f2_p1p2
          
    e1 = ((total * 9) / 16).to_f
    e2 = ((total * 3) / 16).to_f
    e3 = ((total * 3) / 16).to_f
    e4 = ((total) / 16).to_f
                   
          chi_square = (((cross_obj.f2_wild - e1) ** 2) / e1 + ((cross_obj.f2_p1 - e2) ** 2) / e2 + ((cross_obj.f2_p2 - e3) ** 2) / e3 + ((cross_obj.f2_p1p2 - e4) ** 2) / e4).to_f
          
          if (chi_square >= 3.84) # This is the value that tell us if the genes are linked or not
            
            puts "Recording: #{cross_obj.parent1.mutantgene.name} is genetically linked to #{cross_obj.parent2.mutantgene.name} with chisquare score #{chi_square}"
            
            cross_obj.parent1.mutantgene.linked = cross_obj.parent2.mutantgene # We change the 'linked' property of the instance ('linked' = instance of the gene to which is linked)
            cross_obj.parent2.mutantgene.linked = cross_obj.parent1.mutantgene # Do the same thing in the linked gene

          
          end
        
  end
#-----------------------------------------------------



end