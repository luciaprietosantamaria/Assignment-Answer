#----------------------------------------------------- 
# Bioinformatic programming challenges
# Assignment1: SEED STOCK DATABASE Object
# author: Lucía Prieto Santamaría
#----------------------------------------------------- 


require './Gene.rb' #Import of Gene object


class SeedStock
  
  attr_accessor :seedstock_id 
  attr_accessor :mutantgene
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :grams_remaining
  
  @@total_seedstock_objects = Hash.new # Class variable that will save in a hash all the instances of SeedStock created (key: seedstock ID)
  
  
  
#-----------------------------------------------------   
  def initialize (params = {})
    
    @seedstock_id = params.fetch(:seedstock_id, "X000")
    @mutantgene = params.fetch(:mutantgene, "Instance of Gene object not found") # This instance variable will save the Gene object which ID corresponds to 'Mutant_Gene_ID'in the seedstock file
    @last_planted = params.fetch(:last_planted, "00/00/0000")
    @storage = params.fetch(:storage, "cama00")
    @grams_remaining = params.fetch(:grams_remaining, "0")
    
    @@total_seedstock_objects[seedstock_id] = self # Everytime a SeedStock object is initialized we add it to the hash that contains all the instances of this object
    
  end
#-----------------------------------------------------  

  
#-----------------------------------------------------   
  def self.all_stock
    # Class method to get the hash with all the instances of SeedStock object
    
    return @@total_seedstock_objects
  
  end
#----------------------------------------------------- 
  

#-----------------------------------------------------     
  def self.get_seed_stock (ssid)
    # Class method that returns the instance of the SeedStock object identified by a given seedstock ID
    
    if @@total_seedstock_objects.has_key? (ssid) # We first check that the instance corresponding to the given ID has been intialized
      return @@total_seedstock_objects[ssid]
    else
      abort "Error: You have not created an SeedStock object with #{ssid} yet"
    end
    
  end  
#-----------------------------------------------------  


#-----------------------------------------------------  
  def self.load_from_file (filename)
    # Class method to create instances of SeedStock object based on an input file
    
    unless File.exists? (filename) # We check that the filename stands for an existing file
      abort "Error: The file #{filename} does not exist" # If not, we stop the program
    end
    
    fss = File.open(filename, "r")
    $head = fss.readline # We save the file's header in a global variable to reuse it later
  
    
    fss.each_line do |line|
      
      seedstockid, mutgeneid, lastplanted, store, gremain  = line.split("\t")
      # We split the line in each of the columns to have the different properties of the seed stock
     
      # We create instances for each seedstock, giving the different properties in each column
      SeedStock.new(
                      :seedstock_id => seedstockid,
                      :mutantgene => Gene.get_gene(mutgeneid), # We call the Gene method that retrieves the instance by giving the geneID
                      :last_planted => lastplanted,
                      :storage => store,
                      :grams_remaining => gremain.to_i
                      )

    end
    
    fss.close
    
  end
#-----------------------------------------------------  
  

#-----------------------------------------------------
  def plant(grams)
    #Method that extracts a given number of grams from a seed stock
    
    grams = grams.to_i
    
    if grams < @grams_remaining
      @grams_remaining -= grams
    else # There cannot be negative remaining grams
      @grams_remaining = 0 
      puts "WARNING: we have run out of Seed Stock #{@seedstock_id}"
    end
     
  end
#-----------------------------------------------------  
  
  
#-----------------------------------------------------  
  def self.write_database (filename)
    #Class method to output in a new file the update version of the seed stock after planting 7 grams
    
    if File.exists? (filename) 
      File.delete (filename) # We remove the file in case it exits to update it
    end
    
    fndb = File.open(filename, "a+")
    fndb.puts $head.to_s # The first line will be the same header as the previous one
    
    date = Time.now.strftime("%d/%m/%Y") # The current time when we are planting
    
    @@total_seedstock_objects.each do |ssid, ssobject|
      ssobject.plant(7) # We call the method to update the remaining grams of the instances
      ssobject.last_planted = date # We update the property 'last_planted' to the current date
      
      # We print the output in the file (Seed_Stock Mutant_Gene_ID	Last_Planted	Storage	Grams_Remaining)
      fndb.puts "#{ssid}\t#{ssobject.mutantgene.gene_id}\t#{date}\t#{ssobject.storage}\t#{ssobject.grams_remaining}"
    end
    
  end
#-----------------------------------------------------

  
  
end