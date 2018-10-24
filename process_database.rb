#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment1: PROCESS_DATABASE, main script
# author: Lucía Prieto Santamaría
#-----------------------------------------------------------------

# This program process 3 data files:
    #1) seed_stock_data.tsv --> information about seeds in the genebank
    #2) gene_information.tsv --> information about genes
    #3) cross_data.tsv --> information about the crosses made

# The program 2 main tasks are:
    #1) Simulation of planting 7 grams of seeds of each stock, creation of a new file with updated data
    #2) Determination of linked genes 

#------------------------------------------------------------------

puts ""

#Import the diferent modules needed
require './Gene.rb'
require './SeedStock.rb'
require './HybridCross.rb'

# Input arguments
geneinf_file, seedstock_file, cross_file, newstock_file = ARGV


unless ARGV.length == 4 # We check user inputs the filenames correctly
    abort " USAGE: process_database.rb gene_information.tsv seed_stock_data.tsv cross_data.tsv new_stock_file.tsv"
end

#-----------------------------------------------------------------


Gene.load_from_file(geneinf_file) # Load the gene information to Gene objects


SeedStock.load_from_file(seedstock_file) # Load the seed stock DB information to SeedStock objects
SeedStock.write_database(newstock_file) # Plant 7 grams of each and update the DB in a new file


HybridCross.load_from_file(cross_file) # Load the crosses information to HybridCross objects


HybridCross.all_hybridcross.each do |cross|
    HybridCross.linked(cross) # For each HybridCross object we call "linked" method to determine linked genes
end

puts "\nFinal Report:"
puts ""


# To output the linked genes we check the property linked of every Gene object
Gene.all_genes.each do |id, gene_obj|
    if gene_obj.linked # Only if we find that the property 'linked' is not false
        puts "#{gene_obj.name} is linked to #{gene_obj.linked.name}"
    end
end


puts ""








