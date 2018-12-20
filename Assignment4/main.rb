#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment4: MAIN SCRIPT
# author: Lucía Prieto Santamaría
#-----------------------------------------------------------------
# Searching for Orthologues
#-----------------------------------------------------------------

# Input: 2 files that can both have nucleotide or protein sequences
#    - One of the file will contain the sequences of which you want to
#      find orthologues, we will call this file 'search_file'
#    - The other one contains the genome/proteome where the orthologue
#      is going to be searched, we will call it 'target_file'

# Function: BLAST each sequence of the search_file to each sequence of
# the target_file, and viceversa to determine which sequences are
# orthologues -- RECIPROCAL BEST HIT

# Output: 'output_orthologues.txt', that will contain the orthologues
# found, and how many of them there are

#-----------------------------------------------------------------


#-----------------------------------------------------------------
#¡¡¡¡¡!!!!!!!!
# NOTE: genome file is not added to github due to the size
#¡¡¡¡¡!!!!!!!!
#-----------------------------------------------------------------



#Import the different modules needed
require 'bio'
require 'stringio'
require 'io/console'



#-----------------------------------------------------------------
#-----------------------------------------------------------------
# BLAST PARAMETERS AND ORTHOLOGUES SEARCHING
#-----------------------------------------------------------------
#-----------------------------------------------------------------

# References for the parameters used in BLAST:
#     1) https://www.ncbi.nlm.nih.gov/pubmed/18042555
#       "Choosing BLAST options for better detection of orthologs as reciprocal best hits."
#       (Moreno-Hagelsieb G, Latimer K.)
#     2) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4094424/
#       "Quickly Finding Orthologs as Reciprocal Best Hits with BLAT, LAST, and UBLAST: How Much Do We Miss?"
#       (Natalie Ward and Gabriel Moreno-Hagelsieb)
# The chosen threshold for E-VALUE will be 10(-6) and an overlap of 50% between the query and
# its hit will be required.


# Next steps would be:
# 1) Looking at GO terms for the putative orthologue genes: the more terms they
# share, the higher probability to be related
# 2) Using Machine Learning algorithms to infer how some features (the length of
# the sequences, the conserved regions, the protein physicochemical profile…)
# affect the orthologues search
# 3) Algorithms to infer phylogenetic relationships
# 4) Web pages that do orthologues search to contrast the results





#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FUNCTIONS DECLARATION
#-----------------------------------------------------------------
#-----------------------------------------------------------------



#-----------------------------------------------------------------

def open_file(filename)
  # Method that checks whether the given file exists, and if it does,
  # it opens it.
  
  unless File.exists?(filename) # We check the given file exists
    abort "Error: File #{filename} does not exist"
  end
  
  return File.open(filename, "r")
  
end

#-----------------------------------------------------------------

def create_and_open_file(filename)
  # Method that checks whether a given file exists, in which case it deletes it;
  # and opens it.
  
  if File.exists?(filename) 
    File.delete(filename) # We remove the file in case it exits to update it
  end
  
  return (File.open(filename, "w"))
  
end

#-----------------------------------------------------------------

def ask_input_filename(message)
  # Method that given a message, will ask the user to input the name of the
  # file that contains the sequences
  
  puts message
  STDOUT.flush
  file_name = gets.chomp
  
  return file_name
  
end

#-----------------------------------------------------------------

def ask_input_filetype(message)
  # Method that given a message, will ask the user to input which type of
  # sequence the previous file contains
  
  puts message
  
  STDOUT.flush
  filetype = gets.chomp.downcase

  if filetype == 'g' # It contains nucleotides
      type = 'nucl' 
  elsif filetype == 'p' # It contains aminoacids
      type = 'prot'
  else 
      abort "Error: #{filetype} is not correct.\n'G' or 'P' should be typed instead"
  end
  
  return type
  
end




#-----------------------------------------------------------------
#-----------------------------------------------------------------
# MAIN PROGRAM
#-----------------------------------------------------------------
#-----------------------------------------------------------------


# Global constants definition according to the defined parameters
$E_VAL = 10**-6
$OVERLAP = 50



puts
puts "ORTHOLOGUES SERACHING PROGRAM, using BLAST"


# We will first ask the user to input the different names of the files and their types
search_file = ask_input_filename("Input the name of the first file with a genome/proteome:\n")
type_search_file = ask_input_filetype("The given file is a [g] genome or a [p] proteome:\nG/P -> ")
target_file = ask_input_filename("Input the name of the second file with a genome/proteome:")
type_target_file = ask_input_filetype("The given file is a [g] genome or a [p] proteome:\nG/P -> ")



# We create the object of the output file and print its title
out_file = create_and_open_file('output_ortologues.txt')
out_file.puts "ORTHOLOGUES FOUND in files #{search_file} and #{target_file}\n"


# We create the directory where we will store the databases and obtain the databases names 
system("mkdir Databases")


# We format the names that the databases will have
db_search = search_file.to_s + '_db'
db_target = target_file.to_s + '_db'


# We create the databases
system("makeblastdb -in '#{search_file}' -dbtype #{type_search_file} -out ./Databases/#{db_search.to_s}") 
system("makeblastdb -in '#{target_file}' -dbtype #{type_target_file} -out ./Databases/#{db_target.to_s}")



# We create factories depending on the types of the files
if type_search_file == 'nucl' and type_target_file == 'nucl' # Both files contain genomes
    factory_search = Bio::Blast.local('blastn', "./Databases/#{db_search.to_s}")
    factory_target = Bio::Blast.local('blastn', "./Databases/#{db_target.to_s}")

elsif type_search_file == 'nucl' and type_target_file == 'prot' # First file contains a genome and the second one a proteome
    factory_search = Bio::Blast.local('tblastn', "./Databases/#{db_search.to_s}")
    factory_target = Bio::Blast.local('blastx', "./Databases/#{db_target.to_s}")
    
elsif type_search_file == 'prot' and type_target_file == 'nucl' # First file contains a proteome and the second one a genome
    factory_search = Bio::Blast.local('blastx', "./Databases/#{db_search.to_s}")
    factory_target = Bio::Blast.local('tblastn', "./Databases/#{db_target.to_s}")

elsif type_search_file == 'prot' and type_target_file == 'p' # Both files contain proteomes
    factory_search = Bio::Blast.local('blastp', "./Databases/#{db_search.to_s}")
    factory_target = Bio::Blast.local('blastp', "./Databases/#{db_target.to_s}")

end


# We create Bio::FastaFormat objects for each of the files
ff_search = Bio::FastaFormat.open(search_file)
ff_target = Bio::FastaFormat.open(target_file)


# We create a hash to store the target_file sequences and make them accesible by its id.
target_hash = Hash.new 
ff_target.each do |seq_target|
  target_hash[(seq_target.entry_id).to_s] = (seq_target.seq).to_s 
end


# We need a counter to know how many orthologues have been found
count = 1 


ff_search.each do |seq_search| # We iterate over each sequence in the search_file
  
  search_id = (seq_search.entry_id).to_s # We store the ID in search_file to later know if it is a reciprocal best hit
  report_target = factory_target.query(seq_search)
  
  if report_target.hits[0] # Only if there have been hits continue.
    
    target_id = (report_target.hits[0].definition.match(/(\w+\.\w+)|/)).to_s # We get ID that will correspond to target_file ID
    
    if (report_target.hits[0].evalue <= $E_VAL) and (report_target.hits[0].overlap >= $OVERLAP) # We check the stablished parameters

      report_search = factory_search.query(">#{target_id}\n#{target_hash[target_id]}")
      # We look in the hash with the previous ID to get the sequence and query the factory
      
      if report_search.hits[0] # Again, only continue if there have been hits
        
        match = (report_search.hits[0].definition.match(/(\w+\.\w+)|/)).to_s # We get the ID that will match with the ID in the search_file
        
        if (report_search.hits[0].evalue <= $E_VAL) and (report_search.hits[0].overlap >= $OVERLAP) # Check parameters
          
          if search_id == match # If the match and the search_file ID match, it means that this is a reciprocal best hit
            
            out_file.puts "#{search_id}\t\t#{target_id}" # We write it in the output file
            puts "#{search_id}\t\t#{target_id}"
            
            count += 1
            
          end
        
        end
      
      end
    
    end      
  
  end


end

# Print in the file the number of orthologues
out_file.puts "\n\nNumber of orthologues found: #{count}"


puts "DONE!\n\n"
puts "You can browse the output in the file output_ortologues.txt"

out_file.close
system("rm -r Databases") # We remove the databases
