#-----------------------------------------------------------------
# Bioinformatic programming challenges
# Assignment3: MAIN SCRIPT
# author: Lucía Prieto Santamaría
#-----------------------------------------------------------------
# GFF feature files and visualization
#-----------------------------------------------------------------

# Input: file with a gene list where directed mutagenesis will be
# done

# Function: the program iterates over the exons of the given genes
# searching a target.

# Output: 3 files are generated
#   1) GFF3 with positions refering to genes --> genes.gff3
#   2) GFF3 with positions refering to chromosomes --> chromosomes.gff3
#   3) Genes without target in exons --> genes_without_target.txt

#-----------------------------------------------------------------


#Import the different modules needed
require 'net/http'
require 'bio'


#-----------------------------------------------------------------
#-----------------------------------------------------------------
# FUNCTIONS DECLARATION
#-----------------------------------------------------------------
#-----------------------------------------------------------------



 

def uri_fetch(uri_str)
    # # Function that will be called to get data from an URI  

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

def load_from_file(filename)
  # Method to read the genes from the file
     
  fg = File.open(filename, "r")
  genes = Array.new
  
  fg.each_line do |line|
    genes << line.delete("\n") # We remove end of line \n 
  end
  
  fg.close
  return genes.uniq # We apply uniq preventing posible duplications user may have introduced

end

#-----------------------------------------------------------------

def obtain_gene_info (gene_id)
  # Routine that creates a Bio:EMBL object of the given gene
  
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id}"
  response = uri_fetch(address)
  
  return nil unless response
  record = response.body
  entry = Bio::EMBL.new(record) 
  
  bioseq = entry.to_biosequence
  
  return bioseq
  
end

#-----------------------------------------------------------------

def find_target_in_exon (exon_id, target_sequence_matches, len_seq, exon_position, strand)
   # Method that determines if any of the target's matches are included in a given exon
   # The methods takes into acount if the exon is in the strand + or - (in case is -, the
   # matches that are passed to the function are the ones in the reverse strand)
   
   target_in_exon = Hash.new
   # Key: Array containing matched positions in exon [ini, end]
   # Values: [ID of the exon, strand]

   case strand # We will check if we are working will the foward or reverse strand
   
   when '+' # Foward
     
     target_sequence_matches.each do |match_init|
        
        match_end = match_init + $len_target - 1 
        
        if (match_init >= exon_position[0].to_i) && (match_init <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
           # The condition is established to see whether the target is inside the exon
           target_in_exon[[match_init, match_end]] = [exon_id, '+']
        end
        
     end
   
   
   when '-' # Reverse
      
      target_sequence_matches.each do |match_init|
        
        match_end = match_init + $len_target - 1
        
        if (match_init >= exon_position[0].to_i) && (match_init <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
           # The condition is established to see whether the target is inside the exon
           
           # To work will the hipotetical positions that correspond to the foward strand, we need to convert the positions as follows
           m_end = len_seq - match_end 
           m_init = len_seq - match_init
           target_in_exon[[m_end, m_init]] = [exon_id, '-']
        
        end
        
     end 
   
   
   end
   
   if not target_in_exon.empty? # We check there are targets inside the exon
      return target_in_exon
   end
   
end

#-----------------------------------------------------------------

def get_exons_targets (bio_seq_object)
  # Routine that given the Bio:EMBL object returns a hash in which the keys are
  # the coordinates of the target's matches inside exons.
  
  len_bio_seq = bio_seq_object.length() # Length of the nucleotide sequence
  target_positions_in_exon = Hash.new # Hash that will contain the positions targeted inside exons as keys and the strand as values
  
  # We get the target's matches in both foward and reverse strand 
  target_matches_in_seq_foward = bio_seq_object.gsub(/#{$target}/).map{Regexp.last_match.begin(0)}
  target_matches_in_seq_reverse = bio_seq_object.reverse.tr('atgc', 'tacg').gsub(/#{$target}/).map{Regexp.last_match.begin(0)}
  
  bio_seq_object.features.each do |feature| 
    
    position = feature.position
    
    next unless (feature.feature == 'exon' && (not position =~ /[A-Z]/))
    # We look for the feature type "exon" and we ommit tras-splicing
    
    exon_id = feature.qualifiers[0].value.gsub('exon_id=', '') # We format the string

    if position =~ /complement/ # Exon is in reverse strand ---> (-)
        
        position = position.tr('complement()', '').split('..')
        position_reverse = []
        
        # Getting a 2 elements array containg initial and end position, we convert it to the reverse strand
        
        position.each do |pos|
           position_reverse.insert(0, len_bio_seq - pos.to_i) # We use insert to give the correct order of the coordinates
        end
        
        target_pos_in_exon = find_target_in_exon(exon_id, target_matches_in_seq_reverse, len_bio_seq, position_reverse, '-')
        # We call "find_target_in_exon" to determine which matches are inside of the exon.
        # Here, we pass to the function the matches and the positions of the exon both in the reverse strand
        if not target_pos_in_exon.nil? # If we retrieve a response, we add the targets to the hash
         target_positions_in_exon = target_positions_in_exon.merge(target_pos_in_exon)
        end
    
    
    else # Exon is in foward strand ---> (+)
        
        position = position.split('..') # Getting a 2 elements array containg initial and end position
        
        target_pos_in_exon= find_target_in_exon(exon_id, target_matches_in_seq_foward, len_bio_seq, position, '+')
        # We call "find_target_in_exon" to determine which matches are inside of the exon.
        # Here, we pass to the function the matches and the positions of the exon both in the foward strand
        if not target_pos_in_exon.nil? # If we retrieve a response, we add the targets to the hash
         target_positions_in_exon = target_positions_in_exon.merge(target_pos_in_exon)
        end
    
    end
       

  end
  
  return target_positions_in_exon
  # We return the hash
  
end  

#-----------------------------------------------------------------

def add_features(gene_id, targets, bioseq)
  # Method that iterates over the hash with the target's matched in exons
  # to add them as new features to the Bio:EMBL objects.

  exon_features = Array.new

  targets.each do |target, exonid_strand|
     
     feat = Bio::Feature.new("#{$target.upcase}_in_exon", "#{target[0]}..#{target[1]}")
     
     feat.append(Bio::Feature::Qualifier.new('nucleotide_motif', "#{$target.upcase}_in_#{exonid_strand[0]}"))
     # New feature qualifier according to https://www.ebi.ac.uk/ols/ontologies/so/terms/graph?iri=http://purl.obolibrary.org/obo/SO_0000110
     # nucleotide_motif
     # Description: A region of nucleotide sequence corresponding to a known motif.
     # Synonyms: INSDC_note:nucleotide_motif, nucleotide motif, INSDC_feature:misc_feature
     # Short id: SO:0000714 (iri: http://purl.obolibrary.org/obo/SO:0000714)
     # This format will be needed for the GFF3
     
     feat.append(Bio::Feature::Qualifier.new('strand', exonid_strand[1]))
     
     $gff_genes.puts "#{gene_id}\t.\t#{feat.feature}\t#{target[0]}\t#{target[1]}\t.\t#{exonid_strand[1]}\t.\tID=#{exonid_strand[0]}"
     # We print the feature in the GFF3 gene file
     
     exon_features << feat
  end
  
  bioseq.features.concat(exon_features) # We add the new features created to the existing ones

end

#-----------------------------------------------------------------

def get_chromosome (gene_id, bio_seq_object)
  # Routine that given a Bio:Sequence object returns the chromosome
  # and positions to which the sequence belongs
  
  bs_pa = bio_seq_object.primary_accession
  
  return false unless bs_pa
  
  chrom_array = bs_pa.split(":")
  
  $gff_chr.puts "#{chrom_array[2]}\t.\tgene\t#{chrom_array[3]}\t#{chrom_array[4]}\t.\t+\t.\tID=#{gene_id}"
  # This line will print the information of the gene in the GFF, so we can refer to it as the parent
  
  # We return:
  #   - Chromosome number ---> [2]
  #   - Chromosome gene start position ---> [3]
  #   - Chromosome gene end position ---> [4]
  return chrom_array[2], chrom_array[3], chrom_array[4]
  
end

#-----------------------------------------------------------------

def create_open_file(filename)
  # Method that checks whether a given file exists, in which case it deletes it;
  # and opens it.
  
  if File.exists?(filename) 
    File.delete(filename) # We remove the file in case it exits to update it
  end
  
  return File.open(filename, "a+")
  
end

#-----------------------------------------------------------------

def convert_to_chr(gene, targets, chr)
  # Given the gene ID, the hash containing the targets, and the information
  # about the chromosome, this method translates the coordinates to the ones
  # refering to the chromosome. It prints them on the GFF3 chromosome file
  
  
  targets.each do |positions, exon_strand|
    pos_ini_chr = chr[1].to_i + positions[0].to_i
    pos_end_chr = chr[1].to_i + positions[1].to_i
    
    $gff_chr.puts "#{chr[0]}\t.\tnucleotide_motif\t#{pos_ini_chr}\t#{pos_end_chr}\t.\t#{exon_strand[1]}\t.\tID=#{exon_strand[0]};parent=#{gene}"
  end
   
  
end

#-----------------------------------------------------------------





#-----------------------------------------------------------------
#-----------------------------------------------------------------
# MAIN PROGRAM
#-----------------------------------------------------------------
#-----------------------------------------------------------------


# Check the input file

unless ARGV[0] # We check user provides the gene list file
    abort "USAGE: main.rb geneList.txt" 
end

unless File.exists?(ARGV[0]) # We check the given file exists
    abort "Error: File #{ARGV[0]} does not exist"
end

#-----------------------------------------------------------------

# Starting tasks

$target = "cttctt" # This global variable contains the sequence to search in the exons
$len_target = $target.length() # Global variable with the target's length 

puts "Working on the tasks...\n"

# We create global variables with the objects of opening the files correponding to
#   1) GFF3 with positions refering to genes
#   2) GFF3 with positions refering to chromosomes
#   3) Genes without target in exons
$gff_genes = create_open_file("genes.gff3")
$gff_chr = create_open_file("chromosomes.gff3")
$no_targets = create_open_file("genes_without_target.txt")


# We print the headers of each file
$gff_genes.puts "##gff-version 3"
$gff_chr.puts "##gff-version 3"
$no_targets.puts "Genes without #{$target.upcase} in exons\n\n"


genes_ids = load_from_file(ARGV[0]) 

#-----------------------------------------------------------------

# Main part of the program to search the targets in the exons

genes_ids.each do |gene|
  
  seq_obj = obtain_gene_info(gene) # Create the Bio:EMBL object from each gene
  
  unless seq_obj.nil?
    target_hash = get_exons_targets(seq_obj) # We get the targets inside exons of the gene
    # Target hash contains as keys: [initial position of target, end position of target]
    # and as values: [ID of the exon, strand]
    
    if target_hash.empty?
      $no_targets.puts gene
      # If the gene has no targets in exons, we add it to the file genes_without_target.txt
      
    else  
      add_features(gene, target_hash, seq_obj) # We create new features and add them to each seq_obj
      chr = get_chromosome(gene, seq_obj) # We return the chromosome number and postions
      convert_to_chr(gene, target_hash, chr) # We convert the positions to the ones that correspond in the chromosome
    
    end
    
  end
  
end

puts "DONE!\n\n"
puts "You can browse the output in the files: "
puts "\t- genes.gff3"
puts "\t- chromosomes.gff3"
puts "\t- genes_without_target.txt"
