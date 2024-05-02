##Script version March 2023.
##Use to AA percentages in protein sequences.
##AA = ACDEFGHIKLMNPQRSTVWY
##Usage: ruby aa_composition.rb <.faa file> > <sample name>

puts "AA percentages in #{(ARGV[0])}\n\n"

##Arguments
#ARGV[0] = input protein sequence fasta file (faa)
#ARGV[1] = sample name (used for output file names)

##Read input protein sequence fasta file
faafile = File.read(ARGV[0])
faa = faafile.scan(/^>(\S+).*\n([^>]*)/)
header,seq = $1,$2

##Create output files
outf_percenteach = File.new("aa_composition_percenteach_#{(ARGV[1]).to_s}.txt", "w")
outf_counteach = File.new("aa_composition_counteach_#{(ARGV[1]).to_s}.txt", "w")
outf_detailed = File.new("aa_composition_detailed_#{(ARGV[1]).to_s}.txt", "w")
outf_overall = File.new("aa_composition_overall_#{(ARGV[1]).to_s}.txt", "w")
outf_percenteach.puts "Sequence\tSequence_length\tAlanine\tCysteine\tAspartic_acid\tGlutamic_acid\tPhenylalanine\tGlycine\tHistidine\tIsoleucine\tLysine\tLeucine\tMethionine\tAsparagine\tProline\tGlutamine\tArginine\tSerine\tThreonine\tValine\tTryptophan\tTyrosine"
outf_counteach.puts "Sequence\tSequence_length\tAlanine\tCysteine\tAspartic_acid\tGlutamic_acid\tPhenylalanine\tGlycine\tHistidine\tIsoleucine\tLysine\tLeucine\tMethionine\tAsparagine\tProline\tGlutamine\tArginine\tSerine\tThreonine\tValine\tTryptophan\tTyrosine"
outf_detailed.puts "Sequence\tSequence_length\tAA_abbrev_short\tAA_name\tAA_abbrev_long\tNumber_AA\tPercent_AA"

##Create empty array, hash and variable
acount = Array.new
aperc = Array.new
h = Hash.new
totprotseqlen = 0

##Make array of AA's
aa = 'ACDEFGHIKLMNPQRSTVWY'
aasplit = aa.split('')
#puts aasplit.length
#aasplit.each { |x| puts x }

##Assign each AA to hash with value 0
aasplit.each { |x| h[x] = 0 }
#puts h

##Assign AA names to each AA abbreviation
aanames = %w{Alanine__Ala Cysteine__Cys Aspartic_acid__Asp Glutamic_acid__Glu Phenylalanine__Phe Glycine__Gly Histidine__His Isoleucine__Ile Lysine__Lys Leucine__Leu Methionine__Met Asparagine__Asn Proline__Pro Glutamine__Gln Arginine__Arg Serine__Ser Threonine__Thr Valine__Val Tryptophan__Trp Tyrosine__Tyr}
haanames = Hash[aasplit.zip aanames]
#puts haanames

##Iterate over each sequence and calculate AA percentages
faa.each do |header,seq|
    puts header
    ##Length of each protein sequence
    protseqlen = ((seq.to_s).count("ACDEFGHIKLMNPQRSTVWY")).to_i
    puts "Seq len: #{protseqlen}"
    ##Sum protein sequence lengths cumulatively
    totprotseqlen = totprotseqlen + ((seq.to_s).count("ACDEFGHIKLMNPQRSTVWY")).to_i
    puts "Total sequence length: #{totprotseqlen}"
    ##Count each AA in array
    acount.clear && aperc.clear
    haanames.each do |x,name| 
        ##Count AA per sequence
        countx = ((seq.to_s).count(x)).to_i
        puts "AA count for #{x}: #{countx}"
        acount << countx
        ##Sum counts of each AA cumulatively 
        h[x] = (h[x].to_i + countx.to_i)
        puts "AA sum for #{x}: #{h[x]}"
        ##Calculate AA percentage per sequence
        percx = ((countx.to_f/protseqlen.to_f)*100)
        puts "AA percentage for #{x}: #{percx}"
        aperc << percx
        ##Write results to output file 'detailed'
        outf_detailed.puts "#{header}\t#{protseqlen}\t#{x}\t#{name.gsub('__',"\t")}\t#{countx}\t#{percx}"
    end    
    ##Write results to output file 'each'
    outf_percenteach.puts "#{header}\t#{protseqlen}\t#{aperc.to_s.gsub(', ',"\t").gsub(/\[|\]/,'')}"
    outf_counteach.puts "#{header}\t#{protseqlen}\t#{acount.to_s.gsub(', ',"\t").gsub(/\[|\]/,'')}"
end

##Calculate percentages of each AA overall, and write overall results to output file 'overall'
h.each do |k,v| 
    p = ((v/totprotseqlen.to_f)*100)
    outf_overall.puts "AA_percentage\t#{k}\t#{haanames[k].gsub('__',"\t")}\t#{p}"
end

##Note: AA names and abbreviations:
##Alanine__Ala__A, Cysteine__Cys__C, Aspartic_acid__Asp__D, Glutamic_acid__Glu__E, Phenylalanine__Phe__F, Glycine__Gly__G, Histidine__His__H, Isoleucine__Ile__I, Lysine__Lys__K, Leucine__Leu__L, Methionine__Met__M, Asparagine__Asn__N, Proline__Pro__P, Glutamine__Gln__Q, Arginine__Arg__R, Serine__Ser__S, Threonine__Thr__T, Valine__Val__V, Tryptophan__Trp__W, Tyrosine__Tyr__Y

