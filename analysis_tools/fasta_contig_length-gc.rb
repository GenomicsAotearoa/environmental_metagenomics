#Use to put header and contig length.

fastafile = File.read(ARGV[0])

test = Array.new

totcontiglen = 0

fasta = fastafile.scan(/^>(\S+.*)\n([^>]*)/)
header,contig = $1,$2

fasta.each do |header,contig|
    contiglen = ((contig.to_s).count('/[A-Za-z]/')).to_f
    gccount = ((contig.to_s).count("GgCc")).to_f
    gc = (gccount/contiglen)*100
    puts "#{header}\tLength=#{contiglen}\tGC=#{gc}"
end


