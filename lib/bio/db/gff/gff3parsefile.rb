

module Bio

  module GFFbrowser

    class GFF3ParseFile
      attr_reader :records, :sequences

      include FastLineParser
      
      def initialize fn
        @records = []
        @sequences = []
        fh = File.open(fn) 
        fh.each_line do | line |
          s = line.strip
          if s == '##FASTA'
            break
          end
          next if s.length == 0 or s =~ /^#/
          @records.push FastLineRecord.new(parse_line_fast(s))
        end
        fasta = Bio::GFF::FastaReader.new(fh)
        fasta.each do | id, fastarec |
          @sequences.push FastaRecord.new(id,fastarec)
        end
      end

    end
  end

end
