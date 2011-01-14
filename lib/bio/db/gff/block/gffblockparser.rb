
module Bio
  module GFFbrowser

    module Block

      # The block parser simplifies parsing, by assuming GFF3 is 
      # organised into blocks. All relevant information is 
      # resolved a block at a time.
      class GffBlockParser
        def initialize filename, options
          info "Starting block parser"
          @filename = filename
          @options = options
          @iter = Bio::GFF::GFF3::FileIterator.new(@filename, options[:parser])
        end

        def parse
          @sequencelist = []
          @iter.each_rec do | id, rec |
            p rec
          end
          @iter.each_sequence do | id, bioseq |
            @sequencelist[id] = bioseq.to_s
          end
        end
   
        def each_seq(gfftype) 
          parse
        end

        def each_gene_seq
          each_seq('gene')
        end
        def each_mRNA_seq
          each_seq('mrna')
        end
        def each_exon_seq
          each_seq('exon')
        end
        def each_CDS_seq
          each_seq('cds')
        end
      end
    end
  end
end
