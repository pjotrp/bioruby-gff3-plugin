
module Bio
  module GFFbrowser

    module Block

      # The block parser simplifies parsing, by assuming GFF3 is 
      # organised into blocks. All relevant information is 
      # resolved a block at a time.
      class GffBlockParser
        include FastLineParser 

        def initialize filename, options
          info "Starting block parser"
          @filename = filename
          @options = options
          @iter = Bio::GFF::GFF3::FileIterator.new(@filename)
        end

        def parse
          @inseqidlist = {}
          @sequencelist = {}
          seqid = nil
          block = []
          @iter.each_rec do | fpos, line |
            rec = FastLineRecord.new(parse_line_fast(line))
            if seqid != rec.seqid 
              # starting a new block
              if @inseqidlist[rec.seqid]
                # whoah, not a well formed GFF3 file, we need
                # to drop
                error "GFF3 file not sorted, falling back to line parser"
                raise "ERROR, bailing out"
              end
              parse_block(block) if seqid
              block = []
              seqid = rec.seqid
              @inseqidlist[seqid] = true
            end
            block.push rec
          end
          parse_block(block)
          @iter.each_sequence do | id, bioseq |
            @sequencelist[id] = bioseq.to_s
          end
        end

        def parse_block recs
              p "-------------"
              p recs
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
