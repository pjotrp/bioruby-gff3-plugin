
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

        def parse(gfftype)
          @inseqidlist = {}
          # Fetch FASTA first
          @sequencelist = {}
          if @options[:fasta_filename]
            File.open(@options[:fasta_filename]) do | f |
              fasta = Bio::GFF::FastaReader.new(f)
              fasta.each do | id, fastarec |
                # p fastarec
                @sequencelist[id] = fastarec
              end
            end
          else
            # Embedded FASTA
            @iter.each_sequence do | id, bioseq |
              @sequencelist[id] = bioseq.to_s
            end
          end
          seqid = nil
          recs = []
          @iter.each_rec do | fpos, line |
            rec = FastLineRecord.new(parse_line_fast(line))
            if seqid != rec.seqid 
              # starting a new block
              if @inseqidlist[rec.seqid]
                # not a well formed GFF3 file, we need
                # to drop
                error "GFF3 file not sorted, falling back to line parser"
                raise "ERROR, bailing out"
              end
              parse_block(gfftype,recs,@sequencelist[seqid])  { | id, seq | yield id,seq } if seqid
              recs = []
              seqid = rec.seqid
              @inseqidlist[seqid] = true
            end
            recs.push rec
          end
          parse_block(gfftype,recs,@sequencelist[seqid])  { | id, seq | yield id,seq } if seqid
        end

        # Parse sequence objects sharing the same seqid 
        # and yield each +gfftype+ as an iq,seq
        def parse_block gfftype, recs, sequence
          recs.each do | rec |
            if rec.feature_type.downcase == gfftype
              yield rec.id, sequence[rec.start-1..rec.end-1]
            end
          end
        end
   
        def each_seq(gfftype) 
          parse(gfftype) { | id, seq | yield id,seq }

        end

        def each_gene_seq
          each_seq('gene') { | id, seq | yield id,seq }
        end
        def each_mRNA_seq
          each_seq('mrna') { | id, seq | yield id,seq }

        end
        def each_exon_seq
          each_seq('exon') { | id, seq | yield id,seq }

        end
        def each_CDS_seq
          each_seq('cds') { | id, seq | yield id,seq }

        end
      end
    end
  end
end
