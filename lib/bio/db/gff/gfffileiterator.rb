#
# = bio/db/gff/gfffileiterator.rb - Fetch records from a file
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License

module Bio

  class GFF

    class GFF3

      class FileRecord < Record
        attr_accessor :io_seek
        def initialize io_seek, buf
          @io_seek = io_seek
          super
        end
      end

      # GFF3::FileIterator takes a file and yields GFF3 records with their
      # seek position included in the record.
      class FileIterator
        attr_accessor :fh
        attr_reader :fasta_io_seek

        def initialize filename
          @fh = File.open(filename)
        end

        # Iterate over every record in the file, yielding the record ID and
        # (File)Record, which includes the io_seek position in the file
        def each_rec
          fpos = 0
          @fh.each_line do | line |
            if line.strip == "##FASTA"
              @fasta_io_seek = fpos
              break
            end
            if line.strip.size != 0 and line !~ /^#/
              rec = FileRecord.new(fpos, line)
              lastpos = @fh.tell
              yield rec.id, rec
              @fh.seek(lastpos) # reset filepos, just in case it changed
            end
            fpos = @fh.tell
          end
        end

        # Iterate over all contained FASTA sequences, yielding the ID
        # and sequence as a FASTA record. Normally call each_rec first and 
        # you can test for existing FASTA records if fasta_io_seek != nil
        def each_sequence
          if @fasta_io_seek == nil
            # Find the FASTA location first
            @fh.each_line do | line |
              break if line.strip == "##FASTA"
            end
          else
            @fh.seek(@fasta_io_seek)
          end
          # read FASTA records
          seqs   = []
          header = nil
          @fh.each_line do | line |
            line = line.strip
            next if line =~ /^#/
            if line =~ /^>/  # FASTA record header
              yield fasta_rec(header, seqs) if header
              header = line
              seqs   = []
            else
              seqs << line
            end
          end
          yield fasta_rec(header, seqs) if header
        end

      private

        def fasta_rec header, buf
          fst = Bio::FastaFormat.new(header+"\n"+buf.to_s)
          return fst.definition, fst
        end
      end

    end
  end
end

    


