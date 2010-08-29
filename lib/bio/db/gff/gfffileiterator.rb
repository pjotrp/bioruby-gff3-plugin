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

      # FileIterator takes a file and yields GFF3 records with their
      # seek position included in the record.
      class FileIterator
        attr_accessor :fh

        def initialize filename
          @fh = File.open(filename)
        end

        def each_rec
          fpos = 0
          inside_fasta = false
          @fh.each_line do | line |
            if line.strip == "##FASTA"
              inside_fasta = true
              break
            end
            if line.strip.size != 0 and line !~ /^#/
              rec = FileRecord.new(fpos, line)
              lastpos = @fh.tell
              yield rec.id, rec
              @fh.seek(lastpos)
            end
            fpos = @fh.tell
          end
          if inside_fasta
            print "READING FASTA"
          end
        end

      end

    end
  end
end

    


