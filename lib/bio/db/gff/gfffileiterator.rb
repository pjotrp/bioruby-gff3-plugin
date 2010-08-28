#
# = bio/db/gff/gfffileiterator.rb - Fetch records from a file
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License

module Bio

  class GFF

    class GFF3
 
      class FileIterator
        attr_accessor :fh

        def initialize filename
          @fh = File.open(filename)
        end

        def each_rec
          fpos = 0
          @fh.each_line do | line |
            p fpos
            if line.strip.size != 0 and line !~ /^#/
              rec = Bio::GFF::GFF3::Record.new(line)
              rec.instance_variable_set(:@io_seek, fpos) # set fpos at begin line
              lastpos = @fh.tell
              yield rec.id, rec
              @fh.seek(lastpos)
            end
            fpos = @fh.tell
          end
        end

      end

    end
  end
end

    


