
module Bio
  module GFFbrowser

    class FastaWriter
      def initialize translate
        @do_translate = translate
      end

      def put id, seq
        puts '>'+id
        if @do_translate
          ntseq = Bio::Sequence::NA.new(seq)
          puts ntseq.translate
        else
          puts seq
        end
      end

    end
  end
end
