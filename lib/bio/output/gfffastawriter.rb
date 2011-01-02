
module Bio
  module GFFbrowser

    class FastaWriter
      def initialize translate
        @do_translate = translate
      end

      def put id, seq
        puts '>'+id
        puts seq
      end

    end
  end
end
