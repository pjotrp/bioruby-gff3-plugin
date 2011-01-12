
module Bio

  module GFFbrowser

    class FastaRecord

      attr_reader :id, :seq

      def initialize id, seq
        @id = id
        @seq = seq
      end

      def entry_id
        @id
      end

      def to_s
        @seq
      end
    end

  end
end
