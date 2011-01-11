module Bio
  module GFFbrowser
    
    # Base class for a record
    class Record
    end

    # Using the fast line parser
    class FastLineRecord < Record
      def initialize fields
        @fields = fields
      end

      def seqid
        @fields[GFF3_SEQID]
      end

      alias seqname :seqid

      def phase
        @fields[GFF3_PHASE]
      end

      alias frame :phase

      def start
        @fields[GFF3_START]
      end

      def end
        @fields[GFF3_END]
      end

      def score
        @fields[GFF3_SCORE]
      end

      def strand
        @fields[GFF3_STRAND]
      end

      def feature
        @fields[GFF3_TYPE]
      end

      def source
        @fields[GFF3_SOURCE]
      end

      def attributes 
        @fields[GFF3_ATTRIBUTES]
      end

      def id
        attributes['ID']
      end
    end
  end
end

