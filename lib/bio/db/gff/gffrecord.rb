require 'bio/db/gff/gff3parserec'

module Bio
  module GFFbrowser
    
    # Base class for a record
    class Record
    end

    # Using the fast line parser
    class FastLineRecord < Record

      include FastLineParser

      def initialize fields
        @fields = fields
      end

      def comment
        false
      end

      def seqid
        @seqid_ ||= @fields[GFF3_SEQID]
      end

      alias seqname :seqid

      def phase
        @phase_ ||= @fields[GFF3_PHASE].to_i
      end

      alias frame :phase

      def start
        @start_ ||= @fields[GFF3_START].to_i
      end

      def end
        @end_ ||= @fields[GFF3_END].to_i
      end

      def score
        @score_ ||= @fields[GFF3_SCORE].to_f
      end

      def strand
        @fields[GFF3_STRAND]
      end

      def feature
        @feature_ ||= @fields[GFF3_TYPE]
      end

      alias feature_type :feature
      def source
        @fields[GFF3_SOURCE]
      end

      def attributes 
        @attributes_ ||= parse_attributes_fast(@fields[GFF3_ATTRIBUTES])
      end

      def get_attribute name
        attributes[name]
      end 

      def id
        @id_ ||= attributes['ID']
      end

      alias entry_id :id
    end
  end
end

