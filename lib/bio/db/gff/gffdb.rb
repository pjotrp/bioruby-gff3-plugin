#
# = bio/db/gff/gff3.rb - GFF database class
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License

# Create a GFF3 file parser


module Bio
  module GFFbrowser
    class GFF3
      attr_reader :assembler

      include Digest

      # Initialize a GFF parser
      def initialize filename, options = {}
        cache_recs    = options[:cache_records]
        @assembler = 
          case cache_recs
            when :cache_none 
              NoCache.new(filename, options)
            else
              InMemory.new(filename, options)  # default 
          end
      end

    end # GFF3
  end # GFFbrowser
end # Bio

