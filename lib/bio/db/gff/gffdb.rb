#
# = bio/db/gff/gffdb.rb - GFF database class
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License

# Create db from a GFF file

require 'bio'
require 'bio/db/gff/gffassemble'
require 'bio/db/gff/gffinmemory'

module Bio
  module GFFbrowser
    class GFFdb
      attr_reader :assembler

      include Digest

      # Initialize a GFF parser
      def initialize filename, options = {}
        cache_recs    = options[:cache_records]
        @assembler = 
          case cache_recs
            when :cache_none :
              NoCache.new(filename)
            else
              InMemory.new(filename)  # default 
          end
      end

    end # GFFdb
  end # GFFbrowser
end # Bio

