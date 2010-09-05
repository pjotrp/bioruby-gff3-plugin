#
# = bio/db/gff/gffdb.rb - GFF database class
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License

# Create db from a GFF file

require 'bio'
require 'bio/db/gff/gfffileiterator'
require 'bio/db/gff/gfffasta'
require 'bio/db/gff/gffassemble'
require 'bio/db/gff/gffparser'
require 'bio/db/gff/gffinmemory'
require 'bio/db/gff/gffnocache'

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
            when :cache_none 
              NoCache.new(filename, options)
            else
              InMemory.new(filename, options)  # default 
          end
      end

    end # GFFdb
  end # GFFbrowser
end # Bio

