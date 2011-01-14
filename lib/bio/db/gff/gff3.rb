#
# = bio/db/gff/gff3.rb - GFF database class
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License

# Create a GFF3 file parser

require 'bio/db/gff/digest/gffinmemory'
require 'bio/db/gff/digest/gffnocache'
require 'bio/db/gff/digest/gfflrucache'
require 'bio/db/gff/block/gffblockparser'

module Bio
  module GFFbrowser
    class GFF3
      attr_reader :assembler

      include Digest

      # Initialize a GFF parser
      def initialize filename, options = {}
        options[:parser] = :line if options[:parser] == nil
        cache_recs    = options[:cache_records]
        @assembler = 
          if options[:parser] == :block
            GffBlockParser.new(filename, options)
          else 
            case cache_recs
              when :cache_none 
                NoCache.new(filename, options)
              when :cache_lru
                LruCache.new(filename, options)
              else
                InMemory.new(filename, options)  # default 
            end
          end
      end

    end # GFF3
  end # GFFbrowser
end # Bio

