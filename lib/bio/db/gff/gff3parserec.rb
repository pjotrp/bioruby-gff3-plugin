#
# = Record parsing logic for GFF3
#
# Copyright::  Copyright (C) 2011
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#

module Bio
  module GFFbrowser

    GFF3_SEQID      = 0
    GFF3_SOURCE     = 1
    GFF3_TYPE       = 2
    GFF3_START      = 3
    GFF3_END        = 4
    GFF3_SCORE      = 5
    GFF3_STRAND     = 6
    GFF3_PHASE      = 7
    GFF3_ATTRIBUTES = 8

    # The fast line parse takes the least number of actions
    # to parse a GFF3 line (record).
    #
    module FastLineParser

      include Helpers::Logger

      # Returns a (partial) record, assuming it is a valid GFF3
      # format, no validation takes place, other than field counting (!)
      #
      # Start and End are always set to int values. ID and Parent 
      # are returned in features.
      #
      # The options define further parsing of fields. Options are
      #   
      #
      # An array is returned, with appropriate fields 
      #
      def parse_line_fast string, options = {}
        fs = string.split(/\t/)
        if fs.size != 9 
          error "Record should have 9 fields, but has #{fs.size} <#{string}>"
          return nil
        end

        fs
      end

      def parse_attributes_fast attribstring, options = {}
        Hash[attribstring.split(/;/).map { | a |
          a.split(/=/,2)
        }]
      end
    
    end
  end
end
