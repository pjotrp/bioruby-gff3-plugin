#
# = bio/db/gff/gfflrucache.rb - Assemble mRNA and CDS from GFF by LRU cache
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#
# Fetch information from a GFF file without using RAM - using a least 
# recently used cache 'LRU' - also check out the caching edition, 
# which uses limited amounts of RAM.
#
# In effect it is NoCache with parser recs cached in the LRU hash

require 'bio/db/gff/digest/gffparser'
require 'bio/system/lruhash'

module Bio
  module GFFbrowser

    module Digest

      module LruCacheHelpers 

        # Module to fetch a line from GFF3 file and returns a parsed 
        # record
        module SeekRec
          # Fetch a record using fh and file seek position
          def SeekRec::fetch(fh,fpos,parser)
            return nil if fh==nil or fpos==nil
            fh.seek(fpos)
            if parser == :bioruby 
              GFF::GFF3::BioRubyFileRecord.new(fpos, fh.gets)
            else
              GFF::GFF3::FastParserFileRecord.new(fpos, fh.gets)
            end
          end
        end

        module LruRec
          # Fetch a record using fh and file seek position,
          # utilising the LRU cache
          def fetch(fh,fpos,parser)
            return nil if fh==nil or fpos==nil
            rec = @lru[fpos]
            if rec==nil
              rec = SeekRec::fetch(fh,fpos,parser)
              @lru[fpos] = rec
            end
            rec
          end
        end

        # Helper class which gives Hash-like access to the 
        # no-cache GFF3 file
        class SeekRecList 
          include LruRec

          def initialize fh, parser, lru
            @fh = fh
            @h = {}
            @parser = parser
            @lru = lru
          end

          def []= id, rec
            raise "id #{id} occurs twice!" if @h[id]
            fpos = rec.io_seek
            @h[id] = fpos
          end

          def [](id)
            fpos = @h[id]
            fetch(@fh,fpos,@parser)
          end

          def each 
            @h.each do | id,fpos |
              yield id, self[id]
            end
          end
        end

        # List of ids
        class SeekLinkedRecs < Hash
          include Helpers::Error
          def add id, rec
            info "Adding #{rec.feature_type} <#{id}>"
            self[id] = [] if self[id] == nil
            self[id] << rec.io_seek
          end
          # validation is switched off for LruCache
          def validate_seqname
          end
          # validation is switched off for LruCache
          def validate_nonoverlapping
          end
          # validation is switched off for LruCache
          def validate_shared_parent
          end
        end
      end

      class LruCache
        include Parser
        include LruCacheHelpers
        include Gff3Sequence
        include LruRec

        def initialize filename, options
          @filename = filename
          @options = options
          @iter = Bio::GFF::GFF3::FileIterator.new(@filename, options[:parser])
          @lru = LRUHash.new 10000
        end

        # parse the whole file once and store all seek locations, 
        # rather than the records themselves
        def parse
          info "---- Digest DB and store data in mRNA Hash (LruCache)"
          @count_ids          = Counter.new   # Count ids
          @count_seqnames     = Counter.new   # Count seqnames
          @componentlist      = SeekRecList.new(@iter.fh,@options[:parser],@lru) # Store containers, like genes, contigs
          @orflist            = SeekLinkedRecs.new   # Store linked gene records
          @mrnalist           = SeekLinkedRecs.new   # Store linked mRNA records
          @cdslist            = SeekLinkedRecs.new
          @exonlist           = SeekLinkedRecs.new
          @sequencelist       = {}
          @unrecognized_features = {}
          @iter.each_rec do | id, rec |
            store_record(rec)
          end
          @iter.each_sequence do | id, bioseq |
            @sequencelist[id] = bioseq.to_s
          end
          validate_mrnas 
          validate_cdss
          show_unrecognized_features
          @genelist      = @count_ids.keys 
          read_fasta
        end

        def each_item list
          # p list.class
          fh = @iter.fh
          list.each do | id, io_seeklist |
            recs = []
            io_seeklist.each do | fpos |
              recs << fetch(fh,fpos,@options[:parser])
            end
            seqid = recs[0].seqname
            component = find_component(recs[0])
            if @options[:no_assemble]
              recs.each do | rec |
                yield id, [rec], component
              end
            else
              yield id, recs, component
            end
          end
        end
        
      end
    end
  end
end
