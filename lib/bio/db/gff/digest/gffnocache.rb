#
# = bio/db/gff/gffnocache.rb - Assemble mRNA and CDS from GFF by fseek
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#
# Fetch information from a GFF file without using RAM - also check
# out the caching edition, which uses limited amounts of RAM

require 'bio/db/gff/digest/gffparser'

module Bio
  module GFFbrowser

    module Digest

      module NoCacheHelpers 

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

        # Helper class which gives Hash-like access to the 
        # no-cache GFF3 file
        class SeekRecList 
          def initialize fh, parser
            @fh = fh
            @h = {}
            @parser = parser
          end

          def []= id, rec
            raise "id #{id} occurs twice!" if @h[id]
            fpos = rec.io_seek
            @h[id] = fpos
          end

          def [](id)
            fpos = @h[id]
            SeekRec::fetch(@fh,fpos,@parser)
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
          # validation is switched off for NoCache
          def validate_seqname
          end
          # validation is switched off for NoCache
          def validate_nonoverlapping
          end
          # validation is switched off for NoCache
          def validate_shared_parent
          end
        end
      end

      class NoCache
        include Parser
        include NoCacheHelpers
        include Gff3Sequence

        def initialize filename, options
          @filename = filename
          @options = options
          @iter = Bio::GFF::GFF3::FileIterator.new(@filename)
        end

        # parse the whole file once and store all seek locations, 
        # rather than the records themselves
        def parse
          info "---- Digest DB and store data in mRNA Hash (NoCache)"
          @count_ids          = Counter.new   # Count ids
          @count_seqnames     = Counter.new   # Count seqnames
          @componentlist      = SeekRecList.new(@iter.fh,@options[:parser]) # Store containers, like genes, contigs
          @orflist            = SeekLinkedRecs.new   # Store linked gene records
          @mrnalist           = SeekLinkedRecs.new   # Store linked mRNA records
          @cdslist            = SeekLinkedRecs.new
          @exonlist           = SeekLinkedRecs.new
          @sequencelist       = {}
          @unrecognized_features = {}
          @iter.each_rec do |fpos, line|
            rec = case @options[:parser]
              when :bioruby
                Bio::GFF::GFF3::BioRubyFileRecord.new(fpos, line)
              when :line
                Bio::GFF::GFF3::FastParserFileRecord.new(fpos, line)
              else
                raise 'Unknown parser'
            end
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
              recs << SeekRec::fetch(fh,fpos,@options[:parser])
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
