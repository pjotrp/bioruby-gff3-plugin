#
# = bio/db/gff/gffnocache.rb - Assemble mRNA and CDS from GFF by fseek
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#
# Fetch information from a GFF file without using RAM - also check
# out the caching edition, which uses limited amounts of RAM

module Bio
  module GFFbrowser

    module Digest

      module NoCacheHelpers 

        module SeekRec
          # Fetch a record using fh and file seek position
          def SeekRec::fetch(fh,fpos)
            return nil if fh==nil or fpos==nil
            fh.seek(fpos)
            GFF::GFF3::FileRecord.new(fpos, fh.gets)
          end
        end

        # The hardwired to file RecList
        class SeekRecList 
          def initialize fh
            @fh = fh
            @h = {}
          end

          def []= id, rec
            raise "id #{id} occurs twice!" if @h[id]
            fpos = rec.io_seek
            @h[id] = fpos
          end

          def [](id)
            fpos = @h[id]
            SeekRec::fetch(@fh,fpos)
          end

          def each 
            @h.each do | id,fpos |
              yield id, self[id]
            end
          end
        end

        class SeekLinkedRecs < Hash
          include Helpers::Error
          def add id, rec
            info "Adding #{rec.feature_type} <#{id}>"
            self[id] = [] if self[id] == nil
            self[id] << rec.io_seek
          end
          def validate_seqname
          end
          def validate_nonoverlapping
          end
          def validate_shared_parent
          end
        end
      end

      class NoCache
        include Parser
        include NoCacheHelpers
        include Gff3Sequence

        def initialize filename
          @filename = filename
          @iter = Bio::GFF::GFF3::FileIterator.new(@filename)
        end

        # parse the whole file once and store all seek locations, 
        # rather than the records themselves
        def parse
          info "---- Digest DB and store data in mRNA Hash (NoCache)"
          @count_ids          = Counter.new   # Count ids
          @count_seqnames     = Counter.new   # Count seqnames
          @componentlist      = SeekRecList.new(@iter.fh) # Store containers, like genes, contigs
          @mrnalist           = SeekLinkedRecs.new   # Store linked mRNA records
          @cdslist            = SeekLinkedRecs.new
          @exonlist           = SeekLinkedRecs.new
          @sequencelist       = {}
          unrecognized_features = {}
          @iter.each_rec do | id, rec |
            store_record(rec)
          end
          @iter.each_sequence do | id, seq |
            id = seq.entry_id
            @sequencelist[id] = seq
          end
          validate_mrnas 
          validate_cdss
          show_unrecognized_features unrecognized_features
          @genelist      = @count_ids.keys 
        end

        def each_item list
          p list.class
          fh = @iter.fh
          list.each do | id, io_seeklist |
            recs = []
            io_seeklist.each do | fpos |
              fh.seek(fpos)
              recs << GFF::GFF3::FileRecord.new(fpos, fh.gets)
            end
            seqid = recs[0].seqname
            component = find_component(recs[0])
            yield id, recs, component
          end
        end
        
      end
    end
  end
end
