#
# = bio/db/gff/gffinmemory.rb - Assemble mRNA and CDS from GFF in RAM
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#
# Fetch information from a GFF file

module Bio
  module GFFbrowser

    module Digest

      class InMemory
        include Parser
        # include Helpers
        # include Helpers::Error
        include Gff3Component
        include Gff3Features
        include Gff3Sequence

        def initialize filename
          @gff = Bio::GFF::GFF3.new(File.read(filename))
        end

        # Digest mRNA from the GFFdb and store in Hash
        # Next yield(id, seq) from Hash
        def parse 
          gff = @gff
          info "---- Digest DB and store data in mRNA Hash"
          count_ids       = Counter.new   # Count ids
          count_seqnames  = Counter.new   # Count seqnames
          components      = {} # Store containers, like genes, contigs
          mrnas           = LinkedRecs.new   # Store linked mRNA records
          cdss            = LinkedRecs.new
          exons           = LinkedRecs.new
          sequences       = {}
          unrecognized_features = {}
          gff.records.each do | rec |
            next if rec.comment # skip GFF comments
            id = Record::formatID(rec)
            count_ids.add(id)
            count_seqnames.add(rec.seqname)

            if COMPONENT_TYPES.include?(rec.feature_type)
              # check for container ID
              warn("Container <#{rec.feature_type}> has no ID, so using sequence name instead",id) if rec.id == nil
              components[id] = rec
              info "Added #{rec.feature_type} with component ID #{id}"
            else
              case rec.feature_type
                when 'mRNA' || 'SO:0000234' : mrnas.add(id,rec)
                when 'CDS'  || 'SO:0000316' : cdss.add(id,rec)
                when 'exon' || 'SO:0000147' : exons.add(id,rec)
                else
                  if !IGNORE_FEATURES.include?(rec.feature_type)
                    unrecognized_features[rec.feature_type] = true
                  end
              end
            end
          end
          gff.sequences.each do | seq |
            id = seq.entry_id
            sequences[id] = seq
          end
          validate_mrnas mrnas
          validate_cdss cdss
          show_unrecognized_features unrecognized_features
          @genelist      = count_ids.keys 
          @componentlist = components
          @mrnalist      = mrnas
          @cdslist       = cdss
          @exonlist      = exons
          @sequencelist  = sequences
        end

        def each_item list
          list.each do | id, recs |
            seqid = recs[0].seqname
            component = find_component(recs[0])
            yield id, recs, component
          end
        end
        
      end
    end
  end
end
