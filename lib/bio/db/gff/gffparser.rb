#
# = bio/db/gff/gffparser.rb - Parsing logic for GFF3 file
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#

module Bio
  module GFFbrowser
    module Digest

      module Parser

        include Bio::GFFbrowser::Helpers
        include Bio::GFFbrowser::Helpers::Error
        # include NoCacheHelpers
        include Gff3Component
        include Gff3Features
        # include Gff3Sequence

        def store_record rec
            return if rec.comment # skip GFF comments
            id = Record::formatID(rec)
            @count_ids.add(id)
            @count_seqnames.add(rec.seqname)

            if COMPONENT_TYPES.include?(rec.feature_type)
              # check for container ID
              warn("Container <#{rec.feature_type}> has no ID, so using sequence name instead",id) if rec.id == nil
              @componentlist[id] = rec
              info "Added #{rec.feature_type} with component ID #{id}"
            else
              case rec.feature_type
                when 'mRNA' || 'SO:0000234' : @mrnalist.add(id,rec)
                when 'CDS'  || 'SO:0000316' : @cdslist.add(id,rec)
                when 'exon' || 'SO:0000147' : @exonlist.add(id,rec)
                else
                  if !IGNORE_FEATURES.include?(rec.feature_type)
                    unrecognized_features[rec.feature_type] = true
                  end
              end
            end
        end

        def validate_mrnas 
          # validate gene/container/component seqname is shared
          @mrnalist.validate_seqname
          @mrnalist.validate_shared_parent
        end

        def validate_cdss 
          @cdslist.validate_seqname
          # validate CDS sections do not overlap
          @cdslist.validate_nonoverlapping
          # validate sections share the parent
          @cdslist.validate_shared_parent
          # display unhandled features
        end

        def show_unrecognized_features unrecognized_features
          unrecognized_features.keys.each do | k |
            warn "Feature has no match",k if k
          end
        end

        # Yield the id, recs, component and sequence of mRNAs
        def each_mRNA
          parse if !@mrnalist
          each_item(@mrnalist) { |id, recs, component | yield id, recs, component }
        end

        # Yield the id, recs, and component
        def each_CDS
          parse if !@cdslist
          each_item(@cdslist) { |id, recs, component | yield id, recs, component }
        end

        # Yield the id, recs, and component
        def each_exon
          parse if !@exonlist
          each_item(@exonlist) { |id, recs, component | yield id, recs, component }
        end

        def each_mRNA_seq
          each_mRNA do | id, reclist, component |
            if component
              sequence = @sequencelist[component.seqname]
              # p sequence
              if sequence
                yield description(id,component,reclist), assemble(sequence,component.start,reclist)
              else 
                warn "No sequence information for",id
              end
            end
          end
        end

        def each_CDS_seq
          each_CDS do | id, reclist, component |
            if component
              sequence = @sequencelist[component.seqname]
              # p sequence
              if sequence
                seq = assemble(sequence,component.start,reclist)
                if seq.size % 3 != 0
                  p reclist
                  raise "CDS size #{seq.size} is not a multiple of 3! <#{seq}>"
                end
                yield description(id,component,reclist), seq
              else 
                warn "No sequence information for",id
              end
            end
          end
        end

        def each_exon_seq
          each_exon do | id, reclist, component |
            if component
              sequence = @sequencelist[component.seqname]
              if sequence
                seq = assemble(sequence,component.start,reclist)
                yield description(id,component,reclist), seq
              else 
                warn "No sequence information for",id
              end
            end
          end
        end

      end
    end
  end
end
