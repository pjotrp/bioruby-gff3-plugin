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

      class NoCache
        include Helpers
        include Helpers::Error
        include Gff3Component
        include Gff3Features
        include Gff3Sequence

        def initialize filename
          @filename = filename
        end

        # parse the whole file once and store all seek locations, 
        # rather than the records themselves
        def parse
          info "---- Digest DB and store data in mRNA Hash (NoCache)"
          count_ids       = Counter.new   # Count ids
          count_seqnames  = Counter.new   # Count seqnames
          components      = {} # Store containers, like genes, contigs
          mrnas           = LinkedRecs.new   # Store linked mRNA records
          cdss            = LinkedRecs.new
          exons           = LinkedRecs.new
          sequences       = {}
          unrecognized_features = {}

          fh = File.open(@filename)
          gff = nil

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
          # validate gene/container/component seqname is shared
          mrnas.validate_seqname
          cdss.validate_seqname
          # validate CDS sections do not overlap
          cdss.validate_nonoverlapping
          # validate sections share the parent
          mrnas.validate_shared_parent
          cdss.validate_shared_parent
          # display unhandled features
          unrecognized_features.keys.each do | k |
            warn "Feature has no match",k if k
          end
          # Finally we are going to index the sequences
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
