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
        # include Gff3Component
        # include Gff3Features
        # include Gff3Sequence
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
