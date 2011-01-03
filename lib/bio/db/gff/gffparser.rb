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
        include Gff3Component
        include Gff3Features

        def store_record rec
            return if rec.comment # skip GFF comments
            id = Record::formatID(rec)
            @count_ids.add(id)
            @count_seqnames.add(rec.seqname)

            is_component = COMPONENT_TYPES.include?(rec.feature_type)
            if is_component
              # check for container ID
              warn("Container <#{rec.feature_type}> has no ID, so using sequence name instead",id) if rec.id == nil
              @componentlist[id] = rec
              info "Added #{rec.feature_type} with component ID #{id}"
            end 
            case rec.feature_type
              when 'gene' || 'SO:0000704'
                @orflist.add(id,rec)
              when 'mRNA' || 'SO:0000234'
                @mrnalist.add(id,rec)
              when 'CDS'  || 'SO:0000316'
                @cdslist.add(id,rec)
              when 'exon' || 'SO:0000147'
                @exonlist.add(id,rec)
              else
                if !is_component and !IGNORE_FEATURES.include?(rec.feature_type)
                  @unrecognized_features[rec.feature_type] = true
                end
            end
        end

        def validate_mrnas 
          return if not @options[:validate]
          # validate gene/container/component seqname is shared
          @mrnalist.validate_seqname
          @mrnalist.validate_shared_parent
        end

        def validate_cdss 
          return if not @options[:validate]
          @cdslist.validate_seqname
          # validate CDS sections do not overlap
          @cdslist.validate_nonoverlapping
          # validate sections share the parent
          @cdslist.validate_shared_parent
          # display unhandled features
        end

        def show_unrecognized_features 
          @unrecognized_features.keys.each do | k |
            warn "Feature has no match",k if k
          end
        end

        def read_fasta
          if @options[:fasta_filename]
            File.open(@options[:fasta_filename]) do | f |
              fasta = Bio::GFF::FastaReader.new(f)
              fasta.each do | id, fastarec |
                # p fastarec
                @sequencelist[id] = fastarec
              end
            end
          end
          # p :inmemory, @sequencelist
        end

        # Yield the id, recs, containing component and sequence of genes
        def each_gene
          parse if !@orflist
          each_item(@orflist) { |id, recs, component | yield id, recs, component }
        end

        # Yield the id, recs, containing component and sequence of mRNAs
        def each_mRNA
          parse if !@mrnalist
          each_item(@mrnalist) { |id, recs, component | yield id, recs, component }
        end

        # Yield the id, recs, and containing component
        def each_CDS
          parse if !@cdslist
          each_item(@cdslist) { |id, recs, component | yield id, recs, component }
        end

        # Yield the id, recs, and containing component
        def each_exon
          parse if !@exonlist
          each_item(@exonlist) { |id, recs, component | yield id, recs, component }
        end

        # Yield a unique description and the sequence
        def each_gene_seq
          each_gene do | id, reclist, component |
            if component
              sequence = @sequencelist[component.seqname]
              # p sequence
              if sequence
                yield description(id,component,reclist), assemble(sequence,component.start,reclist)
              else 
                error "No sequence information for",id
              end
            end
          end
        end

        # Yield a unique description and the sequence
        def each_mRNA_seq
          each_mRNA do | id, reclist, component |
            if component
              sequence = @sequencelist[component.seqname]
              # p sequence
              if sequence
                yield description(id,component,reclist), assemble(sequence,component.start,reclist)
              else 
                error "No sequence information for",id
              end
            end
          end
        end

        # Yield a unique description and the sequence
        def each_CDS_seq
          each_CDS do | id, reclist, component |
            if component
              sequence = @sequencelist[component.seqname]
              # p sequence
              if sequence
                seq = assemble(sequence,component.start,reclist,:codonize=>true)
                if seq.size % 3 != 0
                  p reclist # leave this in
                  # raise "CDS size #{seq.size} is not a multiple of 3! <#{seq}>"
                  warn "CDS size is not a multiple of 3",id
                end
                yield description(id,component,reclist), seq
              else 
                error "No sequence information for",id
              end
            end
          end
        end

        # Yield a unique description and the sequence
        def each_exon_seq
          each_exon do | id, reclist, component |
            if component
              sequence = @sequencelist[component.seqname]
              if sequence
                seq = assemble(sequence,component.start,reclist)
                yield description(id,component,reclist), seq
              else 
                error "No sequence information for",id
              end
            end
          end
        end
      end
    end
  end
end
