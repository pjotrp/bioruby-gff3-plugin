#
# = bio/db/gff/gffassemble.rb - Assemble mRNA and CDS from GFF
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#
# Fetch information from a GFF file

module Bio
  module GFFbrowser

    module Helpers

      module Error
        def info str
          $stderr.print "Info: "+str+"\n"
        end

        def warn str,id=''
          $stderr.print "Warning: "+str+" <#{id}>\n"
        end
      end

      # Helper class for counting IDs
      class Counter < Hash
        def add id
          self[id] = 0 if self[id] == nil
          self[id] += 1
        end
      end

      # Helper class for storing linked records based on a shared ID
      class LinkedRecs < Hash
        include Error
        def add id, rec
          info "Adding #{rec.feature_type} <#{id}>"
          self[id] = [] if self[id] == nil
          self[id] << rec
        end

        # Validate all lists belong to the same container/component
        def validate_seqname
          each do | id, rec |
            seqname = rec.first.seqname
            rec.each do | section |
              raise "Non-matching seqname #{section.seqname} in #{seqname}" if section.seqname != seqname
            end
          end
        end
     
        # Validate all lists share the same parent (if available). First checks
        # for Parent attribute, next for mRNA attribute
        def validate_shared_parent
          each do | id, rec |
            parent = rec.first.get_attribute('Parent')
            if parent
              rec.each do | section |
                _parent = section.get_attribute('Parent')
                raise "Non-matching parent #{_parent} and #{parent} in #{id}" if _parent != parent
              end
            end
            parent = rec.first.get_attribute('mRNA')
            if parent
              rec.each do | section |
                _parent = section.get_attribute('mRNA')
                raise "Non-matching parent #{_parent} and #{parent} in #{id}" if _parent != parent
              end
            end
          end
        end
     
        # walk all (CDS) lists for every container/component and 
        # validate they do not overlap
        def validate_nonoverlapping
          each do | id, rec |
            sections = Sections::sort(rec)
            sections.each_with_index do | check, i |
              neighbour = sections[i+1]
              if neighbour and check.intersection(neighbour)
                warn "Overlapping sections for ",id
              end
            end
          end
        end
      end

      class Section < Range
        attr_reader :rec
        def initialize rec
          super(rec.start,rec.end)
          @rec = rec
        end
        def intersection(other)
          raise ArgumentError, 'value must be a Range' unless other.kind_of?(Range)
          min, max = first, exclude_end? ? max : last
          other_min, other_max = other.first, other.exclude_end? ? other.max : other.last
          new_min = self === other_min ? other_min : other === min ? min : nil
          new_max = self === other_max ? other_max : other === max ? max : nil
          new_min && new_max ? new_min..new_max : nil
        end
        alias_method :&, :intersection
        def <=> other
          first <=> other.first
        end
      end

      module Sections
        # Return list of sorted Sections
        def Sections::sort rec
          sections = []
          rec.each do | section |
            sections.push Section.new(section)
          end
          sections.sort
        end
      end

      module Record
        # Format a record ID by, first, getting the ID attribute. If that fails
        # the seqname is used with the start/stop positions.
        def Record::formatID rec
          id = rec.id if rec.id
          if !id
            if rec.seqname
              id = rec.seqname+" #{rec.start} #{rec.end}"
            else
              id = 'unknown'
              warn "Record with unknown ID",rec.to_s
            end
          end
          id
        end
      end

      module Component
        # Walk the component list to find a matching component/container for a
        # record. First use the parent ID. If that is missing go by sequence
        # name.
        def find_component rec
          parent = rec.get_attribute('Parent')
          if @componentlist[parent] 
            # nice, there is a match
            return @componentlist[parent]
          end
          search = rec.seqname
          return @componentlist[search] if @componentlist[search]
          @componentlist.each do | componentid, component |
            # dissemble id
            (id, start, stop) = componentid.split(/ /)
            if id==search and rec.start >= start.to_i and rec.end <= stop.to_i
              return component
            end
          end
          # Ah, painful. At this point the record has no matching container, probably
          # because it has no parent ID and the component has an ID. We have to go by
          # ID for every component individually
          @componentlist.each do | componentid, component |
            if component.seqname==search and rec.start >= component.start and rec.end <= component.end
              return component
            end
          end
          warn "Could not find container/component for",Record::formatID(rec)
        end
      end

      module Gff3Sequence
        # Patch a sequence together from a Sequence string and an array
        # of records
        def assemble sequence, startpos, rec
          retval = ""
          Sections::sort(rec).each do | section |
            retval += sequence[(section.rec.start-startpos)..(section.rec.end-startpos)]
          end
          retval
        end

        # Create a description for output
        def description id, component, rec
          sections = Sections::sort(rec)
          id+' Sequence:'+component.seqname+"_#{component.start}:#{component.end} ("+
           sections.map { |s| "#{s.first}:#{s.last}" }.join(', ')  +")"
        end

      end
    end # Helpers

    module MRNA
      include Helpers
      include Helpers::Error
      include Component
      include Gff3Sequence

      attr_reader :genelist, :componentlist, :contiglist, :cdslist, :mrnalist, :sequencelist

      COMPONENT_TYPES = %w{
        gene SO:0000704 contig transcript Component region
      }
      # Ignore the following features (case sensitive?)
      IGNORE_FEATURES = COMPONENT_TYPES + %w{
        transposon Match similarity UTR
        TF_binding_site intronSO:0000188 polyA_sequence SO:0000610
        polyA_site SO:0000553
        five_prime_UTR SO:0000204 three_prime_UTR SO:0000205
        exon SO:0000147
      }

      # Digest mRNA from the GFFdb and store in Hash
      # Next yield(id, seq) from Hash
      def parse gff
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
        parse(@gff) if !@mrnalist
        each_item(@mrnalist) { |id, recs, component | yield id, recs, component }
      end

      # Yield the id, recs, and component
      def each_CDS
        parse(@gff) if !@cdslist
        each_item(@cdslist) { |id, recs, component | yield id, recs, component }
      end

      # Yield the id, recs, and component
      def each_exon
        parse(@gff) if !@exonlist
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
