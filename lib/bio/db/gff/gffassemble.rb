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
        def info str, id=''
          $stderr.print "Info: "+str+" <#{id}>\n"
        end

        def warn str,id=''
          Kernel.warn "Warning: "+str+" <#{id}>"
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
        include Error
        # Format a record ID by, first, getting the ID attribute. If that fails
        # the seqname is used with the start/stop positions.
        def Record::formatID rec  
          id = rec.id if rec.id
          if !id
            if rec.seqname
              id = "#{rec.seqname} #{rec.start} #{rec.end}".strip
            else
              id = 'unknown'
              $stderr.print "Record with unknown ID"+rec.to_s
            end
          end
          id
        end
      end

      module Gff3Component

        include Error

        COMPONENT_TYPES = %w{
          gene SO:0000704 contig transcript Component region
        }
 
        # Walk the component list to find a matching component/container for a
        # record. First use the parent ID. If that is missing go by sequence
        # name.
        def find_component rec
          parent = rec.get_attribute('Parent')
          if @componentlist[parent] 
            # nice, there is a match
            info "find_component: Matched parent", parent
            return @componentlist[parent]
          end
          search = rec.seqname
          if @componentlist[search]
            info "find_component: Matched seqname", search
            return @componentlist[search] 
          end
          @componentlist.each do | componentid, component |
            # dissemble id
            (id, start, stop) = componentid.split(/ /)
            if id==search and rec.start >= start.to_i and rec.end <= stop.to_i
              info "find_component: Matched column 0 and location", componentid
              return component
            end
          end
          # Ah, painful. At this point the record has no matching container, probably
          # because it has no parent ID and the component has an ID. We have to go by
          # ID for every component individually
          @componentlist.each do | componentid, component |
            if component.seqname==search and rec.start >= component.start and rec.end <= component.end
              # p ["----",search,rec]
              # p component
              info "find_component: Matched (long search) column 0 and location", componentid
              return component
            end
          end
          warn "Could not find container/component for",Record::formatID(rec)
        end
      end

      module Gff3Features

        # Ignore the following features (case sensitive?)
        IGNORE_FEATURES = Gff3Component::COMPONENT_TYPES + %w{
          transposon Match similarity UTR
          TF_binding_site intronSO:0000188 polyA_sequence SO:0000610
          polyA_site SO:0000553
          five_prime_UTR SO:0000204 three_prime_UTR SO:0000205
          exon SO:0000147
        }
      end

      module Gff3Sequence
        # Patch a sequence together from a Sequence string and an array
        # of records. Note that rec positions are 1-based coordinates, relative 
        # to the landmark given in column 1 - in this case the sequence as it
        # is passed in. The container startpos is unused.
        def assemble sequence, startpos, rec, options = { :codonize=>true }
          codonize = options[:codonize]
          retval = ""
          Sections::sort(rec).each do | section |
            # p sequence
            if sequence.kind_of?(Bio::FastaFormat)
              sequence = sequence.seq
            end
            rec1 = section.rec
            # For forward strand features, phase is counted from the start
            # field. For reverse strand features, phase is counted from the end
            # field. 
            frame = 0
            frame = rec1.frame if rec1.frame
            seq = sequence[(rec1.start-1)..rec1.end]
            # if strand is negative, reverse
            reversed = (rec1.strand == '-')
            seq = seq.reverse if reversed
            # correct phase and size to multiple of 3
            if codonize
              seq = seq[frame..-1] if frame != 0 # set phase
              # seq = seq[1..-1] if reversed # test bug
              reduce = seq.size % 3
              seq = seq[0..(seq.size-1-reduce)] if reduce
            end
            retval += seq
          end
          retval
        end

        # Patch a sequence together from a Sequence string and an array
        # of records and translate in the correct direction and frame
        def assembleAA sequence, startpos, rec
          seq = assemble(sequence, startpos, rec)
          ntseq = Bio::Sequence::NA.new(seq)
          ntseq.translate
        end

        # Create a description for output
        def description id, component, rec
          sections = Sections::sort(rec)
          id+' Sequence:'+component.seqname+"_#{component.start}:#{component.end} ("+
           sections.map { |s| "#{s.first}:#{s.last}" }.join(', ')  +")"
        end

      end
    end # Helpers

  end
end
