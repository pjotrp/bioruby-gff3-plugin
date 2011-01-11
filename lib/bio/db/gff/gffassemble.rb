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
              log = Bio::Log::LoggerPlus['bio-gff3']
              log.warn "Record with unknown ID"+rec.to_s.chomp
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
        # is passed in. The following options are available:
        #
        #   :reverse      : do reverse if reverse is indicated (default true)
        #   :complement   : do complement if reverse is indicated (default true)
        #   :phase        : do set CDS phase (default false, normally ignore)
        #   :trim         : make sure sequence is multiple of 3 nucleotide bps (default true)
        #
        # special options:
        #
        #   :raw          : raw sequence (all above false)
        #   :codonize     : codon sequence (reverse, complement, and trim are true)
        #   :fix          : fix errors (default false)
        #
        def assemble sequence, startpos, reclist, options = { :phase=>false, :reverse=>true, :trim=>true, :complement=>true, :fix=>false, :debug=>false }
          # default to nil, if not passed in
          do_debug = options[:debug]
          do_phase = options[:phase]
          do_fix        = options[:fix]
          # default to true, if not passed in
          do_reverse = (options[:reverse] == false ? false : true)
          do_trim    = (options[:trim] == false ? false : true)
          do_complement = (options[:complement] == false ? false : true)
          if options[:raw]
            do_phase = false
            do_reverse = false
            do_trim = false
            do_complement = false
          elsif options[:codonize]
            do_phase = false
            do_reverse = true
            do_trim = true
            do_complement = true
          end
          sectionlist = Sections::sort(reclist)
          rec0 = sectionlist.first.rec
          # we assume ORF is always read in the same direction
          orf_reverse = (rec0.strand == '-')
          orf_frame = startpos - 1
          orf_frameshift = orf_frame % 3
          sectionlist = sectionlist.reverse if orf_reverse
          if do_debug
            debug options.to_s
            debug [:reverse,do_reverse].to_s
            debug [:complement,do_complement].to_s
            debug [:trim,do_trim].to_s
            debug [:orf_reverse, orf_reverse, rec0.strand].to_s
          end

          if sequence.kind_of?(Bio::FastaFormat)
            # BioRuby conversion
            sequence = sequence.seq
          end
          # Generate array of sequences
          seq = sectionlist.map { | section |
            rec = section.rec
            s = sequence[(section.begin-1)..(section.end-1)]
            if do_reverse and orf_reverse
              s = s.reverse 
            end
            # Correct for phase. Unfortunately the use of phase is ambiguous.
            # Here we check whether rec.start is in line with orf_frame. If it
            # is, we correct for phase. Otherwise it is ignored.
            if do_phase and rec.phase
              phase = rec.phase.to_i
              # if ((rec.start-startpos) % 3 == 0) 
              s = s[phase..-1]
              # end
            end
            s
          }
          # p seq
          seq = seq.join
          if do_complement and do_reverse and orf_reverse
            ntseq = Bio::Sequence::NA.new(seq)
            seq = ntseq.forward_complement.upcase
          end
          # This is the place to fix sequences (e.g. the Wormbase bug)
          if do_fix or @options[:fix] or @options[:fix_wormbase]
            if @options[:fix_wormbase] and rec0.id.index('gene1')==0
              # Wormbase gene1 only, so ignore rest
            else
              test_frame = 0
              ntseq = Bio::Sequence::NA.new(seq)
              aaseq = ntseq.translate
              if aaseq.count('*') > 1
                test_frame = 1
                seq = seq[1..-1] 
                ntseq = Bio::Sequence::NA.new(seq)
                aaseq = ntseq.translate
                if aaseq.count('*') > 1
                  test_frame = 2
                  seq = seq[1..-1]
                  ntseq = Bio::Sequence::NA.new(seq)
                  aaseq = ntseq.translate
                  raise 'Validation problem '+rec0.id if aaseq.count('*') > 1
                end
              end
              if test_frame > 0
                warn rec0.id,"Frame adjusted to #{test_frame} (fix)"
              end
            end
          end
          if do_trim
            reduce = seq.size % 3
            seq = seq[0..(seq.size-1-reduce)] if reduce != 0
          end
          if @options[:validate]
            ntseq = Bio::Sequence::NA.new(seq)
            aaseq = ntseq.translate
            raise 'Validate translation problem '+rec0.id+"\n"+seq if aaseq.count('*') > 1
          end

          retval = seq
          retval
        end

        # Patch a sequence together from a Sequence string and an array
        # of records and translate in the correct direction and frame. The options 
        # are the same as for +assemble+, except :trim defaults to true.
        def assembleAA sequence, startpos, reclist, options = { :phase=>false, :reverse=>true, :trim=>true, :complement=>true }
          seq = assemble(sequence, startpos, reclist, options)
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
