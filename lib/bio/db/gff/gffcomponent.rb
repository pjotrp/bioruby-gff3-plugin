#
# = bio/db/gff/gffcomponent.rb - Assemble mRNA and CDS from GFF
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#
# Fetch information from a GFF file

require 'set'

module Bio
  module GFFbrowser

    module Helpers

      module Record
        include Logger
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

        include Logger

        COMPONENT_TYPES = Set.new(%w{
          gene SO:0000704 contig transcript Component region
        })
 
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
          $stderr.print @componentlist
          raise "Could not find <#{search}> container/component for #{Record::formatID(rec)}"
        end
      end

      module Gff3Features

        # Ignore the following features (case sensitive?)
        IGNORE_FEATURES = Gff3Component::COMPONENT_TYPES + Set.new(%w{
          transposon Match similarity UTR hsp Peptide match
          TF_binding_site intronSO:0000188 polyA_sequence SO:0000610
          polyA_site SO:0000553
          five_prime_UTR SO:0000204 three_prime_UTR SO:0000205
          exon SO:0000147
        })
      end

    end
  end
end
