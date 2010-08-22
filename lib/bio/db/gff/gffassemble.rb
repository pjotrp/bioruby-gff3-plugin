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
        def add rec
          parent = rec.get_attribute('Parent')
          if parent
            id = rec.id
            puts "Adding <#{id}>"
            self[id] = [] if self[id] == nil
            self[id] << rec
          else
            warn "record has no parent",id
          end
        end

        def validate
          each do | id, rec |
            sections = []
            rec.each do | section |
              sections.push Section.new(section)
            end
            sections.sort.each_with_index do | check, i |
              neighbour = sections[i+1]
              if neighbour and check.intersection(neighbour)
                warn "Overlapping sections for ",id
              end
            end
          end
        end
      end

      class Section < Range
        def initialize rec
          super(rec.start,rec.end)
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
    end # Helpers
     
    module MRNA
      include Helpers
      include Helpers::Error

      IGNORE_FEATURES = %w{TF_binding_site intronSO:0000188 polyA_sequence SO:0000610 polyA_site SO:0000553 five_prime_UTR SO:0000204 three_prime_UTR SO:0000205}

      # Digest mRNA from the GFFdb and store in Hash
      # Next yield(id, seq) from Hash
      def each_mRNA
        puts "---- Digest DB and store data in mRNA Hash"
        ids   = Counter.new   # Count ids
        seqnames = Counter.new   # Count seqnames
        ignored_features = {}
        genes = []   # Store gene ids
        mrnas = LinkedRecs.new   # Store linked mRNA records
        cdss  = LinkedRecs.new
        @gffs.each do | gff |
          gff.records.each do | rec |
            seqnames.add(rec.seqname)
            id = rec.id
            ids.add(id)
            case rec.feature_type
              when 'gene' || 'SO:0000704' : genes << id 
              when 'mRNA' || 'SO:0000234' : mrnas.add(rec) 
              when 'CDS'  || 'SO:0000316' : cdss.add(rec)
              when 'exon' || 'SO:0000147' : false
              else
                if !IGNORE_FEATURES.include?(rec.feature_type)
                  ignored_features[rec.feature_type] = true
                end
            end
          end
        end
        # validate CDS sections do not overlap
        cdss.validate
        puts "---- Yield each mRNA"
        mrnas.each do | k, v |
          yield k,v
        end
        cdss.each do | k, v |
          yield k,v
        end
 
        ignored_features.keys.each do | k |
          warn "Feature has no match",k
        end
      end
    end
  end
end
