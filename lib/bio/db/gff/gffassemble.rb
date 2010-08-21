#
# = bio/db/gff/gffassemble.rb - Assemble mRNA and CDS from GFF
#
# Copyright::  Copyright (C) 2010
#              Pjotr Prins <pjotr.prins@thebird.nl>
# License::    The Ruby License
#
# Fetch information from a GFF file

module Bio
  class GFF

    class Counter < Hash
      def add id
        self[id] = 0 if self[id] == nil
        self[id] += 1
      end
    end
   
    # Helper class for storing linked records based on a shared ID
    class LinkedRecs < Hash
      def add id, rec
        puts "Adding <#{id}>"
        self[id] = [] if self[id] == nil
        self[id] << rec
      end
    end
   
    # mRNA mixin for fetching GFF info
    module MRNA
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
            parent = rec.get_attribute('Parent')
            case rec.feature_type
              when 'gene' || 'SO:0000704' : genes << id 
              when 'mRNA' || 'SO:0000234' : mrnas.add(parent,rec) if parent
              when 'CDS'  || 'SO:0000316' : cdss.add(parent,rec) if parent
              else ignored_features[rec.feature_type] = true
            end
          end
        end
        puts "---- Yield each mRNA"
        genes.each do | id |
          # result = mrnas[id]
          yield id,mrnas[id]
        end
        cdss.each do | k, v |
          yield k,v
        end
        puts "---- Display features that have no match"
        p ignored_features.keys
      end
    end
  end
end
