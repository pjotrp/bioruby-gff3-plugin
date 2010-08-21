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
       ids   = {}   # Count ids
       ignored_features = {}
       genes = []   # Store gene ids
       mrnas = LinkedRecs.new   # Store array of mRNA records
       @gffs.each do | gff |
         gff.records.each do | rec |
           id = rec.id
           if ['gene','mRNA'].index(rec.feature_type) != nil
             ids[id] = 0 if ids[id] == nil
             ids[id] += 1
             genes << id if rec.feature_type =~ /gene/i
             parent = rec.get_attribute('Parent')
             mrnas.add(parent,rec) if parent
           else
             ignored_features[rec.feature_type] = true
           end
         end
       end
       puts "---- Yield each mRNA"
       genes.each do | id |
         yield id,mrnas[id]
       end
       puts "---- Display features that have no match"
       p ignored_features.keys
     end
   end
  end
end
