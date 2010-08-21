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

   # mRNA mixin for fetching GFF info
   module MRNA
     def each_mRNA
     end
   end
  end
end
