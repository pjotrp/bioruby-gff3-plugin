# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   ruby -I ../bioruby/lib/ ~/.gems/bin/spec spec/gff3_assemble_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio/db/gff/gffdb'

include Bio::GFFbrowser

FASTAFILE="test/data/gff/MhA1_Contig1133.fa"  
GFF3FILE="test/data/gff/MhA1_Contig1133.gff3"

describe GFFdb, "Assemble CDS" do
  before :all do 
    gffdb = Bio::GFFbrowser::GFFdb.new(GFF3FILE, :fasta_filename => FASTAFILE)
    @gff = gffdb.assembler
  end

  it "should have the single contig" do 
    @gff.parse
    @gff.sequencelist.size.should == 1
    @gff.sequencelist["MhA1_Contig1133"].should_not == nil
    @gff.sequencelist["MhA1_Contig1133"].size.should == 33905
  end
end

