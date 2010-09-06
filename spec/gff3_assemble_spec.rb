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
    @gff.parse
    @cdslist = {}
    @gff.each_CDS do | id, reclist, container |
      @cdslist[id] = reclist
    end
  end

  it "should have the single contig" do 
    @gff.sequencelist.size.should == 1
    @gff.sequencelist["MhA1_Contig1133"].should_not == nil
    @gff.sequencelist["MhA1_Contig1133"].size.should == 33905
  end

  it "should have CDS 7838:7980" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds0 = recs[0]
    cds0.start.should == 7838
    cds0.end.should == 7980
  end
  it "should translate CDS 7838:7980"
  it "should have CDS 8065:8308" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds1 = recs[1]
    cds1.start.should == 8065
    cds1.end.should == 8308
  end
  it "should translate CDS 8065:8308"
  it "should assemble 3 CDSs for MhA1_Contig1133.frz3.gene4"
end

