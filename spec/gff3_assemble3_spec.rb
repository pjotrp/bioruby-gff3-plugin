# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   ruby -I ../bioruby/lib/ ~/.gems/bin/spec spec/gff3_assemble3_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio/db/gff/gffdb'

include Bio::GFFbrowser

GFF3FILE3="test/data/gff/test-cds.gff3"

describe GFFdb, "Assemble CDS (extra checks)" do
  before :all do 
    gffdb = Bio::GFFbrowser::GFFdb.new(GFF3FILE3)
    @gff = gffdb.assembler
    @gff.parse
  end

  it "should translate gene MhA1_Contig1040.frz3.gene29" do
    @contigsequence = @gff.sequencelist["MhA1_Contig1040"]
    @componentlist = {}
    @cdslist = {}
    @gff.each_CDS do | id, reclist, component |
      @componentlist[id] = component
      @cdslist[id] = reclist
    end
    name = "cds:MhA1_Contig1040.frz3.gene"
    recs = @cdslist[name]
    component = @componentlist[name]
    p recs
    cds0 = recs[0]
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:raw=>true)
    ntseq.should == "TTAATTAATTTGCCTAGAAAAACAAAGGCATAACATGCTTGCAGTCATCATACGGTAAGAGAGAAACCAACGATATGTTAATAATGTTGATGGGGGAATATCCTCATTAGAATTCTTTTTTGGGTGAATTGAAATTGCCATATTATTAGTATTATTAGAAAATATTAAATTTGTTGATAA"
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:codonize=>true)
    ntseq.should == "TTATCAACAAATTTAATATTTTCTAATAATACTAATAATATGGCAATTTCAATTCACCCAAAAAAGAATTCTAATGAGGATATTCCCCCATCAACATTATTAACATATCGTTGGTTTCTCTCTTACCGTATGATGACTGCAAGCATGTTATGCCTTTGTTTTTCTAGGCAAATTAATTAA"
    aaseq = @gff.assembleAA(@contigsequence,component.start,recs)
    aaseq.should == "LSTNLIFSNNTNNMAISIHPKKNSNEDIPPSTLLTYRWFLSYRMMTASMLCLCFSRQIN*"
  end
  it "should translate gene MhA1_Contig2992.frz3.gene1" do
    @contigsequence = @gff.sequencelist["MhA1_Contig2992"]
    @componentlist = {}
    @cdslist = {}
    @gff.each_CDS do | id, reclist, component |
      @componentlist[id] = component
      @cdslist[id] = reclist
    end
    name = "cds:MhA1_Contig2992.frz3.gene1"
    recs = @cdslist[name]
    component = @componentlist[name]
    p recs
    cds0 = recs[0]
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:raw=>true)
    ntseq.should == ""
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:codonize=>true)
    ntseq.should == ""
    aaseq = @gff.assembleAA(@contigsequence,component.start,recs)
    aaseq.should == ""
  end
end



